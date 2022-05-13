#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Wollbrett, created 2018-01-09
# Based on ./launch_diff_analysis_rna_seq.pl written by Sebastien Moretti
# launch the sex differential expression analysis for RNA Seq
# USAGE: perl launch_sex_diff_analysis_rna_seq.pl
#        perl launch_sex_diff_analysis_rna_seq.pl <EXP_ID>
###################################################

use Getopt::Long;
use File::Basename;
use File::Slurp;
use File::Path qw/make_path/;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Data::Dumper;

$| = 1;

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug) = (0);
my ($path_generes, $abundance_file) = ('', '');
my ($path_library_worm, $path_library) = ('', '');
my ($path_target, $path_processed) = ('', '');
my %opts = ('debug'                 => \$debug,            # more verbose
            'bgee=s'                => \$bgee_connector,   # Bgee connector string
            'path_generes=s'        => \$path_generes,     # rna_seq/all_results_bgee_v15/
            'path_library=s'        => \$path_library,     # bgee_pipeline/source_files/RNA_Seq/RNASeqLibrary.tsv
            'path_library_worm=s'   => \$path_library_worm,# bgee_pipeline/source_files/RNA_Seq/RNASeqLibrary_worm.tsv
            'abundance_file=s'      => \$abundance_file,   # abundance_gene_level+new_tpm+new_fpkm+calls.tsv
            'path_target=s'         => \$path_target,      # name files for replicates rna_seq/bioconductor_bgee_v15/targets/sex/
            'path_processed=s'      => \$path_processed,   # final result dir rna_seq/processed_differential_bgee_v15/sex/
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $path_generes eq '' || $abundance_file eq '' || $path_target eq '' 
	|| $path_processed eq '' || $path_library_worm eq '' || $path_library eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -path_generes=\$(RNASEQALLRES) -abundance_file=\$(RNASEQABUNDANCEFILE) -path_library=\$(RNASEQ_LIB_FILEPATH) -path_library_worm=\$(RNASEQ_LIB_FILEPATH_WORM) -path_target=\$(RNASEQBIOCONDUCTORTARG_SEX) -path_processed=\$(RNASEQDIFFEXPRPATH_SEX)  <expId>
\t-bgee                 Bgee    connector string
\t-path_generes         rna_seq/all_results_bgee_v15/                    directory path
\t-path_library         bgee_pipeline/source_files/RNA_Seq/RNASeqLibrary.tsv
\t-path_library_worm    bgee_pipeline/source_files/RNA_Seq/RNASeqLibrary_worm.tsv
\t-abundance_file       abundance_gene_level+new_tpm+new_fpkm+calls.tsv  file name
\t-path_target          rna_seq/bioconductor_bgee_v15/targets/sex        directory path
\t-path_processed       rna_seq/processed_differential_bgee_v15/sex      directory path
\t-debug                More verbose
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

# Remove "life stage" (UBERON:0000104) from analyses because development stage root!
my $to_exclude = 'UBERON:0000104';

# Define that the comparison factor of the differential expression is the sex
my $comparisonFactor = 'sex';

my %experiments;
# Monster query to get info on libraries for diff. analyses.
# It get info for diff.analyses comparing different organs at a same broad stage and a same sex.
print "\tprepare query\n";
my $selExpr = $dbh->prepare( 
                            # Second query, get info for comparing different organs at a same broad stage
                            # ("anatomy" comparison factor)
                            'SELECT DISTINCT t1.rnaSeqLibraryId, t1.rnaSeqExperimentId, t1.rnaSeqPlatformId, '
                            .'t1.libraryType, t7.anatEntityId, t7.sex, '

                                # Sub-query to map the actually annotated stage to a broad developmental stage
                                # (e.g., in human, we don't want to compare different organs
                                # at stage "Carnegie stage 10", because of heterochrony,
                                # but at stage "organogenesis").
                                .'(SELECT t3.stageId FROM stage AS t3 '
                                   # Use taxon constraints to make sure to get a parent stage valid in the related species
                                   .'INNER JOIN stageTaxonConstraint AS t3bis on t3.stageId = t3bis.stageId '
                                   # Get the stage itself or its parent (left bound - right bound),
                                   .'WHERE t3.stageLeftBound <= t2.stageLeftBound AND t3.stageRightBound >= t2.stageRightBound '
                                   # in the proper species,
                                   .'AND (t3bis.speciesId is null or t3bis.speciesId = t2bis.speciesId) '
                                   # that is broad enough to group organs (groupingStage = 1),
                                   # and that is the closest to the annotated stage (left bound desc order limit 1)
                                   .'AND t3.groupingStage = 1 ORDER BY t3.stageLeftBound DESC LIMIT 1) AS fakeStageId, '
																
                                # Get the species related to this library. We use a GROUP_CONCAT,
                                # it allows to check the related species over all genes of the library,
                                # not by checking only one gene.
                                # TODO we might also check that the species of the library
                                # corresponds to the species of the annotated stage.
                                .'(SELECT GROUP_CONCAT(DISTINCT t4.speciesId SEPARATOR ", ") FROM rnaSeqResult AS t5 STRAIGHT_JOIN '
                                   .'gene AS t4 ON t4.bgeeGeneId = t5.bgeeGeneId WHERE '
                                   .'t5.rnaSeqLibraryId = t1.rnaSeqLibraryId) AS speciesIds '

                                .'FROM rnaSeqLibrary AS t1 '
                                # Join to stage and to taxon constraint tables for the stage sub-query
                                .'INNER JOIN cond AS t6 ON t1.conditionId = t6.conditionId '
                                .'INNER JOIN cond AS t7 ON t6.exprMappedConditionId = t7.conditionId '
                                .'INNER JOIN stage AS t2 ON t6.stageId = t2.stageId '
                                .'INNER JOIN stageTaxonConstraint AS t2bis on t2.stageId = t2bis.stageId '
                                # Don't take into account 'not annotated','hermaphrodite', 'mixed', or'NA' sexes
                                .'WHERE t6.sex IN (\'male\', \'female\') '
                                .'AND t6.stageId != \''.$to_exclude.'\'');


$selExpr->execute()  or die $selExpr->errstr;
while ( my @data = $selExpr->fetchrow_array ){
	# $data[0] => rnaSeqLibraryId
	# $data[1] => rnaSeqExperimentId
	# $data[2] => rnaSeqPlatformId
	# $data[3] => libraryType
	# $data[4] => anatEntityId
	# $data[5] => sex
	# $data[6] => fakeStageId
	# $data[7] => speciesId
    # if multiple species per library => Problem
    if ( $data[7] =~ /,/ ){
        warn "Problem with [$data[0]]: multiple species in it [$data[7]]\n";
        next;
    }

    if ( !defined $ARGV[0] || (defined $ARGV[0] && $data[1] eq $ARGV[0]) ){
		#die if no processed file is found
        if ( !(-e $path_generes.$data[0].'/'.$abundance_file && -s $path_generes.$data[0].'/'.$abundance_file )){
             die "Error, no processed file found for expId: [$data[1]] - libId: [$data[0]]\n";
        }
        # sex comparison => same anat. entity, dev stage AND sex
        $experiments{$data[1]}->{$data[7]}->{$data[6]}->{$data[4]}->{$data[5]}->{$data[0]}->{'libraryType'}      = $data[3];
        $experiments{$data[1]}->{$data[7]}->{$data[6]}->{$data[4]}->{$data[5]}->{$data[0]}->{'rnaSeqPlatformId'} = $data[2];
    }
}
$selExpr->finish;
$dbh->disconnect;

my %exp_to_treat;
for my $exp ( keys %experiments ){
    for my $speciesId ( keys %{$experiments{$exp}} ){
        for my $fakeStageId ( keys %{$experiments{$exp}->{$speciesId}} ){
           	for my $anatEntity (keys %{$experiments{$exp}->{$speciesId}->{$fakeStageId}} ){
            	my $countSex = 0;
                for my $sex ( keys %{$experiments{$exp}->{$speciesId}->{$fakeStageId}->{$anatEntity}} ){
                	# if there are some replicates for that condition (organ/stage/sex)
                	if ( scalar keys %{$experiments{$exp}->{$speciesId}->{$fakeStageId}->{$anatEntity}->{$sex}} > 1 ){
                    	$exp_to_treat{$exp}->{$speciesId}->{$fakeStageId}->{$anatEntity}->{$sex} = $experiments{$exp}->{$speciesId}->{$fakeStageId}->{$anatEntity}->{$sex};
                    	$countSex++;
                	}
            	}
            	# if only one sex (male/female) has replicates for that condition (organ/stage)
            	if ($countSex == 1){
                    print "remove anatEntity when only one sex $exp, $speciesId, $fakeStageId, $anatEntity : sexeNumber = $countSex\n" if($debug);
                    delete $exp_to_treat{$exp}->{$speciesId}->{$fakeStageId}->{$anatEntity};
            	}
            }
            # Remove stage that didn't match the criteria
    		if ( scalar keys %{$exp_to_treat{$exp}->{$speciesId}->{$fakeStageId}} == 0 ){
        		delete $exp_to_treat{$exp}->{$speciesId}->{$fakeStageId};
    		}
        }
        # Remove species that didn't match the criteria
        if ( scalar keys %{$exp_to_treat{$exp}->{$speciesId}} == 0 ){
            delete $exp_to_treat{$exp}->{$speciesId};
        }
    }
    # Remove experiments that didn't match the criteria
    if ( scalar keys %{$exp_to_treat{$exp}} == 0 ){
        delete $exp_to_treat{$exp};
    }
}

if ( scalar keys %exp_to_treat == 0 ){
    die "\tProblem! No experiment to analyze!\n"
}

# get all technical replicates (see get_replicates function)
my %replicates = ( get_replicates($path_library), get_replicates($path_library_worm) );

# filter useful technical replicates and store them in an easily usable hash
# @technical_replicates{$expId}->{$anatEntityId}->{$stageId}->{$sex}
my %technical_replicates;
for my $expId ( sort keys %replicates ){
	for my $anatEntityId ( sort keys %{$replicates{$expId}} ){
		for my $stageId ( sort keys %{$replicates{$expId}->{$anatEntityId}} ){
			for my $sex ( sort keys %{$replicates{$expId}->{$anatEntityId}->{$stageId}} ){
				next if ( $sex ne "\"M\"" &&  $sex ne "\"F\"");
				for my $strain ( sort keys %{$replicates{$expId}->{$anatEntityId}->{$stageId}->{$sex}} ){
					for my $platform ( sort keys %{$replicates{$expId}->{$anatEntityId}->{$stageId}->{$sex}->{$strain}} ){
						for my $replicateNumber ( sort keys %{$replicates{$expId}->{$anatEntityId}->{$stageId}->{$sex}->{$strain}->{$platform}} ){
							# if more than one technical replicate
							if ( scalar keys @{$replicates{$expId}->{$anatEntityId}->{$stageId}->{$sex}->{$strain}->{$platform}->{$replicateNumber}} > 1 ){
								if($sex eq "\"M\""){
									$technical_replicates{$expId}->{$anatEntityId}->{$stageId}->{'male'}->{'alreadyUsed'} = 0;
									$technical_replicates{$expId}->{$anatEntityId}->{$stageId}->{'male'}->{'libIds'} = $replicates{$expId}->{$anatEntityId}->{$stageId}->{$sex}->{$strain}->{$platform}->{$replicateNumber};
									
								}else{
									$technical_replicates{$expId}->{$anatEntityId}->{$stageId}->{'female'}->{'alreadyUsed'} = 0;
									$technical_replicates{$expId}->{$anatEntityId}->{$stageId}->{'female'}->{'libIds'} = $replicates{$expId}->{$anatEntityId}->{$stageId}->{$sex}->{$strain}->{$platform}->{$replicateNumber};
								}
							}
						}
					}
				}
			}
		}
	}
}

# small test to be sure that technical replicates are taken into account
# $technical_replicates{'GSE44612'}->{'UBERON:0007023'}->{'UBERON:0000066'}->{'male'}->{'alreadyUsed'} = 0;
# $technical_replicates{'GSE44612'}->{'UBERON:0007023'}->{'UBERON:0000066'}->{'male'}->{'libIds'} = ['SRX246998','SRX246999','SRX018867','SRX018865'];

print Dumper(\%technical_replicates) if($debug);
print Dumper(\%exp_to_treat) if($debug);

								

if ( !-e "$path_target" ){
    make_path("$path_target");
}

EXP:
for my $exp ( sort keys %exp_to_treat ){
    for my $speciesId ( keys %{$exp_to_treat{$exp}} ){
        if ( !-e "$path_processed$exp" ){
            make_path("$path_processed$exp");
        }
        # Create file .targets
        SINGLE:
        for my $stageId ( keys %{$exp_to_treat{$exp}->{$speciesId}} ){
        	my $stageIdForFile = $stageId;
        	$stageIdForFile =~ tr/:/_/;
        	for my $anatEntityId ( keys %{$exp_to_treat{$exp}->{$speciesId}->{$stageId}} ){
        		my $anatEntityIdForFile = $anatEntityId;
        		$anatEntityIdForFile =~ tr/:/_/;
            	printf("\t%-15s %-15s %-15s %s...", $exp, $speciesId, $stageId, $anatEntityId);
            	my $target = "${path_target}${exp}__${speciesId}__${stageIdForFile}__${anatEntityIdForFile}.target";
            	my $path   = '';
            	open(my $TARGET, '>', "$target")  or die 'Cannot open TARGET file';
            	print {$TARGET} "#organ\tstage\tsex\tbgeeRNASeqLibId\n";
            	for my $sex ( keys %{$exp_to_treat{$exp}->{$speciesId}->{$stageId}->{$anatEntityId}} ){
                    for my $bgeeRNASeqLibId ( keys %{$exp_to_treat{$exp}->{$speciesId}->{$stageId}->{$anatEntityId}->{$sex}} ){
                    	if ( exists $technical_replicates{$exp}->{$anatEntityId}->{$stageId}->{$sex} ) {
                    		my @libIds = @{$technical_replicates{$exp}->{$anatEntityId}->{$stageId}->{$sex}->{'libIds'}};
                    		print "@libIds\n";
                    		if ( grep( /^$bgeeRNASeqLibId$/, @libIds ) ) {
                     			print "$bgeeRNASeqLibId\n";
  								if ( $technical_replicates{$exp}->{$anatEntityId}->{$stageId}->{$sex}->{'alreadyUsed'} ){
  									print "library $bgeeRNASeqLibId is a technical replicate that will not be taken into account\n";
  								}else{
  									$technical_replicates{$exp}->{$anatEntityId}->{$stageId}->{$sex}->{'alreadyUsed'} = 1;
	   								print {$TARGET} "$anatEntityId\t$stageId\t$sex\t$bgeeRNASeqLibId\n";
  								}
							}
						}else{
                   			print {$TARGET} "$anatEntityId\t$stageId\t$sex\t$bgeeRNASeqLibId\n";
                   		}
                    }
                }
	           	close $TARGET;

            	# Launch R script
                my $log = "$path_processed$exp/${speciesId}__sex__${anatEntityIdForFile}__${stageIdForFile}.log";
                my $cmd = "R CMD BATCH --vanilla '--args target_file_path=\"$target\" output_folder_path=\"$path_processed$exp\" input_folder_path=\"$path_generes\" abundance_file_name=\"$abundance_file\" speciesID=\"$speciesId\"' $FindBin::Bin/diff_sex_analysis_rna_seq.R  $log";
                system($cmd)==0  or do{warn "\tsystem [$cmd] failed: $?\n"; map { system("mv $_ $_.PROB") } glob("$path_processed$exp/*.out");};
                print "\n";
            }
        }
    }
}

exit 0;

# Read annotation files to retrieve technical replicates libraries. 
# Technical replicates are not supposed to be taken into account in differential epression analyses.
# These information will then be used to remove libraries before the DE analyses.
# @replicates{$expId}->{$anatEntityId}->{$stageId}->{$sex}->{$strain}->{$platform}->{$replicateNumber}
# The array contains all technical replicate libIds 
sub get_replicates {
	my $path_library = $_[0];
	my %replicates;
	for my $line ( read_file($path_library) ){
		next  if ( $line =~ /^#/ || $line =~ /^\"#/);
		my @columns = split(/\t/, $line);
		# test presence of replicate number
		if( $columns[25] && length $columns[25]){
			$columns[1] =~ s/"//g;
			$columns[5] =~ s/"//g;
			$columns[7] =~ s/"//g;
			$columns[0] =~ s/"//g;
			push ( @{$replicates{$columns[1]}->{$columns[5]}->{$columns[7]}->{$columns[19]}->{$columns[20]}->{$columns[2]}->{$columns[25]}}, $columns[0]);
		}
	}
    return %replicates;
}
