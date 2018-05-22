#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Wollbrett, created 25/01/2018
# USAGE: perl generate_files_sex_diff_rnaseq.pl <speciesId>
############################################

use Getopt::Long;
use File::Basename;
use File::Slurp;
use File::Path qw/make_path/;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Data::Dumper;
use IO::Compress::Gzip qw(gzip $GzipError) ;

$| = 1;

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug) = (0);
my ($path_target, $path_processed, $path_download_file) = ('', '', '');
my %opts = ('debug'             	=> \$debug,            		# more verbose
            'bgee=s'            	=> \$bgee_connector,   		# Bgee connector string
            'path_target=s'     	=> \$path_target,      		# path to the .target files
            'path_processed=s'  	=> \$path_processed,   		# path to the .out files
            'path_download_file=s'  => \$path_download_file,   	# path to the repository where download files should be written
           );
           
# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || ($#ARGV ne -1 && $#ARGV ne 1) || $path_target eq '' || $path_processed eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -path_target=\$(RNASEQBIOCONDUCTORTARG_DEVANDANAT) -path_processed=\$(RNASEQDIFFEXPRPATH_DEVANDANAT) -path_download_file=\$(RNASEQDIFFEXPRPATH_DEVANDANAT)  <speciesId>
\t-bgee             	Bgee connector string
\t-path_target      	rna_seq/bioconductor_bgee_v14/targets/sex				directory path
\t-path_processed   	rna_seq/processed_differential_bgee_v14/sex  			directory path
\t-path_download_file	rna_seq/processed_differential_bgee_v14/sex			 	directory path
\t-debug            	More verbose
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

###################################################################
# Load Anat. Entities and Stages
###################################################################

print 'Start loading anatomical entities and developmental stages from database... \n'  if ( $debug );

my (%anatEntities, %devStages, %species);

my $retrieveAnatSql		= $bgee->prepare('SELECT DISTINCT t2.anatEntityId, t3.anatEntityName from rnaSeqLibrary AS t1 INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId INNER JOIN anatEntity AS t3 ON t2.anatEntityId = t3.anatEntityId');
my $retrieveStagesSql   = $bgee->prepare('SELECT DISTINCT t2.stageId, t3.stageName from rnaSeqLibrary AS t1 INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId INNER JOIN stage AS t3 ON t2.stageId = t3.stageId where t3.groupingStage = 1');
my $retrieveSpeciesSql  = $bgee->prepare('SELECT speciesId, genus, species from species');

#retrieve Anat. Entities
$retrieveAnatSql->execute()  or die $retrieveAnatSql->errstr;
while ( my @data = $retrieveAnatSql->fetchrow_array ){
    $anatEntities{$data[0]} = $data[1];
}
$retrieveAnatSql->finish;

#retrieve Anat. Entities
$retrieveStagesSql->execute()  or die $retrieveStagesSql->errstr;
while ( my @data = $retrieveStagesSql->fetchrow_array ){
    $devStages{$data[0]} = $data[1];
}
$retrieveStagesSql->finish;

#retrieve Species
$retrieveSpeciesSql->execute()  or die $retrieveSpeciesSql->errstr;
while ( my @data = $retrieveSpeciesSql->fetchrow_array ){
    $species{$data[0]} = "$data[1]_$data[2]";
}
$retrieveSpeciesSql->finish;

###################################################################
# Read .target and .out files
###################################################################

# Read target files in order to retrieve all species having sex dif. expression
# This step is mandatory to create files by species
my %speciesWithDifExp; 
for my $file ( glob($path_target.'/*.target') ){
    $file = basename($file);
    # @speciesWithDifExp{$speciesId} = @files
    my $isNeededSpecies = 1;
     if ( $file =~ /^.+?__(.+?)__UBERON_.+?__UBERON_.+?.target$/ ){
     	if ( defined $ARGV[0] ){
            if ( $ARGV[0] ne $1 ){
                $isNeededSpecies = 0;
            }
        }
        if ($isNeededSpecies) {
        	${$speciesWithDifExp{$1}->{$file}} = ();
#         	$speciesWithDifExp{$1} = ();
#         	push(@{$speciesWithDifExp{$1}}, $file);
        }
    }else {
        die "File name [$file] does not correspond to the proper file name pattern\n";
    }
}

     	

# parse .target and .out file and store every useful information into %diffExp
my %diffExp;
for my $speciesId ( keys %speciesWithDifExp ){
	# Retrieve all genes of this species with their names
	my %genes;
	my $retrieveGenesSql = $bgee->prepare('SELECT geneId, geneName FROM gene where speciesId = ?');
	$retrieveGenesSql->execute($speciesId)  or die $retrieveAnatSql->errstr;
	while ( my @data = $retrieveGenesSql->fetchrow_array ){
    	$genes{$data[0]} = $data[1];
	}
	$retrieveGenesSql->finish;
	print "target files for species $speciesId:\n";
	for my $file ( keys %{$speciesWithDifExp{$speciesId}} ){
		print "file : $file\n";# if ($debug);
    	# ${exp}__${speciesId}__{$tage}__{$anatEntity}.target
    	# SRP012682__9606__UBERON_0000113__UBERON_0035805.target
    	# $exp_diff{$speciesId}->{$anatEntityId}->{$stageId}->{$exp}->{$gene}->{'logFC'}
    	#                      ->{'name'}       ->{'name'}  ->{'name'}       ->{'pvalue'}
    	#                                                                    ->{'presentCall'}
    	#                                                                    ->{'name'}
		$file =~ /^(.+?)__.+?__(UBERON_.+?)__(UBERON_.+?).target$/ || die "File name [$file] does not correspond to the proper file name pattern\n";
		my ($exp, $stageId, $anatEntityId) = ($1, $2, $3);
		# Load genes, logFC, pvalue, and present calls from .out files
       	my $outFile = "$path_processed$exp/${speciesId}_sex_${anatEntityId}_${stageId}.out";
       	my $outFileWithPotentialProblems = "$path_processed$exp/${speciesId}_sex_${anatEntityId}_${stageId}.out.PROB";
       	if ( -e $outFileWithPotentialProblems && -s $outFileWithPotentialProblems ){
       		$outFile = $outFileWithPotentialProblems;
       		print "The file $outFile.PROB has been found but it could be corrupted";
       	}elsif ( !(-e $outFile && -s $outFile) ){
			die "Error, processed_differential file not found: [$outFile]\n";
		}
		$anatEntityId =~ tr/_/:/;
		$stageId =~ tr/_/:/;
       	${$diffExp{$speciesId}->{'name'}} = $species{$speciesId};
      	${$diffExp{$speciesId}->{$anatEntityId}->{'name'}} = $anatEntities{$anatEntityId};
       	${$diffExp{$speciesId}->{$anatEntityId}->{$stageId}->{'name'}} = $devStages{$stageId};
       	for my $line ( read_file($outFile, chomp=>1) ){
       		next  if ( $line =~ /^#/ );
            my @tmp = split(/\t/, $line);
            ${$diffExp{$speciesId}->{$anatEntityId}->{$stageId}->{$exp}->{$tmp[0]}->{'name'}} = $genes{$tmp[0]};
            ${$diffExp{$speciesId}->{$anatEntityId}->{$stageId}->{$exp}->{$tmp[0]}->{'pvalue'}} = $tmp[1];
            ${$diffExp{$speciesId}->{$anatEntityId}->{$stageId}->{$exp}->{$tmp[0]}->{'logFC'}} = $tmp[2];
            ${$diffExp{$speciesId}->{$anatEntityId}->{$stageId}->{$exp}->{$tmp[0]}->{'presentCall'}} = $tmp[3];
       	}
    }
}

$bgee->disconnect;

###################################################################
# Write differential expression download files
###################################################################

if ( !-e "$path_download_file" ){
    make_path("$path_download_file");
}

print Dumper(\%diffExp) if ($debug);

for my $speciesId ( keys %diffExp ){
	print "species : $speciesId\n";
	my $speciesName = ${$diffExp{$speciesId}->{'name'}};
	my $downloadFile = "${path_download_file}${speciesName}_diffexpr-sex-single.tsv";
	open(my $DOWNLOADFILE, '>', "$downloadFile")  or die 'Cannot open DOWNLOADFILE file';
	print {$DOWNLOADFILE} "Experiment_ID\tGene_ID\tGene_name\tAnatomical_entity_ID\tAnatomical_entity_name\tDevelopmental_stage_ID\tDevelopmental_stage_name\tlogFC_Male-Female\tpvalue\n";
	for my $anatEntityId ( keys %{$diffExp{$speciesId}} ){
		if( $anatEntityId ne "name" ){
			print "AE : $anatEntityId\n";		
			my $anatEntityName = ${$diffExp{$speciesId}->{$anatEntityId}->{'name'}};
			for my $stageId ( keys %{$diffExp{$speciesId}->{$anatEntityId}} ){
				if( $stageId ne "name"){
					print "stage : $stageId\n";
					my $stageName = ${$diffExp{$speciesId}->{$anatEntityId}->{$stageId}->{'name'}};
					my $numberOfExp = (keys %{ $diffExp{$speciesId}->{$anatEntityId}->{$stageId}}) - 1;
					print "number of exp : $numberOfExp\n";
					if( $numberOfExp > 1){
						warn "We didn't currently implement cases where same species/anatEntity/stage are found in different RNA-Seq experiments ($speciesId, $anatEntityId, $stageId)\n";
					}else{
						for my $expId ( keys %{$diffExp{$speciesId}->{$anatEntityId}->{$stageId}} ){
							if( $expId ne "name" ){
								print "experiment : $expId\n";
								for my $geneId ( keys %{$diffExp{$speciesId}->{$anatEntityId}->{$stageId}->{$expId}} ){
									my $logFC = ${$diffExp{$speciesId}->{$anatEntityId}->{$stageId}->{$expId}->{$geneId}->{'logFC'}};
									my $geneName = ${$diffExp{$speciesId}->{$anatEntityId}->{$stageId}->{$expId}->{$geneId}->{'name'}};
									my $pvalue = ${$diffExp{$speciesId}->{$anatEntityId}->{$stageId}->{$expId}->{$geneId}->{'pvalue'}};
									my $presentCall = ${$diffExp{$speciesId}->{$anatEntityId}->{$stageId}->{$expId}->{$geneId}->{'presentCall'}};
									if ( $presentCall ne 'FALSE' ){
       						            print {$DOWNLOADFILE} "$expId\t$geneId\t$geneName\t$anatEntityId\t$anatEntityName\t$stageId\t$stageName\t$logFC\t$pvalue\n";
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
		
exit 0;	


