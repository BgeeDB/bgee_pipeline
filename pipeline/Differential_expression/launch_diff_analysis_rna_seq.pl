#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Sebastien Moretti, created 2014-10-29
# launch the differential expression analysis for RNA Seq
# USAGE: perl launch_diff_analysis_rna_seq.pl
#        perl launch_diff_analysis_rna_seq.pl <EXP_ID>
###################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

$| = 1;

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug) = (0);
my ($path_generes, $ensRelease) = ('', 0);
my ($path_target, $path_processed) = ('', '');
my %opts = ('debug'             => \$debug,            # more verbose
            'bgee=s'            => \$bgee_connector,   # Bgee connector string
            'path_generes=s'    => \$path_generes,     # rna_seq/gene_results/
            'ensRelease=i'      => \$ensRelease,       # Ensembl version used for RNA Seq mapping
            'path_target=s'     => \$path_target,      # name files for replicates rna_seq/bioconductor/targets/
            'path_processed=s'  => \$path_processed,   # final result dir rna_seq/processed_differential/
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $path_generes eq '' || $ensRelease == 0 || $path_target eq '' || $path_processed eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -path_generes=\$(RNASEQGENERES) -ensRelease=\$(ENSRELEASE) -path_target=\$(RNASEQBIOCONDUCTORTARG) -path_processed=\$(RNASEQDIFFEXPRPATH)  <expId>
\t-bgee             Bgee    connector string
\t-path_generes     rna_seq/gene_results/                     directory path
\t-ensRelease       Ensembl version used for RNA Seq mapping
\t-path_target      rna_seq/bioconductor/targets/             directory path
\t-path_processed   rna_seq/processed_differential/           directory path
\t-debug            More verbose
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


my %experiments;
# Monster query to get info on libraries for diff. analyses.
# It makes the union of two queries: one to get info for diff.analyses comparing one organ at different stages,
# one for analyses comparing different organs at a same broad stage.
# TODO DRY, large parts of the query repeated, also in script launch_diff_analysis-Affy.pl
# XXX: actually, the root of the ontology "life stage" (UBERON:0000104) is filtered afterwards,
# it could be filtered right away in the query
my $selExpr = $dbh->prepare(
                            # First query, get info for comparing a same organ at different stages
                            # ("development" comparison factor)
                            '(SELECT DISTINCT rnaSeqLibraryId, rnaSeqExperimentId, rnaSeqPlatformId, libraryType, anatEntityId, '
                                # Specify the comparison factor
                                .'"development" AS comparisonFactor '

                                # Sub-query to map the actually annotated stage to a not too granular stage
                                # (e.g., in human, we don't want to compare organs at stage "83-year old",
                                # but to consider all info at stage "80 year-old and over human stage").
                                .'(SELECT t3.stageId FROM stage AS t3 '
                                   # Use taxon constraints to make sure to get a parent stage valid in the related species
                                   .'INNER JOIN stageTaxonConstraint AS t3bis on t3.stageId = t3bis.stageId '
                                   # Get the stage itself or its parent (left bound - right bound),
                                   .'WHERE t3.stageLeftBound <= t2.stageLeftBound AND t3.stageRightBound >= t2.stageRightBound '
                                   # in the proper species,
                                   .'AND (t3bis.speciesId is null or t3bis.speciesId = t2bis.speciesId) '
                                   # that is not too granular (tooGranular = 0),
                                   # and that is the closest to the annotated stage (left bound desc order limit 1)
                                   .'AND t3.tooGranular = 0 ORDER BY t3.stageLeftBound DESC LIMIT 1) AS fakeStageId, '

                                # Get the species related to this library. We use a GROUP_CONCAT,
                                # it allows to check the related species over all genes of the library,
                                # not by checking only one gene.
                                # TODO we might also check that the species of the library
                                # corresponds to the species of the annotated stage.
                                .'(SELECT GROUP_CONCAT(DISTINCT t4.speciesId SEPARATOR ", ") FROM gene AS t4 INNER
                                   JOIN rnaSeqResult AS t5 ON t4.geneId = t5.geneId WHERE
                                   t5.rnaSeqLibraryId = t1.rnaSeqLibraryId) AS speciesIds

                                FROM rnaSeqLibrary AS t1 '
                                # Join to stage and to taxon constraint tables for the stage sub-query
                                .'INNER JOIN stage AS t2 ON t1.stageId = t2.stageId
                                INNER JOIN stageTaxonConstraint AS t2bis on t2.stageId = t2bis.stageId
                             )

                             UNION ALL '

                            # Second query, get info for comparing different organs at a same broad stage
                            # ("anatomy" comparison factor)
                            .'(SELECT DISTINCT rnaSeqLibraryId, rnaSeqExperimentId, rnaSeqPlatformId, libraryType, anatEntityId, '
                                # Specify the comparison factor
                                .'"anatomy" AS comparisonFactor, '

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
                                .'(SELECT GROUP_CONCAT(DISTINCT t4.speciesId SEPARATOR ", ") FROM gene AS t4 INNER
                                   JOIN rnaSeqResult AS t5 ON t4.geneId = t5.geneId WHERE
                                   t5.rnaSeqLibraryId = t1.rnaSeqLibraryId) AS speciesIds

                                FROM rnaSeqLibrary AS t1 '
                                # Join to stage and to taxon constraint tables for the stage sub-query
                                .'INNER JOIN stage AS t2 ON t1.stageId = t2.stageId
                                INNER JOIN stageTaxonConstraint AS t2bis on t2.stageId = t2bis.stageId
                             )');
$selExpr->execute()  or die $selExpr->errstr;
while ( my @data = $selExpr->fetchrow_array ){
    # Multiple species per library => Problem
    if ( $data[6] =~ /,/ ){
        warn "Problem with [$data[1]]: multiple species in it [$data[6]]\n";
        next;
    }

    if ( !defined $ARGV[0] || (defined $ARGV[0] && $data[1] eq $ARGV[0]) ){
        my $fileName;
        if ( -e $path_generes.'/'.$data[1].'/'.$data[0].'_out_'.$ensRelease && -s $path_generes.'/'.$data[1].'/'.$data[0].'_out_'.$ensRelease ){
            $fileName = $data[0].'_out_'.$ensRelease;
        }
        else {
            die "Error, no processed file found for expId: [$data[1]] - libId: [$data[0]]\n";
        }
        if ( $data[7] eq 'development' ){
            # Development comparison => Single organ comparison
            $experiments{$data[1]}->{$data[6]}->{'a_'.$data[4]}->{$data[5]}->{$fileName}->{'path'}             = $path_generes;
            $experiments{$data[1]}->{$data[6]}->{'a_'.$data[4]}->{$data[5]}->{$fileName}->{'bgeeRNASeqLibId'}  = $data[0];
            $experiments{$data[1]}->{$data[6]}->{'a_'.$data[4]}->{$data[5]}->{$fileName}->{'libraryType'}      = $data[3];
            $experiments{$data[1]}->{$data[6]}->{'a_'.$data[4]}->{$data[5]}->{$fileName}->{'rnaSeqPlatformId'} = $data[2];
        }
        elsif ( $data[7] eq 'anatomy' ){
            # Anatomy comparison => Single stage comparison
            $experiments{$data[1]}->{$data[6]}->{'s_'.$data[5]}->{$data[4]}->{$fileName}->{'path'}             = $path_generes;
            $experiments{$data[1]}->{$data[6]}->{'s_'.$data[5]}->{$data[4]}->{$fileName}->{'bgeeRNASeqLibId'}  = $data[0];
            $experiments{$data[1]}->{$data[6]}->{'s_'.$data[5]}->{$data[4]}->{$fileName}->{'libraryType'}      = $data[3];
            $experiments{$data[1]}->{$data[6]}->{'s_'.$data[5]}->{$data[4]}->{$fileName}->{'rnaSeqPlatformId'} = $data[2];
        }
    }
}
$selExpr->finish;
$dbh->disconnect;

# Remove "life stage" (UBERON:0000104) from analyses because development stage root!
my $to_exclude = 'UBERON:0000104';

my %exp_to_treat;
for my $exp ( keys %experiments ){
    for my $speciesId ( keys %{$experiments{$exp}} ){
        for my $single ( keys %{$experiments{$exp}->{$speciesId}} ){
            if ( $single eq "s_$to_exclude" ){ # Remove "life stage" (UBERON:0000104) from analyses because development stage root!
                delete $exp_to_treat{$exp}->{$speciesId}->{$single};
                next;
            }
            my $count = 0;
            for my $condition ( keys %{$experiments{$exp}->{$speciesId}->{$single}} ){
                # Remove "life stage" (UBERON:0000104) from analyses because development stage root!
                next  if ( $condition eq $to_exclude );
                # if there are some replicates for that condition (organ/stage)
                if ( scalar keys %{$experiments{$exp}->{$speciesId}->{$single}->{$condition}} > 1 ){
                    $exp_to_treat{$exp}->{$speciesId}->{$single}->{$condition} = $experiments{$exp}->{$speciesId}->{$single}->{$condition};
                    $count++;
                }
            }
            # If at least 3 conditions had replicates
            if ( $count < 3 ){
                delete $exp_to_treat{$exp}->{$speciesId}->{$single};
            }
        }
        # Remove speciesId that didn't match the criteria
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


#FIXME for this release we mix read_type
EXP:
for my $exp ( sort keys %exp_to_treat ){
    for my $speciesId ( keys %{$exp_to_treat{$exp}} ){
        if ( !-e $path_processed.'/'.$exp ){
            mkdir "$path_processed/$exp";
        }

        # Create file .targets
        SINGLE:
        for my $single ( keys %{$exp_to_treat{$exp}->{$speciesId}} ){
            my $nbr_conditions = scalar keys %{$exp_to_treat{$exp}->{$speciesId}->{$single}};
            # Skip case with several repeats but low number of conditions
            next SINGLE  if ( $nbr_conditions < 3 );

            my $comparisonFactor = $single =~ /^a_/ ? 'development' : 'anatomy';
            printf("\t%-15s %-15s %-15s %s...", $exp, $speciesId, $comparisonFactor, $single);
            # Target file $exp__$speciesId must contain 1 stage and several organs OR 1 organ for several stages
            # path/exp__speciesId__singleorgan  OR  path/exp__speciesId__singlestage
            my $target = "$path_target/${exp}___${speciesId}__$single.target";
            my $path   = '';
            open(my $TARGET, '>', "$target")  or die 'Cannot open TARGET file';
            print {$TARGET} "#organ\tstage\tfile\tpath\tbgeeRNASeqLibId\tspeciesId\tcomparisonFactor\n";
            for my $condition ( keys %{$exp_to_treat{$exp}->{$speciesId}->{$single}} ){
                for my $file ( keys %{$exp_to_treat{$exp}->{$speciesId}->{$single}->{$condition}} ){
                    my ($organ, $stage) = ('', '');
                    if ( $single =~ /^a_(.+)$/ ){
                        $organ = $1;
                        $stage = $condition;
                    }
                    elsif ( $single =~ /^s_(.+)$/ ){
                        $stage = $1;
                        $organ = $condition;
                    }
                    else {
                        die "Problem with organ/stage for [$exp][$speciesId][$single][$condition]";
                    }
                    print {$TARGET} "$organ\t$stage\t$file\t".
                       $exp_to_treat{$exp}->{$speciesId}->{$single}->{$condition}->{$file}->{'path'}."\t".
                       $exp_to_treat{$exp}->{$speciesId}->{$single}->{$condition}->{$file}->{'bgeeRNASeqLibId'}."\t".
                       "$speciesId\t$comparisonFactor\n";
                }
            }
            close $TARGET;

            # Launch R script
            my $log = "$path_processed/$exp/${exp}___${speciesId}__$single.log";
            my $cmd = "R CMD BATCH --vanilla '--args target_file_path=\"$target\" output_folder_path=\"$path_processed/$exp\" input_folder_path=\"$path_generes/$exp\" speciesID=\"$speciesId\"' $FindBin::Bin/diff_analysis_rna_seq.R  $log";
            system($cmd)==0  or do{warn "\tsystem [$cmd] failed: $?\n"; map { system("mv $_ $_.PROB") } glob("$path_processed/$exp/*.out");};
            print "\n";
        }
    }
}

exit 0;

