#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Roux, created 05/05/09
# launch the differential expression analysis for Affymetrix
# USAGE: perl launch_diff_analysis_affy.pl
#        perl launch_diff_analysis_affy.pl <EXP_ID>
###################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

$| = 1;

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug) = (0);
my ($path_mas5, $path_schuster ) = ('', '');
my ($path_target, $path_processed) = ('', '');
my %opts = ('debug'             => \$debug,            # more verbose
            'bgee=s'            => \$bgee_connector,   # Bgee connector string
            'path_mas5=s'       => \$path_mas5,        # Affymetrix/processed_mas5/
            'path_schuster=s'   => \$path_schuster,    # Affymetrix/processed_schuster/
            'path_target=s'     => \$path_target,      # name files for replicates Affymetrix/bioconductor/targets/
            'path_processed=s'  => \$path_processed,   # final result dir Affymetrix/processed_differential/
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $path_mas5 eq '' || $path_schuster eq '' || $path_target eq '' || $path_processed eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -path_mas5=\$(MAS5PATH) -path_schuster=\$(SCHUSTERPATH) -path_target=\$(BIOCONDUCTORTARG) -path_processed=\$(DIFFEXPRPATH)  <expId>
\t-bgee             Bgee    connector string
\t-path_mas5        Affymetrix/processed_mas5/              directory path
\t-path_schuster    Affymetrix/processed_schuster/          directory path
\t-path_target      Affymetrix/bioconductor/targets/        directory path
\t-path_processed   Affymetrix/processed_differential/      directory path
\t-debug            More verbose
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


my %experiments;
# Monster query to get info on chips for diff. analyses.
# It makes the union of two queries: one to get info for diff.analyses comparing one organ at different stages,
# one for analyses comparing different organs at a same broad stage.
# TODO DRY, large parts of the query repeated, also in script launch_diff_analysis-RNASeq.pl
# XXX: actually, the root of the ontology "life stage" (UBERON:0000104) is filtered afterwards,
# it could be filtered right away in the query
my $selExpr = $dbh->prepare(
                            # First query, get info for comparing a same organ at different stages
                            # ("development" comparison factor)
                            '(SELECT DISTINCT affymetrixChipId, microarrayExperimentId, chipTypeId, normalizationType,
                                anatEntityId, bgeeAffymetrixChipId, '
                                # Specify the comparison factor
                                .'"development" AS comparisonFactor, '

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

                                # Get the species related to this chip. We use a GROUP_CONCAT,
                                # it allows to check the related species over all genes of the chip,
                                # not by checking only one gene.
                                # TODO we might also check that the species of the chip
                                # corresponds to the species of the annotated stage.
                                .'(SELECT GROUP_CONCAT(DISTINCT t4.speciesId SEPARATOR ", ") FROM gene AS t4 INNER
                                   JOIN affymetrixProbeset AS t5 ON t4.geneId = t5.geneId WHERE
                                   t5.bgeeAffymetrixChipId = t1.bgeeAffymetrixChipId) AS speciesIds

                                FROM affymetrixChip AS t1 '
                                # Join to stage and to taxon constraint tables for the stage sub-query
                                .'INNER JOIN stage AS t2 ON t1.stageId = t2.stageId
                                INNER JOIN stageTaxonConstraint AS t2bis on t2.stageId = t2bis.stageId
                             )

                             UNION ALL '

                            # Second query, get info for comparing different organs at a same broad stage
                            # ("anatomy" comparison factor)
                            .'(SELECT DISTINCT affymetrixChipId, microarrayExperimentId, chipTypeId, normalizationType,
                                anatEntityId, bgeeAffymetrixChipId, '
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

                                # Get the species related to this chip. We use a GROUP_CONCAT,
                                # it allows to check the related species over all genes of the chip,
                                # not by checking only one gene.
                                # TODO we might also check that the species of the chip
                                # corresponds to the species of the annotated stage.
                                .'(SELECT GROUP_CONCAT(DISTINCT t4.speciesId SEPARATOR ", ") FROM gene AS t4 INNER
                                   JOIN affymetrixProbeset AS t5 ON t4.geneId = t5.geneId WHERE
                                   t5.bgeeAffymetrixChipId = t1.bgeeAffymetrixChipId) AS speciesIds

                                FROM affymetrixChip AS t1 '
                                # Join to stage and to taxon constraint tables for the stage sub-query
                                .'INNER JOIN stage AS t2 ON t1.stageId = t2.stageId
                                INNER JOIN stageTaxonConstraint AS t2bis on t2.stageId = t2bis.stageId)');

$selExpr->execute()  or die $selExpr->errstr;
while ( my @data = $selExpr->fetchrow_array ){
    # Multiple species per chip => Problem
    if ( $data[7] =~ /,/ ){
        warn "Problem with [$data[1]]: multiple species in it [$data[7]]\n";
        next;
    }

    if ( !defined $ARGV[0] || (defined $ARGV[0] && $data[1] eq $ARGV[0]) ){
        if ( $data[3] eq 'MAS5' ){    # normalization by MAS5
            if ( $data[8] eq 'anatomy' ){
                # Anatomy comparison => Single stage comparison
                $experiments{$data[1]}->{$data[2]}->{'s_'.$data[5]}->{$data[4]}->{$data[0]}->{'path'}                 = $path_mas5;
                $experiments{$data[1]}->{$data[2]}->{'s_'.$data[5]}->{$data[4]}->{$data[0]}->{'bgeeAffymetrixChipId'} = $data[6];
                $experiments{$data[1]}->{$data[2]}->{'s_'.$data[5]}->{$data[4]}->{$data[0]}->{'speciesId'}            = $data[7];
            }
            elsif ( $data[8] eq 'development' ){
                # Development comparison => Single organ comparison
                $experiments{$data[1]}->{$data[2]}->{'a_'.$data[4]}->{$data[5]}->{$data[0]}->{'path'}                 = $path_mas5;
                $experiments{$data[1]}->{$data[2]}->{'a_'.$data[4]}->{$data[5]}->{$data[0]}->{'bgeeAffymetrixChipId'} = $data[6];
                $experiments{$data[1]}->{$data[2]}->{'a_'.$data[4]}->{$data[5]}->{$data[0]}->{'speciesId'}            = $data[7];
            }
        }
        elsif( $data[3] eq 'gcRMA' ){  # normalization by gcRNA
            # There are different possibilites for the file name,
            # based on the affymetrixChipId: with an added suffix .cel, .CEL, or no added suffix
            # and then, always .out at the end
            my $fileName;
            if ( -e $path_schuster.'/'.$data[1].'/'.$data[0].'.cel.out' ){
                $fileName = $data[0].'.cel.out';
            }
            elsif ( -e $path_schuster.'/'.$data[1].'/'.$data[0].'.CEL.out' ){
                $fileName = $data[0].'.CEL.out';
            }
            elsif ( -e $path_schuster.'/'.$data[1].'/'.$data[0].'.out' ){
                $fileName = $data[0].'.out';
            }
            else {
                die "Error, no processed_schuster file found for expId: [$data[1]] - chipId: [$data[0]]\n";
            }

            if ( $data[8] eq 'anatomy' ){
                # Anatomy comparison => Single stage comparison
                $experiments{$data[1]}->{$data[2]}->{'s_'.$data[5]}->{$data[4]}->{$fileName}->{'path'}                 = $path_schuster;
                $experiments{$data[1]}->{$data[2]}->{'s_'.$data[5]}->{$data[4]}->{$fileName}->{'bgeeAffymetrixChipId'} = $data[6];
                $experiments{$data[1]}->{$data[2]}->{'s_'.$data[5]}->{$data[4]}->{$fileName}->{'speciesId'}            = $data[7];
            }
            elsif ( $data[8] eq 'development' ){
                # Development comparison => Single organ comparison
                $experiments{$data[1]}->{$data[2]}->{'a_'.$data[4]}->{$data[5]}->{$fileName}->{'path'}                 = $path_schuster;
                $experiments{$data[1]}->{$data[2]}->{'a_'.$data[4]}->{$data[5]}->{$fileName}->{'bgeeAffymetrixChipId'} = $data[6];
                $experiments{$data[1]}->{$data[2]}->{'a_'.$data[4]}->{$data[5]}->{$fileName}->{'speciesId'}            = $data[7];
            }
        }
    }
}
$selExpr->finish;
$dbh->disconnect;

# Remove "life stage" (UBERON:0000104) from analyses because development stage root!
my $to_exclude = 'UBERON:0000104';

my %exp_to_treat;
for my $exp ( keys %experiments ){
    for my $array ( keys %{$experiments{$exp}} ){
        for my $single ( keys %{$experiments{$exp}->{$array}} ){
            if ( $single eq "s_$to_exclude" ){ # Remove "life stage" (UBERON:0000104) from analyses because development stage root!
                delete $exp_to_treat{$exp}->{$array}->{$single};
                next;
            }
            my $count = 0;
            for my $condition ( keys %{$experiments{$exp}->{$array}->{$single}} ){
                # Remove "life stage" (UBERON:0000104) from analyses because development stage root!
                next  if ( $condition eq $to_exclude );
                # if there are some replicates for that condition (organ/stage)
                if ( scalar keys %{$experiments{$exp}->{$array}->{$single}->{$condition}} > 1 ){
                    $exp_to_treat{$exp}->{$array}->{$single}->{$condition} = $experiments{$exp}->{$array}->{$single}->{$condition};
                    $count++;
                }
            }

            # If at least 3 conditions had replicates
            if ( $count < 3 ){
                delete $exp_to_treat{$exp}->{$array}->{$single};
            }
        }
        # Remove array/chiptype that didn't match the criteria
        if ( scalar keys %{$exp_to_treat{$exp}->{$array}} == 0 ){
            delete $exp_to_treat{$exp}->{$array};
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

for my $exp ( sort keys %exp_to_treat ){
    for my $array ( keys %{$exp_to_treat{$exp}} ){
        if ( !-e $path_processed.'/'.$exp ){
            mkdir "$path_processed/$exp";
        }
        # If mix of mas5 & Schuster, simply remove MAS5 because less valuable
        # R will crash with both input_type and if we create 2 different target files
        # we may have issues with out files and how they will be found by next scripts
        my %mix;
        for my $single ( keys %{$exp_to_treat{$exp}->{$array}} ){
            for my $condition ( keys %{$exp_to_treat{$exp}->{$array}->{$single}} ){
                for my $file ( keys %{$exp_to_treat{$exp}->{$array}->{$single}->{$condition}} ){
                    $mix{ $exp_to_treat{$exp}->{$array}->{$single}->{$condition}->{$file}->{'path'} }++;
                }
            }
        }

        # Create file .targets
        SINGLE:
        for my $single ( keys %{$exp_to_treat{$exp}->{$array}} ){
            my $nbr_conditions = scalar keys %{$exp_to_treat{$exp}->{$array}->{$single}};
            # Skip case with several repeats but low number of conditions
            next SINGLE  if ( $nbr_conditions < 3 );

            my $comparisonFactor = $single =~ /^a_/ ? 'development' : 'anatomy';
            printf("\t%-15s %-15s %-15s %s...", $exp, $array, $comparisonFactor, $single);
            # Target file $exp__$array must contain 1 stage and several organs OR 1 organ for several stages
            # path/exp__array__singleorgan  OR  path/exp__array__singlestage
            my $target = "$path_target/${exp}___${array}__$single.target";
            my $path   = '';
            open(my $TARGET, '>', "$target")  or die 'Cannot open TARGET file';
            print {$TARGET} "#organ\tstage\tfile\tpath\tbgeeAffymetrixChipId\tspeciesId\tcomparisonFactor\n";
            for my $condition ( keys %{$exp_to_treat{$exp}->{$array}->{$single}} ){
                for my $file ( keys %{$exp_to_treat{$exp}->{$array}->{$single}->{$condition}} ){
                    # Can get $path from the last loop/loop/loop iteration because we don't mix MAS5 and Shuster in target files !
                    $path = $exp_to_treat{$exp}->{$array}->{$single}->{$condition}->{$file}->{'path'};
                    # Mix of MAS5 & Schuster, so skip MAS5
                    if ( scalar keys %mix > 1 ){
                        next  if ( $path =~ /processed_mas5/ );
                    }
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
                        die "Problem with organ/stage for [$exp][$array][$single][$comparisonFactor][$condition]";
                    }
                    print {$TARGET} "$organ\t$stage\t$file\t".
                        $exp_to_treat{$exp}->{$array}->{$single}->{$condition}->{$file}->{'path'}."\t".
                        $exp_to_treat{$exp}->{$array}->{$single}->{$condition}->{$file}->{'bgeeAffymetrixChipId'}."\t".
                        $exp_to_treat{$exp}->{$array}->{$single}->{$condition}->{$file}->{'speciesId'}."\t".
                        "$comparisonFactor\n";
                }
            }
            close $TARGET;
            if ( scalar keys %mix > 1 ){
                warn "\nMix of MAS & Schuster in [$exp][$array]\n";
                $path = $path_schuster;
            }

            # Launch R script
            my $log = "$path_processed/$exp/${exp}___${array}__$single.log";
            my $cmd = '';
            if ( $path eq $path_schuster ){
                $cmd = "R CMD BATCH --vanilla '--args target_file_path=\"$target\" output_folder_path=\"$path_processed/$exp\" input_folder_path=\"$path_schuster/$exp\" input_type=\"Schuster\" array_type=\"$array\"' $FindBin::Bin/diff_analysis_affy.R  $log";
            }
            elsif ( $path eq $path_mas5 ){
                $cmd = "R CMD BATCH --vanilla '--args target_file_path=\"$target\" output_folder_path=\"$path_processed/$exp\" input_folder_path=\"$path_mas5/$exp\"     input_type=\"MAS5\"     array_type=\"$array\"' $FindBin::Bin/diff_analysis_affy.R  $log";
            }
            system($cmd)==0  or do{warn "\tsystem [$cmd] failed: $?\n"; map { system("mv $_ $_.PROB") } glob("$path_processed/$exp/*.out");};
            print "\n";
        }
    }
}

exit 0;

