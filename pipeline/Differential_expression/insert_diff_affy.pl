#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Roux, created 03/07/09
# USAGE: perl insert_diff_affy.pl <experimentId chipTypeId>
############################################

use Getopt::Long;
use File::Basename;
use File::Slurp;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

$| = 1;

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug) = (0);
my ($dry_run) = (0);
my ($path_target, $path_processed) = ('', '');
my %opts = ('debug'             => \$debug,            # more verbose
            'dry_run'           => \$dry_run,          # run scripts with no actual insertions/updates performed
            'bgee=s'            => \$bgee_connector,   # Bgee connector string
            'path_target=s'     => \$path_target,      # name files for replicates Affymetrix/bioconductor/targets/
            'path_processed=s'  => \$path_processed,   # final result dir Affymetrix/processed_differential/
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || ($#ARGV ne -1 && $#ARGV ne 1) || $path_target eq '' || $path_processed eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -path_target=\$(BIOCONDUCTORTARG) -path_processed=\$(DIFFEXPRPATH)   <expId> <chipTypeId>
\t-bgee             Bgee    connector string
\t-path_target      Affymetrix/bioconductor/targets/        directory path
\t-path_processed   Affymetrix/processed_differential/      directory path
\t-debug            More verbose
\t-dry_run          run scripts with no actual insertions/updates performed
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);



###################################################################
# Read the Experiment, chiptype and conditions and fill dea tables
###################################################################
print 'Reading informations on differential expression analyses... '  if ( $debug );

my %experiments;
# all the experiment/array that have been analyzed
for my $file ( glob($path_target.'/*.target') ){
    $file = basename($file);
    # ${exp}___${array}__$single.target
    # GSE19667___A-AFFY-44__a_CL:0002328.target
    # $exp_to_treat{$exp}->{$array}->{$single}->{$condition}
    if ( $file =~ /^(.+?)___(.+?)__(.+?).target$/ ){
        if ( defined $ARGV[0] && defined $ARGV[1] ){
            if ( $ARGV[0] eq $1 && $ARGV[1] eq $2 ){
                $experiments{$1}->{$2}->{$3} = ();
            }
        }
        else {
            $experiments{$1}->{$2}->{$3} = ();
        }
    }
    else {
        die "File name [$file] does not correspond to the proper file name pattern\n";
    }
}

# remove experiments that are not anymore in the annotation
my $selMicExp = $bgee->prepare('SELECT microarrayExperimentId FROM microarrayExperiment WHERE microarrayExperimentId=?');
for my $exp ( keys %experiments ){
    $selMicExp->execute($exp)  or die $selMicExp->errstr;
    my $inserted = $selMicExp->fetchrow_array;
    if ( !defined $inserted ){
        delete $experiments{$exp};
        die "Error, microarrayExperimentId not found in the database: [$exp]\n";
    }
}
$selMicExp->finish;

for my $exp ( keys %experiments ){
    for my $chiptype ( keys %{$experiments{$exp}} ){
        for my $single ( keys %{$experiments{$exp}->{$chiptype}} ){
            # Get condition count to properly retrieve .out file
            my $field      = $single =~ /^a_/ ? 2 : 1;
            my @conditions = `grep -v '^#' $path_target/${exp}___${chiptype}__$single.target | cut -f$field | uniq`;
            my $type       = $single =~ /^a_/ ? '1_'.scalar(@conditions) : scalar(@conditions).'_1';

            # Retrieve the chips used in each analysis
            my $previousOrganId = undef;
            my $previousStageId = undef;
            for my $line ( read_file("$path_target/${exp}___${chiptype}__$single.target", chomp=>1) ){
                next  if ( $line =~ /^#/ );
                my @tmp = split(/\t/, $line);

                #organ  stage  file  path  bgeeAffymetrixChipId  speciesId  comparisonFactor
                my $organId          = $tmp[0];
                my $stageId          = $tmp[1];
                my $affymetrixChipId = $tmp[4];
                my $speciesId        = $tmp[5];
                my $condition = $single =~ /^a_/ ? $stageId : $organId;

                if ( !defined $previousOrganId || $previousOrganId ne $organId || $previousStageId ne $stageId ){
                    my $tempOrganId = $organId;
                    my $tempStageId = $stageId;
                    # to generate file names, in diff_analysis.R, : are replaced by _ in organIds and stageIds
                    $tempOrganId =~ s{:}{_}g;
                    $tempStageId =~ s{:}{_}g;
                    my $file = "$path_processed/$exp/${chiptype}_${tempOrganId}_${tempStageId}_$type.out";
                    if ( -e $file && -s $file ){
                        $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'file'} = $file;
                    }
                    else {
                        die "Error, processed_differential file not found: [$file]\n";
                    }
                }
                $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'speciesId'} = $speciesId;
                $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'chips'}->{$affymetrixChipId}++;

                $previousOrganId = $organId;
                $previousStageId = $stageId;
            }
        }
    }
}
print "Done\n"  if ( $debug );


##########################################################################
# /!\ should check that experiments and chips are still in affymetrixChip
# The annotation shoud not have changed
##########################################################################
print "Inserting informations on differential expression analyses...\n"  if ( $debug );

my $insDea             = $bgee->prepare('INSERT INTO differentialExpressionAnalysis (detectionType, microarrayExperimentId, comparisonFactor) VALUES (?, ?, ?)');
my $insDeaChipGrp      = $bgee->prepare('INSERT INTO deaSampleGroup (deaId, anatEntityId, stageId) VALUES (?, ?, ?)');
my $insDeaChipGrp2Affy = $bgee->prepare('INSERT INTO deaSampleGroupToAffymetrixChip (deaSampleGroupId, bgeeAffymetrixChipId) VALUES (?, ?)');

for my $exp ( keys %experiments ){
    for my $chiptype ( keys %{$experiments{$exp}} ){
        for my $single ( keys %{$experiments{$exp}->{$chiptype}} ){
            printf("Inserting analysis for \t%-15s %-15s %s\n", $exp, $chiptype, $single);

            my $comparisonFactor = $single =~ /^a_/ ? 'development' : 'anatomy';
            # limma - MCM
            # one id per exp/chip
            my $dea_id = 0;
            if (!$dry_run) {
                $insDea->execute('Limma - MCM', $exp, $comparisonFactor)  or warn "Parameters: [Limma - MCM] - [$exp] - [$comparisonFactor]\n";
                $dea_id = $bgee->{'mysql_insertid'};
            } else {
                printf("INSERT INTO differentialExpressionAnalysis: %s - %s - %s\n", 'Limma - MCM', $exp, $comparisonFactor);
            }

            for my $condition ( keys %{$experiments{$exp}->{$chiptype}->{$single}} ){
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
                    die "Problem with organ/stage for [$exp][$chiptype][$single][$comparisonFactor][$condition]";
                }

                $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'} = 0;
                if (!$dry_run) {
                    $insDeaChipGrp->execute($dea_id, $organ, $stage)  or warn "Parameters: [$dea_id] - [$organ] - [$stage]\n";
                    $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'} = $bgee->{'mysql_insertid'};
                } else {
                    printf("INSERT INTO deaSampleGroup: %s - %s - %s\n", $dea_id, $organ, $stage);
                }

                # fill table deaSampleGroupToAffymetrixChip
                for my $chip ( keys %{$experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'chips'}} ){
                    #print "$experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}\t$chip\n";
                    if (!$dry_run) {
                        $insDeaChipGrp2Affy->execute($experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}, $chip)
                        or warn "Parameters: [$experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}] - [$chip]\n";
                    } else {
                        printf("INSERT INTO deaSampleGroupToAffymetrixChip: %s - %s",
                            $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}, $chip)
                    }
                }
            }
        }
    }
}
$insDea->finish;
$insDeaChipGrp->finish;
$insDeaChipGrp2Affy->finish;
print "Done\n"  if ( $debug );


###########################################
# Retrieve chips mapping (pbsets to genes)
###########################################
my $retrieveProbeset = $bgee->prepare('SELECT STRAIGHT_JOIN DISTINCT t2.affymetrixProbesetId, t2.geneId FROM affymetrixChip AS t1
                                       INNER JOIN affymetrixProbeset AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                                       WHERE t1.chipTypeId = ? AND t2.reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"');

# retrieve annotation
my %annotation;
for my $exp ( keys %experiments ){
    for my $chiptype ( keys %{$experiments{$exp}} ){
        if ( !exists $annotation{$chiptype} ){
            print "Retrieving probeset-gene mappings for [$chiptype]... "  if ( $debug );
            my %annot_pbsets;

            # Retrieve the mapping from affymetrixProbeset: this way,
            # we have the guarantee to only use probesets not removed, filtered, etc,
            # and we already have the mapping
            $retrieveProbeset->execute($chiptype)  or die $retrieveProbeset->errstr;
            while ( my @data = $retrieveProbeset->fetchrow_array ){
                # If we already have a mapping, there is an ambiguity.
                # Should never happen...
                if (exists $annot_pbsets{$data[0]}) {
                    $annot_pbsets{$data[0]} = '';
                    warn "Ambiguous mapping: ".$data[0]." - ".$data[1]."\n";
                } else {
                    $annot_pbsets{$data[0]} = $data[1];
                }
            }

            # keep the mapping to use it later
            $annotation{$chiptype} = \%annot_pbsets;
            print "Done\n"  if ( $debug );
        }
    }
}
$retrieveProbeset->finish;


##########################################
# fill deaAffymetrixProbesetSummary table
##########################################
print "Inserting on deaAffymetrixProbesetSummary table...\n"  if ( $debug );
my $insDeaSum = $bgee->prepare('INSERT INTO deaAffymetrixProbesetSummary (deaAffymetrixProbesetSummaryId, deaSampleGroupId, geneId, foldChange, differentialExpressionAffymetrixData, deaRawPValue) VALUES (?, ?, ?, ?, ?, ?)');

for my $exp ( keys %experiments ){
    for my $chiptype ( keys %{$experiments{$exp}} ){
        for my $single ( keys %{$experiments{$exp}->{$chiptype}} ){
            printf("Inserting deaAffymetrixProbesetSummary for \t%-15s %-15s %s\n", $exp, $chiptype, $single);

            for my $condition ( keys %{$experiments{$exp}->{$chiptype}->{$single}} ){
                printf("Condition: %s\n", $condition) if ($debug);
                # read the result file
                #TODO add a proper header to read this file with Utils::read_spreadsheet()
                printf("Reading file %s\n", $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'file'}) if ($debug);
                for my $line ( read_file($experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'file'}, chomp=>1) ){
                    next  if ( $line =~ /^#/ );
                    my @tmp = split(/\t/, $line);
                    #probe_set_ID   p_value   logFC   present_calls
                    printf("Read line: %s\n", $line) if ($debug);

                    # if the probeset is mapped to an ensembl gene
                    if ( exists $annotation{$chiptype}->{$tmp[0]} && $annotation{$chiptype}->{$tmp[0]} ne ''){
                        # 'not expressed'
                        if ( $tmp[3] eq 'FALSE' ){
                            if (!$dry_run) {
                                $insDeaSum->execute($tmp[0], $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}, $annotation{$chiptype}->{$tmp[0]}, $tmp[2], 'not expressed', $tmp[1])
                                or warn "Parameters: [$tmp[0]] - [$experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}] - [$annotation{$chiptype}->{$tmp[0]}] - [$tmp[2]] - [not expressed] - [$tmp[1]]\n";
                            } else {
                                printf("INSERT INTO deaAffymetrixProbesetSummary: %s - %s - %s - %s - %s - %s",
                                    $tmp[0], $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}, $annotation{$chiptype}->{$tmp[0]}, $tmp[2], 'not expressed', $tmp[1]);
                            }
                        }
                        # p_value threshold for 'high quality'
                        elsif ( $tmp[1] < 0.01 ){
                            if (!$dry_run) {
                                $insDeaSum->execute($tmp[0], $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}, $annotation{$chiptype}->{$tmp[0]}, $tmp[2], 'high quality', $tmp[1])
                                or warn "Parameters: [$tmp[0]] - [$experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}] - [$annotation{$chiptype}->{$tmp[0]}] - [$tmp[2]] - [high quality] - [$tmp[1]]\n";
                            } else {
                                printf("INSERT INTO deaAffymetrixProbesetSummary: %s - %s - %s - %s - %s - %s",
                                    $tmp[0], $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}, $annotation{$chiptype}->{$tmp[0]}, $tmp[2], 'high quality', $tmp[1]);
                            }
                        }
                        # p_value threshold for 'poor quality'
                        elsif ( $tmp[1] < 0.05 ){
                            if (!$dry_run) {
                                $insDeaSum->execute($tmp[0], $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}, $annotation{$chiptype}->{$tmp[0]}, $tmp[2], 'poor quality', $tmp[1])
                                or warn "Parameters: [$tmp[0]] - [$experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}] - [$annotation{$chiptype}->{$tmp[0]}] - [$tmp[2]] - [poor quality] - [$tmp[1]]\n";
                            } else {
                                printf("INSERT INTO deaAffymetrixProbesetSummary: %s - %s - %s - %s - %s - %s",
                                    $tmp[0], $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}, $annotation{$chiptype}->{$tmp[0]}, $tmp[2], 'poor quality', $tmp[1]);
                            }
                        }
                        # p_value threshold for 'no diff expression'
                        else {
                            if (!$dry_run) {
                                $insDeaSum->execute($tmp[0], $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}, $annotation{$chiptype}->{$tmp[0]}, $tmp[2], 'no diff expression', $tmp[1])
                                or warn "Parameters: [$tmp[0]] - [$experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}] - [$annotation{$chiptype}->{$tmp[0]}] - [$tmp[2]] - [no diff expression] - [$tmp[1]]\n";
                            } else {
                                printf("INSERT INTO deaAffymetrixProbesetSummary: %s - %s - %s - %s - %s - %s",
                                    $tmp[0], $experiments{$exp}->{$chiptype}->{$single}->{$condition}->{'group'}, $annotation{$chiptype}->{$tmp[0]}, $tmp[2], 'no diff expression', $tmp[1]);
                            }
                        }
                    } else {
                        if ( !exists $annotation{$chiptype}->{$tmp[0]} ) {
                            printf("No mapping for probeset ID %s in chip type ID %s\n", $tmp[0], $chiptype) if ($debug);
                        } else {
                            printf("Ambiguous mapping for probeset ID %s in chip type ID %s\n", $tmp[0], $chiptype) if ($debug);
                        }
                    }
                }
            }
        }
    }
}
$insDeaSum->finish;
$bgee->disconnect;
exit 0;

