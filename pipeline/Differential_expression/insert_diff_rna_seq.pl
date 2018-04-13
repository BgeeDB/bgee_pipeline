#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Roux, created 03/07/09, adapted for RNA-seq by smoretti (2014-11-04)
# USAGE: perl insert_diff_rna_seq.pl <experimentId speciesId>
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
my ($path_target, $path_processed) = ('', '');
my %opts = ('debug'             => \$debug,            # more verbose
            'bgee=s'            => \$bgee_connector,   # Bgee connector string
            'path_target=s'     => \$path_target,      # name files for replicates rna_seq/bioconductor/targets/
            'path_processed=s'  => \$path_processed,   # final result dir rna_seq/processed_differential/
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || ($#ARGV ne -1 && $#ARGV ne 1) || $path_target eq '' || $path_processed eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -path_target=\$(RNASEQBIOCONDUCTORTARG) -path_processed=\$(RNASEQDIFFEXPRPATH)  <expId> <speciesId>
\t-bgee             Bgee    connector string
\t-path_target      rna_seq/bioconductor/targets/        directory path
\t-path_processed   rna_seq/processed_differential/      directory path
\t-debug            More verbose
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);



###################################################################
# Read the Experiment, speciesId and conditions and fill dea tables
###################################################################
print 'Reading informations on differential expression analyses... '  if ( $debug );

my %experiments;
# all the experiment/array that have been analyzed
for my $file ( glob($path_target.'/*.target') ){
    $file = basename($file);
    # ${exp}___${speciesId}__$single.target
    # GSE30352___10090__s_UBERON:0000113.target
    # $exp_to_treat{$exp}->{$speciesId}->{$single}->{$condition}
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
my $selMicExp = $bgee->prepare('SELECT rnaSeqExperimentId FROM rnaSeqExperiment WHERE rnaSeqExperimentId=?');
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
    for my $speciesId ( keys %{$experiments{$exp}} ){
        for my $single ( keys %{$experiments{$exp}->{$speciesId}} ){
            # Get condition count to properly retrieve .out file
            my $field      = $single =~ /^a_/ ? 2 : 1;
            my @conditions = `grep -v '^#' $path_target/${exp}___${speciesId}__$single.target | cut -f$field | uniq`;
            my $type       = $single =~ /^a_/ ? '1_'.scalar(@conditions) : scalar(@conditions).'_1';

            # Retrieve the chips used in each analysis
            my $previousOrganId = undef;
            my $previousStageId = undef;
            for my $line ( read_file("$path_target/${exp}___${speciesId}__$single.target", chomp=>1) ){
                next  if ( $line =~ /^#/ );
                my @tmp = split(/\t/, $line);

                #organ  stage  file  path  bgeeRNASeqLibId  speciesId  comparisonFactor
                my $organId          = $tmp[0];
                my $stageId          = $tmp[1];
                my $rnaSeqLibraryId  = $tmp[4];
                my $condition = $single =~ /^a_/ ? $stageId : $organId;

                if ( !defined $previousOrganId || $previousOrganId ne $organId || $previousStageId ne $stageId ){
                    my $tempOrganId = $organId;
                    my $tempStageId = $stageId;
                    # to generate file names, in diff_analysis.R, : are replaced by _ in organIds and stageIds
                    $tempOrganId =~ s{:}{_}g;
                    $tempStageId =~ s{:}{_}g;
                    my $file = "$path_processed/$exp/${speciesId}_${tempOrganId}_${tempStageId}_$type.out";
                    if ( -e $file && -s $file ){
                        $experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'file'} = $file;
                    }
                    else {
                        die "Error, processed_differential file not found: [$file]\n";
                    }
                }
                $experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'chips'}->{$rnaSeqLibraryId}++;

                $previousOrganId = $organId;
                $previousStageId = $stageId;
            }
        }
    }
}
print "Done\n"  if ( $debug );


##########################################################################
# /!\ should check that experiments and lib are still in rnaSeqLibrary
# The annotation shoud not have changed
##########################################################################
print "Inserting informations on differential expression analyses...\n"  if ( $debug );

my $insDea             = $bgee->prepare('INSERT INTO differentialExpressionAnalysis (detectionType, rnaSeqExperimentId, comparisonFactor) VALUES (?, ?, ?)');
my $insDeaChipGrp      = $bgee->prepare('INSERT INTO deaSampleGroup (deaId, anatEntityId, stageId) VALUES (?, ?, ?)');
my $insDeaChipGrp2Affy = $bgee->prepare('INSERT INTO deaSampleGroupToRnaSeqLibrary (deaSampleGroupId, rnaSeqLibraryId) VALUES (?, ?)');

for my $exp ( keys %experiments ){
    for my $speciesId ( keys %{$experiments{$exp}} ){
        for my $single ( keys %{$experiments{$exp}->{$speciesId}} ){
            printf("\t%-15s %-15s %s\n", $exp, $speciesId, $single);

            my $comparisonFactor = $single =~ /^a_/ ? 'development' : 'anatomy';
            # limma - MCM
            # one id per exp/chip
            $insDea->execute('Limma - MCM', $exp, $comparisonFactor)  or warn "Parameters: [Limma - MCM] - [$exp] - [$comparisonFactor]\n";
            my $dea_id = $bgee->{'mysql_insertid'};

            for my $condition ( keys %{$experiments{$exp}->{$speciesId}->{$single}} ){
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
                    die "Problem with organ/stage for [$exp][$speciesId][$single][$comparisonFactor][$condition]";
                }

                $insDeaChipGrp->execute($dea_id, $organ, $stage)  or warn "Parameters: [$dea_id] - [$organ] - [$stage]\n";
                $experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'} = $bgee->{'mysql_insertid'};

                # fill table deaSampleGroupToRnaSeqLibrary
                for my $chip ( keys %{$experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'chips'}} ){
                    #print "$experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'}\t$chip\n";
                    $insDeaChipGrp2Affy->execute($experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'}, $chip)
                        or warn "Parameters: [$experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'}] - [$chip]\n";
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
# Retrieve chips mapping (libs to genes)
###########################################
my $retrieveRNASeqLib = $bgee->prepare('SELECT DISTINCT t1.geneId FROM gene AS t1 WHERE t1.speciesId = ?
                                        AND EXISTS (SELECT 1 FROM rnaSeqResult AS t2
                                        WHERE t2.geneId = t1.geneId AND t2.reasonForExclusion = "'
                                            .$Utils::CALL_NOT_EXCLUDED.'" limit 1)');

# retrieve annotation
my %annotation;
for my $exp ( keys %experiments ){
    for my $speciesId ( keys %{$experiments{$exp}} ){
        if ( !exists $annotation{$speciesId} ){
            print "Retrieving mapping of [$speciesId]... "  if ( $debug );

            my %annot_libs;

            # Retrieve all genes not excluded in at least on RNA-seq library for the species considered
            $retrieveRNASeqLib->execute($speciesId)  or die $retrieveRNASeqLib->errstr;
            while ( my @data = $retrieveRNASeqLib->fetchrow_array ){
                $annot_libs{$data[0]} = 1;
            }

            # keep the mapping to use it later
            $annotation{$speciesId} = \%annot_libs;
            print "Done\n"  if ( $debug );
        }
    }
}
$retrieveRNASeqLib->finish;


##########################################
# fill deaRNASeqSummary table
##########################################
print "Inserting on deaRNASeqSummary table...\n"  if ( $debug );
my $insDeaSum = $bgee->prepare('INSERT INTO deaRNASeqSummary (geneSummaryId, deaSampleGroupId, foldChange, differentialExpressionRNASeqData, deaRawPValue) VALUES (?, ?, ?, ?, ?)');

for my $exp ( keys %experiments ){
    for my $speciesId ( keys %{$experiments{$exp}} ){
        for my $single ( keys %{$experiments{$exp}->{$speciesId}} ){
            printf("\t%-15s %-15s %s\n", $exp, $speciesId, $single);

            for my $condition ( keys %{$experiments{$exp}->{$speciesId}->{$single}} ){
                # read the result file
                #TODO add a proper header to read this file with Utils::read_spreadsheet()
                for my $line ( read_file($experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'file'}, chomp=>1) ){
                    next  if ( $line =~ /^#/ );
                    my @tmp = split(/\t/, $line);
                    #gene_names   p_value   logFC   present_calls

                    # if the gene exists and was seen as not excluded at least once
                    if ( exists $annotation{$speciesId}->{$tmp[0]} ){
                        # 'not expressed'
                        if ( $tmp[3] eq 'FALSE' ){
                             $insDeaSum->execute($tmp[0], $experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'}, $tmp[2], 'not expressed', $tmp[1])
                                 or warn "Parameters: [$tmp[0]] - [$experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'}] - [$tmp[2]] - [not expressed] - [$tmp[1]]\n";
                        }
                        # p_value threshold for 'high quality'
                        elsif ( $tmp[1] < 0.01 ){
                            $insDeaSum->execute($tmp[0], $experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'}, $tmp[2], 'high quality', $tmp[1])
                                or warn "Parameters: [$tmp[0]] - [$experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'}] - [$tmp[2]] - [high quality] - [$tmp[1]]\n";
                        }
                        # p_value threshold for 'poor quality'
                        elsif ( $tmp[1] < 0.05 ){
                            $insDeaSum->execute($tmp[0], $experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'}, $tmp[2], 'poor quality', $tmp[1])
                                or warn "Parameters: [$tmp[0]] - [$experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'}] - [$tmp[2]] - [poor quality] - [$tmp[1]]\n";
                        }
                        # p_value threshold for 'no diff expression'
                        else {
                            $insDeaSum->execute($tmp[0], $experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'}, $tmp[2], 'no diff expression', $tmp[1])
                                or warn "Parameters: [$tmp[0]] - [$experiments{$exp}->{$speciesId}->{$single}->{$condition}->{'group'}] - [$tmp[2]] - [no diff expression] - [$tmp[1]]\n";
                        }
                    } else {
                        printf("Gene %s does not exist or was never seen as not excluded in RNA-seq libraries\n", $tmp[0]) if ($debug);
                    }
                }
            }
        }
    }
}
$insDeaSum->finish;
$bgee->disconnect;
exit 0;

