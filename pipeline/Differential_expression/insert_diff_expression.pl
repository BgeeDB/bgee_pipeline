#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Roux, created 06/07/09
# USAGE: perl insert_diff_expression_affy.pl (organ stage)
##########################################################

use Getopt::Long;
use List::MoreUtils qw(uniq);

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

$| = 1;

# Define arguments & their default value
my ($bgee_connector, $datatype) = ('', '');
my ($debug) = (0);
my %opts = ('debug'     => \$debug,            # more verbose
            'bgee=s'    => \$bgee_connector,   # Bgee connector string
            'type=s'    => \$datatype,         # Affymetrix or RNASeq
           );
my %allowed_datatype = ('affy'   => 1,
                        'rnaseq' => 1,
                       );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || !exists $allowed_datatype{lc $datatype} || ($#ARGV ne -1 && $#ARGV ne 1) ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -type=affy  [<anatEntityId> <stageId>]
\t-bgee             Bgee connector string
\t-type             Datatype: affy or rnaseq
\t-debug            More verbose
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);



## Should we use the probes which are always absent or present/low quality (pre-filtering step performed in insert_affy.pl)?
## if the probesets do not work we should not use them
## but if they don't work they will not be detected as differentially expressed between conditions
## -> OK to keep them!


###########################
# retrieve all sampleGroups
###########################
my $selDeaId;
my $sql = '';
if ( lc $datatype eq 'affy' ){
    $sql = 'SELECT g.deaSampleGroupId, g.anatEntityId, g.stageId, a.comparisonFactor, s.geneId,        s.foldChange, s.deaRawPValue, s.reasonForExclusion, s.deaAffymetrixProbesetSummaryId
            FROM deaSampleGroup g, differentialExpressionAnalysis a, deaAffymetrixProbesetSummary s
            WHERE g.deaId=a.deaId AND s.deaSampleGroupId=g.deaSampleGroupId';
}
elsif ( lc $datatype eq 'rnaseq' ){
    $sql = 'SELECT g.deaSampleGroupId, g.anatEntityId, g.stageId, a.comparisonFactor, s.geneSummaryId, s.foldChange, s.deaRawPValue, s.reasonForExclusion
            FROM deaSampleGroup g, differentialExpressionAnalysis a, deaRNASeqSummary s
            WHERE g.deaId=a.deaId AND s.deaSampleGroupId=g.deaSampleGroupId';
}
else {
    die "Invalid datatype\n";
}

# All experiments
if ( $#ARGV eq -1 ){
    $selDeaId = $bgee->prepare($sql);
    $selDeaId->execute()  or die $selDeaId->errstr;
}
# or only one organ/stage
else {
    $sql .= ' AND g.anatEntityId=? AND g.stageId=?';
    $selDeaId = $bgee->prepare($sql);
    $selDeaId->execute($ARGV[0], $ARGV[1])  or die $selDeaId->errstr;
}


# Store only not excluded results
my %chipsgroups;
# store association between comparisonFactor -> stage -> organ -> gene and all sample groups,
# even those with excluded results.
my %allChipsGroups;
while ( my @data = $selDeaId->fetchrow_array ){
    # comparisonFactor -> stage -> organ -> gene

    # For a same condition, analyses can sometimes generate contradicting differential
    # expression calls ('over-expressed', 'under-expressed'). In that case, we take
    # the call with the best p-value. But p-values can sometimes be equal, or of level
    # of similar significance (e.g., p=0.001 and p=0.0011); in that case, it is not possible
    # to simply use the best p-value, and we need to use a voting system.
    #
    # The principle of the voting system is to cluster p-values by level of significance
    # (using the log10 of p-values rounded to no decimal), and, among the cluster of best p-values,
    # to choose the winning call based on the total number of unique conditions
    # studied, in the analyses supporting this call.
    #
    # So, we need to store p-values, rounded log10 p-values, and fold changes,
    # associated to their sampleGroupId, to be able to get back to the analysis that
    # the sample is part of, to retrieve the conditions that were compared in the analysis.

    # store for analysis ONLY not excluded results
    if ( $data[7] eq $Utils::CALL_NOT_EXCLUDED ){
        # default value in case a p-value = 0 (no log possible)
        my $roundedLog10 = -100000;
        if ( $data[6] != 0 ){
            $roundedLog10 = sprintf("%.0f", log($data[6])/log(10));
        }
        my $geneOrProbeId = $data[4];
        # With affymetrix data, there can be several probesets for a same gene in a same sample group,
        # and we want to store all of them
        if ( defined $data[8] ){
            $geneOrProbeId = $data[8];
        }
        $chipsgroups{$data[3]}->{$data[2]}->{$data[1]}->{$data[4]}->{'sampleGroupId'}->{$data[0]}->{$geneOrProbeId}->{'log 10 p-value'} = $roundedLog10;
        $chipsgroups{$data[3]}->{$data[2]}->{$data[1]}->{$data[4]}->{'sampleGroupId'}->{$data[0]}->{$geneOrProbeId}->{'p-value'}        = $data[6];
        $chipsgroups{$data[3]}->{$data[2]}->{$data[1]}->{$data[4]}->{'sampleGroupId'}->{$data[0]}->{$geneOrProbeId}->{'fold change'}    = $data[5];
    }

    # store all sample IDs for update queries
    push @{ $allChipsGroups{$data[3]}->{$data[2]}->{$data[1]}->{$data[4]}->{'sampleGroupId'} }, $data[0];
}
$selDeaId->finish;

# Store the relation between sampleGroupId and analysisId,
# and the relation between analysisId and conditions compared, for the voting system.
my %sampleToAnalysis;
my %analysisToCondition;
my $sampleAnalysisStmt = $bgee->prepare('SELECT deaSampleGroupId, deaId, anatEntityId, stageId FROM deaSampleGroup');
$sampleAnalysisStmt->execute()  or die $sampleAnalysisStmt->errstr;
while ( my @data = $sampleAnalysisStmt->fetchrow_array ){
    # sampleToAnalysis: sampleGroupeId = analysisId
    $sampleToAnalysis{$data[0]} = $data[1];
    # analysisToCondition: analysisId -> anatEntityId -> stageId
    $analysisToCondition{$data[1]}->{$data[2]}->{$data[3]} = 1;
}
$sampleAnalysisStmt->finish;


########################
# Fill expression table
########################
print "Filling differentialExpression...\n"  if ( $debug );

# Insert differential expression
my $count   = 0;
my $sqlI = '';
my $upDeaSum = '';
if ( lc $datatype eq 'affy' ){
    $sqlI     = 'INSERT INTO differentialExpression (geneId, anatEntityId, stageId, diffExprCallAffymetrix, diffExprAffymetrixData, comparisonFactor, bestPValueAffymetrix) VALUES (?, ?, ?, ?, ?, ?, ?)
                 ON DUPLICATE KEY UPDATE diffExprCallAffymetrix = GREATEST(diffExprCallAffymetrix + 0, VALUES(diffExprCallAffymetrix)),
                                         diffExprAffymetrixData = GREATEST(diffExprAffymetrixData + 0, VALUES(diffExprAffymetrixData)),
                                         bestPValueAffymetrix   = LEAST(bestPValueAffymetrix,          VALUES(bestPValueAffymetrix)),
                                         differentialExpressionId=LAST_INSERT_ID(differentialExpressionId)'; # MySQL VALUES can deal with non-indice
    $upDeaSum = 'UPDATE deaAffymetrixProbesetSummary SET differentialExpressionId="__deaid__" WHERE geneId="__geneid__"        AND deaSampleGroupId IN (__grpids__)';
}
elsif ( lc $datatype eq 'rnaseq' ){
    $sqlI     = 'INSERT INTO differentialExpression (geneId, anatEntityId, stageId, diffExprCallRNASeq,     diffExprRNASeqData,     comparisonFactor, bestPValueRNASeq) VALUES (?, ?, ?, ?, ?, ?, ?)
                 ON DUPLICATE KEY UPDATE diffExprCallRNASeq = GREATEST(diffExprCallRNASeq + 0, VALUES(diffExprCallRNASeq)),
                                         diffExprRNASeqData = GREATEST(diffExprRNASeqData + 0, VALUES(diffExprRNASeqData)),
                                         bestPValueRNASeq   = LEAST(bestPValueRNASeq,          VALUES(bestPValueRNASeq)),
                                         differentialExpressionId=LAST_INSERT_ID(differentialExpressionId)'; # MySQL VALUES can deal with non-indice
    $upDeaSum = 'UPDATE deaRNASeqSummary             SET differentialExpressionId="__deaid__" WHERE geneSummaryId="__geneid__" AND deaSampleGroupId IN (__grpids__)';
}
else {
    die "Invalid datatype\n";
}

my $insDea = $bgee->prepare($sqlI);
for my $comparisonFactor ( sort keys %chipsgroups ){
    for my $stage ( sort keys %{$chipsgroups{$comparisonFactor}} ){
        for my $organ ( sort keys %{$chipsgroups{$comparisonFactor}->{$stage}} ){
            printf("\t%-15s %-15s %s\n", $comparisonFactor, $organ, $stage);
            for my $gene ( sort keys %{$chipsgroups{$comparisonFactor}->{$stage}->{$organ}} ){
                my ($quality, $direction, $lowestPVal) = get_diff_expression_status($gene, $organ, $stage, $comparisonFactor);
                # Insert diff expression
                $insDea->execute($gene, $organ, $stage, $direction, $quality, $comparisonFactor, $lowestPVal)  or die $insDea->errstr;
                my $diff_expression_id = $bgee->{'mysql_insertid'};

                # Update affymetrixProbeset
                my $sampleGroupId_list = join(', ', sort { $a <=> $b } uniq @{ $allChipsGroups{$comparisonFactor}->{$stage}->{$organ}->{$gene}->{'sampleGroupId'} });
                my $sql  = $upDeaSum; # To not substitute in $upDeaSum
                $sql =~ s{__deaid__}{$diff_expression_id};
                $sql =~ s{__geneid__}{$gene};
                $sql =~ s{__grpids__}{$sampleGroupId_list};
                my $upId = $bgee->prepare($sql);
                $upId->execute()  or die $upId->errstr;
                $upId->finish;
                $count++;
            }
        }
    }
}
$insDea->finish;

print "Insertion $datatype done: $count\n";
$bgee->disconnect;
exit 0;


sub get_diff_expression_status {
    my ($gene, $organ, $stage, $comparisonFactor) = @_;
    my $info = $chipsgroups{$comparisonFactor}->{$stage}->{$organ}->{$gene};
    # $info->{'sampleGroupId'}->{sample group id}->{gene or probeset ID}->{'log 10 p-value'} = rounded log10 p-value
    # $info->{'sampleGroupId'}->{sample group id}->{gene or probeset ID}->{'p-value'} = real p-value
    # $info->{'sampleGroupId'}->{sample group id}->{gene or probeset ID}->{'fold change'} = fold change

    # Default values
    my $low_threshold    = 0.05;
    my $high_threshold   = 0.01;
    my $quality          = '';
    my $direction        = '';
    # As we use a voting system, the best winning p-value might not be the overall lowest p-value,
    # so we store the best winning p-value separately.
    my $winningLowestPVal = 1;

    # First, we identify the lowest p-value and rounded log10 p-value for this call
    my $lowestPVal      = 1;
    my $lowestLog10PVal = 1;
    for my $sampleGroupId ( sort keys %{$info->{'sampleGroupId'}} ){
        for my $geneOrProbesetId ( sort keys %{$info->{'sampleGroupId'}->{$sampleGroupId}} ){
            my $log10PVal = $info->{'sampleGroupId'}->{$sampleGroupId}->{$geneOrProbesetId}->{'log 10 p-value'};
            if ( $log10PVal < $lowestLog10PVal ){
                $lowestLog10PVal = $log10PVal;
            }
            my $pVal = $info->{'sampleGroupId'}->{$sampleGroupId}->{$geneOrProbesetId}->{'p-value'};
            if ( $pVal < $lowestPVal ){
                $lowestPVal = $pVal;
            }
        }
    }

    # Then, we check whether there are conflicts of any kind (e.g., 'over-expression'
    # vs. 'no diff expression'), and, more specifically, whether there are over-expressed
    # and under-expressed calls that could be in conflict (with both the best rounded log10 p-value).
    my $conflict = 0;
    my $hasBestOverExpression = 0;
    my $hasBestUnderExpression = 0;
    my $lastCall = '';
    for my $sampleGroupId ( sort keys %{$info->{'sampleGroupId'}} ){
        for my $geneOrProbesetId ( sort keys %{$info->{'sampleGroupId'}->{$sampleGroupId}} ){
            my $pVal       = $info->{'sampleGroupId'}->{$sampleGroupId}->{$geneOrProbesetId}->{'p-value'};
            my $log10PVal  = $info->{'sampleGroupId'}->{$sampleGroupId}->{$geneOrProbesetId}->{'log 10 p-value'};
            my $foldChange = $info->{'sampleGroupId'}->{$sampleGroupId}->{$geneOrProbesetId}->{'fold change'};

            my $iterateCall = '';
            if ( $pVal >= $low_threshold ){
                $iterateCall = 'no diff expression';
            }
            else {
                if ( $foldChange >= 0 ){
                    $iterateCall = 'over-expression';
                    if ( $log10PVal == $lowestLog10PVal ){
                        $hasBestOverExpression = 1;
                    }
                }
                else {
                    $iterateCall = 'under-expression';
                    if ( $log10PVal == $lowestLog10PVal ){
                        $hasBestUnderExpression = 1;
                    }
                }
            }
            if ( $lastCall ne '' && $lastCall ne $iterateCall ){
                $conflict = 1;
            }
            $lastCall = $iterateCall;
        }
    }

    # Now, we try to identify the most likely differential expression call.
    # If the lowest p-value is not significant, it's easy...
    if ( $lowestPVal >= $low_threshold ){
        $direction = 'no diff expression';
        $quality   = 'high quality';
        $winningLowestPVal = $lowestPVal;
    }
    else {
        # It happens sometimes that p-values of similar level of significance lead to
        # conflicting diff expression calls; for such cases, we use a voting system.
        # It compares the total number of unique conditions studied in the analyses
        # that produce the calls with the best p-values.

        # How to assert in perl?
        # assert $hasBestOverExpression || $hasBestUnderExpression

        if ( $hasBestOverExpression && !$hasBestUnderExpression ){
            # This one is easy, over-expression, best p-value wins
            $direction         = 'over-expression';
            $winningLowestPVal = $lowestPVal;
        }
        elsif ( !$hasBestOverExpression && $hasBestUnderExpression ){
            # This one is also easy, under-expression, best p-value wins
            $direction         = 'under-expression';
            $winningLowestPVal = $lowestPVal;
        }
        else {
            # Conflict over/under expressiobn found. We use the voting system.

            # store the unique conditions studied in the analyses that produced
            # 'over-expression' calls, among the calls in the group of lowest p-values.
            # anatEntityId -> stageId
            my %overConditions;
            # lowest real p-value supporting over-expression calls.
            my $lowestPValOver  = 1;
            # store the unique conditions studied in the analyses that produced
            # 'under-expression' calls, among the calls in the group of lowest p-values.
            # anatEntityId -> stageId
            my %underConditions;
            # lowest real p-value supporting under-expression calls.
            my $lowestPValUnder = 1;

            for my $sampleGroupId ( sort keys %{$info->{'sampleGroupId'}} ){
                for my $geneOrProbesetId ( sort keys %{$info->{'sampleGroupId'}->{$sampleGroupId}} ){
                    #use only results having a call with a p-value among the best ones.
                    if ( $info->{'sampleGroupId'}->{$sampleGroupId}->{$geneOrProbesetId}->{'log 10 p-value'} != $lowestLog10PVal ){
                        next;
                    }
                    # check that p-value is really significant (maybe the cluster of best p-values
                    # include both significant and non-significant p-values).
                    my $pVal = $info->{'sampleGroupId'}->{$sampleGroupId}->{$geneOrProbesetId}->{'p-value'};
                    if ( $pVal >= $low_threshold ){
                        next;
                    }

                    my $foldChange = $info->{'sampleGroupId'}->{$sampleGroupId}->{$geneOrProbesetId}->{'fold change'};
                    if ( $foldChange >= 0 && $pVal < $lowestPValOver ){
                        $lowestPValOver = $pVal;
                    }
                    elsif ($foldChange < 0 && $pVal < $lowestPValUnder ){
                        $lowestPValUnder = $pVal;
                    }

                    # find and store the conditions studied in the analysis this sample group is part of.
                    my $analysisId = $sampleToAnalysis{$sampleGroupId};
                    for my $anatEntityId ( sort keys %{$analysisToCondition{$analysisId}} ){
                        for my $stageId ( sort keys %{$analysisToCondition{$analysisId}->{$anatEntityId}} ){
                            if ( $foldChange >= 0 ){
                                $overConditions{$anatEntityId}->{$stageId} = 1;
                            }
                            else {
                                $underConditions{$anatEntityId}->{$stageId} = 1;
                            }
                        }
                    }
                }
            }

            # I don't know how to assert in perl:
            # assert keys %overConditions || keys %underConditions
            # assert $lowestPValOver < 1 || $lowestPValUnder < 1

            # find out who won
            my $overConditionCount  = 0;
            for my $anatEntityId ( keys %overConditions ){
                $overConditionCount += scalar keys %{$overConditions{$anatEntityId}};
            }
            my $underConditionCount = 0;
            for my $anatEntityId ( keys %underConditions ){
                $underConditionCount += scalar keys %{$underConditions{$anatEntityId}};
            }

            if ( $overConditionCount > $underConditionCount ||
                 # Ultimately, if we can't still decide, we simply use the best p-value...
                 ($overConditionCount == $underConditionCount && $lowestPValOver <= $lowestPValUnder) ){
                $direction         = 'over-expression';
                $winningLowestPVal = $lowestPValOver;
            }
            else {
                $direction         = 'under-expression';
                $winningLowestPVal = $lowestPValUnder;
            }
        }

        # get the quality
        if ( $winningLowestPVal >= $high_threshold ){ # && p-val < $low_threshold
            $quality   = 'poor quality';
        }
        else{
            $quality   = 'high quality';
        }
    }

    # Lower the quality in case of conflict
    if ( $conflict ){
        $quality = 'poor quality';
    }

    return ($quality, $direction, $winningLowestPVal);
}

