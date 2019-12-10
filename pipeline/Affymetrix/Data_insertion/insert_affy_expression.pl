#!/usr/bin/env perl
# Julien Roux, created 24/07/09
# modified 15/08/16, JR: revamped whole script for bgee v14 with condition Ids
# modified Nov 2016, Frederic Bastian: split former create_files_affy.pl script between insert_affy.pl and insert_affy_expression.pl.
# last modified Dec. 2016 for use of new microarrayExperimentExpression table

# Perl core modules
use strict;
use warnings;
use diagnostics;
use File::Slurp;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

$| = 1; # no buffering of output


# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug)          = (0);
my %opts = ('bgee=s'     => \$bgee_connector,   # Bgee connector string
            'debug'      => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\t-debug            printing the update/insert SQL queries, not executing them
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);


##########################################
# GET CONDITIONS STUDIED WITH AFFYMETRIX #
##########################################
print "Retrieving conditions...\n";

my $queryConditions = $bgee->prepare('SELECT DISTINCT t2.exprMappedConditionId FROM affymetrixChip AS t1
                                      INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId');
$queryConditions->execute()  or die $queryConditions->errstr;
my @exprMappedConditions = ();
while ( my @data = $queryConditions->fetchrow_array ){
    push(@exprMappedConditions, $data[0]);
}

print "Done, ", scalar(@exprMappedConditions), " conditions retrieved.\n";

#####################################################
# EXAMINE EXPRESSION CALLS ACROSS ALL SAMPLES       #
# TO KNOW WHICH PROBESETS ARE NEVER SEEN AS PRESENT #
#####################################################
print "Examining all results to update probesets never seen as 'present'...\n";
my $findPresentGenes = $bgee->prepare('CREATE TEMPORARY TABLE tempPresentAffy (PRIMARY KEY (affymetrixProbesetId, chipTypeId)) ENGINE=InnoDB
                                        AS (
                                            SELECT DISTINCT t1.affymetrixProbesetId, t2.chipTypeId
                                            FROM affymetrixProbeset AS t1
                                            INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                                            WHERE t1.detectionFlag = "'.$Utils::PRESENT_CALL.'"
                                        )');
$findPresentGenes->execute()  or die $findPresentGenes->errstr;

# The "not excluded" status has to be set first, for the next query to properly take into account "undefined" status
my $presentUp = $bgee->prepare('UPDATE affymetrixProbeset AS t1
                                INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                                INNER JOIN tempPresentAffy AS t3
                                ON t1.affymetrixProbesetId = t3.affymetrixProbesetId AND t2.chipTypeId = t3.chipTypeId
                                SET reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"');
$presentUp->execute()  or die $presentUp->errstr;
my $preFilteringUp = $bgee->prepare('UPDATE affymetrixProbeset AS t1
                                     INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                                     LEFT OUTER JOIN tempPresentAffy AS t3
                                     ON t1.affymetrixProbesetId = t3.affymetrixProbesetId AND t2.chipTypeId = t3.chipTypeId
                                     SET reasonForExclusion = "'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"
                                     WHERE t3.affymetrixProbesetId IS NULL AND t1.reasonForExclusion != "'.$Utils::EXCLUDED_FOR_UNDEFINED.'"');
$preFilteringUp->execute()  or die $preFilteringUp->errstr;

my $dropTempPresent = $bgee->prepare('DROP TABLE tempPresentAffy');
$dropTempPresent->execute()  or die $dropTempPresent->errstr;

print "Done\n";

##########################################
# PREPARE QUERIES                        #
##########################################

# insert/update expression
my $insUpExpr   = $bgee->prepare('INSERT INTO expression (bgeeGeneId, conditionId) VALUES (?, ?)
                                  ON DUPLICATE KEY UPDATE expressionId=LAST_INSERT_ID(expressionId)');

# Query to update affymetrixProbeset with the expressionId and reasonForExclusion
my $updResult = $bgee->prepare('UPDATE affymetrixProbeset AS t1
                                INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                                INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                SET expressionId = ?, reasonForExclusion = ?
                                WHERE t1.bgeeGeneId = ? AND t3.exprMappedConditionId = ?
                                    AND t1.reasonForExclusion != "'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"');

# Insertion into microarrayExperimentExpression
my $insExpSummary = $bgee->prepare('INSERT INTO microarrayExperimentExpression
                                    (expressionId, microarrayExperimentId,
                                        presentHighMicroarrayChipCount, presentLowMicroarrayChipCount,
                                        absentHighMicroarrayChipCount, absentLowMicroarrayChipCount,
                                        microarrayExperimentCallDirection, microarrayExperimentCallQuality)
                                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)');

# query to get all the Affymetrix probesets for a condition
my $queryResults = $bgee->prepare("SELECT t1.bgeeGeneId, t2.microarrayExperimentId, t1.bgeeAffymetrixChipId,
                                          t1.affymetrixProbesetId, t1.detectionFlag, t1.affymetrixData
                                   FROM affymetrixProbeset AS t1
                                   INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                                   INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                   WHERE t3.exprMappedConditionId = ? AND t1.reasonForExclusion = '".$Utils::CALL_NOT_EXCLUDED."'");

##########################################
# SUBROUTINES                            #
##########################################
# subroutine to insert/update in expression table
# return the relevant expressionId
sub add_expression {
    my ($bgeeGeneId, $exprMappedConditionId) = @_;

    $insUpExpr->execute($bgeeGeneId, $exprMappedConditionId)  or die $insUpExpr->errstr;
    return $bgee->{'mysql_insertid'}; # expr_id
}

##########################################
# ITERATING CONDITIONS TO INSERT DATA    #
##########################################

print "Processing conditions...\n";
for my $exprMappedConditionId ( @exprMappedConditions ){
    print "\tconditionId: $exprMappedConditionId\n";

    # retrieve genes results for this condition
    print "\t\tRetrieving genes results...\n";
    $queryResults->execute($exprMappedConditionId)  or die $queryResults->errstr;
    my %results = ();
    while ( my @data = $queryResults->fetchrow_array ){
        #data[0] = bgeeGeneId, data[1] = experimentId, data[2] = chipId, data[3] = probesetId, data[4] = detectionFlag, data[5] = affymetrixData
        #important to retrieve the probeset IDs for first per-chip reconciliation
        $results{$data[0]}->{$data[1]}->{$data[2]}->{$data[3]}->{'call'}    = $data[4];
        $results{$data[0]}->{$data[1]}->{$data[2]}->{$data[3]}->{'quality'} = $data[5];
    }
    print "\t\tDone, ", scalar(keys %results), ' genes retrieved. ',
          "Generating expression summary...\n";

    # now iterating the genes to insert expression data
    # (one row for a gene-condition)
    for my $geneId ( keys %results ){

        # aggregate probesets from all chips for this gene/condition to process quality summary
        my %pbsetsTemp = ();
        for my $expId ( keys %{$results{$geneId}} ) {
            my $countTemp  = 0;
            for my $chip ( keys %{$results{$geneId}->{$expId}} ) {
                my $pst = 0;
                my $pst_high_qual = 0;
                my $abs_high_qual = 0;
                # First, summarize expression call and quality for all probesets of the gene
                # pst high > pst low > abs high > abs low
                for my $pbset ( keys %{$results{$geneId}->{$expId}->{$chip}} ) {
                    # present wins in all cases: at least one probeset present leads to summarize call as present
                    if ($results{$geneId}->{$expId}->{$chip}->{$pbset}->{'call'} eq $Utils::PRESENT_CALL ||
                            $results{$geneId}->{$expId}->{$chip}->{$pbset}->{'call'} eq $Utils::MARGINAL_CALL){
                        $pst = 1;
                        # is there at least one present high quality?
                        if ($results{$geneId}->{$expId}->{$chip}->{$pbset}->{'quality'} eq $Utils::HIGH_QUAL){
                            $pst_high_qual = 1;
                        }
                    }
                    if ($results{$geneId}->{$expId}->{$chip}->{$pbset}->{'call'} eq $Utils::ABSENT_CALL
                            and $results{$geneId}->{$expId}->{$chip}->{$pbset}->{'quality'} eq $Utils::HIGH_QUAL){
                        $abs_high_qual = 1;
                    }
                }
                # add the summarize call by chip to the hash of all calls seen for this gene/experiment/condition
                if ($pst_high_qual eq 1){
                    $pbsetsTemp{$expId}->{$countTemp}->{'call'}    = $Utils::PRESENT_CALL;
                    $pbsetsTemp{$expId}->{$countTemp}->{'quality'} = $Utils::HIGH_QUAL;
                }
                elsif ($pst eq 1){
                    $pbsetsTemp{$expId}->{$countTemp}->{'call'}    = $Utils::PRESENT_CALL;
                    $pbsetsTemp{$expId}->{$countTemp}->{'quality'} = $Utils::LOW_QUAL;
                }
                elsif ($abs_high_qual eq 1){
                    $pbsetsTemp{$expId}->{$countTemp}->{'call'}    = $Utils::ABSENT_CALL;
                    $pbsetsTemp{$expId}->{$countTemp}->{'quality'} = $Utils::HIGH_QUAL;
                }
                else {
                    $pbsetsTemp{$expId}->{$countTemp}->{'call'}    = $Utils::ABSENT_CALL;
                    $pbsetsTemp{$expId}->{$countTemp}->{'quality'} = $Utils::LOW_QUAL;
                }
                $countTemp++;
                #TODO: the above code could be added as a function in insert_expression_utils.pl, receiving as argument the %{$pbsets_summary{$gene}->{'pbsets'}} hash.
            }
        }

        # Now, for each gene and condition, we calculate the final quality scores
        my ($reasonForExclusion, $summary) = Utils::summarizeExperimentCallAndQuality(\%pbsetsTemp);

        my $expressionId   = undef;

        # insert or update the expression table if not only undefined calls
        if ( $reasonForExclusion eq $Utils::CALL_NOT_EXCLUDED ){
            if ( $debug ){
                print "INSERT INTO expression (bgeeGeneId, conditionId) VALUES ($geneId, $exprMappedConditionId)...\n";
            } else {
                $expressionId = add_expression($geneId, $exprMappedConditionId);
            }
        } else {
            warn "No calls available for geneId $geneId in exprMappedConditionId $exprMappedConditionId\n";
        }

        # Now update the related affymetrixProbesets
        if ( $debug ){
            my $printExprId = 'notRetrievedForLog';
            print "UPDATE affymetrixProbeset AS t1 ",
                  "INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId ",
                  "INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId ",
                  "SET expressionId = $printExprId, reasonForExclusion = $reasonForExclusion ",
                  "WHERE t1.bgeeGeneId = $geneId AND t3.exprMappedConditionId = $exprMappedConditionId ",
                  "and t1.reasonForExclusion != '".$Utils::EXCLUDED_FOR_PRE_FILTERED."'\n";
            if ( $reasonForExclusion eq $Utils::CALL_NOT_EXCLUDED ) {
                for my $expId ( keys %{ $summary } ) {
                    print "INSERT INTO microarrayExperimentExpression ",
                               "(expressionId, microarrayExperimentId, ",
                               "presentHighMicroarrayChipCount, presentLowMicroarrayChipCount, ",
                               "absentHighMicroarrayChipCount, absentLowMicroarrayChipCount, ",
                               "microarrayExperimentCallDirection, microarrayExperimentCallQuality) ",
                           "VALUES ($printExprId, $expId, ",
                               "$summary->{$expId}->{'pstHighEvidenceCount'}, ",
                               "$summary->{$expId}->{'pstLowEvidenceCount'}, ",
                               "$summary->{$expId}->{'absHighEvidenceCount'}, ",
                               "$summary->{$expId}->{'absLowEvidenceCount'}, ",
                               "$summary->{$expId}->{'expCall'}, $summary->{$expId}->{'expCallQuality'})\n";
                }
            }
        } else {
            $updResult->execute($expressionId, $reasonForExclusion, $geneId, $exprMappedConditionId)
                or die $updResult->errstr;
            if ( $reasonForExclusion eq $Utils::CALL_NOT_EXCLUDED ) {
                for my $expId ( keys %{ $summary } ) {
                    $insExpSummary->execute($expressionId, $expId,
                        $summary->{$expId}->{'pstHighEvidenceCount'},
                        $summary->{$expId}->{'pstLowEvidenceCount'},
                        $summary->{$expId}->{'absHighEvidenceCount'},
                        $summary->{$expId}->{'absLowEvidenceCount'},
                        $summary->{$expId}->{'expCall'}, $summary->{$expId}->{'expCallQuality'})
                    or die $insExpSummary->errstr;
                }
            }
        }
    }
}

$bgee->disconnect;
print "Done\n"  if ( $debug );

exit 0;
