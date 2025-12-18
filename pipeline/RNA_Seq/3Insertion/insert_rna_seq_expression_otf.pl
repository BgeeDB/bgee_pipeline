#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Parallel::ForkManager;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../..";
use Utils;


# Frederic Bastian, created November 2012, last update Dec. 2016

# USAGE: perl insert_rna_seq_expression.pl -bgee=connection_string <OPTIONAL: -debug>
# After the insertion of bulk RNA-Seq data, this script inserts data
# into the expression table and update the rnaSeqLibraryAnnotatedSampleGeneResult table.
# -debug: if provided, run in verbose mode (print the update/insert SQL queries, not executing them)
#TODO: for Bgee 16 we have to homogeneise how expression IDs are created among RNASeq datatypes and
# refactor the code to run that step only once for all RNASeq (bulk, full-length, droplet)

# Julien Wollbrett, updated Nov. 2025

# Added new columns for bulk RNA-Seq expression data insertion:
# now insert meanScore, meanPValue, weight and numberObs (number of observation) for each observed expression call
# this implementation has been done as part of the on the fly propagation revamp.
# The new columns will be used to compute global expression calls/score without querying the datatype specific tables.
# USAGE: perl insert_rna_seq_expression_otf.pl -bgee=connection_string -number_threads=N <OPTIONAL: -debug>

# steps:
# retrieve all processed expression values for the selected datatype (here bulk RNA-Seq but could be also used for single cell)
# insert if not already exist expression data for all unique combination of conditionId/geneId 
# calculate meanRank, meanPValue, weight and numberObs for each combination conditionId/geneId
# update expression table with these new values (do not check if values already present, just update them)
# update rnaSeqLibraryAnnotatedSampleGeneResult table with the corresponding expressionId

# Define arguments & their default value
my $bgee_connector = '';
my $debug          = 0;
my $number_threads = 1;
my %opts = (
    'bgee=s'           => \$bgee_connector,
    'number_threads=i' => \$number_threads,
    'debug'            => \$debug,
);

my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\t-number_threads   Number of threads used to insert expression (default 1)
\t-debug            printing the update/insert SQL queries, not executing them
\n";
    exit 1;
}

# Validate number of threads
if ($number_threads < 1) {
    warn "number_threads < 1; defaulting to 1\n";
    $number_threads = 1;
}

# Main (single) DB connection only used to fetch the list of conditions and do the absent-call update
my $bgee = Utils::connect_bgee_db($bgee_connector);
$bgee->{RaiseError} = 1;

$| = 1;

print "Retrieving conditions...\n";

my $queryConditions = $bgee->prepare(
    'SELECT DISTINCT t2.exprMappedConditionId'.
    ' FROM rnaSeqLibraryAnnotatedSample AS t1'.
    ' INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId'.
    ' INNER JOIN rnaSeqLibrary AS t3 ON t3.rnaSeqLibraryId = t1.rnaSeqLibraryId'.
    ' WHERE t3.rnaSeqTechnologyIsSingleCell = 0'.
    ' AND NOT EXISTS (SELECT 1 FROM rnaSeqLibraryAnnotatedSample AS t4'.
    ' INNER JOIN rnaSeqLibraryAnnotatedSampleGeneResult AS t5'.
    ' ON t5.rnaSeqLibraryAnnotatedSampleId = t4.rnaSeqLibraryAnnotatedSampleId'.
    ' WHERE t1.rnaSeqLibraryId = t4.rnaSeqLibraryId'.
    ' AND t5.expressionId IS NOT NULL)'
);
$queryConditions->execute() or die $queryConditions->errstr;
my @exprMappedConditions;
while ( my @data = $queryConditions->fetchrow_array ){
    push @exprMappedConditions, $data[0];
}
$queryConditions->finish;
print 'Done, ', scalar(@exprMappedConditions), " conditions retrieved.\n";

# Absent-call exclusion update (run on main handle)
my $updExclSql = 
    'UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1'.
    ' INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId'.
    ' INNER JOIN rnaSeqLibrary AS t5 ON t2.rnaSeqLibraryId = t5.rnaSeqLibraryId'.
    ' INNER JOIN gene as t3 on t1.bgeeGeneId = t3.bgeeGeneId'.
    ' INNER JOIN rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls AS t4'.
    ' ON t5.rnaSeqPopulationCaptureId = t4.rnaSeqPopulationCaptureId'.
    ' AND t3.geneBioTypeId = t4.geneBioTypeId'.
    ' SET t1.reasonForExclusion = "'.$Utils::EXCLUDED_FOR_ABSENT_CALLS.'"'.
    ' WHERE t1.pValue > 0.05'.
    ' AND t1.expressionId IS NULL'.
    ' AND t5.rnaSeqTechnologyIsSingleCell = 0 '.
    ' AND t1.reasonForExclusion !="'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"'.
    ' AND t1.reasonForExclusion !="'.$Utils::EXCLUDED_FOR_ABSENT_CALLS.'"';

if ($debug) {
    print $updExclSql, "\n";
} else {
    my $updExclResult = $bgee->prepare($updExclSql);
    $updExclResult->execute() or die $updExclResult->errstr;
    $updExclResult->finish;
}
print "Done excluding absent calls based on biotype\n";

# Disconnect the main handle (children will open their own connections)
$bgee->disconnect;

##########################################
# PREPARE SQL (used by child processes)  #
##########################################
my $insUpExprQuery = 'INSERT INTO expression (bgeeGeneId, conditionId) VALUES (?, ?) '.
                     'ON DUPLICATE KEY UPDATE expressionId=LAST_INSERT_ID(expressionId)';

my $updResultQuery = 'UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1 '.
                     ' INNER JOIN rnaSeqLibraryAnnotatedSample AS t2'.
                     ' ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId'.
                     ' INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId'.
                     ' SET expressionId = ?'.
                     ' WHERE bgeeGeneId = ? and t3.exprMappedConditionId = ?'.
                     ' AND t1.reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"'.
                     ' AND t1.expressionId IS NULL';

my $queryResultsQuery = 'SELECT distinct t1.bgeeGeneId, t1.rnaSeqLibraryAnnotatedSampleId, t1.rawRank, t1.pValue, t4.mappedReadsCount'.
                        ' t2.rnaSeqLibraryAnnotatedSampleMaxRank, t2.rnaSeqLibraryAnnotatedSampleDistinctRankCount'
                        ' FROM rnaSeqLibraryAnnotatedSampleGeneResult AS t1'.
                        ' INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId'.
                        ' INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId'.
                        ' INNER JOIN rnaSeqLibrary AS t4 ON t2.rnaSeqLibraryId = t4.rnaSeqLibraryId'.
                        ' WHERE t3.exprMappedConditionId = ? AND t1.reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"'.
                        ' AND t4.rnaSeqTechnologyIsSingleCell = 0';

my $updExprMetricsQuery = 'UPDATE expression SET bulkScore = ?, bulkPValue = ?, bulkWeight = ?, bulkNumberObs = ? WHERE expressionId = ?';

##########################################
# PROCESS CONDITIONS IN PARALLEL         #
##########################################
my $pm = Parallel::ForkManager->new($number_threads);
print "Processing conditions with $number_threads thread(s)...\n";

CONDITION:
for my $exprMappedConditionId (@exprMappedConditions) {
    my $pid = $pm->start and next CONDITION;  # parent continues loop

    # CHILD PROCESS START
    my $bgee_thread = Utils::connect_bgee_db($bgee_connector);
    # ensure manual commits in child
    eval {
        $bgee_thread->{AutoCommit} = 0;

        # prepare child statements
        my $queryResults = $bgee_thread->prepare($queryResultsQuery);
        my $updResult = $bgee_thread->prepare($updResultQuery);
        my $insUpExpr = $bgee_thread->prepare($insUpExprQuery);
        my $updExprMetrics = $bgee_thread->prepare($updExprMetricsQuery);

        print "\t[child $$] conditionId: $exprMappedConditionId\n" unless $debug;

        # retrieve genes results for this condition
        print "\t\t[child $$] Retrieving genes results...\n" unless $debug;
        $queryResults->execute($exprMappedConditionId) or die $queryResults->errstr;
        my %results;

        while ( my @data = $queryResults->fetchrow_array ){
            $results{$data[0]}{$data[1]}{'rawRank'}       = $data[2];
            $results{$data[0]}{$data[1]}{'pValue'}        = $data[3];
            $results{$data[0]}{$data[1]}{'mappedReadsCount'} = $data[4];
            $results{$data[0]}{$data[1]}{'maxRank'} = $data[5];
            $results{$data[0]}{$data[1]}{'distinctRankCount'} = $data[6];
        }
        $queryResults->finish;
        print "\t\t[child $$] Done, ", scalar(keys %results), " genes retrieved. Generating expression summary...\n" unless $debug;

        for my $geneId ( keys %results ) {
            # Insert new expression ID or retrieve already existing one
            my $expressionId;
            if ($debug) {
                print "INSERT INTO expression (bgeeGeneId, conditionId) VALUES ($geneId, $exprMappedConditionId)...\n";
            } else {
                $insUpExpr->execute($geneId, $exprMappedConditionId) or die $insUpExpr->errstr;
                # DBD::mysql populates mysql_insertid when using LAST_INSERT_ID() technique
                $expressionId = $bgee_thread->{mysql_insertid};
                if (!defined $expressionId || $expressionId == 0) {
                    die "Failed to retrieve expressionId for gene $geneId, condition $exprMappedConditionId";
                }
            }

            my %weightedExprMetrics = calculate_expression_summary(\%{ $results{$geneId} }, $geneId, $exprMappedConditionId);

            if ($debug) {
                print "UPDATE expression SET meanScore = ".$weightedExprMetrics{'meanScore'}.
                      ", meanPValue = ".$weightedExprMetrics{'meanPValue'}.
                      ", weight = ".$weightedExprMetrics{'mappedReadsCount'}.
                      ", numberObs = ".$weightedExprMetrics{'numberObs'}.
                      " WHERE expressionId = $expressionId\n";
                print "UPDATE rnaSeqLibraryAnnotatedSampleGeneResult ... SET expressionId = notRetrievedForLog WHERE bgeeGeneId = $geneId AND exprMappedConditionId = $exprMappedConditionId ...\n";
            } else {
                $updExprMetrics->execute(
                    $weightedExprMetrics{'meanRank'},
                    $weightedExprMetrics{'meanPValue'},
                    $weightedExprMetrics{'mappedReadsCount'},
                    $weightedExprMetrics{'numberObs'},
                    $expressionId
                ) or die $updExprMetrics->errstr;

                $updResult->execute($expressionId, $geneId, $exprMappedConditionId) or die $updResult->errstr;
            }
        }

        # commit child transaction
        if (!$debug) {
            $bgee_thread->commit or die "Commit failed in child $$: $!";
        }

        # finish prepared statements
        $insUpExpr->finish;
        $updResult->finish;
        $updExprMetrics->finish;

        $bgee_thread->disconnect;
        print "\t[child $$] Done processing conditionId: $exprMappedConditionId\n" unless $debug;

        1;
    } or do {
        my $err = $@ || 'Unknown error';
        warn "\t[child $$] Error processing condition $exprMappedConditionId: $err\n";
        # try to rollback if possible
        eval { $bgee_thread->rollback if $bgee_thread && $bgee_thread->{AutoCommit} == 0 };
        eval { $bgee_thread->disconnect if $bgee_thread };
        # exit child with non-zero so parent can detect failure via run_on_finish if desired
        $pm->finish(1);
    };

    # normal child exit
    $pm->finish(0);
}

$pm->wait_all_children;

print "Done inserting expression IDs\n";

exit 0;

##########################################
# SUBROUTINES                            #
##########################################
sub calculate_expression_summary {
    my ($geneData, $geneId, $conditionId) = @_;
    my %summary;

    my $score = 0;
    my $pValue = 0;
    my $weight = 0;
    my $obsCount = 0;

    # for now we decided to use the distinct rank count as the weighting factor but at some point we descussed about using
    # the mapped reads count instead
    foreach my $sampleId (keys %$geneData) {
        my $data = $geneData->{$sampleId};

        my $rawRank = $data->{'rawRank'};
        my $distinctRanks = $data->{'distinctRankCount'};
        my $maxRank = $data->{'maxRank'};
        my $obsPValue = $data->{'pValue'};

        die "rawRank est null" if !defined $rawRank;
        die "distinctRanks est null" if !defined $distinctRanks;
        die "maxRank est null" if !defined $maxRank;
        die "obsPValue est null" if !defined $obsPValue;

        my $range = $maxRank - 1;
        my $adjustedRank = $rawRank - 1;

        # calculate the score for that sample
        my $expressionScore = 100 - ($adjustedRank * 99 / $range);

        # weight the score by the distinct rank count
        $weight += $distinctRanks;
        $score += $expressionScore * $weight;
        $pValue += $obsPValue;
        $obsCount += 1;


    }

    if ($obsCount == 0) {
        die "No observations found for geneId $geneId in conditionId $conditionId";
    }
    $summary{'score'} = $score / $weight;
    #XXX: the pValue has to be multiplied by 2 at some points. For now we decided not to multiply by 2 in the datatbase but to do it on the fly when retrieving the value for propagation
    #TODO: confirm that the multiplication will be done on the fly
    $summary{'pvalue'} = $pValue / $obsCount;

    $summary{'weight'} = $weight;
    $summary{'numberObs'} = $obsCount;

    return %summary;
}