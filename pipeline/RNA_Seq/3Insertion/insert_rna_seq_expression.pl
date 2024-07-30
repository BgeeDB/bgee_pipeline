#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Parallel::ForkManager;


# Frederic Bastian, created November 2012, last update Dec. 2016

# USAGE: perl insert_rna_seq_expression.pl -bgee=connection_string <OPTIONAL: -debug>
# After the insertion of bulk RNA-Seq data, this script inserts data
# into the expression table and update the rnaSeqLibraryAnnotatedSampleGeneResult table.
# -debug: if provided, run in verbose mode (print the update/insert SQL queries, not executing them)
#TODO: for Bgee 16 we have to homogeneise how expression IDs are created among RNASeq datatypes and
# refactor the code to run that step only once for all RNASeq (bulk, full-length, droplet)

#############################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my $bgee_connector = '';
my $debug          = 0;
my $number_threads = 0;
my %opts = ('bgee=s'                => \$bgee_connector,   # Bgee connector string
            'number_threads=s'      => \$number_threads,   # number of threads defining the number of parallel updates in the database
            'debug'                 => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\t-number_threads   Number of threads used to insert expression
\t-debug            printing the update/insert SQL queries, not executing them
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

$| = 1;


##########################################
# GET CONDITIONS STUDIED WITH BULK RNA-SEQ    #
##########################################
print "Retrieving conditions...\n";

# retrieve conditions only for bulk RNASeq for which at least one library does not have any expressionId
my $queryConditions = $bgee->prepare('SELECT DISTINCT t2.exprMappedConditionId'.
                                     ' FROM rnaSeqLibraryAnnotatedSample AS t1'.
                                     ' INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId'.
                                     ' INNER JOIN rnaSeqLibrary AS t3 ON t3.rnaSeqLibraryId = t1.rnaSeqLibraryId'.
                                     ' WHERE t3.rnaSeqTechnologyIsSingleCell = 0'.
                                     ' AND NOT EXISTS (SELECT 1 FROM rnaSeqLibraryAnnotatedSample AS t4'.
                                     ' INNER JOIN rnaSeqLibraryAnnotatedSampleGeneResult AS t5'.
                                     ' ON t5.rnaSeqLibraryAnnotatedSampleId = t4.rnaSeqLibraryAnnotatedSampleId'.
                                     ' WHERE t1.rnaSeqLibraryId = t4.rnaSeqLibraryId'.
                                     ' AND t5.expressionId IS NOT NULL)');
$queryConditions->execute()  or die $queryConditions->errstr;
my @exprMappedConditions = ();
while ( my @data = $queryConditions->fetchrow_array ){
    push(@exprMappedConditions, $data[0]);
}
$queryConditions->finish;
print 'Done, ', scalar(@exprMappedConditions), " conditions retrieved.\n";

# Absent calls are not considered for all biotypes depending on the protocol used to generate a library.
# biotypes not considered for absent calls depending on the protocols are described in the table rnaSeqProtocolToBiotypeExcludedAbsentCalls.
# This query update reason of exclusion for such calls if they are not already exluded for prefiltering (gene never present)
my $updExclResult = $bgee->prepare( 'UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1'.
                                    ' INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId'.
                                    ' INNER JOIN rnaSeqLibrary AS t5 ON t2.rnaSeqLibraryId = t5.rnaSeqLibraryId'.
                                    ' INNER JOIN gene as t3 on t1.bgeeGeneId = t3.bgeeGeneId'.
                                    ' INNER JOIN rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls AS t4'.
                                    ' ON t5.rnaSeqPopulationCaptureId = t4.rnaSeqPopulationCaptureId'.
                                    ' AND t3.geneBioTypeId = t4.geneBioTypeId'.
                                    ' SET t1.reasonForExclusion = "'.$Utils::EXCLUDED_FOR_ABSENT_CALLS.'"'.
                                    ' WHERE t1.pValue > 0.05'.
                                    ' AND t1.expressionId IS NULL'.
                                    # condition to remove once the script is used for all RNASeq datatypes
                                    ' AND t5.rnaSeqTechnologyIsSingleCell = 0 '.
                                    #TODO: have to keep this condition as we do not generate calls for Bgee 15.2 and
                                    # some genes were excluded for prefiltering for Bgee 15.0. No genes prefiltered anymore
                                    # so this condition can be removed for the next major release.
                                    ' AND t1.reasonForExclusion !="'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"'.
                                    #for single cell all absent calls have already been inserted with a reason for exclusion
                                    # equal to "absent calls not reliatble". This filtering avoid updating again those data
                                    ' AND t1.reasonForExclusion !="'.$Utils::EXCLUDED_FOR_ABSENT_CALLS.'"');
if ( $debug ){
    print                           'UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1'.
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
} else {
    $updExclResult->execute()  or die $updExclResult->errstr;
}
$updExclResult->finish;
$bgee->disconnect;

print "Done excluding absent calls based on biotype\n";

##########################################
# PREPARE QUERIES                        #
##########################################

# Insert/update expression
my $insUpExprQuery =        'INSERT INTO expression (bgeeGeneId, conditionId) VALUES (?, ?) '.
                            'ON DUPLICATE KEY UPDATE expressionId=LAST_INSERT_ID(expressionId)';

# Query to update rnaSeqResult with the expressionId and reasonForExclusion for present calls
my $updResultQuery =        'UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1 '.
                            ' INNER JOIN rnaSeqLibraryAnnotatedSample AS t2'.
                            ' ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId'.
                            ' INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId'.
                            ' SET expressionId = ?'.
                            ' WHERE bgeeGeneId = ? and t3.exprMappedConditionId = ?'.
                            ' AND t1.reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"'.
                            ' AND t1.expressionId IS NULL';

# query to get all the RNA-Seq results for a condition
my $queryResultsQuery =     'SELECT distinct t1.bgeeGeneId FROM rnaSeqLibraryAnnotatedSampleGeneResult AS t1'.
                            ' INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId'.
                            ' INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId'.
                            ' WHERE t3.exprMappedConditionId = ? AND t1.reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"';

##########################################
# ITERATING CONDITIONS TO INSERT DATA    #
##########################################

my $pm = new Parallel::ForkManager($number_threads);
print "Processing conditions...\n";
for my $exprMappedConditionId ( @exprMappedConditions ){
    #start parallelization
    my $pid = $pm->start and next;

    # start thread specific connection to the database
    my $bgee_thread = Utils::connect_bgee_db($bgee_connector);

    # prepare queries
    my $queryResults = $bgee_thread->prepare($queryResultsQuery);
    my $updResult = $bgee_thread->prepare($updResultQuery);
    my $insUpExpr   = $bgee_thread->prepare($insUpExprQuery);
    print "\tconditionId: $exprMappedConditionId\n";

    # retrieve genes results for this condition
    print "\t\tRetrieving genes results...\n";
    $queryResults->execute($exprMappedConditionId)  or die $queryResults->errstr;
    my %results = ();
    while ( my @data = $queryResults->fetchrow_array ){
        # $data[0] = bgeeGeneId
        $results{$data[0]} = 1;
    }
    print "\t\tDone, ", scalar(keys %results), ' genes retrieved. ',
          "Generating expression summary...\n";

    # now iterating the genes to insert expression data
    # (one row for a gene-condition)
    for my $geneId ( keys %results ){

        my $expressionId = undef;

        # insert or update the expression table if not only undefined calls
        if ( $debug ){
            print "INSERT INTO expression (bgeeGeneId, conditionId) VALUES ($geneId, $exprMappedConditionId)...\n";
        } else {
            $insUpExpr->execute($geneId, $exprMappedConditionId)  or die $insUpExpr->errstr;
            $expressionId = $bgee_thread->{'mysql_insertid'};
        }

        # Now update the related rnaSeqResults
        if ( $debug ){
            my $printExprId = 'notRetrievedForLog';
            print "UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1 ",
                  "INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ",
                  "ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId ",
                  "INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId ",
                  "SET expressionId = $printExprId ",
                  "WHERE bgeeGeneId = $geneId AND t3.exprMappedConditionId = $exprMappedConditionId ",
                  "AND t1.reasonForExclusion = \"".$Utils::CALL_NOT_EXCLUDED."\"".
                  "AND t1.expressionId IS NULL\n";

        } else {
            $updResult->execute($expressionId, $geneId, $exprMappedConditionId)
                or die $updResult->errstr;
        }
    }
    $insUpExpr->finish;
    $updResult->finish;
    $queryResults->finish;
    $bgee_thread->disconnect;
    $pm->finish;
}
$pm->wait_all_children;

print "Done inserting expression IDs\n"  if ( $debug );

exit 0;

