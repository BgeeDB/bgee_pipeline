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


# J. Wollbrett - January 2026
# - One single script for all RNA-Seq datatypes (bulk, full-length, droplet). Possibility to target one or several datatypes
#   with the options -allDatatypes, -bulkOnly, -fullLengthOnly, -dropletOnly
# - Possibility to remove pvalues for calls based on biotype with the -removePvaluesNotTargeted option
# - Parallelization of the insertion process with the numberThreads option
# - Update the reason for exclusion for genes with biotype not targeted by library preparation as "biotype not targeted"
# - this script only insert expressionId, bgeeGeneId, conditionId in the expression table. the datatype specific fields like
#   score, pvalue, weight and numberOfObs will be inserted in a next step once the ranks per library have been computed.

#############################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my $bgee_connector = '';
my $debug          = 0;
my $number_threads = 0;
my $remove_not_targeted = 0;
my $all_datatypes = 0;
my $bulk_only = 0;
my $full_length_only = 0;
my $droplet_only = 0;

my %opts = ('bgee=s'                => \$bgee_connector,   # Bgee connector string
            'numberThreads=s'       => \$number_threads,   # number of threads defining the number of parallel updates in the database
            'debug'                 => \$debug,
            'removePvaluesNotTargeted'     => \$remove_not_targeted,
            'allDatatypes'          => \$all_datatypes,
            'bulkOnly'              => \$bulk_only,
            'fullLengthOnly'        => \$full_length_only,
            'dropletOnly'           => \$droplet_only,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t-bgee                 Bgee connector string
\t-numberThreads       Number of threads used to insert expression
\t-debug                Printing the update/insert SQL queries, not executing them
\t-removePvaluesNotTargeted    Remove pValues for calls based on biotype not targeted by library preparation
\t-allDatatypes         Process all RNA-Seq datatypes (bulk, full-length, droplet)
\t-bulkOnly             Process only bulk RNA-Seq datatype
\t-fullLengthOnly       Process only full-length RNA-Seq datatype
\t-dropletOnly          Process only droplet RNA-Seq datatype
\n";
    exit 1;
}

# first check that only one datatype option is provided
my $datatype_options_provided = 0;
$datatype_options_provided++  if ( $all_datatypes );
$datatype_options_provided++  if ( $bulk_only );
$datatype_options_provided++  if ( $full_length_only );
$datatype_options_provided++  if ( $droplet_only );
if ( $datatype_options_provided > 1 ){
    print "\n\tError: only one datatype option can be provided among -allDatatypes, -bulkOnly, -fullLengthOnly, -dropletOnly\n\n";
    exit 1;
}

$| = 1;

if ($all_datatypes || $bulk_only ){
    print "Processing bulk RNA-Seq datatype...\n";
    insert_expression_one_datatype(0, 0, $remove_not_targeted);
}
if ($all_datatypes || $full_length_only ){
    print "Processing full-length single-cell RNA-Seq datatype...\n";
    insert_expression_one_datatype(1, 0, $remove_not_targeted);
}
if ($all_datatypes || $droplet_only ){
    print "Processing droplet single-cell RNA-Seq datatype...\n";
    insert_expression_one_datatype(1, 1, $remove_not_targeted);
}

#XXX in the future we could imagine always inserting for all datatypes at the same time as it is always supposed to be done that way
#XXX we keep it like that for now in case calls filtering per datatype would be needed in the future.
sub insert_expression_one_datatype {
    my ($is_single_cell, $is_sample_multiplexing, $remove_not_targeted) = @_;

    # Bgee db connection
    my $bgee = Utils::connect_bgee_db($bgee_connector);

    print "Retrieving conditions for " . ($is_single_cell ? ($is_sample_multiplexing ? "droplet single-cell" :
            "full-length single-cell") : "bulk") . " RNA-Seq datatype...\n";

    # retrieve conditions only for bulk RNASeq for which at least one library does not have any expressionId
    my $queryConditions = $bgee->prepare('SELECT DISTINCT c.exprMappedConditionId'.
                                        ' FROM rnaSeqLibrary AS l'.
                                        ' JOIN rnaSeqLibraryAnnotatedSample AS s'.
                                        '   ON s.rnaSeqLibraryId = l.rnaSeqLibraryId'.
                                        ' JOIN cond AS c'.
                                        '   ON c.conditionId = s.conditionId'.
                                        ' LEFT JOIN ('.
                                        '     SELECT DISTINCT s2.rnaSeqLibraryId'.
                                        '     FROM rnaSeqLibraryAnnotatedSample AS s2'.
                                        '     JOIN rnaSeqLibraryAnnotatedSampleGeneResult AS gr'.
                                        '       ON gr.rnaSeqLibraryAnnotatedSampleId = s2.rnaSeqLibraryAnnotatedSampleId'.
                                        '      AND gr.expressionId IS NOT NULL'.
                                        ' ) AS has_expr'.
                                        '   ON has_expr.rnaSeqLibraryId = l.rnaSeqLibraryId'.
                                        ' WHERE l.rnaSeqTechnologyIsSingleCell = '.$is_single_cell.
                                        '   AND l.sampleMultiplexing = '.$is_sample_multiplexing.
                                        '   AND has_expr.rnaSeqLibraryId IS NULL');
    $queryConditions->execute()  or die $queryConditions->errstr;
    my @exprMappedConditions = ();
    while ( my @data = $queryConditions->fetchrow_array ){
        push(@exprMappedConditions, $data[0]);
    }
    $queryConditions->finish;
    print 'Done, ', scalar(@exprMappedConditions), " conditions retrieved.\n";

    # calls are not considered for all biotypes depending on the protocol used to generate a library.
    # biotypes not considered depending on the protocols are described in the table rnaSeqProtocolToBiotypeExcludedAbsentCalls.
    # This query update reason of exclusion for such calls
    my $updExclResult = $bgee->prepare( 'UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1'.
                                        ' INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId'.
                                        ' INNER JOIN rnaSeqLibrary AS t5 ON t2.rnaSeqLibraryId = t5.rnaSeqLibraryId'.
                                        ' INNER JOIN gene as t3 on t1.bgeeGeneId = t3.bgeeGeneId'.
                                        ' INNER JOIN rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls AS t4'.
                                        ' ON t5.rnaSeqPopulationCaptureId = t4.rnaSeqPopulationCaptureId'.
                                        ' AND t3.geneBioTypeId = t4.geneBioTypeId'.
                                        ' SET t1.reasonForExclusion = "'.$Utils::BIOTYPE_NOT_TARGETED.'"'.
                                        if ( $remove_not_targeted ) {
                                            ' , t1.pValue = NULL '
                                        }.
                                        ' WHERE t1.expressionId IS NULL'.
                                        ' AND t5.rnaSeqTechnologyIsSingleCell = '.$is_single_cell.''.
                                        ' AND t5.sampleMultiplexing = '.$is_sample_multiplexing.''.
                                        ' AND t1.reasonForExclusion ="'.$Utils::CALL_NOT_EXCLUDED.'"');
    if ( $debug ){
        print                           'UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1'.
                                        ' INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId'.
                                        ' INNER JOIN rnaSeqLibrary AS t5 ON t2.rnaSeqLibraryId = t5.rnaSeqLibraryId'.
                                        ' INNER JOIN gene as t3 on t1.bgeeGeneId = t3.bgeeGeneId'.
                                        ' INNER JOIN rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls AS t4'.
                                        ' ON t5.rnaSeqPopulationCaptureId = t4.rnaSeqPopulationCaptureId'.
                                        ' AND t3.geneBioTypeId = t4.geneBioTypeId'.
                                        ' SET t1.reasonForExclusion = "'.$Utils::BIOTYPE_NOT_TARGETED.'"'.
                                        if ( $remove_not_targeted ) {
                                            ' , t1.pValue = NULL '
                                        }.
                                        ' WHERE t1.expressionId IS NULL'.
                                        ' AND t5.rnaSeqTechnologyIsSingleCell = '.$is_single_cell.
                                        ' AND t5.sampleMultiplexing = '.$is_sample_multiplexing.
                                        ' AND t1.reasonForExclusion ="'.$Utils::CALL_NOT_EXCLUDED.'"' . "\n";
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

    # Query to update rnaSeqLibraryAnnotatedSampleGeneResult with the expressionId and reasonForExclusion calls
    my $updResultQuery =        'UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1 '.
                                ' INNER JOIN rnaSeqLibraryAnnotatedSample AS t2'.
                                ' ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId'.
                                ' INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId'.
                                ' INNER JOIN rnaSeqLibrary AS t4 ON t2.rnaSeqLibraryId = t4.rnaSeqLibraryId'.
                               ' SET expressionId = ?'.
                                ' WHERE bgeeGeneId = ? and t3.exprMappedConditionId = ?'.
                                ' AND t1.reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"'.
                                ' AND t4.rnaSeqTechnologyIsSingleCell = '.$is_single_cell.
                                ' AND t4.sampleMultiplexing = '.$is_sample_multiplexing.
                                ' AND t1.expressionId IS NULL';

    # query to get all the RNA-Seq results for a condition
    my $queryResultsQuery =     'SELECT distinct t1.bgeeGeneId FROM rnaSeqLibraryAnnotatedSampleGeneResult AS t1'.
                                ' INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId'.
                                ' INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId'.
                                ' INNER JOIN rnaSeqLibrary AS t4 ON t2.rnaSeqLibraryId = t4.rnaSeqLibraryId'.
                                ' WHERE t3.exprMappedConditionId = ? AND t1.reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"'.
                                ' AND t4.rnaSeqTechnologyIsSingleCell = '.$is_single_cell.
                                ' AND t4.sampleMultiplexing = '.$is_sample_multiplexing;

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

        # retrieve genes results for this condition and the given datatype
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

            # Now update the related rnaSeqLibraryAnnotatedSampleGeneResult
            if ( $debug ){
                my $printExprId = 'notRetrievedForLog';
                print "UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1 \n".
                                " INNER JOIN rnaSeqLibraryAnnotatedSample AS t2\n".
                                " ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId\n".
                                " INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId\n".
                                " INNER JOIN rnaSeqLibrary AS t4 ON t2.rnaSeqLibraryId = t4.rnaSeqLibraryId\n".
                                " SET expressionId = ?\n".
                                " WHERE bgeeGeneId = ? and t3.exprMappedConditionId = ?\n".
                                " AND t1.reasonForExclusion = \"".$Utils::CALL_NOT_EXCLUDED."\"".
                                " AND t4.rnaSeqTechnologyIsSingleCell = ".$is_single_cell."\n".
                                " AND t4.sampleMultiplexing = ".$is_sample_multiplexing."\n".
                                " AND t1.expressionId IS NULL\n;";

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

    print "Done inserting expression IDs for $is_single_cell and $is_sample_multiplexing\n"  if ( $debug );

}

exit 0;

