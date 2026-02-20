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


# J. Wollbrett, January 2026
# - One single script for all RNA-Seq datatypes (bulk, full-length, droplet). Possibility to target one or several datatypes
#   with the options -allDatatypes, -bulk, -fullLength, -droplet
# - Parallelization of the insertion process with the numberThreads option
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
my $all_datatypes = 0;
my $bulk = 0;
my $full_length = 0;
my $droplet = 0;

my %opts = ('bgee=s'                => \$bgee_connector,   # Bgee connector string
            'numberThreads=s'       => \$number_threads,   # number of threads defining the number of parallel updates in the database
            'debug'                 => \$debug,
            'allDatatypes'          => \$all_datatypes,
            'bulk'                  => \$bulk,
            'fullLength'            => \$full_length,
            'droplet'               => \$droplet,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t-bgee                 Bgee connector string
\t-numberThreads        Number of threads used to insert expression (default: 1)
\t-debug                Printing the update/insert SQL queries, not executing them
\t-allDatatypes         Process all RNA-Seq datatypes (bulk, full-length, droplet)
\t-bulk                 Process only bulk RNA-Seq datatype
\t-fullLength           Process only full-length RNA-Seq datatype
\t-droplet              Process only droplet RNA-Seq datatype
\n";
    exit 1;
}

# Validate numberThreads
if ( $number_threads < 0 || $number_threads =~ /\D/ ){
    print "\tError: numberThreads must be a non-negative integer\n";
    exit 1;
}
$number_threads = 1 if $number_threads == 0;

# first check that only one datatype option is provided
my $datatype_options_provided = 0;
$datatype_options_provided = 3  if ( $all_datatypes );
$datatype_options_provided++  if ( $bulk );
$datatype_options_provided++  if ( $full_length );
$datatype_options_provided++  if ( $droplet );
if ( $datatype_options_provided == 0 ){
    print "\tError: at least one datatype option must be provided among -allDatatypes, -bulk, -fullLength, -droplet\n";
    exit 1;
}
if ( $datatype_options_provided > 3 ){
    print "\tError: if -allDatatypes is selected then  -bulk, -fullLength, -droplet should not be selected\n";
    exit 1;
}

$| = 1;

if ($all_datatypes || $bulk ){
    print "Processing bulk RNA-Seq...\n";
    insert_expression_one_datatype(0, 0);
}
if ($all_datatypes || $full_length ){
    print "Processing full-length single-cell RNA-Seq datatype...\n";
    insert_expression_one_datatype(1, 0);
}
if ($all_datatypes || $droplet ){
    print "Processing droplet single-cell RNA-Seq datatype...\n";
    insert_expression_one_datatype(1, 1);
}

#XXX in the future we could imagine always inserting for all datatypes at the same time as it is always supposed to be done that way
#XXX we keep it like that for now in case calls filtering per datatype would be needed in the future.
sub insert_expression_one_datatype {
    my ($is_single_cell, $is_sample_multiplexing) = @_;

    # Bgee db connection
    my $bgee = Utils::connect_bgee_db($bgee_connector);

    print "Retrieving conditions for " . ($is_single_cell ? ($is_sample_multiplexing ? "droplet single-cell" :
            "full-length single-cell") : "bulk") . " RNA-Seq datatype...\n";

    # retrieve conditions for which at least one library does not have any expressionId
    my $queryConditions = $bgee->prepare('SELECT DISTINCT c.exprMappedConditionId'.
                                        ' FROM rnaSeqLibrary AS lib'.
                                        ' JOIN rnaSeqLibraryAnnotatedSample AS las'.
                                        '   ON las.rnaSeqLibraryId = lib.rnaSeqLibraryId'.
                                        ' JOIN cond AS c'.
                                        '   ON c.conditionId = las.conditionId'.
                                        ' LEFT JOIN ('.
                                        '     SELECT DISTINCT as2.rnaSeqLibraryId'.
                                        '     FROM rnaSeqLibraryAnnotatedSample AS as2'.
                                        '     JOIN rnaSeqLibraryAnnotatedSampleGeneResult AS gr'.
                                        '       ON gr.rnaSeqLibraryAnnotatedSampleId = as2.rnaSeqLibraryAnnotatedSampleId'.
                                        '      AND gr.expressionId IS NOT NULL'.
                                        ' ) AS has_expr'.
                                        '   ON has_expr.rnaSeqLibraryId = lib.rnaSeqLibraryId'.
                                        ' WHERE lib.rnaSeqTechnologyIsSingleCell = '.$is_single_cell.
                                        '   AND lib.sampleMultiplexing = '.$is_sample_multiplexing.
                                        '   AND has_expr.rnaSeqLibraryId IS NULL');
    $queryConditions->execute()  or die $queryConditions->errstr;
    my @exprMappedConditions = ();
    while ( my @data = $queryConditions->fetchrow_array ){
        push(@exprMappedConditions, $data[0]);
    }
    $queryConditions->finish;
    print 'Done, ', scalar(@exprMappedConditions), " conditions retrieved.\n";

    $bgee->disconnect;

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

        # Disable autocommit to do a single commit at the end of the condition
        $bgee_thread->{'AutoCommit'} = 0  if ( !$debug );

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
        eval {
            for my $geneId ( keys %results ){

                my $expressionId = undef;

                # insert or update the expression table
                if ( $debug ){
                    print "INSERT INTO expression (bgeeGeneId, conditionId) VALUES ($geneId, $exprMappedConditionId)...\n";
                } else {
                    $insUpExpr->execute($geneId, $exprMappedConditionId)  or die $insUpExpr->errstr;
                    $expressionId = $bgee_thread->{'mysql_insertid'};
                }

                # Now update the related rnaSeqLibraryAnnotatedSampleGeneResult
                if ( $debug ){
                    print "UPDATE rnaSeqLibraryAnnotatedSampleGeneResult...".
                                    " SET expressionId = ?".
                                    " WHERE bgeeGeneId = $geneId and exprMappedConditionId = $exprMappedConditionId\n";
                } else {
                    $updResult->execute($expressionId, $geneId, $exprMappedConditionId)
                        or die $updResult->errstr;
                }
            }
            1; # eval success
        } or do {
            my $error = $@ || 'Unknown error';
            warn "Error processing condition $exprMappedConditionId: $error\n";
            # Rollback on error
            if ( !$debug ){
                $bgee_thread->rollback();
            }
        };

        # Commit the transaction for this condition
        if ( !$debug ){
            $bgee_thread->commit()  or die $bgee_thread->errstr;
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

