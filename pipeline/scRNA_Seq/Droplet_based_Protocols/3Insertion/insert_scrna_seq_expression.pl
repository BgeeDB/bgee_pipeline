#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Parallel::ForkManager;
use Data::Dumper;

# Julien Wollbrett, created April 2023

# USAGE: perl insert_scrna_seq_expression.pl -bgee=connection_string <OPTIONAL: -debug>
# After the insertion of full length single cell RNA-Seq data, this script inserts the data
# into the expression table and update the scRnaSeqFullLengthResult table.
# -debug: if provided, run in verbose mode (print the update/insert SQL queries, not executing them)

#############################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my ($bgee_connector) = ('');
my $number_threads   = '';
my ($debug)          = (0);
my %opts = ('bgee=s'            => \$bgee_connector,   # Bgee connector string
            'number_threads=s'  => \$number_threads,   # Number of threads to run in parallel
            'debug'             => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $number_threads eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\t-number_threads   number of insertion/update threads to run in parallel
\t-debug            printing the update/insert SQL queries, not executing them
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

$| = 1;


##########################################
# GET CONDITIONS STUDIED WITH RNA-SEQ    #
##########################################
print "Retrieving conditions...\n";

# multipleLibraryIndividualSample is equal to 1 only for target based RNA-Seq
my $queryConditions = $bgee->prepare('SELECT DISTINCT t2.exprMappedConditionId FROM rnaSeqLibraryAnnotatedSample AS t1
                                      INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId where 
                                      t1.multipleLibraryIndividualSample = 1');
$queryConditions->execute()  or die $queryConditions->errstr;
my @exprMappedConditions = ();
while ( my @data = $queryConditions->fetchrow_array ){
    push(@exprMappedConditions, $data[0]);
}
$queryConditions->finish;
print 'Done, ', scalar(@exprMappedConditions), " conditions retrieved.\n";

# retrieve excluded biotypes for population captures
my $queryPopCaptureExcludedAbsentCalls = 
    $bgee->prepare('SELECT rnaSeqPopulationCaptureId, geneBiotypeId FROM rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls');
$queryPopCaptureExcludedAbsentCalls->execute()  or die $queryPopCaptureExcludedAbsentCalls->errstr;
my %PopCaptureExcludedAbsentCalls = ();
while ( my @data = $queryPopCaptureExcludedAbsentCalls->fetchrow_array ){
    $PopCaptureExcludedAbsentCalls{$data[0]}{$data[1]} = 1;
}
my @polyAExcludedAbsentCalls = keys %{$PopCaptureExcludedAbsentCalls{'polyA'}};
$queryPopCaptureExcludedAbsentCalls->finish;
$bgee->disconnect;

##########################################
# SUBROUTINES                            #
##########################################
sub insert_expression {
    my ($insUpExpr, $bgee_thread, $geneId, $exprMappedConditionId) = @_;
    # insert or update the expression table
    if ( $debug ){
        print "INSERT INTO expression (bgeeGeneId, conditionId) VALUES ($geneId, $exprMappedConditionId)...\n";
        return undef;
    } else {
        $insUpExpr->execute($geneId, $exprMappedConditionId)  or die $insUpExpr->errstr;
        return $bgee_thread->{'mysql_insertid'}; # expr_id
    }
}

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

    ## PREPARE QUERIES
    # Insert/update expression
    my $insUpExpr   = $bgee_thread->prepare('INSERT INTO expression (bgeeGeneId, conditionId) VALUES (?, ?)
                                ON DUPLICATE KEY UPDATE expressionId=LAST_INSERT_ID(expressionId)');
    # Query to update rnaSeqLibraryAnnotatedSampleGeneResult with the expressionId for target based
    # to select only target based we filter on multipleLibraryIndividualSample = 1
    my $updResult = $bgee_thread->prepare("UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1
                                INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId
                                INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                SET expressionId = ?
                                WHERE bgeeGeneId = ? and t3.exprMappedConditionId = ?
                                AND t2.multipleLibraryIndividualSample = 1");
    # when absent calls are not reliable for the biotype of the geneId we first check that present calls exists for this condition/gene
    # before creating an expressionId. It avoids having expressionIds never used in the table rnaSeqLibraryAnnotatedSampleGeneResults
    my $countCallsNotExcluded = $bgee_thread->prepare("SELECT COUNT(*) from rnaSeqLibraryAnnotatedSampleGeneResult as t1
                                INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId
                                INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                WHERE bgeeGeneId = ? and t3.exprMappedConditionId = ?
                                AND t1.reasonForExclusion = '".$Utils::CALL_NOT_EXCLUDED."' AND t2.multipleLibraryIndividualSample = 1");
    my $updPresentResult = $bgee_thread->prepare("UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1
                                INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId
                                INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                SET expressionId = ?
                                WHERE bgeeGeneId = ? and t3.exprMappedConditionId = ?
                                AND t1.reasonForExclusion = '".$Utils::CALL_NOT_EXCLUDED."' AND t2.multipleLibraryIndividualSample = 1");


    # query to get all the RNA-Seq results for a condition
    my $queryResults = $bgee_thread->prepare("SELECT t1.bgeeGeneId, t4.geneBioTypeId
                                FROM rnaSeqLibraryAnnotatedSampleGeneResult AS t1
                                INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId
                                INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                INNER JOIN gene AS t4 ON t1.bgeeGeneId = t4.bgeeGeneId
                                WHERE t3.exprMappedConditionId = ?");

    print "\tconditionId: $exprMappedConditionId\n";
    # retrieve genes results for this condition
    print "\t\tRetrieving genes results...\n";
    $queryResults->execute($exprMappedConditionId)  or die $queryResults->errstr;
    my %results = ();
    while ( my @data = $queryResults->fetchrow_array ){
        $results{$data[0]}{'geneBiotypeId'} = $data[1];
    }
    print "\t\tDone, ", scalar(keys %results), ' genes retrieved. ',
          "Generating expression summary...\n";
    # now iterating the genes to insert expression data
    # (one row for a gene-condition)
    for my $geneId ( keys %results ){

        # As for Bgee 15.1 all target-based protocols have polyA as populationCapture
        #TODO: the script will have to be updated if other population captures are considered
        if (!grep(/^$results{$geneId}->{'geneBiotypeId'}$/, @polyAExcludedAbsentCalls)) {
            # insert or update the expression table
            my $expressionId = insert_expression($insUpExpr, $bgee_thread, $geneId, $exprMappedConditionId);

            # Now update the related rnaSeqLibraryAnnotatedSampleGeneResult
            if ( $debug ){
                my $printExprId = 'notRetrievedForLog';
                print "UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1",
                      "INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId",
                      "INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId ",
                      "SET expressionId=$printExprId ",
                      "WHERE bgeeGeneId=$geneId and t3.exprMappedConditionId = $exprMappedConditionId ",
                      "AND t2.multipleLibraryIndividualSample = 1\n";
            } else {
                $updResult->execute($expressionId, $geneId, $exprMappedConditionId)
                    or die $updResult->errstr;
            }
        } else {
            # count the number of calls for which expressionId have to be created
            $countCallsNotExcluded->execute($geneId, $exprMappedConditionId) or die $countCallsNotExcluded->errstr;
            my $countNotExcludedCalls = ();
            while ( my @data = $countCallsNotExcluded->fetchrow_array ){
                $countNotExcludedCalls =$data[0];
            }
            if ($countNotExcludedCalls > 0) {

                # insert or update the expression table
                my $expressionId = insert_expression($insUpExpr, $bgee_thread, $geneId, $exprMappedConditionId);

                # Now update the related rnaSeqLibraryAnnotatedSampleGeneResult
                
                if ( $debug ){
                    my $printExprId = 'notRetrievedForLog';
                    print "UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1",
                          "INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId",
                          "INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId ",
                          "SET expressionId=$printExprId ",
                          "WHERE bgeeGeneId=$geneId and t3.exprMappedConditionId = $exprMappedConditionId ",
                          "AND t1.reasonForExclusion = '".$Utils::CALL_NOT_EXCLUDED."' AND t2.multipleLibraryIndividualSample = 1\n";
                } else {
                    $updPresentResult->execute($expressionId, $geneId, $exprMappedConditionId)
                        or die $updResult->errstr;
                }
            }
        }
    }
    $insUpExpr->finish;
    $updResult->finish;
    $updPresentResult->finish;
    $countCallsNotExcluded->finish;
    $queryResults->finish;
    $bgee_thread->disconnect;
    $pm->finish;
}
$pm->wait_all_children;

print "Done\n"  if ( $debug );

exit 0;

