#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Parallel::ForkManager;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/..";
use Utils;


# -----------------------------------------------------------------------------
# Script: insert_rna_seq_expression_otf.pl
# Author: Julien Wollbrett
# Created: Mar. 2026
#
# Description:
#   This script updates the 'expression' table in the Bgee database with summarized
#   RNA-Seq expression metrics (p-value, score, weight, numberObs) for each geneId/conditionId
#   combination and for each RNA-Seq datatype (bulk, full-length single-cell, droplet-based single-cell).
#
#   It must be run after the 'expression' table is already populated with all bgeeGeneId/conditionId
#   combinations for all RNA-Seq datatypes.
#
# Usage Example:
#   perl insert_rna_seq_expression_otf.pl -bgee=<BGEECMD> -numberThreads=4 -bulk
#
#   Options:
#     -bgee             Bgee connector string (required)
#     -numberThreads    Number of threads to use for parallel processing (default: 1)
#     -debug            Print SQL queries instead of executing them
#     -bulk             Process bulk RNA-Seq data
#     -fullLength       Process full-length single-cell RNA-Seq data
#     -droplet          Process droplet-based single-cell RNA-Seq data
#     -allDatatypes     Process all RNA-Seq datatypes (bulk, full-length sc, droplet-based sc)
#     -processAll       Reset the columns of the selected datatypes to NULL before processing, in order to
#                       reprocess all conditions. Without it, only conditions with missing scores are processed
#
# Steps performed by the script:
#   1. Retrieve all relevant conditionIds from the 'expression' table (exprMappedConditionId).
#   2. For each conditionId, retrieve associated rnaSeqLibraryAnnotatedSampleIds.
#   3. Summarize expression info (score, pValue, weight, numberObs) for each geneId based on all
#      rnaSeqLibraryAnnotatedSampleIds linked to the conditionId.
#   4. Update the 'expression' table with the summarized metrics for each geneId/conditionId/datatype.
#
# Dependencies:
#   - Perl modules: strict, warnings, diagnostics, Parallel::ForkManager, Getopt::Long, FindBin
#   - Custom module: Utils (must be in the parent directory)
#   - Database: Bgee schema with required tables and columns
# -----------------------------------------------------------------------------


###############################################################################
# Argument parsing and default values
###############################################################################

my $bgee_connector = '';
my $debug          = 0;
my $number_threads = 1;
my $all_datatypes  = 0;
my $bulk           = 0;
my $full_length    = 0;
my $droplet        = 0;
my $process_all    = 0;

# Column prefixes for different RNA-Seq datatypes
my $BULK_COLUMN_PREFIX        = "bulk";
my $FULL_LENGTH_COLUMN_PREFIX = "fullLength";
my $DROPLET_COLUMN_PREFIX     = "droplet";


my %opts = (
    'bgee=s'           => \$bgee_connector,
    'numberThreads=i'  => \$number_threads,
    'debug'            => \$debug,
    'allDatatypes'     => \$all_datatypes,
    'bulk'             => \$bulk,
    'fullLength'       => \$full_length,
    'droplet'          => \$droplet,
    'processAll'       => \$process_all,
);

my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
	e.g. $0  -bgee=\$(BGEECMD)
	-bgee             Bgee connector string (required)
	-numberThreads    Number of threads used to insert expression (default 1)
	-debug            Print SQL queries, do not execute them
	-bulk             Process bulk RNA-Seq data
	-fullLength       Process full-length single-cell RNA-Seq data
	-droplet          Process droplet-based single-cell RNA-Seq data
	-allDatatypes     Process all RNA-Seq datatypes (bulk, full-length sc, droplet-based sc)
    -processAll       Reset pValue, score, weight and numberObs to NULL for the selected datatypes before processing,
                      so that all conditions are reprocessed. Without it, only conditions with missing scores are processed
                      (relaunch without -processAll to resume an interrupted run).
\n";
    exit 1;
}


# Validate number of threads
if ($number_threads < 1) {
    warn "numberThreads < 1; defaulting to 1\n";
    $number_threads = 1;
}

# Ensure only one datatype option is provided
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
    print "\tError: if -allDatatypes is selected then -bulk, -fullLength, -droplet should not be selected\n";
    exit 1;
}

# Datatypes to process: [column prefix, rnaSeqTechnologyIsSingleCell, sampleMultiplexing]
my @datatypesToProcess = ();
push @datatypesToProcess, [$BULK_COLUMN_PREFIX,        0, 0]  if ( $all_datatypes || $bulk );
push @datatypesToProcess, [$FULL_LENGTH_COLUMN_PREFIX, 1, 0]  if ( $all_datatypes || $full_length );
push @datatypesToProcess, [$DROPLET_COLUMN_PREFIX,     1, 1]  if ( $all_datatypes || $droplet );

# -----------------------------------------------------------------------------
# Database connection and data retrieval
# -----------------------------------------------------------------------------

# Main (single) DB connection only used to fetch the list of conditions and do the absent-call update
my $bgee = Utils::connect_bgee_db($bgee_connector);
$bgee->{RaiseError} = 1;

$| = 1;

# Retrieve max rank values for each species/datatype combination
print "Retrieving species/datatype max ranks...\n";
my %speciesMaxRanks = ();
my $queryMaxRanks = $bgee->prepare(
    'SELECT speciesId, rnaSeqPopulationCaptureId, rnaSeqTechnologyIsSingleCell, sampleMultiplexing, maxRank '.
    'FROM rnaSeqPopulationCaptureSpeciesMaxRank'
);
$queryMaxRanks->execute() or die $queryMaxRanks->errstr;
while ( my @data = $queryMaxRanks->fetchrow_array ){
    $speciesMaxRanks{$data[0]}{$data[1]}{$data[2]}{$data[3]} = $data[4];
}
$queryMaxRanks->finish;

# -----------------------------------------------------------------------------
# With -processAll, all rows of the datatypes to process are first reset to NULL.
# The list of conditions to process is then always retrieved with the same
# "IS NULL" filter, whether -processAll was provided or not. It means that if the
# run is interrupted (e.g. some queries end with an out of time error because too
# many threads overloaded the database), it can be relaunched without -processAll
# and will continue from where it stopped.
# -----------------------------------------------------------------------------
if ( $process_all ){
    for my $datatype (@datatypesToProcess) {
        my $prefix = $datatype->[0];
        my $resetQuery = 'UPDATE expression SET '.
            $prefix.'PValue = NULL, '.$prefix.'Score = NULL, '.
            $prefix.'Weight = NULL, '.$prefix.'NumberObs = NULL '.
            'WHERE '.$prefix.'PValue IS NOT NULL OR '.$prefix.'Score IS NOT NULL OR '.
            $prefix.'Weight IS NOT NULL OR '.$prefix.'NumberObs IS NOT NULL';
        if ( $debug ){
            print "[DEBUG] $resetQuery\n";
        } else {
            print "Resetting $prefix columns to NULL...\n";
            my $resetRows = $bgee->do($resetQuery);
            print "Done, $resetRows rows reset.\n";
        }
    }
}

# Retrieve all conditionIds that need to be processed
print "Retrieving conditions...\n";

my $sqlConditionFilter = 'WHERE '.join(' OR ', map { $_->[0].'Score IS NULL' } @datatypesToProcess);
my $queryConditions = $bgee->prepare(
    'SELECT DISTINCT conditionId FROM expression ' . $sqlConditionFilter
);
$queryConditions->execute() or die $queryConditions->errstr;
my @exprMappedConditions;
while ( my @data = $queryConditions->fetchrow_array ){
    push @exprMappedConditions, $data[0];
}
$queryConditions->finish;
print 'Done, ', scalar(@exprMappedConditions), " conditions retrieved.\n";

# SQL query to retrieve gene/sample results for a given conditionId and datatype
my $sqlQuery = 'SELECT t4.bgeeGeneId, t3.speciesId, t1.rnaSeqLibraryAnnotatedSampleMaxRank, '.
                't1.rnaSeqLibraryAnnotatedSampleDistinctRankCount, t2.rnaSeqPopulationCaptureId, t4.pValue, t4.rawRank '.
                'FROM rnaSeqLibraryAnnotatedSample AS t1 '.
                'INNER JOIN rnaSeqLibrary AS t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId '.
                'INNER JOIN cond AS t3 ON t1.conditionId = t3.conditionId '.
                'INNER JOIN rnaSeqLibraryAnnotatedSampleGeneResult AS t4  '.
                'ON t1.rnaSeqLibraryAnnotatedSampleId = t4.rnaSeqLibraryAnnotatedSampleId '.
                'WHERE t3.exprMappedConditionId = ? AND t4.expressionId IS NOT NULL '.
                'AND t2.rnaSeqTechnologyIsSingleCell = ? AND t2.sampleMultiplexing = ? '.
                'ORDER BY t4.bgeeGeneId';

# -----------------------------------------------------------------------------
# Parallel processing of conditions
# -----------------------------------------------------------------------------

my $pm = Parallel::ForkManager->new($number_threads);

for my $exprMappedConditionId (@exprMappedConditions) {
    print "Processing conditionId: $exprMappedConditionId\n";
    # Forks and returns the pid for the child
    my $pid = $pm->start and next;
    my $bgee_thread = Utils::connect_bgee_db($bgee_connector);
    # For each datatype, call the insert_expression_scores subroutine
    for my $datatype (@datatypesToProcess) {
        insert_expression_scores($bgee_thread, $exprMappedConditionId, $datatype->[1], $datatype->[2],
                                 $sqlQuery, \%speciesMaxRanks);
    }
    $bgee_thread->disconnect;
    $pm->finish; # Terminates the child process
}
$pm->wait_all_children;

# -----------------------------------------------------------------------------
# Subroutine: insert_expression_scores
#   For a given conditionId and datatype, summarizes gene expression metrics and
#   updates the 'expression' table.
# -----------------------------------------------------------------------------
sub insert_expression_scores {
    my ($bgee_thread, $exprMappedConditionId, $isSingleCell, $isSampleMultiplexing, $sqlQuery, $speciesMaxRanks) = @_;

    my $sth = $bgee_thread->prepare($sqlQuery);
    $sth->execute($exprMappedConditionId, $isSingleCell, $isSampleMultiplexing) or die $sth->errstr;
    my $previousGeneId = 0;
    my %results;
    # Aggregate results for each geneId
    while ( my @data = $sth->fetchrow_array ){
        my $geneId = $data[0];
        my $speciesId = $data[1];
        my $sampleMaxRank = $data[2];
        my $sampleDistinctRankCount = $data[3];
        my $populationCaptureId = $data[4];
        my $pValue = $data[5];
        my $rawRank = $data[6];
        my $maxRank = $speciesMaxRanks->{$speciesId}->{$populationCaptureId}->{$isSingleCell}->{$isSampleMultiplexing};
        if (!defined $maxRank) {
            warn "Undefined maxRank for populationCaptureId=$populationCaptureId, isSingleCell=$isSingleCell, isSampleMultiplexing=$isSampleMultiplexing (geneId=$geneId, conditionId=$exprMappedConditionId) no score will be computed for this geneId/conditionId/datatype combination\n";
            next;
        }
        # For each geneId, accumulate weighted sums and counts
        if ($geneId != $previousGeneId) {
            $previousGeneId = $geneId;
            $results{$geneId}{'pValue'} = $pValue * $sampleDistinctRankCount;
            $results{$geneId}{'weight'} = $sampleDistinctRankCount;
            $results{$geneId}{'numberObs'} = 1;
            $results{$geneId}{'score'} = calculate_score($rawRank, $sampleMaxRank, $maxRank) * $sampleDistinctRankCount;
        } else {
            $results{$geneId}{'pValue'} += $pValue * $sampleDistinctRankCount;
            $results{$geneId}{'weight'} += $sampleDistinctRankCount;
            $results{$geneId}{'numberObs'} += 1;
            $results{$geneId}{'score'} += calculate_score($rawRank, $sampleMaxRank, $maxRank) * $sampleDistinctRankCount;
        }
    }
    $sth->finish;
    # Update the expression table with summarized metrics for each geneId/conditionId
    my $updateQuery = 'UPDATE expression SET ';
    my $datatypePrefix = $BULK_COLUMN_PREFIX;
    if ($isSingleCell) {
        if ($isSampleMultiplexing) {
            $datatypePrefix = $DROPLET_COLUMN_PREFIX;
        } else {
            $datatypePrefix = $FULL_LENGTH_COLUMN_PREFIX;
        }
    }
    $updateQuery .= $datatypePrefix."PValue = ?, ";
    $updateQuery .= $datatypePrefix."Score = ?, ";
    $updateQuery .= $datatypePrefix."Weight = ?, ";
    $updateQuery .= $datatypePrefix."NumberObs = ? ";
    $updateQuery .= 'WHERE conditionId = ? AND bgeeGeneId = ?';
    my $sthUpdate = $bgee_thread->prepare($updateQuery);
    for my $geneId (keys %results) {
        my $pValue = $results{$geneId}{'pValue'};
        my $score = $results{$geneId}{'score'};
        my $weight = $results{$geneId}{'weight'};
        my $numberObs = $results{$geneId}{'numberObs'};
        my $weightedPValue = $pValue / $weight;
        my $weightedScore = $score / $weight;
        if ($debug) {
            print "[DEBUG] update expression set $datatypePrefix"."PValue = $weightedPValue, $datatypePrefix"."Score = $weightedScore, $datatypePrefix"."Weight = $weight, $datatypePrefix"."NumberObs = $numberObs where conditionId = $exprMappedConditionId and bgeeGeneId = $geneId\n";
        } else {
            $sthUpdate->execute($weightedPValue, $weightedScore, $weight, $numberObs, $exprMappedConditionId, $geneId) or die $sthUpdate->errstr;
        }
    }
    $sthUpdate->finish;
}

# -----------------------------------------------------------------------------
# Subroutine: calculate_score
#   Normalizes the raw rank to a 0-100 scale based on sample and species max ranks.
#   Returns the normalized score, or undef if input is invalid.
# -----------------------------------------------------------------------------
sub calculate_score {
    my ($rawRank, $sampleMaxRank, $speciesMaxRank) = @_;
    if (defined $rawRank && $rawRank >= 1 && defined $sampleMaxRank && $sampleMaxRank >= 1 && defined $speciesMaxRank && $speciesMaxRank >= 1 && $rawRank <= $sampleMaxRank && $sampleMaxRank <= $speciesMaxRank) {
        my $normalizedRank = $rawRank / $sampleMaxRank * $speciesMaxRank;
        my $range = $speciesMaxRank - 1;
        my $adjustedRank = $normalizedRank - 1;
        return 100 - ($adjustedRank * 99 / $range);
    } else {
        warn "Invalid input for score calculation: rawRank=$rawRank, sampleMaxRank=$sampleMaxRank, speciesMaxRank=$speciesMaxRank\n";
        return undef;
    }
}
