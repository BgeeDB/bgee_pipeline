#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw( time );
use diagnostics;

# Compute inSituScore, inSituPValue, inSituWeight, inSituNumberObs for the expression table.
#
# Logic (mirrors ranks_global.pl for in situ):
# * For each condition (batch), build a temp table mapping exprMappedConditionId -> inSituSpotId.
# * For each condition, create an inSituRanking temp table:
#     - Compute scoreSum per gene using the same scoring as ranks_global.pl:
#         present/high  = +1,  present/low  = +0.5
#         absent/low    = -0.5, absent/high  = -1
#     - Also accumulate, per gene, a per-spot evidence weight (2 for a HIGH_QUAL spot,
#       1 for a LOW_QUAL spot), used both to weight the average p-value and as the
#       stored inSituWeight. This mirrors how sampleDistinctRankCount weighs more
#       informative RNA-Seq libraries more heavily, and is specific to the gene
#       (unlike a condition-wide constant): a gene seen on more/better spots gets
#       more weight than one seen on a single low-quality spot.
#     - Dense-rank genes by scoreSum (descending) within the condition.
# * Normalize the dense rank to a 0-100 score:
#     rank 1 (highest scoreSum) -> 100,  rank maxRank -> ~1
#     formula: score = 100 - (rawRank - 1) * 99 / (maxRank - 1)    [maxRank > 1]
#              score = 100                                           [maxRank = 1]
# * inSituPValue    = weightedPValueSum / weightSum (weighted average pValue of spots,
#                    weighted the same way as inSituWeight, mirroring the RNA-Seq formula)
# * inSituWeight    = weightSum (per-gene evidence weight, see above)
# * inSituNumberObs = spotCount (total number of spots observed)
# * Update expression table with these values.
#
# TODO (weighting change, needs review before running on production data):
#   inSituWeight used to be maxRank, a value constant across all genes of a condition
#   (the number of distinct scoreSum groups in that condition), unrelated to how much
#   evidence backs any specific gene. It is now weightSum, a per-gene weight (2 points per
#   HIGH_QUAL spot, 1 point per LOW_QUAL spot observed for that gene), mirroring how
#   sampleDistinctRankCount weighs RNA-Seq libraries. inSituPValue is now this same
#   per-spot weight applied to the average (weightedPValueSum / weightSum), instead of the
#   previous unweighted average (pValueSum / spotCount).
#   Please verify: (1) the 2:1 HIGH_QUAL:LOW_QUAL ratio is the intended relative
#   confidence between spot qualities, (2) weighting inSituPValue itself (not just
#   inSituWeight) by spot quality is desired, and (3) this stays an arithmetic mean on
#   purpose (unlike bulk/single-cell RNA-Seq, in situ is not expected to suffer from
#   dropout, so no harmonic mean is used here for now).
#   Since this changes previously computed values, conditions already processed will need
#   to be reprocessed (e.g. via -processAll, or by resetting the relevant columns to NULL)
#   for the new formula to apply.

use Parallel::ForkManager;

use FindBin;
use lib "$FindBin::Bin/..";
use Utils;
use Getopt::Long;
use POSIX qw/ceil/;

$|=1;

my ($bgee_connector) = ('');
my ($parallel_jobs)  = (10);
my ($conds_per_job)  = (100);
my (@cond_ids)       = ();
my ($cond_offset)    = (0);
my ($cond_count)     = (0);
my ($process_all)    = (0);

my %opts = ('bgee=s'          => \$bgee_connector,
            'parallel_jobs=i' => \$parallel_jobs,
            'conds_per_job=i' => \$conds_per_job,
            'cond_ids=s'      => \@cond_ids,
            'cond_offset=i'   => \$cond_offset,
            'cond_count=i'    => \$cond_count,
            'processAll'      => \$process_all,
           );

my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD)
\t-bgee           Bgee connector string
\t-parallel_jobs  Number of parallel threads (default: 10)
\t-conds_per_job  Number of conditions per thread (default: 100)
\t-cond_ids       Comma-separated list of exprMappedConditionIds to process
\t-cond_offset    Offset for condition retrieval (requires -cond_count)
\t-cond_count     Row count for condition retrieval
\t-processAll     Process all conditions, not only those with inSituNumberObs = 0
\n";
    exit 1;
}
if ($cond_offset < 0 || $cond_count < 0) {
    die('cond_offset and cond_count cannot be negative');
}
if ($cond_offset > 0 && $cond_count == 0) {
    die('cond_count must be provided if cond_offset is provided');
}
if ($parallel_jobs <= 0 || $conds_per_job <= 0) {
    die("Invalid argument parallel_jobs/conds_per_job\n");
}

@cond_ids = split(/,/, join(',', @cond_ids));
if ($cond_offset > 0 && @cond_ids) {
    die("Not possible to provide condition IDs and offset parameters at the same time\n");
}

# Retrieve conditions to process
my @conditions = @cond_ids;
if (!@cond_ids) {
    my $dbh = Utils::connect_bgee_db($bgee_connector);

    my $condSql = 'SELECT DISTINCT c.exprMappedConditionId '.
                  'FROM inSituSpot AS s '.
                  'INNER JOIN cond AS c ON s.conditionId = c.conditionId '.
                  'INNER JOIN expression AS e '.
                  '    ON e.conditionId = c.exprMappedConditionId AND e.bgeeGeneId = s.bgeeGeneId '.
                  "WHERE s.reasonForExclusion = '".$Utils::CALL_NOT_EXCLUDED."' ".
                  'AND s.expressionId IS NOT NULL';
    if (!$process_all) {
        $condSql .= ' AND e.inSituNumberObs IS NULL';
    }
    if ($cond_count > 0) {
        $condSql .= ' ORDER BY c.exprMappedConditionId LIMIT '.$cond_offset.', '.$cond_count;
    }

    my $queryConditions = $dbh->prepare($condSql);
    $queryConditions->execute()  or die $queryConditions->errstr;
    @conditions = map { $_->[0] } @{$queryConditions->fetchall_arrayref};
    $dbh->disconnect();
}

my $conditionCount = @conditions;
if ($conditionCount == 0) {
    print "No conditions to process.\n";
    exit 0;
}
if ($conditionCount < $conds_per_job) {
    $conds_per_job = $conditionCount;
}
print "In situ score computations per condition, $conditionCount conditions retrieved...\n";

my $iterationCount = ceil($conditionCount/$conds_per_job);
my $parallel = $parallel_jobs;
if ($iterationCount < $parallel_jobs) {
    $parallel = $iterationCount;
}

if ($parallel == 1) {
    while ( my @next_conds = splice(@conditions, 0, $conds_per_job) ) {
        print("\nStart batch of ".scalar(@next_conds)." conditions, process ID $$...\n");
        compute_update_insitu_scores(\@next_conds);
        print("\nDone batch, process ID $$.\n");
    }
} else {
    print("In situ score computations with $parallel threads and $conds_per_job conditions per thread...\n");
    my $pm = new Parallel::ForkManager($parallel);
    while ( my @next_conds = splice(@conditions, 0, $conds_per_job) ) {
        my $pid = $pm->start and next;
        print("\nStart batch of ".scalar(@next_conds)." conditions, process ID $$...\n");
        compute_update_insitu_scores(\@next_conds);
        print("\nDone batch, process ID $$.\n");
        $pm->finish;
    }
    $pm->wait_all_children;
}
print("In situ score computations done\n");

exit 0;


sub compute_update_insitu_scores {
    my ($batchRef) = @_;
    my $batchLength = scalar @{ $batchRef };

    # Each thread gets its own DB connection
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
    my $dbh_thread = Utils::connect_bgee_db($bgee_connector);

    #****************************************
    # PREPARE QUERIES
    #****************************************

    # Temp table mapping exprMappedConditionId -> inSituSpotId for the whole batch.
    # Mirrors inSituToGlobalCond in ranks_global.pl, but works on raw conditions:
    # we join cond directly on exprMappedConditionId instead of going through globalCondToCond.
    my $inSituToExprCondTableName = 'inSituToExprCond';
    my $inSituToExprCondSql =
        'CREATE TEMPORARY TABLE '.$inSituToExprCondTableName.' '.
        '(PRIMARY KEY(exprMappedConditionId, inSituSpotId)) '.
        'SELECT DISTINCT c.exprMappedConditionId, s.inSituSpotId '.
        'FROM cond AS c '.
        'INNER JOIN inSituSpot AS s ON c.conditionId = s.conditionId '.
        "WHERE s.reasonForExclusion = '".$Utils::CALL_NOT_EXCLUDED."' ".
        'AND s.expressionId IS NOT NULL '.
        'AND c.exprMappedConditionId';
    if ($batchLength == 1) {
        $inSituToExprCondSql .= ' = ?';
    } else {
        $inSituToExprCondSql .= ' IN (';
        $inSituToExprCondSql .= join(', ', ('?') x $batchLength);
        $inSituToExprCondSql .= ')';
    }
    my $inSituToExprCondStmt = $dbh_thread->prepare($inSituToExprCondSql);

    my $isInSituDataInCondStmt = $dbh_thread->prepare(
        'SELECT 1 FROM '.$inSituToExprCondTableName.' WHERE exprMappedConditionId = ? LIMIT 1');

    # Create inSituRanking temp table for one condition.
    # Mirrors inSituRankingStmt in ranks_global.pl, with the addition of weightedPValueSum,
    # weightSum and spotCount (used to compute inSituPValue, inSituWeight and inSituNumberObs
    # after ranking).
    my $inSituRankingStmt = $dbh_thread->prepare(
        "CREATE TEMPORARY TABLE inSituRanking (PRIMARY KEY(bgeeGeneId))
         SELECT STRAIGHT_JOIN s.bgeeGeneId,
             SUM(
                 IF(s.detectionFlag = '$Utils::PRESENT_CALL' AND s.inSituData = '$Utils::HIGH_QUAL',  1,
                 IF(s.detectionFlag = '$Utils::PRESENT_CALL',                                         0.5,
                 IF(s.detectionFlag = '$Utils::ABSENT_CALL'  AND s.inSituData = '$Utils::HIGH_QUAL', -1,
                 IF(s.detectionFlag = '$Utils::ABSENT_CALL',                                         -0.5, 0))))
             ) AS scoreSum,
             -- Per-spot evidence weight: 2 for a HIGH_QUAL spot, 1 for a LOW_QUAL spot.
             -- Kept as integers (rather than 1/0.5) so the sum fits the 'bigint unsigned'
             -- inSituWeight column without a schema migration; the 2:1 ratio is
             -- mathematically equivalent to 1:0.5 for a weighted average or a weighted
             -- harmonic mean (a uniform scale factor on all weights cancels out).
             SUM(s.pValue * IF(s.inSituData = '$Utils::HIGH_QUAL', 2, 1)) AS weightedPValueSum,
             SUM(IF(s.inSituData = '$Utils::HIGH_QUAL', 2, 1))            AS weightSum,
             COUNT(s.inSituSpotId) AS spotCount,
             0000000.00           AS rawRank
         FROM $inSituToExprCondTableName
         INNER JOIN inSituSpot AS s ON $inSituToExprCondTableName.inSituSpotId = s.inSituSpotId
         WHERE $inSituToExprCondTableName.exprMappedConditionId = ?
           AND s.expressionId IS NOT NULL
         GROUP BY s.bgeeGeneId");
    my $dropInSituRankingStmt = $dbh_thread->prepare('DROP TABLE inSituRanking');

    # Dense ranking (same approach as ranks_global.pl)
    my $inSituRanksForCondStmt = $dbh_thread->prepare(
        'SELECT bgeeGeneId, scoreSum FROM inSituRanking ORDER BY scoreSum DESC');
    my $inSituRankUpdateStart   = 'UPDATE inSituRanking SET rawRank = ? WHERE bgeeGeneId ';
    my $inSituUpdateRankingStmt = $dbh_thread->prepare($inSituRankUpdateStart.'= ?');

    # Retrieve max rank (for score normalization) and per-gene data after ranking
    my $getMaxRankStmt = $dbh_thread->prepare('SELECT MAX(rawRank) FROM inSituRanking');
    my $getResultsStmt = $dbh_thread->prepare(
        'SELECT bgeeGeneId, rawRank, weightedPValueSum, weightSum, spotCount FROM inSituRanking');

    # Update the expression table
    my $updateExprStmt = $dbh_thread->prepare(
        'UPDATE expression '.
        'SET inSituScore = ?, inSituPValue = ?, inSituWeight = ?, inSituNumberObs = ? '.
        'WHERE conditionId = ? AND bgeeGeneId = ?');

    #****************************************
    # EXECUTE BATCH QUERY FOR THE WHOLE CONDITION BATCH
    #****************************************
    $inSituToExprCondStmt->execute(@{ $batchRef }) or die $inSituToExprCondStmt->errstr;

    #****************************************
    # COMPUTE SCORES PER CONDITION
    #****************************************
    for my $k ( 0..$batchLength-1 ) {
        Utils::start_transaction($dbh_thread);
        my $condId = ${$batchRef}[$k];

        $isInSituDataInCondStmt->execute($condId) or die $isInSituDataInCondStmt->errstr;
        if ($isInSituDataInCondStmt->fetchrow_arrayref()) {

            # Build the ranking table for this condition
            $inSituRankingStmt->execute($condId) or die $inSituRankingStmt->errstr;

            # Dense-rank genes by scoreSum (same subroutine as ranks_global.pl)
            $inSituRanksForCondStmt->execute() or die $inSituRanksForCondStmt->errstr;
            my @inSituResults = map { {"id" => $_->[0], "val" => $_->[1]} }
                    @{$inSituRanksForCondStmt->fetchall_arrayref};
            compute_dense_ranking_update_tmp_ranks($dbh_thread, \@inSituResults,
                    $inSituUpdateRankingStmt, $inSituRankUpdateStart);

            # Get maxRank for 0-100 normalization
            $getMaxRankStmt->execute() or die $getMaxRankStmt->errstr;
            my ($maxRank) = $getMaxRankStmt->fetchrow_array();

            # Compute score/pValue/weight/numberObs for each gene and update expression table.
            # weight = weightSum, the per-gene sum of per-spot quality weights (2 for HIGH_QUAL,
            # 1 for LOW_QUAL), analogous to sampleDistinctRankCount in RNA-Seq: unlike the
            # previous maxRank-based weight (constant for every gene in the condition), this
            # scales with how many spots observed this gene and how reliable those spots are.
            # pValue is the weighted average of per-spot pValues, using the same weights.
            # numberObs = spotCount (total number of spots observed).
            $getResultsStmt->execute() or die $getResultsStmt->errstr;
            while ( my @row = $getResultsStmt->fetchrow_array() ) {
                my ($geneId, $rawRank, $weightedPValueSum, $weightSum, $spotCount) = @row;
                my $score     = calculate_score($rawRank, $maxRank);
                my $pValue    = $weightedPValueSum / $weightSum;
                my $weight    = $weightSum;
                my $numberObs = $spotCount;
                $updateExprStmt->execute($score, $pValue, $weight, $numberObs, $condId, $geneId)
                    or die $updateExprStmt->errstr;
            }

            $dropInSituRankingStmt->execute() or die $dropInSituRankingStmt->errstr;
        }

        $dbh_thread->commit() or die('Failed commit');
        printf("conditionId: %s - PID: %s - %d/%d\n", $condId, $$, $k+1, $batchLength);
    }

    $inSituToExprCondStmt->finish    or die('Failed finish');
    $isInSituDataInCondStmt->finish  or die('Failed finish');
    $inSituRankingStmt->finish       or die('Failed finish');
    $dropInSituRankingStmt->finish   or die('Failed finish');
    $inSituRanksForCondStmt->finish  or die('Failed finish');
    $inSituUpdateRankingStmt->finish or die('Failed finish');
    $getMaxRankStmt->finish          or die('Failed finish');
    $getResultsStmt->finish          or die('Failed finish');
    $updateExprStmt->finish          or die('Failed finish');

    $dbh_thread->disconnect();
}


# Normalize dense rank to a 0-100 score (same direction as RNA-Seq):
#   rank 1        -> 100.00  (highest scoreSum = most expressed)
#   rank maxRank  ->   1.00  (lowest scoreSum  = least expressed)
#   maxRank = 1   -> 100.00  (single gene in condition)
sub calculate_score {
    my ($rawRank, $maxRank) = @_;
    unless (defined $rawRank && defined $maxRank && $maxRank >= 1 && $rawRank >= 1) {
        warn "Invalid input for score calculation: rawRank=$rawRank, maxRank=$maxRank\n";
        return undef;
    }
    return 100.00 if $maxRank == 1;
    return sprintf("%.2f", 100.0 - ($rawRank - 1) * 99.0 / ($maxRank - 1));
}


# Identical to compute_dense_ranking_update_tmp_ranks in ranks_global.pl
sub compute_dense_ranking_update_tmp_ranks {
    my ($dbh, $resultArrRef, $updateRankingStmt, $updateRankingStmtStart) = @_;

    my %sorted = Utils::dense_ranking(@{ $resultArrRef });
    my %reverseHash = Utils::revhash(%sorted);

    for my $rank ( keys %reverseHash ){
        my $geneIds_arrRef = $reverseHash{$rank};
        my @geneIds_arr    = @$geneIds_arrRef;
        my $exprCount      = scalar @geneIds_arr;
        if ( $exprCount == 1 ){
            $updateRankingStmt->execute($rank, $geneIds_arr[0]) or die $updateRankingStmt->errstr;
        } else {
            my $query = $updateRankingStmtStart.'IN (';
            for ( my $i = 0; $i < $exprCount; $i++ ){
                $query .= ', ' if $i > 0;
                $query .= '?';
            }
            $query .= ')';
            my $rankMultiUpdateStmt = $dbh->prepare($query);
            $rankMultiUpdateStmt->execute($rank, @geneIds_arr) or die $rankMultiUpdateStmt->errstr;
        }
    }
}
