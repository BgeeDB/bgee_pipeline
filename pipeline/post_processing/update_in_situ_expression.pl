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
# * Normalize the dense rank between datatypes, then map it to a 1-100 score:
#     - condMaxRank    = max dense rank in the condition (MAX(rawRank) of inSituRanking)
#     - speciesMaxRank = max rank of the species over all datatypes. In practice the max
#                        of rnaSeqPopulationCaptureSpeciesMaxRank for that species, the
#                        RNA-Seq max ranks being far above the in situ ones.
#     - rankNorm = (rawRank + rawRank * speciesMaxRank / condMaxRank) / 2
#                  (same formula as normalize_ranks.pl and ranks_global.pl)
#     - score    = 100 - (rankNorm - 1) * 99 / (speciesMaxRank - 1)
#   The score is no longer stretched over the whole 1-100 range inside each condition:
#   a condition with few distinct scoreSum groups has a low ranking resolution, so its
#   genes stay in the upper part of the scale (the worst-ranked gene of a condition tends
#   to ~50 as soon as condMaxRank << speciesMaxRank, and rank 1 tends to 100 as the
#   condition gets finer). This is what makes in situ scores comparable between
#   conditions, and comparable with the RNA-Seq scores.
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
#   to be reprocessed (e.g. via -processAll) for the new formula to apply.

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
my ($process_all)    = (0);

my %opts = ('bgee=s'          => \$bgee_connector,
            'parallel_jobs=i' => \$parallel_jobs,
            'conds_per_job=i' => \$conds_per_job,
            'processAll'      => \$process_all,
           );

my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD)
\t-bgee           Bgee connector string
\t-parallel_jobs  Number of parallel threads (default: 10)
\t-conds_per_job  Number of conditions per thread (default: 100)
\t-processAll     Reset inSituScore, inSituPValue, inSituWeight and inSituNumberObs to NULL before
\t                processing, so that all conditions are reprocessed. Without it, only conditions
\t                with missing scores are processed (relaunch without -processAll to resume an
\t                interrupted run)
\n";
    exit 1;
}
if ($parallel_jobs <= 0 || $conds_per_job <= 0) {
    die("Invalid argument parallel_jobs/conds_per_job\n");
}

# -----------------------------------------------------------------------------
# With -processAll, the in situ columns are first reset to NULL. The conditions to
# process are then always retrieved with the same "IS NULL" filter, whether
# -processAll was provided or not. It means that if the run is interrupted (e.g.
# some queries end with an out of time error because too many threads overloaded
# the database), it can be relaunched without -processAll and will continue from
# where it stopped.
# -----------------------------------------------------------------------------
if ($process_all) {
    my $dbh = Utils::connect_bgee_db($bgee_connector);
    my $resetSql = 'UPDATE expression SET inSituScore = NULL, inSituPValue = NULL, '.
                   'inSituWeight = NULL, inSituNumberObs = NULL '.
                   'WHERE inSituScore IS NOT NULL OR inSituPValue IS NOT NULL '.
                   'OR inSituWeight IS NOT NULL OR inSituNumberObs IS NOT NULL';
    print "Resetting in situ columns to NULL...\n";
    my $resetRows = $dbh->do($resetSql)  or die $dbh->errstr;
    print "Done, $resetRows rows reset.\n";
    $dbh->disconnect();
}

# Retrieve conditions to process
my @conditions = ();
# speciesId of each condition to process, needed to normalize its ranks
my %conditionSpecies = ();
my $dbh = Utils::connect_bgee_db($bgee_connector);
my $condSql = 'SELECT DISTINCT c.exprMappedConditionId, c.speciesId '.
              'FROM inSituSpot AS s '.
              'INNER JOIN cond AS c ON s.conditionId = c.conditionId '.
              'INNER JOIN expression AS e '.
              '    ON e.conditionId = c.exprMappedConditionId AND e.bgeeGeneId = s.bgeeGeneId '.
              "WHERE s.reasonForExclusion = '".$Utils::CALL_NOT_EXCLUDED."' ".
              'AND s.expressionId IS NOT NULL '.
              'AND e.inSituNumberObs IS NULL';
my $queryConditions = $dbh->prepare($condSql);
$queryConditions->execute()  or die $queryConditions->errstr;
for my $row ( @{$queryConditions->fetchall_arrayref} ){
    push @conditions, $row->[0];
    $conditionSpecies{$row->[0]} = $row->[1];
}

# Max rank of each species over all datatypes, used as the common scale onto which the
# ranks of each condition are normalized. Only RNA-Seq stores per species max ranks, but
# they are far above the in situ ones (dense ranks over scoreSum), so the max over all
# RNA-Seq protocols of a species is its max rank over all datatypes.
my %speciesMaxRanks = ();
my $queryMaxRanks = $dbh->prepare(
    'SELECT speciesId, MAX(maxRank) FROM rnaSeqPopulationCaptureSpeciesMaxRank GROUP BY speciesId');
$queryMaxRanks->execute()  or die $queryMaxRanks->errstr;
while ( my @data = $queryMaxRanks->fetchrow_array ){
    $speciesMaxRanks{$data[0]} = $data[1];
}
$queryMaxRanks->finish;

# Check all species up front, to fail before starting rather than in the middle of a run
my %speciesToProcess = map { $_ => 1 } values %conditionSpecies;
for my $speciesId ( sort keys %speciesToProcess ){
    if ( !defined $speciesMaxRanks{$speciesId} || $speciesMaxRanks{$speciesId} <= 1 ){
        die "No usable max rank for species $speciesId in rnaSeqPopulationCaptureSpeciesMaxRank, ".
            "in situ ranks can not be normalized\n";
    }
}
$dbh->disconnect();


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

            # Max rank of the condition and max rank of the species over all datatypes,
            # both needed to normalize the ranks between datatypes
            $getMaxRankStmt->execute() or die $getMaxRankStmt->errstr;
            my ($condMaxRank) = $getMaxRankStmt->fetchrow_array();
            my $speciesMaxRank = $speciesMaxRanks{ $conditionSpecies{$condId} };

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
                my $score     = calculate_score($rawRank, $condMaxRank, $speciesMaxRank);
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


# Normalize the dense rank of a condition between datatypes, then map it to a 1-100 score
# (same direction and same formula as the RNA-Seq scores):
#   rankNorm = (rawRank + rawRank * speciesMaxRank / condMaxRank) / 2
#   score    = 100 - (rankNorm - 1) * 99 / (speciesMaxRank - 1)
# As rawRank <= condMaxRank <= speciesMaxRank, rankNorm stays in [1, speciesMaxRank] and
# the score in [1, 100]. Only a condition ranking as many genes as the species max rank
# can reach the bottom of the scale: the lowest rank of a coarser condition tends to ~50,
# which is the point of the normalization (a gene ranked last among a handful of distinct
# scoreSum groups is far less informative than one ranked last among thousands).
# A condition with a single rank group (condMaxRank = 1) therefore scores ~50.5, and no
# longer 100 as when each condition was normalized on its own max rank.
# Dies if input is invalid.
sub calculate_score {
    my ($rawRank, $condMaxRank, $speciesMaxRank) = @_;
    unless (defined $rawRank && $rawRank >= 1 && defined $condMaxRank && $condMaxRank >= 1
            && defined $speciesMaxRank && $speciesMaxRank > 1
            && $rawRank <= $condMaxRank && $condMaxRank <= $speciesMaxRank) {
        die "Invalid input for score calculation: rawRank=$rawRank, condMaxRank=$condMaxRank, ".
            "speciesMaxRank=$speciesMaxRank\n";
    }
    my $rankNorm = ($rawRank + ($rawRank * $speciesMaxRank / $condMaxRank)) / 2;
    return sprintf("%.2f", 100.0 - ($rankNorm - 1) * 99.0 / ($speciesMaxRank - 1));
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
