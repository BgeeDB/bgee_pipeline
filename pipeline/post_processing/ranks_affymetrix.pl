#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw( time );
use diagnostics;

# Updates the ranks in affymetrixProbeset table and affymetrix mean rank in the expression table.
# Philippe Moret, created Oct 2015.
# Frederic Bastian, last updated June 2016.
# Frederic Bastian, last updated Feb. 2017: adapt to new conditions and new schema in Bgee 14

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Getopt::Long;

$|=1;

# Define arguments and their default value
my ($bgee_connector) = ('');
my ($ranks_computed) = (0);
my %opts = ('bgee=s'         => \$bgee_connector, # Bgee connector string
            'ranks_computed' => \$ranks_computed,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\t-ranks_computed   Skip generation of raw ranks per chip
\n";
    exit 1;
}

my $dbh = Utils::connect_bgee_db($bgee_connector);

#Set to 0 in order to disable autocommit to optimize speed
my $auto = 0;

if ($auto == 0) {
    $dbh->{AutoCommit} = 0;
}

# NOTE: this script assumes that we don't use in Bgee chip types allowing to hybridize samples
# from two different species at the same time.
#
# Reasonning of the computations:
# 1) compute fractional ranks of genes in table affymetrixProbeset, for each affymetrix chip,
# based on the highest signal intensity from all the probesets mapped to a given gene.
# 2) normalize ranks between chips in a same condition-species of the expression table (mapped condition).
# First, for each chip type, compute its max rank over all conditions. Then, normalize ranks of each chip,
# based on the max rank of this type of chip, as compared to the max of the max ranks of other chip types
# present in the same mapped condition. The idea is to correct for the different genomic coverage
# of different chip types. We do not normalize simply based on the max rank in a given condition,
# to not penalize conditions with a lower number of expressed genes, and thus higher number
# of ex-aequo genes, and lower fractional max ranks.
# 3) compute weighted mean of normalized ranks per gene and mapped condition,
# weighted by the number of distinct ranks in each chip: we assume that chips with higher
# number of distinct ranks have a higher power for ranking genes. The higher power
# at ranking genes of a chip is taken into account by weighting the mean
# by the number of distinct ranks in the chip, not by "normalizing" away chips with lower max ranks,
# again, for not penalizing conditions with lower number of expressed genes;
# 5) also, insert max of max ranks of chip types represented in each mapped condition, in condition table
# (will allow to normalize ranks between conditions and data types), and sum of numbers of distinct ranks
# per gene and condition, in expression table (used to compute weigthed mean over all data types in a condition)


##############################################
# COMPUTE RANKS PER CHIP                     #
##############################################
if ( !$ranks_computed ) {

#    # Clean potentially already computed ranks
#    my $cleanProbeset = $dbh->prepare("UPDATE affymetrixProbeset SET rank = null");
#    my $cleanChip     = $dbh->prepare("UPDATE affymetrixChip SET chipMaxRank = null,
#                                                                 chipDistinctRankCount = null");
#    my $cleanChipType = $dbh->prepare("UPDATE chipType SET chipTypeMaxRank = null");
#    printf("Cleaning existing data: ");
#    $cleanProbeset->execute() or die $cleanProbeset->errstr;
#    $cleanChip->execute() or die $cleanProbeset->errstr;
#    $cleanChipType->execute() or die $cleanProbeset->errstr;
#    if ($auto == 0) {
#        $dbh->commit() or die("Failed commit");
#    }
#    printf("Done\n");

    # Queries to compute gene ranks per chip. We rank over all not-excluded genes, not only over expressed genes
    my $affymChipStmt            = $dbh->prepare("SELECT t1.bgeeAffymetrixChipId FROM affymetrixChip AS t1
                                                 WHERE EXISTS (SELECT 1 FROM affymetrixProbeset AS t2
                                                    WHERE t2.expressionId IS NOT NULL
                                                    AND t2.bgeeAffymetrixChipId = t1.bgeeAffymetrixChipId) ");
    my $affymProbeSetStmt        = $dbh->prepare("SELECT bgeeGeneId, MAX(normalizedSignalIntensity) AS maxIntensity
                                                  FROM affymetrixProbeset
                                                  WHERE bgeeAffymetrixChipId = ?
                                                  AND reasonForExclusion NOT IN ('$Utils::EXCLUDED_FOR_PRE_FILTERED', '$Utils::EXCLUDED_FOR_UNDEFINED')
                                                  GROUP BY bgeeGeneId
                                                  ORDER BY maxIntensity DESC");
    # if several genes at a same rank, we'll update them at once with a 'bgeeGeneId IN (?,?, ...)' clause.
    # If only one gene at a given rank, updated with the prepared statement below.
    my $rankUpdateStart   = "UPDATE affymetrixProbeset SET rank = ?
                             WHERE bgeeAffymetrixChipId = ? AND bgeeGeneId ";
    my $affymeProbeSetUpdateStmt = $dbh->prepare($rankUpdateStart."= ?");



    $affymChipStmt->execute() or die $affymChipStmt->errstr;

    my @libs = map { $_->[0] } @{$affymChipStmt->fetchall_arrayref};
    my $l = @libs;

    printf("Found %d Affymetrix chips\n", $l);
    my $done = 0;
    my $start = time();
    my $timeref = $start;
    for my $k (0..$l-1) {
        my $affymetrixChipId = $libs[$k];
        printf("[%d/%d] Chips: %s\tIdx:%d", $done+1, $l , $affymetrixChipId, $k);
        #The rank of a gene is based on his best signal intensity
        $affymProbeSetStmt->execute($affymetrixChipId);

        my $mid = time();

        my @results = map { {"id" => $_->[0], "val" => $_->[1]} } @{$affymProbeSetStmt->fetchall_arrayref};

        my %sorted = Utils::fractionnal_ranking(@results);
        # we get ranks as keys, with reference to an array of gene IDs with that rank as value
        my %reverseHash = Utils::revhash(%sorted);

        for my $rank(keys(%reverseHash)) {
            my $geneIds_arrRef = $reverseHash{$rank};
            my @geneIds_arr = @$geneIds_arrRef;
            my $geneCount = scalar @geneIds_arr;
            if ($geneCount == 1) {
                my $geneId = $geneIds_arr[0];
                $affymeProbeSetUpdateStmt->execute($rank, $affymetrixChipId, $geneId) or die $affymeProbeSetUpdateStmt->errstr;
            } else {
                my $query = $rankUpdateStart."IN (";
                for (my $i = 0; $i < $geneCount; $i++) {
                    if ($i > 0) {
                        $query .= ", ";
                    }
                    $query .= "?";
                }
                $query .= ")";
                my $affyRankMultiUpdateStmt = $dbh->prepare($query);
                $affyRankMultiUpdateStmt->execute($rank, $affymetrixChipId, @geneIds_arr) or die $affyRankMultiUpdateStmt->errstr;
            }
        }

        $done += 1;

        if ($auto == 0) {
            $dbh->commit() or die("Failed commit");
        }


        my $end = time();
        my $rem = ($end-$start)/$done * ($l-$done);
        printf("\t%.2fs (%.2fs)\tElasped:%.2fs \tRemaining: %.2fs\n", $end-$timeref, $mid-$timeref,$end - $start, $rem);
        $timeref = $end;
    }

    my $end = time();
    printf("%.2fs To update ranks in affymetrix probeset\n", $end - $start);

    # ##############
    # Store max rank and number of distinct ranks per chip, and max rank per chip type
    my $sql =
    "UPDATE affymetrixChip AS t0
     INNER JOIN (
         SELECT t1.bgeeAffymetrixChipId, MAX(t1.rank) AS maxRank, COUNT(DISTINCT t1.rank) AS distinctRankCount
         FROM affymetrixProbeset AS t1
         WHERE t1.reasonForExclusion NOT IN ('$Utils::EXCLUDED_FOR_PRE_FILTERED', '$Utils::EXCLUDED_FOR_UNDEFINED')
         GROUP BY t1.bgeeAffymetrixChipId
     ) AS ranks ON t0.bgeeAffymetrixChipId = ranks.bgeeAffymetrixChipId
     SET t0.chipMaxRank = ranks.maxRank, t0.chipDistinctRankCount = ranks.distinctRankCount
     WHERE EXISTS (
         SELECT 1 FROM affymetrixProbeset AS t2
         WHERE t2.expressionId IS NOT NULL AND t2.bgeeAffymetrixChipId = t0.bgeeAffymetrixChipId
     )";

    my $t0 = time();
    printf("Inserting max ranks and distinct rank counts in affymetrixChip table: ");
    my $maxRankChipStmt = $dbh->prepare($sql);
    $maxRankChipStmt->execute() or die $maxRankChipStmt->errstr;
    printf("Done in %.2fs\n", (time() - $t0));


    # Store max ranks per chip type over all conditions,
    # for later normalization
    my $affyChipTypeMaxStmt = $dbh->prepare(
    "UPDATE chipType AS t1
     INNER JOIN (
         SELECT chipTypeId, MAX(chipMaxRank) AS maxRank
         FROM affymetrixChip
         GROUP BY chipTypeId
     ) AS chipMaxRanks ON t1.chipTypeId = chipMaxRanks.chipTypeId
     SET t1.chipTypeMaxRank = chipMaxRanks.maxRank");

    $t0 = time();
    printf("Inserting max ranks per chip type: ");
    $affyChipTypeMaxStmt->execute() or die $affyChipTypeMaxStmt->errstr;
    printf("Done in %.2fs\n", (time() - $t0));


    if ($auto == 0) {
        $dbh->commit() or die("Failed commit");
    }
}


####################################################################################
# NORMALIZE RANKS PER CONDITION AND COMPUTE WEIGHTED MEAN RANKS PER GENE-CONDITION #
####################################################################################

# Retrieve codition parameter combinations to compute Affymetrix ranks for.
my $condParamCombinationsArrRef = Utils::get_cond_param_combinations($dbh, $Utils::AFFY_DATA_TYPE);
print "All condition parameter combinations to compute for Affymetrix data:\n";
print "@$_\n" for @{$condParamCombinationsArrRef};

# connection will be opened/closed at each combination iteration to delete all temp tables
$dbh->disconnect();

for my $condParamCombArrRef ( @{$condParamCombinationsArrRef} ){
    my @condParamComb = @{$condParamCombArrRef};
    print "***** combination @condParamComb *****\n";

    $dbh = Utils::connect_bgee_db($bgee_connector);
    if ($auto == 0) {
        $dbh->{AutoCommit} = 0;
    }
    # Queries to first clean data
#    my $cleanExpr = $dbh->prepare("UPDATE globalExpression AS t1
#                                   INNER JOIN globalCond AS t2 ON t1.globalConditionId = t2.globalConditionId
#                                   SET affymetrixMeanRank              = null,
#                                       affymetrixMeanRankNorm          = null,
#                                       affymetrixDistinctRankSum       = null,
#                                       affymetrixGlobalMeanRank        = null,
#                                       affymetrixGlobalMeanRankNorm    = null,
#                                       affymetrixGlobalDistinctRankSum = null
#                                   WHERE ".Utils::get_cond_param_comb_sql_clause($condParamCombArrRef, "t2"));
#    my $cleanCond = $dbh->prepare("UPDATE globalCond SET affymetrixMaxRank = null,
#                                                     affymetrixGlobalMaxRank = null
#                                   WHERE ".Utils::get_cond_param_comb_sql_clause($condParamCombArrRef, "globalCond"));
#
#    printf("Cleaning existing data: ");
#    $cleanExpr->execute() or die $cleanExpr->errstr;
#    $cleanCond->execute() or die $cleanCond->errstr;
#    if ($auto == 0) {
#        $dbh->commit() or die("Failed commit");
#    }
#    printf("Done\n");


    # connection will be opened/closed at each combination iteration to delete all temp tables
    $dbh->disconnect();

    # we will compute two different rank information for each combination: one taking into account
    # all rank info mapped to a given condition ('self'), or all rank info mappd to a given condition
    # plus all its sub-conditions ('global').
    my $selfRanks = 0;
    my $globalRanks = 0;

    while (!$selfRanks || !$globalRanks) {
        if (!$selfRanks) {
            print "** Computation of self ranks\n"
        } else {
            print "** Computation of global ranks\n"
        }
        $dbh = Utils::connect_bgee_db($bgee_connector);
        if ($auto == 0) {
            $dbh->{AutoCommit} = 0;
        }

        # Store an association between each globalCondition and the chips considered in it
        my $sql =
        "CREATE TEMPORARY TABLE globalCondToChip (
             PRIMARY KEY(bgeeAffymetrixChipId, globalConditionId), INDEX(globalConditionId))
             SELECT DISTINCT t1.globalConditionId, t4.bgeeAffymetrixChipId ".
             # Retrieve the valid raw conditions mapped to each globalCondition
             "FROM globalCond AS t1
              INNER JOIN globalCondToCond AS t2 ON t1.globalConditionId = t2.globalConditionId ";
        if (!$selfRanks) {
            $sql .= "AND t2.conditionRelationOrigin = 'self' ";
        } else {
            $sql .= "AND t2.conditionRelationOrigin IN ('self', 'descendant') ";
        }
        $sql .= "INNER JOIN cond AS t3 ON t3.exprMappedConditionId = t2.conditionId ".
                # retrieve the chips present in this globalCondition
                "INNER JOIN affymetrixChip AS t4 ON t3.conditionId = t4.conditionId".
                # Use only globalConditions for the requested condition parameter combination
                " WHERE ".Utils::get_cond_param_comb_sql_clause($condParamCombArrRef, "t1");

        my $t0 = time();
        printf("Creating temp table mapping Affymetrix chips to globalConditions: ");
        my $chipToGlobalCondStmt = $dbh->prepare($sql);
        $chipToGlobalCondStmt->execute() or die $chipToGlobalCondStmt->errstr;
        printf("Done in %.2fs\n", (time() - $t0));


        # Store the max rank of the chip type with the highest max rank
        # present in each valid global condition (this is different than retrieving the max rank
        # of the *chips* present in the condition)
        $sql =
        "CREATE TEMPORARY TABLE maxForCond (PRIMARY KEY(globalConditionId))
             SELECT t1.globalConditionId, MAX(t3.chipTypeMaxRank) AS maxRank ".
             # Retrieve all chips present in each globalCondition
             "FROM globalCondToChip AS t1 ".
             # and the max rank from all chip types present in this condition
             "INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
              INNER JOIN chipType AS t3 ON t2.chipTypeId = t3.chipTypeId
              GROUP BY t1.globalConditionId;";

        $t0 = time();
        printf("Creating temp table for chip type max rank per condition: ");
        my $affyCondMaxStmt = $dbh->prepare($sql);
        $affyCondMaxStmt->execute() or die $affyCondMaxStmt->errstr;
        printf("OK in %.2fs\n", (time() - $t0));


        # 1) We normalized the ranks between chips in a same condition of the expression table (mapped condition),
        # based on the max rank that can get their respective chip type to, and the max of the max ranks
        # of chip types represented in the condition.
        # 2) We compute a weighted mean of the normalized gene ranks from all chips in a condition,
        # weighted by the number of distinct ranks in each chip. We also sum the numbers of distinct ranks
        # in chips where a gene is expressed, to be able to later compute a weighted mean per gene-condition
        # over all data types.

        # Now, the problem is that the ranks from individual chips are normalized differently
        # depending on the condition they are considered in. If we computed these normalized ranks
        # at once, this would result as having potentially a significant portion of the affymetrixProbeset
        # table duplicated, leading in using potentially too much disk space. For this reason,
        # we normalize and average ranks per condition
        # ###################
        # prepare queries to be run for each globalCondition
        # ###################

        # Normalize ranks between chip types for each valid globalCondition
        $sql = "CREATE TEMPORARY TABLE chipNormalizedRank
           SELECT STRAIGHT_JOIN
                  affymetrixProbeset.bgeeGeneId, affymetrixProbeset.bgeeAffymetrixChipId, ".
           # within data-type normalization:
           # (we use the best normalized rank of genes from each Affymetrix chip,
           # to not count genes multiple times because of multiple probesets)
           "      (MIN(affymetrixProbeset.rank) +
                        (MIN(affymetrixProbeset.rank) * (maxForCond.maxRank / chipType.chipTypeMaxRank)) )
                   /2 AS rank
           FROM globalCondToChip ".
           # get the max rank in the global condition
           "INNER JOIN maxForCond ON globalCondToChip.globalConditionId = maxForCond.globalConditionId ".
           # get the max rank of the chip type corresponding to each chip
           "INNER JOIN affymetrixChip ON globalCondToChip.bgeeAffymetrixChipId = affymetrixChip.bgeeAffymetrixChipId
            INNER JOIN chipType ON affymetrixChip.chipTypeId = chipType.chipTypeId ".
           # get for each gene and chip the min rank from the probesets
           "INNER JOIN affymetrixProbeset ON affymetrixChip.bgeeAffymetrixChipId = affymetrixProbeset.bgeeAffymetrixChipId ".
           # where and group by.
           "WHERE affymetrixProbeset.expressionId IS NOT NULL AND globalCondToChip.globalConditionId = ?
            GROUP BY affymetrixProbeset.bgeeGeneId, affymetrixProbeset.bgeeAffymetrixChipId";
        my $affyChipNormRankStmt = $dbh->prepare($sql);
        my $keyAffyChipNormRankStmt = $dbh->prepare("ALTER TABLE chipNormalizedRank
                                                     ADD PRIMARY KEY(bgeeGeneId, bgeeAffymetrixChipId)");
        my $dropAffyChipNormRank = $dbh->prepare("DROP TABLE chipNormalizedRank");

        # compute weighted mean normalized ranks, and sum of numbers of distinct ranks
        $sql = "CREATE TEMPORARY TABLE weightedMeanRank
                SELECT STRAIGHT_JOIN
                      chipNormalizedRank.bgeeGeneId,
                      SUM(chipNormalizedRank.rank * affymetrixChip.chipDistinctRankCount)
                          /SUM(affymetrixChip.chipDistinctRankCount) AS meanRank,
                      SUM(affymetrixChip.chipDistinctRankCount) AS distinctRankCountSum
                FROM chipNormalizedRank
                INNER JOIN affymetrixChip ON affymetrixChip.bgeeAffymetrixChipId = chipNormalizedRank.bgeeAffymetrixChipId
                GROUP BY chipNormalizedRank.bgeeGeneId";
        my $affyWeightedMeanStmt = $dbh->prepare($sql);
        my $dropAffyWeightedMeanStmt = $dbh->prepare("DROP TABLE weightedMeanRank");

        #update the expression table
        $sql = "UPDATE weightedMeanRank
                STRAIGHT_JOIN globalExpression ".
                # we build the quer this way in order to benefit from the clustered index
                # on (bgeeGeneId, globalConditionId) of the globalExpression table
               "ON globalExpression.bgeeGeneId = weightedMeanRank.bgeeGeneId
                    AND globalExpression.globalConditionId = ? ";
        if (!$selfRanks) {
            $sql .= "SET globalExpression.affymetrixMeanRank = weightedMeanRank.meanRank,
                globalExpression.affymetrixDistinctRankSum = weightedMeanRank.distinctRankCountSum";
        } else {
            $sql .= "SET globalExpression.affymetrixGlobalMeanRank = weightedMeanRank.meanRank,
                globalExpression.affymetrixGlobalDistinctRankSum = weightedMeanRank.distinctRankCountSum";
        }
        my $expressionUpdateMeanRank = $dbh->prepare($sql);


        # ###################
        # Run computations per globalCondition
        # ###################
        my $queryConditions = $dbh->prepare('SELECT DISTINCT globalConditionId FROM globalCondToChip');
        $queryConditions->execute()  or die $queryConditions->errstr;
        my @conditions = ();
        while ( my @data = $queryConditions->fetchrow_array ){
            push(@conditions, $data[0]);
        }
        my $conditionCount = scalar(@conditions);
        print "Computing ranks per global condition, $conditionCount conditions retrieved.\n";

        my $i = 0;
        for my $globalConditionId ( @conditions ){
            $i++;

            $t0 = time();
#            printf("Creating temp table for within-data-type normalized ranks per gene-chip-globalCond: ");
            $affyChipNormRankStmt->execute($globalConditionId) or die $affyChipNormRankStmt->errstr;
            $keyAffyChipNormRankStmt->execute() or die $keyAffyChipNormRankStmt->errstr;
#            printf("OK in %.2fs\n", (time() - $t0));

            $t0 = time();
#            printf("Creating temp table for weighted mean ranks: ");
            $affyWeightedMeanStmt->execute() or die $affyWeightedMeanStmt->errstr;
#            printf("OK in %.2fs\n", (time() - $t0));

            $t0 = time();
#            printf("Updating expression table with normalized mean ranks and sum of numbers of distinct ranks: ");
            $expressionUpdateMeanRank->execute($globalConditionId) or die $expressionUpdateMeanRank->errstr;
#            printf("OK in %.2fs\n", (time() - $t0));

            $dropAffyChipNormRank->execute() or die $dropAffyChipNormRank->errstr;
            $dropAffyWeightedMeanStmt->execute() or die $dropAffyWeightedMeanStmt->errstr;

            if (($i / 100 - int($i / 100)) == 0) {
                printf("$i conditions done.\n");
            }
        }

        # update the globalCond table
        $sql = "UPDATE globalCond ".
               # get the max ranks per condition
               "INNER JOIN maxForCond ON globalCond.globalConditionId = maxForCond.globalConditionId ";
        if (!$selfRanks) {
            $sql .= "SET globalCond.affymetrixMaxRank = maxForCond.maxRank";
        } else {
            $sql .= "SET globalCond.affymetrixGlobalMaxRank = maxForCond.maxRank";
        }
        $t0 = time();
        printf("Updating condition table with max ranks: ");
        my $updateExpressionMaxRank = $dbh->prepare($sql);
        $updateExpressionMaxRank->execute() or die $updateExpressionMaxRank->errstr;
        printf("OK in %.2fs\n", (time() - $t0));

        if ($auto == 0) {
            $dbh->commit() or die("Failed commit");
        }

        # closing the connection will destroy the temporary tables
        $dbh->disconnect();

        if (!$selfRanks) {
            $selfRanks = 1;
        } else {
            $globalRanks = 1;
        }
    }
}

exit 0;

