#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw( time );
use diagnostics;

# Philippe Moret, created Oct 2015.
# Frederic Bastian, updated June 2016: Updates the In Situ mean ranks in the expression table
# Frederic Bastian, last updated Feb. 2017: adapt to new conditions and new schema in Bgee 14

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Getopt::Long;

$|=1;

# Define arguments and their default value
my ($bgee_connector) = ('');
my %opts = ('bgee=s'         => \$bgee_connector, # Bgee connector string
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\n";
    exit 1;
}

#Set to 0 in order to disable autocommit to optimize speed
my $auto = 0;

# =============
# Reasonning of the computations:
# 1) First, compute a score for each gene-condition, based on the detectionFlag
# of in situ evidence (present, absent) and the quality level (high quality, poor quality).
# Stored in temp table inSituRanking.
# 2) "dense ranking" of the genes in each condition, based on the score computed.
# Ranks stored in the temp table inSituRanking.
# All experiments are all considered together (no ranking first performed per experiment,
# as for other data types). This is because each in situ experiment usually study
# a limited number of genes.
# 3) Insert ranks into expression table
# 4) for each condition, retrieve max rank
# 5) insert max rank per condition into condition table (will allow to compute weigthed mean
# between data types in a condition, and to later normalize ranks between conditions and data types)
# =============


##############################################
# COMPUTE RANKS PER CONDITION                #
##############################################
my $dbh = Utils::connect_bgee_db($bgee_connector);
if ($auto == 0) {
    $dbh->{AutoCommit} = 0;
}
#    # Queries to first clean data
#    my $cleanExpr = $dbh->prepare("UPDATE globalExpression SET inSituRank = NULL,
#                                                               inSituRankNorm = NULL,
#                                                               inSituGlobalRank = NULL,
#                                                               inSituGlobalRankNorm = NULL");
#    my $cleanCond = $dbh->prepare("UPDATE globalCond SET inSituMaxRank = NULL,
#                                                         inSituGlobalMaxRank = NULL");
#
#    printf("Cleaning existing data: ");
#    $cleanExpr->execute() or die $cleanExpr->errstr;
#    $cleanCond->execute() or die $cleanCond->errstr;
#    if ($auto == 0) {
#        $dbh->commit() or die("Failed commit");
#    }
#    printf("Done\n");

# Retrieve codition parameter combinations to compute in situ ranks for.
my $condParamCombinationsArrRef = Utils::get_cond_param_combinations($dbh, $Utils::IN_SITU_DATA_TYPE);
print "All condition parameter combinations to compute for in situ data:\n";
print "@$_\n" for @{$condParamCombinationsArrRef};

# connection will be opened/closed at each combination iteration to delete all temp tables
$dbh->disconnect();

for my $condParamCombArrRef ( @{$condParamCombinationsArrRef} ){
    my @condParamComb = @{$condParamCombArrRef};
    print "***** combination @condParamComb *****\n";


    # connection will be opened/closed at each combination iteration to delete all temp tables

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

        # Store an association between each globalCondition and the inSituSpots considered in it
        my $sql =
        "CREATE TEMPORARY TABLE globalCondToSpot (
             PRIMARY KEY(inSituSpotId, globalConditionId), INDEX(globalConditionId))
             SELECT DISTINCT t1.globalConditionId, t4.inSituSpotId ".
             # Retrieve the valid raw conditions mapped to each globalCondition
             "FROM globalCond AS t1
              INNER JOIN globalCondToCond AS t2 ON t1.globalConditionId = t2.globalConditionId ";
        if (!$selfRanks) {
            $sql .= "AND t2.conditionRelationOrigin = 'self' ";
        } else {
            $sql .= "AND t2.conditionRelationOrigin IN ('self', 'descendant') ";
        }
        $sql .= "INNER JOIN cond AS t3 ON t3.exprMappedConditionId = t2.conditionId ".
                # retrieve the spots present in this globalCondition
                "INNER JOIN inSituSpot AS t4 ON t3.conditionId = t4.conditionId".
                # Use only globalConditions for the requested condition parameter combination
                " WHERE ".Utils::get_cond_param_comb_sql_clause($condParamCombArrRef, "t1");

        my $t0 = time();
        printf("Creating temp table mapping in situ spots to globalConditions: ");
        my $spotToGlobalCondStmt = $dbh->prepare($sql);
        $spotToGlobalCondStmt->execute() or die $spotToGlobalCondStmt->errstr;
        printf("Done in %.2fs\n", (time() - $t0));

        # give a score to spots, depending on their detectionFlag and quality:
        # present - high quality = 1
        # present - low quality = 0.5
        # absent - low quality = -0.5
        # absent - high quality = -1
        # We sum these sores for each gene in a given mapped condition of the expression table,
        # and we will rank genes in each mapped condition based on this score.
        $sql = "CREATE TEMPORARY TABLE inSituRanking (
            PRIMARY KEY(bgeeGeneId), INDEX(rank))
            SELECT STRAIGHT_JOIN inSituSpot.bgeeGeneId,

            SUM(
                IF(inSituSpot.detectionFlag = '$Utils::PRESENT_CALL' AND inSituData = '$Utils::HIGH_QUAL', 1,
                    IF(inSituSpot.detectionFlag = '$Utils::PRESENT_CALL', 0.5,
                        IF(inSituSpot.detectionFlag = '$Utils::ABSENT_CALL' AND inSituData = '$Utils::HIGH_QUAL', -1,
                            IF(inSituSpot.detectionFlag = '$Utils::ABSENT_CALL', -0.5, 0))))
            ) AS scoreSum,
            0 AS rank

            FROM globalCondToSpot
            INNER JOIN inSituSpot ON globalCondToSpot.inSituSpotId = inSituSpot.inSituSpotId
            WHERE globalCondToSpot.globalConditionId = ?
            AND inSituSpot.detectionFlag != '$Utils::UNDEFINED_CALL' AND inSituSpot.inSituData != 'no data'
            AND inSituSpot.reasonForExclusion NOT IN ('$Utils::EXCLUDED_FOR_PRE_FILTERED', '$Utils::EXCLUDED_FOR_UNDEFINED')
            GROUP BY inSituSpot.bgeeGeneId;";

        my $rankingStmt = $dbh->prepare($sql);


        # Queries to compute gene ranks per condition
        my $ranksForCondStmt = $dbh->prepare("
            SELECT bgeeGeneId, scoreSum FROM inSituRanking
            ORDER BY scoreSum DESC;
        ");
        # if several genes at a same rank, we'll update them at once with a 'geneId IN (?,?, ...)' clause.
        # If only one gene at a given rank, updated with the prepared statement below.
        my $rankUpdateStart   = "UPDATE inSituRanking SET rank = ?
                                 WHERE bgeeGeneId ";
        my $updateRankingStmt = $dbh->prepare($rankUpdateStart."= ?");
        my $dropSpotRankStmt = $dbh->prepare("DROP TABLE inSituRanking");

        #update the expression table
        $sql = "UPDATE inSituRanking
                STRAIGHT_JOIN globalExpression ".
                # we build the query this way in order to benefit from the clustered index
                # on (bgeeGeneId, globalConditionId) of the globalExpression table
               "ON globalExpression.bgeeGeneId = inSituRanking.bgeeGeneId
                    AND globalExpression.globalConditionId = ? ";
        if (!$selfRanks) {
            $sql .= "SET globalExpression.inSituRank = inSituRanking.rank";
        } else {
            $sql .= "SET globalExpression.inSituGlobalRank = inSituRanking.rank";
        }
        my $expressionUpdateMeanRank = $dbh->prepare($sql);


        $sql = "UPDATE globalCond ";
        if (!$selfRanks) {
            $sql .= "SET globalCond.inSituMaxRank ";
        } else {
            $sql .= "SET globalCond.inSituGlobalMaxRank ";
        }
        $sql .= "= (SELECT MAX(rank) FROM inSituRanking)
                 WHERE globalConditionId = ?";
        my $updateExpressionMaxRank = $dbh->prepare($sql);


        # ###################
        # Run computations per globalCondition
        # ###################
        my $queryConditions = $dbh->prepare('SELECT DISTINCT globalConditionId FROM globalCondToSpot');
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
#            printf("Creating temp table for in situ ranking: ");
            $rankingStmt->execute($globalConditionId) or die $rankingStmt->errstr;
#            printf("OK in %.2fs\n", (time() - $t0));

            $ranksForCondStmt->execute() or die $ranksForCondStmt->errstr;

            my @results = map { {"id" => $_->[0], "val" => $_->[1]} } @{$ranksForCondStmt->fetchall_arrayref};

            my %sorted = Utils::dense_ranking(@results);
            # we get ranks as keys, with reference to an array of expression IDs with that rank as value
            my %reverseHash = Utils::revhash(%sorted);

            for my $rank(keys(%reverseHash)) {
                my $geneIds_arrRef = $reverseHash{$rank};
                my @geneIds_arr = @$geneIds_arrRef;
                my $exprCount = scalar @geneIds_arr;
                if ($exprCount == 1) {
                    my $geneId = $geneIds_arr[0];
                    $updateRankingStmt->execute($rank,  $geneId) or die $updateRankingStmt->errstr;
                } else {
                    my $query = $rankUpdateStart."IN (";
                    for (my $i = 0; $i < $exprCount; $i++) {
                        if ($i > 0) {
                            $query .= ", ";
                        }
                        $query .= "?";
                    }
                    $query .= ")";
                    my $rankMultiUpdateStmt = $dbh->prepare($query);
                    $rankMultiUpdateStmt->execute($rank,  @geneIds_arr) or die $rankMultiUpdateStmt->errstr;
                }
            }

            $t0 = time();
#            printf("Updating expression table with in situ ranks: ");
            $expressionUpdateMeanRank->execute($globalConditionId) or die $expressionUpdateMeanRank->errstr;
#            printf("OK in %.2fs\n", (time() - $t0));

            $t0 = time();
#            printf("Updating globalCond table with in situ max rank: ");
            $updateExpressionMaxRank->execute($globalConditionId) or die $updateExpressionMaxRank->errstr;
#            printf("OK in %.2fs\n", (time() - $t0));

            $dropSpotRankStmt->execute() or die $dropSpotRankStmt->errstr;

            if (($i / 100 - int($i / 100)) == 0) {
                printf("$i conditions done.\n");
            }
        }

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
