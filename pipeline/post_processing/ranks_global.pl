#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw( time );
use diagnostics;

# Compute the ranks for globalExpression calls in each global conditions
# Frederic Bastian, Apr. 2021

use Parallel::ForkManager;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Getopt::Long;
use POSIX qw/ceil/;

$|=1;

# Define arguments and their default value
my ($bgee_connector) = ('');
my ($parallel_jobs)  = (15); # default 15 parallel threads used to compute ranks
my ($conds_per_job)  = (100); # default 100 conditions per thread
my (@cond_ids)       = ();
my ($cond_offset)    = (0);
my ($cond_count)     = (0);
my %opts = ('bgee=s'          => \$bgee_connector, # Bgee connector string
            'parallel_jobs=i' => \$parallel_jobs,
            'conds_per_job=i' => \$conds_per_job,
            'cond_ids=s'      => \@cond_ids,
            'cond_offset=i'   => \$cond_offset,
            'cond_count=i'    => \$cond_count,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD)
\t-bgee           Bgee connector string
\t-parallel_jobs  Number of threads used to compute ranks
\t-conds_per_job  Number of conditions per thread
\t-cond_ids       a comma-separated list of global condition IDs to treat, instead of providing cond_offset and cond_count
\t-cond_offset    The offset parameter to retrieve conditions to compute ranks for
\t-cond_count     The row_count parameter to retrieve conditions to compute ranks for
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

# Reasoning of the computations:
# First, ranks for each RNA-Seq library and Affymetrix chips are computed
# (see related scripts in this folder).
# Then:
# For RNA-Seq data:
# * update global expression table: compute weighted mean of ranks per gene and condition,
# weighted by the number of distinct ranks in each library: we assume that libraries with higher
# number of distinct ranks have a higher power for ranking genes.
# For RNA-Seq data, we used to apply the same max rank to all conditions of a same species,
# because we considered that we always had access to all expression information
# for the same set of genes in all libraries. This is not true anymore (and it actually never was),
# as we do now distinguish between different RNA-Seq protocols targeting different biotypes.
#
# So we have a similar approach to Affymetrix data, where we store
# for each condition th max rank among the *chip types* present in a condition
# (and not the max rank among the actual chips present in a condition): we store
# the max rank among the *protocols* present in a condition. We will then do
# a first "between-protocol" normalization, as we do a first "between-chip-type" normalization
# in Affymetrix data.
# The higher power at ranking genes of a library (for instance, thanks to a higher number of mapped reads)
# is taken into account by weighting the mean by the number of distinct ranks in the library,
# not by "normalizing" away libraries with lower max ranks; this would penalize conditions
# with a lower number of expressed genes, and thus with more ex-aequo ranked genes, corresponding
# to genes receiving 0 read.
# * also, insert max ranks per mapped condition in condition table
# (will allow to normalize ranks between conditions and data types), and sum of distinct rank counts
# per gene and condition, in expression table (used to compute weigthed mean over all data types in a condition)
#
# For Affymetrix data:
# * normalize ranks between chips in a same condition-species of the global expression table (global condition).
# First, for each chip type, compute its max rank over all conditions. Then, normalize ranks of each chip,
# based on the max rank of this type of chip, as compared to the max of the max ranks of other chip types
# present in the same mapped condition. The idea is to correct for the different genomic coverage
# of different chip types. We do not normalize simply based on the max rank in a given condition,
# to not penalize conditions with a lower number of expressed genes, and thus higher number
# of ex-aequo genes, and lower fractional max ranks.
# * compute weighted mean of normalized ranks per gene and mapped condition,
# weighted by the number of distinct ranks in each chip: we assume that chips with higher
# number of distinct ranks have a higher power for ranking genes. The higher power
# at ranking genes of a chip is taken into account by weighting the mean
# by the number of distinct ranks in the chip, not by "normalizing" away chips with lower max ranks,
# again, for not penalizing conditions with lower number of expressed genes;
# * also, insert max of max ranks of chip types represented in each mapped condition, in condition table
# (will allow to normalize ranks between conditions and data types), and sum of numbers of distinct ranks
# per gene and condition, in expression table (used to compute weigthed mean over all data types in a condition)
#
# For EST data:
# * "dense ranking" of the genes in each condition, based on the EST counts.
# Ranks stored in the temp table estRanking. All libraries in a condition are all considered
# together (no ranking first performed per library, as for other data types). This is because
# each library usually provides information about a very low number of genes.
# * Insert ranks into expression table.
# * for each condition, retrieve max rank.
# * insert max rank per condition into expression table (will allow to compute weigthed mean
# between data types in a condition, and to later normalize ranks between conditions and data types)
#
# For in situ hybridization data:
# * First, compute a score for each gene-condition, based on the detectionFlag
# of in situ evidence (present, absent) and the quality level (high quality, poor quality).
# Stored in temp table inSituRanking.
# * "dense ranking" of the genes in each condition, based on the score computed.
# Ranks stored in the temp table inSituRanking.
# All experiments are all considered together (no ranking first performed per experiment,
# as for other data types). This is because each in situ experiment usually study
# a limited number of genes.
# * Insert ranks into expression table
# * for each condition, retrieve max rank
# * insert max rank per condition into condition table (will allow to compute weigthed mean
# between data types in a condition, and to later normalize ranks between conditions and data types)
#
#
# A message that I wrote when transitioning to Bgee 15. Maybe it'll help to understand our reaosoning
# for the within-data-type normalization. It's not formally code documentation, but I wanted
# to save that message:
#
# # Affymetrix
# For Affymetrix data, we perform a first within-data-type normalization, before
# the between-data-types normalization: on different chip types we have different genomic coverages,
# and we need to normalize for that.
# So for each chip type, we retrieve the maximum rank among all chips of this chip type.
# And then for each chip, we normalize the ranks by using the max of all chip type max ranks,
# as compared to the max rank of the chip type of the chip being normalized:
# normalizedGeneChipRank = geneChipRank * (1 + MaxOfChipTypeRankMax/ChipTypeRankMax) / 2
# (and actually what we use is the max of the chip types used *in the condition* we're normalizing
# the ranks for)
#
# # RNA-Seq
# For RNA-Seq data, we don't do that. Our assumption (which was wrong, but even more so now)
# was that we had access to the same set of genes to be ranked in all libraries, the genomic coverage
# was the same in all libraries, and we wouldn't need within-data-type normalization.
# But now we distinguish between different RNA-Seq protocols: the genomic coverage is officially
# not the same anymore.
# So, wouldn't we need a normalization between different protocols, as we normalize between
# different Affymetrix chip types? Based on the maximum rank for each protocol among the libraries
# generated using this protocol. It would be:
# normalizedGeneLibraryRank = geneLibraryRank * (1 + MaxOfProtocolRankMax/ProtocolRankMax) / 2
#
# # Between-data-types normalization
# Then we normalize between conditions and data types. It is notably based on the maximum rank
# in each condition for each data type. For Affymetrix data, the max rank considered for a condition
# is the max rank of the chip types used in the condition (and not the max rank of the chips used
# in the condition, same logic as described above).
# For RNA-Seq data, we used to apply the same max rank to all the conditions of a species
# (same idea as above, genomic coverage being the same).
# So again, shouldn't we apply the same principle to RNA-Seq data? The max rank considered
# for a condition would be the max rank of the protocols used in that condition?
#
# # scRna-Seq full-length
# OK, I guess I consider that we have access to only one protocol. We don't store protocol information.
# And so, no within-data-type normalization, same max rank considered in all conditions of a species?


sub compute_update_global_ranks {
    my ($batchRef, $selfRanks) = @_;
    my $batchLength = scalar @{ $batchRef };

    # Connection to database in the parent process must have been closed
    # before calling this sub, otherwise it will generate errors
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289).
    # Get a new database connection for each thread.
    my $dbh_thread = Utils::connect_bgee_db($bgee_connector);

    #****************************************
    # PREPARE QUERIES
    #****************************************
    my $estRankField = 'estRank';
    my $estMaxRankField = 'estMaxRank';
    my $inSituRankField = 'inSituRank';
    my $inSituMaxRankField = 'inSituMaxRank';
    my $affyRankField = 'affymetrixMeanRank';
    my $affyDistinctRankSum = 'affymetrixDistinctRankSum';
    my $affyMaxRankField = 'affymetrixMaxRank';
    my $rnaSeqRankField = 'rnaSeqMeanRank';
    my $rnaSeqDistinctRankSum = 'rnaSeqDistinctRankSum';
    my $rnaSeqMaxRankField = 'rnaSeqMaxRank';
    my $scRnaSeqRankField = 'scRnaSeqFullLengthMeanRank';
    my $scRnaSeqDistinctRankSum = 'scRnaSeqFullLengthDistinctRankSum';
    my $scRnaSeqMaxRankField = 'scRnaSeqFullLengthMaxRank';
    if (!$selfRanks) {
        $estRankField = 'estGlobalRank';
        $estMaxRankField = 'estGlobalMaxRank';
        $inSituRankField = 'inSituGlobalRank';
        $inSituMaxRankField = 'inSituGlobalMaxRank';
        $affyRankField = 'affymetrixGlobalMeanRank';
        $affyDistinctRankSum = 'affymetrixGlobalDistinctRankSum';
        $affyMaxRankField = 'affymetrixGlobalMaxRank';
        $rnaSeqRankField = 'rnaSeqGlobalMeanRank';
        $rnaSeqDistinctRankSum = 'rnaSeqGlobalDistinctRankSum';
        $rnaSeqMaxRankField = 'rnaSeqGlobalMaxRank';
        $scRnaSeqRankField = 'scRnaSeqFullLengthGlobalMeanRank';
        $scRnaSeqDistinctRankSum = 'scRnaSeqFullLengthGlobalDistinctRankSum';
        $scRnaSeqMaxRankField = 'scRnaSeqFullLengthGlobalMaxRank';
    }

    # ***
    # *** Queries to update all ranks for all data types at the end ***
    # ***
    # Queries to create temp tables
    my $exprResultsTempTable = 'tempExprResults';
    # Note: for now we do the computation one globalConditionId at a time,
    # but we let open the possibility to compute several at once,
    # in which case the PRIMARY KEY and the SELECT clause should contain globalConditionId
    my $createExprResultsTempTableStmt  = $dbh_thread->prepare(
            'CREATE TEMPORARY TABLE '.$exprResultsTempTable.' '.
            # Super important to have this primary key to trigger the clauses 'ON DUPLICATE KEY UPDATE'
            '(PRIMARY KEY(bgeeGeneId))
             SELECT t1.bgeeGeneId,
                    t1.'.$estRankField.', t1.'.$inSituRankField.', t1.'.$affyRankField.',
                    t1.'.$rnaSeqRankField.', t1.'.$scRnaSeqRankField.',
                    t1.'.$affyDistinctRankSum.', t1.'.$rnaSeqDistinctRankSum.',
                    t1.'.$scRnaSeqDistinctRankSum.'
             FROM globalExpression AS t1 LIMIT 0;');
    my $dropExprResultsTempTableStmt = $dbh_thread->prepare('DROP TABLE '.$exprResultsTempTable);
    my $condResultsTempTable = 'tempCondResults';
    # Note: for now we do the computation one globalConditionId at a time,
    # but we let open the possibility to compute several at once.
    # This is why we keep this temp table mechanism for globalCond, while there will be only
    # one row inserted into this temp table per global condition.
    my $createCondResultsTempTableStmt  = $dbh_thread->prepare(
            'CREATE TEMPORARY TABLE '.$condResultsTempTable.' '.
            # Super important to have this primary key to trigger the clauses 'ON DUPLICATE KEY UPDATE'
            '(PRIMARY KEY(globalConditionId)) '.
            'SELECT t2.globalConditionId, t2.'.$estMaxRankField.', t2.'.$inSituMaxRankField.',
                    t2.'.$affyMaxRankField.', t2.'.$rnaSeqMaxRankField.', t2.'.$scRnaSeqMaxRankField.'
             FROM globalCond AS t2 LIMIT 0;');
    my $dropCondResultsTempTableStmt = $dbh_thread->prepare('DROP TABLE '.$condResultsTempTable);

    # Queries to finally update the globalExpression and globalCond tables.
    # Note: if we were to do the computations for several globalConditionIds at a time,
    # the ON clause of the INNER JOIN should be adapted accordingly.
    my $finalExprUpdateStmt = $dbh_thread->prepare(
            'UPDATE '.$exprResultsTempTable.' AS t2 STRAIGHT_JOIN globalExpression AS t1
                 ON t1.bgeeGeneId = t2.bgeeGeneId AND t1.globalConditionId = ?
             SET t1.'.$estRankField.' = t2.'.$estRankField.',
                 t1.'.$inSituRankField.' = t2.'.$inSituRankField.',
                 t1.'.$affyRankField.' = t2.'.$affyRankField.',
                 t1.'.$rnaSeqRankField.' = t2.'.$rnaSeqRankField.',
                 t1.'.$scRnaSeqRankField.' = t2.'.$scRnaSeqRankField.',
                 t1.'.$affyDistinctRankSum.' = t2.'.$affyDistinctRankSum.',
                 t1.'.$rnaSeqDistinctRankSum.' = t2.'.$rnaSeqDistinctRankSum.',
                 t1.'.$scRnaSeqDistinctRankSum.' = t2.'.$scRnaSeqDistinctRankSum);
    my $finalCondUpdateStmt = $dbh_thread->prepare(
            'UPDATE '.$condResultsTempTable.' AS t2 STRAIGHT_JOIN globalCond AS t1
                 ON t1.globalConditionId = t2.globalConditionId
             SET t1.'.$estMaxRankField.' = t2.'.$estMaxRankField.',
                 t1.'.$inSituMaxRankField.' = t2.'.$inSituMaxRankField.',
                 t1.'.$affyMaxRankField.' = t2.'.$affyMaxRankField.',
                 t1.'.$rnaSeqMaxRankField.' = t2.'.$rnaSeqMaxRankField.',
                 t1.'.$scRnaSeqMaxRankField.' = t2.'.$scRnaSeqMaxRankField);

    # ***
    # *** For EST data ***
    # ***
    # Query to store an association between each globalCondition and the EST libraries considered in it
    my $estToGlobalCondTableName = 'estToGlobalCond';
    my $estToGlobalCondStmt = $dbh_thread->prepare(cond_to_data_query_string(
        $estToGlobalCondTableName, 'estLibrary', 'estLibraryId', $batchLength, $selfRanks));
    my $isEstDataInGlobalCondStmt = $dbh_thread->prepare('SELECT 1 FROM '
            .$estToGlobalCondTableName.' WHERE globalConditionId = ? LIMIT 1');
    # query to rank genes by EST data
    # buid this table by using the fact that there is always an expressionId associated to each EST
    my $estRankingStmt = $dbh_thread->prepare('CREATE TEMPORARY TABLE estRanking (
            PRIMARY KEY(bgeeGeneId))
            SELECT STRAIGHT_JOIN expressedSequenceTag.bgeeGeneId, COUNT(1) AS estCount, 0000000.00 AS rawRank
            FROM '.$estToGlobalCondTableName.'
            INNER JOIN expressedSequenceTag ON '.$estToGlobalCondTableName.'.estLibraryId = expressedSequenceTag.estLibraryId
            WHERE '.$estToGlobalCondTableName.'.globalConditionId = ?
            GROUP BY expressedSequenceTag.bgeeGeneId');
    my $estDropLibRankStmt   = $dbh_thread->prepare('DROP TABLE estRanking');

    # Queries to compute gene ranks per condition from EST data
    my $estRanksForCondStmt = $dbh_thread->prepare('SELECT bgeeGeneId, estCount FROM estRanking
            ORDER BY estCount DESC');
    # if several genes at a same rank, we'll update them at once with a 'bgeeGeneId IN (?,?, ...)' clause.
    # If only one gene at a given rank, updated with the prepared statement below.
    my $estRankUpdateStart   = 'UPDATE estRanking SET rawRank = ?
                             WHERE bgeeGeneId ';
    my $estUpdateRankingStmt = $dbh_thread->prepare($estRankUpdateStart.'= ?');

    # Queries to insert into temp tables
    my $insertESTExprResults = $dbh_thread->prepare(
            'INSERT INTO '.$exprResultsTempTable.'
                 (bgeeGeneId, '.$estRankField.')
             SELECT bgeeGeneId, rawRank FROM estRanking
             ON DUPLICATE KEY UPDATE '.$estRankField.' = VALUES('.$estRankField.')');
    my $insertESTCondResults = $dbh_thread->prepare(
            'INSERT INTO '.$condResultsTempTable.'
                 (globalConditionId, '.$estMaxRankField.')
             VALUES(?, (SELECT MAX(rawRank) FROM estRanking))
             ON DUPLICATE KEY UPDATE '.$estMaxRankField.' = VALUES('.$estMaxRankField.')');


    # ***
    # *** For in situ data ***
    # ***
    # Query to store an association between each globalCondition and the in situ spots considered in it
    my $inSituToGlobalCondTableName = 'inSituToGlobalCond';
    my $inSituToGlobalCondStmt = $dbh_thread->prepare(cond_to_data_query_string(
        $inSituToGlobalCondTableName, 'inSituSpot', 'inSituSpotId', $batchLength, $selfRanks));
    my $isInSituDataInGlobalCondStmt = $dbh_thread->prepare('SELECT 1 FROM '
            .$inSituToGlobalCondTableName.' WHERE globalConditionId = ? LIMIT 1');
    # give a score to spots, depending on their detectionFlag and quality:
    # present - high quality = 1
    # present - low quality = 0.5
    # absent - low quality = -0.5
    # absent - high quality = -1
    # We sum these sores for each gene in a given mapped condition of the expression table,
    # and we will rank genes in each mapped condition based on this score.
    my $inSituRankingStmt = $dbh_thread->prepare("CREATE TEMPORARY TABLE inSituRanking (
            PRIMARY KEY(bgeeGeneId))
            SELECT STRAIGHT_JOIN inSituSpot.bgeeGeneId,

            SUM(
                IF(inSituSpot.detectionFlag = '$Utils::PRESENT_CALL' AND inSituData = '$Utils::HIGH_QUAL', 1,
                    IF(inSituSpot.detectionFlag = '$Utils::PRESENT_CALL', 0.5,
                        IF(inSituSpot.detectionFlag = '$Utils::ABSENT_CALL' AND inSituData = '$Utils::HIGH_QUAL', -1,
                            IF(inSituSpot.detectionFlag = '$Utils::ABSENT_CALL', -0.5, 0))))
            ) AS scoreSum,
            0000000.00 AS rawRank

            FROM $inSituToGlobalCondTableName
            INNER JOIN inSituSpot ON $inSituToGlobalCondTableName.inSituSpotId = inSituSpot.inSituSpotId
            WHERE $inSituToGlobalCondTableName.globalConditionId = ?
            AND inSituSpot.expressionId IS NOT NULL
            GROUP BY inSituSpot.bgeeGeneId");
    my $inSituDropSpotRankStmt  = $dbh_thread->prepare('DROP TABLE inSituRanking');

    # Queries to compute gene ranks per condition
    my $inSituRanksForCondStmt = $dbh_thread->prepare('SELECT bgeeGeneId, scoreSum FROM inSituRanking
            ORDER BY scoreSum DESC');
    # if several genes at a same rank, we'll update them at once with a 'geneId IN (?,?, ...)' clause.
    # If only one gene at a given rank, updated with the prepared statement below.
    my $inSituRankUpdateStart   = 'UPDATE inSituRanking SET rawRank = ?
                             WHERE bgeeGeneId ';
    my $inSituUpdateRankingStmt = $dbh_thread->prepare($inSituRankUpdateStart.'= ?');

    # Queries to insert into temp tables
    my $insertInSituExprResults = $dbh_thread->prepare(
            'INSERT INTO '.$exprResultsTempTable.'
                 (bgeeGeneId, '.$inSituRankField.')
             SELECT bgeeGeneId, rawRank FROM inSituRanking
             ON DUPLICATE KEY UPDATE '.$inSituRankField.' = VALUES('.$inSituRankField.')');
    my $insertInSituCondResults = $dbh_thread->prepare(
            'INSERT INTO '.$condResultsTempTable.'
                 (globalConditionId, '.$inSituMaxRankField.')
             VALUES(?, (SELECT MAX(rawRank) FROM inSituRanking))
             ON DUPLICATE KEY UPDATE '.$inSituMaxRankField.' = VALUES('.$inSituMaxRankField.')');


    # ***
    # *** For Affymetrix data ***
    # ***
    # Query to store an association between each globalCondition and the Affymetrix chips considered in it
    my $affyToGlobalCondTableName = 'affyToGlobalCond';
    my $affyToGlobalCondStmt = $dbh_thread->prepare(cond_to_data_query_string(
        $affyToGlobalCondTableName, 'affymetrixChip', 'bgeeAffymetrixChipId', $batchLength, $selfRanks));
    my $isAffyDataInGlobalCondStmt = $dbh_thread->prepare('SELECT 1 FROM '
            .$affyToGlobalCondTableName.' WHERE globalConditionId = ? LIMIT 1');

    # For each chip type, we have computed its max rank over any condition.
    # Now we retrieve the chip types present in the requested global condition,
    # and store the max of their max ranks (this is different than retrieving the max rank
    # of the *chips* present in the condition; the idea is to correct for the different
    # genetic coverages of the chip types, while not penalizing conditions with lower numbers
    # of expressed genes, and thus more ex-aequo genes with no signal of expression detected)
    my $affyCondMaxStmt = $dbh_thread->prepare(
            'SELECT STRAIGHT_JOIN MAX(t3.chipTypeMaxRank) AS maxRank '.
            # Retrieve all chips present in each globalCondition
            'FROM '.$affyToGlobalCondTableName.' AS t1 '.
            # and the max rank from all chip types present in this condition
            'INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
             INNER JOIN chipType AS t3 ON t2.chipTypeId = t3.chipTypeId
             WHERE t1.globalConditionId = ?');

    # Normalize ranks between chip types for the globalCondition
    my $affyChipNormRankStmt = $dbh_thread->prepare(
            'CREATE TEMPORARY TABLE chipNormalizedRank (
                PRIMARY KEY(bgeeGeneId, bgeeAffymetrixChipId))
             SELECT STRAIGHT_JOIN
                 affymetrixProbeset.bgeeGeneId, affymetrixProbeset.bgeeAffymetrixChipId, '.
             # between-chip-types normalization:
             # (we use the best normalized rank of genes from each Affymetrix chip,
             # to not count genes multiple times because of multiple probesets)
          '      (MIN(affymetrixProbeset.rawRank) +
                       (MIN(affymetrixProbeset.rawRank) * (? / chipType.chipTypeMaxRank)) )
                  /2 AS rawRank
             FROM '.$affyToGlobalCondTableName.' '.
             # get the max rank of the chip type corresponding to each chip
             'INNER JOIN affymetrixChip ON '.$affyToGlobalCondTableName.'.bgeeAffymetrixChipId = affymetrixChip.bgeeAffymetrixChipId
              INNER JOIN chipType ON affymetrixChip.chipTypeId = chipType.chipTypeId '.
             # get for each gene and chip the min rank from the probesets
             'INNER JOIN affymetrixProbeset ON affymetrixChip.bgeeAffymetrixChipId = affymetrixProbeset.bgeeAffymetrixChipId '.
             # where and group by.
             'WHERE affymetrixProbeset.expressionId IS NOT NULL AND '.$affyToGlobalCondTableName.'.globalConditionId = ?
              GROUP BY affymetrixProbeset.bgeeGeneId, affymetrixProbeset.bgeeAffymetrixChipId');
    my $dropAffyChipNormRank = $dbh_thread->prepare('DROP TABLE chipNormalizedRank');

    # Queries to insert into temp tables
    my $insertAffyExprResults = $dbh_thread->prepare(
            'INSERT INTO '.$exprResultsTempTable.'
                 (bgeeGeneId, '.$affyRankField.', '.$affyDistinctRankSum.')
             SELECT STRAIGHT_JOIN
                   chipNormalizedRank.bgeeGeneId,
                   SUM(chipNormalizedRank.rawRank * affymetrixChip.chipDistinctRankCount)
                       /SUM(affymetrixChip.chipDistinctRankCount) AS meanRank,
                   SUM(affymetrixChip.chipDistinctRankCount) AS distinctRankCountSum
             FROM chipNormalizedRank
             INNER JOIN affymetrixChip ON affymetrixChip.bgeeAffymetrixChipId = chipNormalizedRank.bgeeAffymetrixChipId
             GROUP BY chipNormalizedRank.bgeeGeneId
             ON DUPLICATE KEY UPDATE '.$affyRankField.' = VALUES('.$affyRankField.'),
                 '.$affyDistinctRankSum.' = VALUES('.$affyDistinctRankSum.')');
    my $insertAffyCondResults = $dbh_thread->prepare(
            'INSERT INTO '.$condResultsTempTable.'
                 (globalConditionId, '.$affyMaxRankField.')
             VALUES(?, ?)
             ON DUPLICATE KEY UPDATE '.$affyMaxRankField.' = VALUES('.$affyMaxRankField.')');


    # ***
    # *** For RNA-Seq data and scRNA-Seq data ***
    # ***
    # Query to store an association between each globalCondition and the RNA-Seq libs considered in it
    my $rnaSeqToGlobalCondTableName = 'rnaSeqToGlobalCond';
    my $scRnaSeqToGlobalCondTableName = 'scRnaSeqToGlobalCond';
    my $rnaSeqToGlobalCondStmt = getRnaSeqToGlobalCondStmt($dbh_thread,
            $rnaSeqToGlobalCondTableName, $batchLength, $selfRanks, 0);
    my $scRnaSeqToGlobalCondStmt = getRnaSeqToGlobalCondStmt($dbh_thread,
            $scRnaSeqToGlobalCondTableName, $batchLength, $selfRanks, 1);

    my $isRnaSeqDataInGlobalCondStmt = $dbh_thread->prepare('SELECT 1 FROM '
            .$rnaSeqToGlobalCondTableName.' WHERE globalConditionId = ? LIMIT 1');
    my $isScRnaSeqDataInGlobalCondStmt = $dbh_thread->prepare('SELECT 1 FROM '
            .$scRnaSeqToGlobalCondTableName.' WHERE globalConditionId = ? LIMIT 1');
    # For RNA-Seq data, we used to apply the same max rank to all conditions of a same species,
    # because we considered that we always had access to all expression information
    # for the same set of genes in all libraries. This is not true anymore (and it actually never was),
    # as we do now distinguish between different RNA-Seq protocols targeting different biotypes.
    #
    # So we have a similar approach to Affymetrix data, where we store
    # for each condition the max rank among the *chip types* present in a condition
    # (and not the max rank among the actual chips present in a condition): we store
    # the max rank among the *protocols* present in a condition. We will then do
    # a first "between-protocol" normalization, as we do a first "between-chip-type" normalization
    # in Affymetrix data.
    my $rnaSeqCondMaxStmt = getRnaSeqCondMaxStmt($dbh_thread, $rnaSeqToGlobalCondTableName);
    my $scRnaSeqCondMaxStmt = getRnaSeqCondMaxStmt($dbh_thread, $scRnaSeqToGlobalCondTableName);

    # Normalize ranks between protocols for the globalCondition
    my $rnaSeqAnnotSampleNormRankStmt = getRnaSeqAnnotSampleNormRankStmt($dbh_thread, $rnaSeqToGlobalCondTableName,
            'rnaSeqAnnotSampleNormalizedRank');
    my $scRnaSeqAnnotSampleNormRankStmt = getRnaSeqAnnotSampleNormRankStmt($dbh_thread, $scRnaSeqToGlobalCondTableName,
            'scRnaSeqAnnotSampleNormalizedRank');
    my $dropRNASeqLibNormRank = $dbh_thread->prepare('DROP TABLE rnaSeqAnnotSampleNormalizedRank');
    my $dropScRNASeqLibNormRank = $dbh_thread->prepare('DROP TABLE scRnaSeqAnnotSampleNormalizedRank');

    # Queries to insert into temp tables
    my $insertRnaSeqExprResults = getInsertRnaSeqExprResultsStmt($dbh_thread, $exprResultsTempTable,
            'rnaSeqAnnotSampleNormalizedRank', $rnaSeqRankField, $rnaSeqDistinctRankSum);
    my $insertScRnaSeqExprResults = getInsertRnaSeqExprResultsStmt($dbh_thread, $exprResultsTempTable,
            'scRnaSeqAnnotSampleNormalizedRank', $scRnaSeqRankField, $scRnaSeqDistinctRankSum);

    my $insertRnaSeqCondResults = $dbh_thread->prepare(
            'INSERT INTO '.$condResultsTempTable.'
                 (globalConditionId, '.$rnaSeqMaxRankField.')
             VALUES(?, ?)
             ON DUPLICATE KEY UPDATE '.$rnaSeqMaxRankField.' = VALUES('.$rnaSeqMaxRankField.')');
    my $insertScRnaSeqCondResults = $dbh_thread->prepare(
            'INSERT INTO '.$condResultsTempTable.'
                 (globalConditionId, '.$scRnaSeqMaxRankField.')
             VALUES(?, ?)
             ON DUPLICATE KEY UPDATE '.$scRnaSeqMaxRankField.' = VALUES('.$scRnaSeqMaxRankField.')');



    #****************************************
    # EXECUTE GENERAL QUERIES FOR THE WHOLE CONDITION BATCH
    #****************************************
    $estToGlobalCondStmt->execute(@{ $batchRef }) or die $estToGlobalCondStmt->errstr;
    $inSituToGlobalCondStmt->execute(@{ $batchRef }) or die $inSituToGlobalCondStmt->errstr;
    $affyToGlobalCondStmt->execute(@{ $batchRef }) or die $affyToGlobalCondStmt->errstr;
    $rnaSeqToGlobalCondStmt->execute(@{ $batchRef }) or die $rnaSeqToGlobalCondStmt->errstr;
    $scRnaSeqToGlobalCondStmt->execute(@{ $batchRef }) or die $scRnaSeqToGlobalCondStmt->errstr;


    #****************************************
    # COMPUTE RANKS PER GLOBAL CONDITION
    #****************************************
    for my $k ( 0..$batchLength-1 ) {
        Utils::start_transaction($dbh_thread);
        my $globalCondId = ${$batchRef}[$k];
        my $hasEstData = 0;
        my $hasInSituData = 0;
        my $hasAffyData = 0;
        my $hasRnaSeqData = 0;
        my $hasScRnaSeqData = 0;

        # ***
        # *** Create the temp results table ***
        # ***
        $createExprResultsTempTableStmt->execute() or die $createExprResultsTempTableStmt->errstr;
        $createCondResultsTempTableStmt->execute() or die $createCondResultsTempTableStmt->errstr;

        # ***
        # *** For EST data ***
        # ***
        $isEstDataInGlobalCondStmt->execute($globalCondId) or die $isEstDataInGlobalCondStmt->errstr;
        if ($isEstDataInGlobalCondStmt->fetchrow_arrayref()) {
            $hasEstData = 1;
            $estRankingStmt->execute($globalCondId) or die $estRankingStmt->errstr;
            $estRanksForCondStmt->execute() or die $estRanksForCondStmt->errstr;
            my @estResults = map { {"id" => $_->[0], "val" => $_->[1]} }
                    @{$estRanksForCondStmt->fetchall_arrayref};
            compute_dense_ranking_update_tmp_ranks($dbh_thread, \@estResults,
                    $estUpdateRankingStmt, $estRankUpdateStart);
            # Insert the results in the temp result tables
            $insertESTExprResults->execute() or die $insertESTExprResults->errstr;
            $insertESTCondResults->execute($globalCondId) or die $insertESTCondResults->errstr;
            # Drop the temp ranking table
            $estDropLibRankStmt->execute() or die $estDropLibRankStmt->errstr;
        }

        # ***
        # *** For in situ data ***
        # ***
        $isInSituDataInGlobalCondStmt->execute($globalCondId) or die $isInSituDataInGlobalCondStmt->errstr;
        if ($isInSituDataInGlobalCondStmt->fetchrow_arrayref()) {
            $hasInSituData = 1;
            $inSituRankingStmt->execute($globalCondId) or die $inSituRankingStmt->errstr;
            $inSituRanksForCondStmt->execute() or die $inSituRanksForCondStmt->errstr;
            my @inSituResults = map { {"id" => $_->[0], "val" => $_->[1]} }
                    @{$inSituRanksForCondStmt->fetchall_arrayref};
            compute_dense_ranking_update_tmp_ranks($dbh_thread, \@inSituResults,
                    $inSituUpdateRankingStmt, $inSituRankUpdateStart);
            # Insert the results in the temp result tables
            $insertInSituExprResults->execute() or die $insertInSituExprResults->errstr;
            $insertInSituCondResults->execute($globalCondId) or die $insertInSituCondResults->errstr;
            # Drop the temp ranking table
            $inSituDropSpotRankStmt->execute() or die $inSituDropSpotRankStmt->errstr;
        }

        # ***
        # *** For Affymetrix data ***
        # ***
        $isAffyDataInGlobalCondStmt->execute($globalCondId) or die $isAffyDataInGlobalCondStmt->errstr;
        if ($isAffyDataInGlobalCondStmt->fetchrow_arrayref()) {
            $hasAffyData = 1;
            $affyCondMaxStmt->execute($globalCondId) or die $affyCondMaxStmt->errstr;
            my $chipMax = undef;
            if (my @chipMaxRes = $affyCondMaxStmt->fetchrow_array()) {
                $chipMax = $chipMaxRes[0];
            }
            # Since we checked if there was Affy data in the condition we should always have a max rank
            if (!$chipMax) {
                die("Could not find chipMax for globalConditionId: $globalCondId\n");
            }
            $affyChipNormRankStmt->execute($chipMax, $globalCondId) or die $affyChipNormRankStmt->errstr;
            # Insert the results in the temp result tables
            $insertAffyExprResults->execute() or die $insertAffyExprResults->errstr;
            $insertAffyCondResults->execute($globalCondId, $chipMax) or die $insertAffyCondResults->errstr;
            # Drop the temp ranking table
            $dropAffyChipNormRank->execute() or die $dropAffyChipNormRank->errstr;
        }

        # ***
        # *** For RNA-Seq data ***
        # ***
        $isRnaSeqDataInGlobalCondStmt->execute($globalCondId) or die $isRnaSeqDataInGlobalCondStmt->errstr;
        if ($isRnaSeqDataInGlobalCondStmt->fetchrow_arrayref()) {
            $hasRnaSeqData = 1;
            $rnaSeqCondMaxStmt->execute($globalCondId) or die $rnaSeqCondMaxStmt->errstr;
            my $rnaSeqLibMax = undef;
            if (my @rnaSeqLibMaxRes = $rnaSeqCondMaxStmt->fetchrow_array()) {
                $rnaSeqLibMax = $rnaSeqLibMaxRes[0];
            }
            # Since we checked if there was RNA-Seq data in the condition we should always have a max rank
            if (!$rnaSeqLibMax) {
                die("Could not find bulk rnaSeqLibMax for globalConditionId: $globalCondId\n");
            }
            $rnaSeqAnnotSampleNormRankStmt->execute($rnaSeqLibMax, $globalCondId) or die $rnaSeqAnnotSampleNormRankStmt->errstr;
            # Insert the results in the temp result tables
            $insertRnaSeqExprResults->execute() or die $insertRnaSeqExprResults->errstr;
            $insertRnaSeqCondResults->execute($globalCondId, $rnaSeqLibMax) or die $insertRnaSeqCondResults->errstr;
            # Drop the temp ranking table
            $dropRNASeqLibNormRank->execute() or die $dropRNASeqLibNormRank->errstr;
        }

        # ***
        # *** For scRNA-Seq ***
        # ***
        $isScRnaSeqDataInGlobalCondStmt->execute($globalCondId) or die $isScRnaSeqDataInGlobalCondStmt->errstr;
        if ($isScRnaSeqDataInGlobalCondStmt->fetchrow_arrayref()) {
            $hasScRnaSeqData = 1;
            $scRnaSeqCondMaxStmt->execute($globalCondId) or die $scRnaSeqCondMaxStmt->errstr;
            my $rnaSeqLibMax = undef;
            if (my @rnaSeqLibMaxRes = $scRnaSeqCondMaxStmt->fetchrow_array()) {
                $rnaSeqLibMax = $rnaSeqLibMaxRes[0];
            }
            # Since we checked if there was RNA-Seq data in the condition we should always have a max rank
            if (!$rnaSeqLibMax) {
                die("Could not find single-cell rnaSeqLibMax for globalConditionId: $globalCondId\n");
            }
            $scRnaSeqAnnotSampleNormRankStmt->execute($rnaSeqLibMax, $globalCondId) or die $scRnaSeqAnnotSampleNormRankStmt->errstr;
            # Insert the results in the temp result tables
            $insertScRnaSeqExprResults->execute() or die $insertScRnaSeqExprResults->errstr;
            $insertScRnaSeqCondResults->execute($globalCondId, $rnaSeqLibMax) or die $insertScRnaSeqCondResults->errstr;
            # Drop the temp ranking table
            $dropScRNASeqLibNormRank->execute() or die $dropScRNASeqLibNormRank->errstr;
        }

        # ***
        # *** Update globalExpression and globalCond tables and drop temp result tables
        # ***
        if ($hasEstData || $hasInSituData || $hasAffyData || $hasRnaSeqData || $hasScRnaSeqData) {
            $finalExprUpdateStmt->execute($globalCondId) or die $finalExprUpdateStmt->errstr;
            $finalCondUpdateStmt->execute() or die $finalCondUpdateStmt->errstr;
        }
        $dropExprResultsTempTableStmt->execute() or die $dropExprResultsTempTableStmt->errstr;
        $dropCondResultsTempTableStmt->execute() or die $dropCondResultsTempTableStmt->errstr;

        # ***
        # *** Commit
        # ***
        $dbh_thread->commit() or die('Failed commit');
        printf("globalConditionId: %s - PID: %s - %d/%d\n", $globalCondId, $$, $k+1, $batchLength);
    }

    $createExprResultsTempTableStmt->finish or die('Failed finish');
    $dropExprResultsTempTableStmt->finish or die('Failed finish');
    $createCondResultsTempTableStmt->finish or die('Failed finish');
    $dropCondResultsTempTableStmt->finish or die('Failed finish');
    $finalExprUpdateStmt->finish or die('Failed finish');
    $finalCondUpdateStmt->finish or die('Failed finish');

    $estToGlobalCondStmt->finish or die('Failed finish');
    $isEstDataInGlobalCondStmt->finish or die('Failed finish');
    $estRankingStmt->finish or die('Failed finish');
    $estDropLibRankStmt->finish or die('Failed finish');
    $estRanksForCondStmt->finish or die('Failed finish');
    $estUpdateRankingStmt->finish or die('Failed finish');
    $insertESTExprResults->finish or die('Failed finish');
    $insertESTCondResults->finish or die('Failed finish');

    $inSituToGlobalCondStmt->finish or die('Failed finish');
    $isInSituDataInGlobalCondStmt->finish or die('Failed finish');
    $inSituRankingStmt->finish or die('Failed finish');
    $inSituDropSpotRankStmt->finish or die('Failed finish');
    $inSituRanksForCondStmt->finish or die('Failed finish');
    $inSituUpdateRankingStmt->finish or die('Failed finish');
    $insertInSituExprResults->finish or die('Failed finish');
    $insertInSituCondResults->finish or die('Failed finish');

    $affyToGlobalCondStmt->finish or die('Failed finish');
    $isAffyDataInGlobalCondStmt->finish or die('Failed finish');
    $affyCondMaxStmt->finish or die('Failed finish');
    $affyChipNormRankStmt->finish or die('Failed finish');
    $dropAffyChipNormRank->finish or die('Failed finish');
    $insertAffyExprResults->finish or die('Failed finish');
    $insertAffyCondResults->finish or die('Failed finish');

    $rnaSeqToGlobalCondStmt->finish or die('Failed finish');
    $isRnaSeqDataInGlobalCondStmt->finish or die('Failed finish');
    $rnaSeqCondMaxStmt->finish or die('Failed finish');
    $rnaSeqAnnotSampleNormRankStmt->finish or die('Failed finish');
    $dropRNASeqLibNormRank->finish or die('Failed finish');
    $insertRnaSeqExprResults->finish or die('Failed finish');
    $insertRnaSeqCondResults->finish or die('Failed finish');

    $scRnaSeqToGlobalCondStmt->finish or die('Failed finish');
    $isScRnaSeqDataInGlobalCondStmt->finish or die('Failed finish');
    $scRnaSeqCondMaxStmt->finish or die('Failed finish');
    $scRnaSeqAnnotSampleNormRankStmt->finish or die('Failed finish');
    $dropScRNASeqLibNormRank->finish or die('Failed finish');
    $insertScRnaSeqExprResults->finish or die('Failed finish');
    $insertScRnaSeqCondResults->finish or die('Failed finish');

    $dbh_thread->disconnect();
}

sub compute_dense_ranking_update_tmp_ranks {
    my ($dbh, $resultArrRef, $updateRankingStmt, $updateRankingStmtStart) = @_;

    my %sorted = Utils::dense_ranking(@{ $resultArrRef });
    # we get ranks as keys, with reference to an array of expression IDs with that rank as value
    my %reverseHash = Utils::revhash(%sorted);

    for my $rank ( keys %reverseHash ){
        my $geneIds_arrRef = $reverseHash{$rank};
        my @geneIds_arr = @$geneIds_arrRef;
        my $exprCount = scalar @geneIds_arr;
        if ( $exprCount == 1 ){
            my $geneId = $geneIds_arr[0];
            $updateRankingStmt->execute($rank, $geneId)  or die $updateRankingStmt->errstr;
        } else {
            my $query = $updateRankingStmtStart.'IN (';
            for ( my $i = 0; $i < $exprCount; $i++ ){
                if ( $i > 0 ){
                    $query .= ', ';
                }
                $query .= '?';
            }
            $query .= ')';
            my $rankMultiUpdateStmt = $dbh->prepare($query);
            $rankMultiUpdateStmt->execute($rank, @geneIds_arr)  or die $rankMultiUpdateStmt->errstr;
        }
    }
}

sub cond_to_data_query_string {
    my ($tempTableName, $sourceDataTableName, $sourceDataIdName, $condCount, $selfRanks, $isRnaSeq, $isSingleCell) = @_;

    my $sql = 'CREATE TEMPORARY TABLE '.$tempTableName.' (
                   PRIMARY KEY(globalConditionId, '.$sourceDataIdName.'))
               SELECT DISTINCT t1.globalConditionId, t4.'.$sourceDataIdName.' '.
               # Retrieve the valid raw conditions mapped to each globalCondition
               'FROM globalCond AS t1
               INNER JOIN globalCondToCond AS t2 ON t1.globalConditionId = t2.globalConditionId ';
    if ( $selfRanks ){
        $sql .= "AND t2.conditionRelationOrigin = 'self' ";
    } else {
        $sql .= "AND t2.conditionRelationOrigin IN ('self', 'descendant') ";
    }
    $sql .= 'INNER JOIN cond AS t3 ON t3.exprMappedConditionId = t2.conditionId '.
            # retrieve the data present in this globalCondition
            'INNER JOIN '.$sourceDataTableName.' AS t4 ON t3.conditionId = t4.conditionId ';
    if ($isRnaSeq) {
        $sql .= 'INNER JOIN rnaSeqLibrary AS t5 ON t4.rnaSeqLibraryId = t5.rnaSeqLibraryId '.
                'WHERE t5.rnaSeqTechnologyIsSingleCell = '.($isSingleCell ? '1' : '0').' AND ';
    } else {
        $sql .= 'WHERE ';
    }
    $sql .= 't1.globalConditionId';
    if ($condCount == 1) {
        $sql .= ' = ?';
    } else {
        $sql .= ' IN (';
        for (my $i = 0; $i < $condCount; $i++) {
            if ($i > 0) {
                $sql .= ', ';
            }
            $sql .= '?';
        }
        $sql .= ')'
    }

    return $sql;
}

sub getRnaSeqToGlobalCondStmt {
    my ($dbh, $rnaSeqToGlobalCondTableName, $batchLength, $selfRanks, $isSingleCell) = @_;
    return $dbh->prepare(cond_to_data_query_string(
        $rnaSeqToGlobalCondTableName, 'rnaSeqLibraryAnnotatedSample', 'rnaSeqLibraryAnnotatedSampleId',
        $batchLength, $selfRanks, 1, $isSingleCell));
}

sub getRnaSeqCondMaxStmt {
    my ($dbh, $rnaSeqToGlobalCondTableName) = @_;
    return $dbh->prepare(
            'SELECT STRAIGHT_JOIN MAX(t4.maxRank) AS maxRank '.
            # Retrieve all libraries present in the globalCondition
            'FROM '.$rnaSeqToGlobalCondTableName.' AS t1 '.
            # and the max rank from all rna-seq protocols present in this condition
            'INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId
             INNER JOIN rnaSeqLibrary AS t3 ON t2.rnaSeqLibraryId = t3.rnaSeqLibraryId
             INNER JOIN globalCond AS t3bis ON t1.globalConditionId = t3bis.globalConditionId
             INNER JOIN rnaSeqPopulationCaptureSpeciesMaxRank AS t4 ON t3.rnaSeqPopulationCaptureId = t4.rnaSeqPopulationCaptureId
                 AND t3bis.speciesId = t4.speciesId
             WHERE t1.globalConditionId = ?');
}

sub getRnaSeqAnnotSampleNormRankStmt {
    my ($dbh, $rnaSeqToGlobalCondTableName, $rnaSeqAnnotSampleNormalizedRankTableName) = @_;
    return $dbh->prepare(
            'CREATE TEMPORARY TABLE '.$rnaSeqAnnotSampleNormalizedRankTableName.' (
                PRIMARY KEY(bgeeGeneId, rnaSeqLibraryAnnotatedSampleId))
             SELECT STRAIGHT_JOIN
                 rnaSeqLibraryAnnotatedSampleGeneResult.bgeeGeneId, rnaSeqLibraryAnnotatedSampleGeneResult.rnaSeqLibraryAnnotatedSampleId, '.
             # between-protocols normalization:
          '      (rnaSeqLibraryAnnotatedSampleGeneResult.rawRank +
                       (rnaSeqLibraryAnnotatedSampleGeneResult.rawRank * (? / rnaSeqPopulationCaptureSpeciesMaxRank.maxRank)) )
                  /2 AS rawRank
             FROM '.$rnaSeqToGlobalCondTableName.' '.
             # get the max rank of the protocol corresponding to each library
             'INNER JOIN rnaSeqLibraryAnnotatedSample ON '.$rnaSeqToGlobalCondTableName.'.rnaSeqLibraryAnnotatedSampleId = rnaSeqLibraryAnnotatedSample.rnaSeqLibraryAnnotatedSampleId
              INNER JOIN rnaSeqLibrary ON rnaSeqLibraryAnnotatedSample.rnaSeqLibraryId = rnaSeqLibrary.rnaSeqLibraryId
              INNER JOIN globalCond ON '.$rnaSeqToGlobalCondTableName.'.globalConditionId = globalCond.globalConditionId
              INNER JOIN rnaSeqPopulationCaptureSpeciesMaxRank ON rnaSeqLibrary.rnaSeqPopulationCaptureId = rnaSeqPopulationCaptureSpeciesMaxRank.rnaSeqPopulationCaptureId
                  AND globalCond.speciesId = rnaSeqPopulationCaptureSpeciesMaxRank.speciesId '.
             # get for each gene and library the rank from the libraries
             'INNER JOIN rnaSeqLibraryAnnotatedSampleGeneResult ON rnaSeqLibraryAnnotatedSample.rnaSeqLibraryAnnotatedSampleId = rnaSeqLibraryAnnotatedSampleGeneResult.rnaSeqLibraryAnnotatedSampleId '.
             'WHERE rnaSeqLibraryAnnotatedSampleGeneResult.expressionId IS NOT NULL AND '.$rnaSeqToGlobalCondTableName.'.globalConditionId = ?');
}

sub getInsertRnaSeqExprResultsStmt {
    my ($dbh, $exprResultsTempTable, $rnaSeqAnnotSampleNormalizedRankTableName,
        $rnaSeqRankField, $rnaSeqDistinctRankSum) = @_;
    return $dbh->prepare(
            'INSERT INTO '.$exprResultsTempTable.'
                 (bgeeGeneId, '.$rnaSeqRankField.', '.$rnaSeqDistinctRankSum.')
             SELECT STRAIGHT_JOIN '.
                   $rnaSeqAnnotSampleNormalizedRankTableName.'.bgeeGeneId,
                   SUM('.$rnaSeqAnnotSampleNormalizedRankTableName.'.rawRank * rnaSeqLibraryAnnotatedSample.rnaSeqLibraryAnnotatedSampleDistinctRankCount)
                       /SUM(rnaSeqLibraryAnnotatedSample.rnaSeqLibraryAnnotatedSampleDistinctRankCount) AS meanRank,
                   SUM(rnaSeqLibraryAnnotatedSample.rnaSeqLibraryAnnotatedSampleDistinctRankCount) AS distinctRankCountSum
             FROM '.$rnaSeqAnnotSampleNormalizedRankTableName.
            ' INNER JOIN rnaSeqLibraryAnnotatedSample ON rnaSeqLibraryAnnotatedSample.rnaSeqLibraryAnnotatedSampleId = '.$rnaSeqAnnotSampleNormalizedRankTableName.'.rnaSeqLibraryAnnotatedSampleId
             GROUP BY '.$rnaSeqAnnotSampleNormalizedRankTableName.'.bgeeGeneId
             ON DUPLICATE KEY UPDATE '.$rnaSeqRankField.' = VALUES('.$rnaSeqRankField.'),
                 '.$rnaSeqDistinctRankSum.' = VALUES('.$rnaSeqDistinctRankSum.')');
}

# As of Bgee 15.0, we only compute global ranks. If we wanted to change that,
# we would need a loop
my $selfRanks = 0;

my @conditions = @cond_ids;
if (!@cond_ids) {
    my $dbh = Utils::connect_bgee_db($bgee_connector);

#    # Retrieve codition parameter combinations to compute ranks for.
#    my $condParamCombinationsArrRef = Utils::get_cond_param_combinations($dbh);
#    print "All condition parameter combinations to compute:\n";
#    print "@$_\n"  for @{$condParamCombinationsArrRef};


    my $condSql = 'SELECT DISTINCT globalConditionId FROM globalCond ';
    if ($selfRanks) {
        $condSql .= 'WHERE estMaxRank IS NULL AND inSituMaxRank IS NULL AND affymetrixMaxRank IS NULL
                     AND rnaSeqMaxRank IS NULL AND scRnaSeqFullLengthMaxRank IS NULL';
    } else {
        $condSql .= 'WHERE estGlobalMaxRank IS NULL AND inSituGlobalMaxRank IS NULL
                     AND affymetrixGlobalMaxRank IS NULL AND rnaSeqGlobalMaxRank IS NULL
                     AND scRnaSeqFullLengthGlobalMaxRank IS NULL';
    }
    if ($cond_count > 0) {
        $condSql .= ' ORDER BY globalConditionId
                      LIMIT '.$cond_offset.', '.$cond_count;
    }
    my $queryConditions = $dbh->prepare($condSql);
    $queryConditions->execute()  or die $queryConditions->errstr;
    @conditions = map { $_->[0] } @{$queryConditions->fetchall_arrayref};

    # Disconnect the DBI connection open in parent process, otherwise it will generate errors
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
    $dbh->disconnect();
}
my $conditionCount = @conditions;
if ($conditionCount < $conds_per_job) {
    $conds_per_job = $conditionCount;
}
print "Rank computations per global condition, $conditionCount conditions retrieved...\n";
# We are going to compute/store ranks using parallelization,
# each child process will be responsible to compute/update ranks for a batch of global conditions
my $iterationCount = ceil($conditionCount/$conds_per_job);
my $parallel = $parallel_jobs;
if ($iterationCount < $parallel_jobs) {
    $parallel = $iterationCount;
}


# To save cores needed, we run the computation in the main thread when only one job is requested.
# And we do that in different loops to not need to use ForManager if not needed.
if ($parallel == 1) {
    while ( my @next_conds = splice(@conditions, 0, $conds_per_job) ) {
        print("\nStart batch of $conds_per_job conditions, process ID $$...\n");
        compute_update_global_ranks(\@next_conds, $selfRanks);
        print("\nDone batch of $conds_per_job conditions, process ID $$.\n");
    }
} else {
    print("Rank computations per condition with $parallel threads and $conds_per_job conditions per thread...\n");
    my $pm = new Parallel::ForkManager($parallel);
    while ( my @next_conds = splice(@conditions, 0, $conds_per_job) ) {
        # Forks and returns the pid for the child
        my $pid = $pm->start and next;
        print("\nStart batch of $conds_per_job conditions, process ID $$...\n");
        compute_update_global_ranks(\@next_conds, $selfRanks);
        print("\nDone batch of $conds_per_job conditions, process ID $$.\n");
        $pm->finish;
    }
    $pm->wait_all_children;
}
print("Rank computations per global condition done\n");

exit 0;
