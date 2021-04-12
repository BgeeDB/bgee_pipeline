#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw( time );
use diagnostics;

# Updates the ranks in affymetrixProbeset table and affymetrix mean rank in the expression table.
# Philippe Moret, created Oct 2015.
# Frederic Bastian, last updated June 2016.
# Frederic Bastian, last updated Feb. 2017: adapt to new conditions and new schema in Bgee 14
# Frederic Bastian, last updated Apr. 2021: allow parallelization as for RNA-Seq data

use Parallel::ForkManager;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Getopt::Long;
use POSIX qw/ceil/;

$|=1;

# Define arguments and their default value
my ($bgee_connector) = ('');
my ($parallel_jobs) = (15); # default 15 parallel threads used to compute ranks
my ($sample_offset) = (0);
my ($sample_count) = (0);
my %opts = ('bgee=s'          => \$bgee_connector, # Bgee connector string
            'parallel_jobs=i' => \$parallel_jobs,
            'sample_offset=i' => \$sample_offset,
            'sample_count=i'  => \$sample_count,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\t-parallel_jobs    Number of threads used to compute ranks
\t-sample_offset    The offset parameter to retrieve chips to compute ranks for
\t-sample_count     The row_count parameter to retrieve chips to compute ranks for
\n";
    exit 1;
}
if ($sample_offset < 0 || $sample_count < 0) {
    die('sample_offset and sample_count cannot be negative');
}
if ($sample_offset > 0 && $sample_count == 0) {
    die('sample_count must be provided if sample_offset is provided');
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


my $dbh = Utils::connect_bgee_db($bgee_connector);

# We set autocommit to 1 so that we can define transaction isolation level
# *before* starting the next transaction, see https://www.perlmonks.org/?node_id=1074673
$dbh->{'AutoCommit'} = 1;

##############################################
# COMPUTE RANKS PER CHIP                     #
##############################################
sub compute_update_rank_batch {
    my ($batchRef, $pid) = @_;
    my $batchLength = scalar @{ $batchRef };

    # Connection to database in the parent process must have been closed
    # before calling this sub, otherwise it will generate errors
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289).
    # Get a new database connection for each thread.
    my $dbh_thread = Utils::connect_bgee_db($bgee_connector);
    Utils::start_transaction($dbh_thread);

    my $affymProbesetStmt = $dbh_thread->prepare("SELECT bgeeGeneId, MAX(normalizedSignalIntensity) AS maxIntensity
                                           FROM affymetrixProbeset
                                           WHERE bgeeAffymetrixChipId = ?
                                           AND reasonForExclusion = '$Utils::CALL_NOT_EXCLUDED'
                                           GROUP BY bgeeGeneId
                                           ORDER BY maxIntensity DESC");

    # if several genes at a same rank, we'll update them at once with a 'bgeeGeneId IN (?,?, ...)' clause.
    # If only one gene at a given rank, updated with the prepared statement below.
    my $rankUpdateStart   = 'UPDATE affymetrixProbeset SET rank = ?
                             WHERE bgeeAffymetrixChipId = ? AND bgeeGeneId ';
    my $affymeProbesetUpdateStmt = $dbh_thread->prepare($rankUpdateStart.'= ?');
                                              
    for my $k ( 0..$batchLength-1 ) {

        # ======= Compute ranks ========
        my $sampleId = ${$batchRef}[$k];
        $affymProbesetStmt->execute($sampleId) or die $affymProbesetStmt->errstr;

        my @results = map { {'id' => $_->[0], 'val' => $_->[1]} } @{$affymProbesetStmt->fetchall_arrayref};
        my %sorted = Utils::fractionnal_ranking(@results);
        # we get ranks as keys, with reference to an array of gene IDs with that rank as value
        my %reverseHash = Utils::revhash(%sorted);

        # ======= Update ranks ========
        for my $rank ( keys %reverseHash ){
            my $geneIds_arrRef = $reverseHash{$rank};
            my @geneIds_arr = @$geneIds_arrRef;
            my $geneCount = scalar @geneIds_arr;
            # with the multi-update query, sometimes there are thousands of genes
            # with equal ranks and the query is extremely slow. Apparently it's better
            # to alway use the single-gene update query
#            for ( my $i = 0; $i < $geneCount; $i++ ){
#                $affymeProbesetUpdateStmt->execute($rank, $sampleId, $geneIds_arr[$i])
#                    or die $affymeProbesetUpdateStmt->errstr;
#            }
            if ( $geneCount == 1 ){
                my $geneId = $geneIds_arr[0];
                $affymeProbesetUpdateStmt->execute($rank, $sampleId, $geneId)
                    or die $affymeProbesetUpdateStmt->errstr;
            } else {
                my $query = $rankUpdateStart.'IN (';
                for ( my $i = 0; $i < $geneCount; $i++ ){
                    if ( $i > 0 ){
                        $query .= ', ';
                    }
                    $query .= '?';
                }
                $query .= ')';
                my $affyRankMultiUpdateStmt = $dbh_thread->prepare($query);
                $affyRankMultiUpdateStmt->execute($rank, $sampleId, @geneIds_arr)
                    or die $affyRankMultiUpdateStmt->errstr;
            }
        }

        #print status
        printf("Chip: %s - PID: %s - %d/%d\n", $sampleId, $pid, $k+1, $batchLength);
    }

    $dbh_thread->commit() or die('Failed commit');
    $dbh_thread->disconnect();
}

## Clean potentially already computed ranks
#my $cleanProbeset = $dbh->prepare("UPDATE affymetrixProbeset SET rank = null");
#my $cleanChip     = $dbh->prepare("UPDATE affymetrixChip SET chipMaxRank = null,
#                                                             chipDistinctRankCount = null");
#my $cleanChipType = $dbh->prepare("UPDATE chipType SET chipTypeMaxRank = null");
#printf("Cleaning existing data: ");
#$cleanProbeset->execute() or die $cleanProbeset->errstr;
#$cleanChip->execute() or die $cleanProbeset->errstr;
#$cleanChipType->execute() or die $cleanProbeset->errstr;
#if ($auto == 0) {
#    $dbh->commit() or die("Failed commit");
#}
#printf("Done\n");

# Queries to compute gene ranks per chip. We rank over all not-excluded genes, not only over expressed genes
my $sqlChip = "SELECT t1.bgeeAffymetrixChipId FROM affymetrixChip AS t1
               WHERE EXISTS (SELECT 1 FROM affymetrixProbeset AS t2
                   WHERE t2.reasonForExclusion = '$Utils::CALL_NOT_EXCLUDED'
                   AND t2.bgeeAffymetrixChipId = t1.bgeeAffymetrixChipId
               )";
if ($sample_count > 0) {
    $sqlChip .= ' ORDER BY t1.bgeeAffymetrixChipId
                  LIMIT '.$sample_offset.', '.$sample_count;
}
my $affymChipStmt = $dbh->prepare($sqlChip);
$affymChipStmt->execute()  or die $affymChipStmt->errstr;

my @libs = map { $_->[0] } @{$affymChipStmt->fetchall_arrayref};
my $l = @libs;

printf("Found %d Affymetrix chips\n", $l);
# We are going to compute/store ranks using parallelization,
# each child process will be responsible to compute/update ranks for a batch of chips
my $batchSize = 100;
my $iterationCount = ceil($l/$batchSize);
my $parallel = $parallel_jobs;
if ($iterationCount < $parallel_jobs) {
    $parallel = $iterationCount;
}

print("Rank computation per chip...\n");
# Disconnect the DBI connection open in parent process, otherwise it will generate errors
# (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
$dbh->disconnect();
my $pm = new Parallel::ForkManager($parallel);
while ( my @next_libs = splice(@libs, 0, $batchSize) ) {
    # Forks and returns the pid for the child
    my $pid = $pm->start and next;
    print("\nStart batch of $batchSize chips, process ID $pid...\n");
    compute_update_rank_batch(\@next_libs, $pid);
    print("\nDone batch of $batchSize chips, process ID $pid.\n");
    $pm->finish;
}
$pm->wait_all_children;

print("Rank computation per chip done\n");

# ##############
# Store max rank and number of distinct ranks per chip, and max rank per chip type
# NOTE: as of Bgee 15, since we can run this script on the cluster for parallelization,
# these updates should be done afterwards (for intance, in the makefile)
# after launching all the jobs
#$dbh = Utils::connect_bgee_db($bgee_connector);
#$dbh->{'AutoCommit'} = 1;
#my $sql =
#"UPDATE affymetrixChip AS t0
# INNER JOIN (
#     SELECT t1.bgeeAffymetrixChipId, MAX(t1.rank) AS maxRank, COUNT(DISTINCT t1.rank) AS distinctRankCount
#     FROM affymetrixProbeset AS t1
#     WHERE t1.reasonForExclusion NOT IN ('$Utils::EXCLUDED_FOR_PRE_FILTERED', '$Utils::EXCLUDED_FOR_UNDEFINED')
#     GROUP BY t1.bgeeAffymetrixChipId
# ) AS ranks ON t0.bgeeAffymetrixChipId = ranks.bgeeAffymetrixChipId
# SET t0.chipMaxRank = ranks.maxRank, t0.chipDistinctRankCount = ranks.distinctRankCount
# WHERE EXISTS (
#     SELECT 1 FROM affymetrixProbeset AS t2
#     WHERE t2.expressionId IS NOT NULL AND t2.bgeeAffymetrixChipId = t0.bgeeAffymetrixChipId
# )";
#
#my $t0 = time();
#printf('Inserting max ranks and distinct rank counts in affymetrixChip table: ');
#my $maxRankChipStmt = $dbh->prepare($sql);
#$maxRankChipStmt->execute() or die $maxRankChipStmt->errstr;
#printf("Done in %.2fs\n", (time() - $t0));
#
#
## Store max ranks per chip type over all conditions,
## for later normalization
#my $affyChipTypeMaxStmt = $dbh->prepare(
#"UPDATE chipType AS t1
# INNER JOIN (
#     SELECT chipTypeId, MAX(chipMaxRank) AS maxRank
#     FROM affymetrixChip
#     GROUP BY chipTypeId
# ) AS chipMaxRanks ON t1.chipTypeId = chipMaxRanks.chipTypeId
# SET t1.chipTypeMaxRank = chipMaxRanks.maxRank");
#
#$t0 = time();
#printf('Inserting max ranks per chip type: ');
#$affyChipTypeMaxStmt->execute() or die $affyChipTypeMaxStmt->errstr;
#printf("Done in %.2fs\n", (time() - $t0));

exit 0;