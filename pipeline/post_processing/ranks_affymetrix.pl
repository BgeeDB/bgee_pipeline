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
my ($chips_per_job)  = (100); # default 100 chips per thread
my (@chip_ids) = ();
my ($sample_offset) = (0);
my ($sample_count) = (0);
my %opts = ('bgee=s'          => \$bgee_connector, # Bgee connector string
            'parallel_jobs=i' => \$parallel_jobs,
            'chips_per_job=i' => \$chips_per_job,
            'chip_ids=s'      => \@chip_ids,
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
\t-chips_per_job    Number of chips per thread
\t-chip_ids         A comma-separated list of bgee internal chip IDs to treat, instead of providing sample_offset and sample_count
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
if ($parallel_jobs <= 0 || $chips_per_job <= 0) {
    die("Invalid argument parallel_jobs/libs_per_job\n");
}

@chip_ids = split(/,/, join(',', @chip_ids));
if ($sample_offset > 0 && @chip_ids) {
    die("Not possible to provide chip IDs and offset parameters at the same time\n");
}

# NOTE: this script assumes that we don't use in Bgee chip types allowing to hybridize samples
# from two different species at the same time.
#
# Reasonning of the computations:
# 1) compute fractional ranks of genes in table affymetrixProbeset, for each affymetrix chip,
# based on the highest signal intensity from all the probesets mapped to a given gene.


##############################################
# COMPUTE RANKS PER CHIP                     #
##############################################
sub compute_update_rank_batch {
    my ($batchRef) = @_;
    my $batchLength = scalar @{ $batchRef };

    for my $k ( 0..$batchLength-1 ) {
        my $sampleId = ${$batchRef}[$k];

        # Connection to database in the parent process must have been closed
        # before calling this sub, otherwise it will generate errors
        # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289).
        # Get a new database connection for each thread.
        my $dbh_thread = Utils::connect_bgee_db($bgee_connector);
        Utils::start_transaction($dbh_thread);

        # ======= Compute ranks ========
        my $affymProbesetStmt = $dbh_thread->prepare(
                "SELECT bgeeGeneId, MAX(normalizedSignalIntensity) AS maxIntensity
                 FROM affymetrixProbeset
                 WHERE bgeeAffymetrixChipId = ?
                 AND expressionId IS NOT NULL
                 GROUP BY bgeeGeneId
                 ORDER BY maxIntensity DESC");
        $affymProbesetStmt->execute($sampleId) or die $affymProbesetStmt->errstr;

        my @results = map { {'id' => $_->[0], 'val' => $_->[1]} } @{$affymProbesetStmt->fetchall_arrayref};
        my %sorted = Utils::fractionnal_ranking(@results);
        # we get ranks as keys, with reference to an array of gene IDs with that rank as value
        my %reverseHash = Utils::revhash(%sorted);

        # ======= Update ranks ========
        my $affyTmpTableStmt = $dbh_thread->prepare(
                'CREATE TEMPORARY TABLE affyChipRanking
                 SELECT bgeeGeneId, bgeeAffymetrixChipId, rank
                 FROM affymetrixProbeset
                 LIMIT 0');
        $affyTmpTableStmt->execute() or die $affyTmpTableStmt->errstr;
        my $sqlInsertTmpRanks = 'INSERT INTO affyChipRanking (bgeeGeneId, bgeeAffymetrixChipId, rank)
                VALUES ';
        my @insertTmpRanksValues = ();
        my $valueCount = 0;
        # We use this reverseHash because we used to use UPDATE statements with IN(...) clause
        for my $rank ( keys %reverseHash ){
            my $geneIds_arrRef = $reverseHash{$rank};
            my @geneIds_arr = @$geneIds_arrRef;
            my $geneCount = scalar @geneIds_arr;
            for ( my $i = 0; $i < $geneCount; $i++ ){
                if ($valueCount > 0) {
                    $sqlInsertTmpRanks .= ', ';
                }
                $sqlInsertTmpRanks .= '(?, ?, ?)';
                push (@insertTmpRanksValues, ($geneIds_arr[$i], $sampleId, $rank));
                $valueCount++;
            }
        }
        my $insertTmpRanksStmt = $dbh_thread->prepare($sqlInsertTmpRanks);
        $insertTmpRanksStmt->execute(@insertTmpRanksValues) or die $insertTmpRanksStmt->errstr;

        my $updateAffyChipsStmt = $dbh_thread->prepare('UPDATE affymetrixProbeset AS t1
                INNER JOIN affyChipRanking AS t2
                    ON t1.bgeeGeneId = t2.bgeeGeneId AND t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId '.
                    # To be able to use the index on affymetrixProbeset, plus this constraint is correct
                   'AND t1.expressionId IS NOT NULL
                SET t1.rank = t2.rank');
        $updateAffyChipsStmt->execute() or die $updateAffyChipsStmt->errstr;

        my $dropAffyTmpTableStmt = $dbh_thread->prepare('DROP TABLE affyChipRanking');
        $dropAffyTmpTableStmt->execute() or die $dropAffyTmpTableStmt->errstr;

        # ======= Commit ========
        $dbh_thread->commit() or die('Failed commit');
        $dbh_thread->disconnect();

        #print status
        printf("Chip: %s - PID: %s - %d/%d\n", $sampleId, $$, $k+1, $batchLength);
    }
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

my @chips = @chip_ids;
if (!@chip_ids) {
    my $dbh = Utils::connect_bgee_db($bgee_connector);
    my $sqlChip = "SELECT t1.bgeeAffymetrixChipId FROM affymetrixChip AS t1
                   WHERE EXISTS (SELECT 1 FROM affymetrixProbeset AS t2
                       WHERE t2.expressionId IS NOT NULL
                       AND t2.bgeeAffymetrixChipId = t1.bgeeAffymetrixChipId
                   ) AND NOT EXISTS (SELECT 1 FROM affymetrixProbeset AS t2
                      WHERE t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                      AND t2.rank IS NOT NULL)";
    if ($sample_count > 0) {
        $sqlChip .= ' ORDER BY t1.bgeeAffymetrixChipId
                      LIMIT '.$sample_offset.', '.$sample_count;
    }
    my $affymChipStmt = $dbh->prepare($sqlChip);
    $affymChipStmt->execute()  or die $affymChipStmt->errstr;

    @chips = map { $_->[0] } @{$affymChipStmt->fetchall_arrayref};

    # Disconnect the DBI connection open in parent process, otherwise it will generate errors
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
    $dbh->disconnect();
}
my $l = @chips;
if (!$l) {
	print "Nothing to be done, exiting.\n";
	exit 0;
}
if ($l < $chips_per_job) {
    $chips_per_job = $l;
}

printf("Found %d Affymetrix chips\n", $l);
# We are going to compute/store ranks using parallelization,
# each child process will be responsible to compute/update ranks for a batch of chips
my $iterationCount = ceil($l/$chips_per_job);
my $parallel = $parallel_jobs;
if ($iterationCount < $parallel_jobs) {
    $parallel = $iterationCount;
}

# To save cores needed, we run the computation in the main thread when only one job is requested.
# And we do that in different loops to not need to use ForManager if not needed.
if ($parallel == 1) {
    while ( my @next_chips = splice(@chips, 0, $chips_per_job) ) {
        print("\nStart batch of $chips_per_job chips, process ID $$...\n");
        compute_update_rank_batch(\@next_chips);
        print("\nDone batch of $chips_per_job chips, process ID $$.\n");
    }
} else {
    print("Rank computation per chip with $parallel threads and $chips_per_job chips per thread......\n");
    my $pm = new Parallel::ForkManager($parallel);
    while ( my @next_chips = splice(@chips, 0, $chips_per_job) ) {
        # Forks and returns the pid for the child
        # See https://stackoverflow.com/a/1673011/1768736 about PID
        my $pid = $pm->start and next;
        print("\nStart batch of $chips_per_job chips, process ID $$...\n");
        compute_update_rank_batch(\@next_chips);
        print("\nDone batch of $chips_per_job chips, process ID $$.\n");
        $pm->finish;
    }
    $pm->wait_all_children;
}

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