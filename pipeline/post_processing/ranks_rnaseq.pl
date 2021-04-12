#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw( time );
use diagnostics;

# Updates the ranks in rnaSeqResult table and RNA-Seq mean rank in the expression table.
# Philippe Moret, created Oct 2015.
# Frederic Bastian, updated June 2016.
# Frederic Bastian, updated Feb. 2017: adapt to new conditions and new schema in Bgee 14
# Frederic Bastian, updated Jan. 2020: parallelize rank computations
# Frederic Bastian, updated Apr. 2021: this script is now responsible only for computing
# ranks for RNA-Seq library; improve parallelization possibilities.

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
my %opts = ('bgee=s'        => \$bgee_connector, # Bgee connector string
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
\t-sample_offset    The offset parameter to retrieve libraries to compute ranks for
\t-sample_count     The row_count parameter to retrieve libraries to compute ranks for
\n";
    exit 1;
}
if ($sample_offset < 0 || $sample_count < 0) {
    die('sample_offset and sample_count cannot be negative');
}
if ($sample_offset > 0 && $sample_count == 0) {
    die('sample_count must be provided if sample_offset is provided');
}

# Reasonning of the computations:
# 1) identify the valid set of genes that should be considered for ranking in all libraries:
# the set of all genes that received at least one read over all libraries in Bgee.
# 2) compute gene fractional ranks in table rnaSeqResult, for each RNA-Seq library, based on TPM values
# 3) update expression table: compute weighted mean of ranks per gene and condition,
# weighted by the number of distinct ranks in each library: we assume that libraries with higher
# number of distinct ranks have a higher power for ranking genes.
# Note that we do not "normalize" ranks between samples before computing the mean, as for Affymetrix data:
# all libraries are used to produce ranking over always the same set of genes in a given species,
# so the genomic coverage is always the same, and no "normalization" is required. The higher power
# at ranking genes of a library (for instance, thanks to a higher number of mapped reads)
# is taken into account by weighting the mean by the number of distinct ranks in the library,
# not by "normalizing" away libraries with lower max ranks; this would penalize conditions
# with a lower number of expressed genes, and thus with more ex-aequo ranked genes, corresponding
# to genes receiving 0 read.
# 4) also, insert max ranks per mapped condition in condition table
# (will allow to normalize ranks between conditions and data types), and sum of distinct rank counts
# per gene and condition, in expression table (used to compute weigthed mean over all data types in a condition)



my $dbh = Utils::connect_bgee_db($bgee_connector);

# We set autocommit to 1 so that we can define transaction isolation level
# *before* starting the next transaction, see https://www.perlmonks.org/?node_id=1074673
$dbh->{'AutoCommit'} = 1;

##############################################
# IDENTIFY VALID GENES                       #
##############################################
# We rank all genes that have received at least one read in any condition.
# So we always rank the same set of gene in a given species over all libraries.
# We don't use a temp table to be able to close/open the connection for each condition parameter combination.
# So we have to drop the table at the end.
# NOTE: as of Bgee 15, since we can run this script on the cluster for parallelization,
# this table should be produced beforehand (for intance, by the makefile)
# before launching all the jobs
#my $dropValidGenesStmt = $dbh->prepare('DROP TABLE IF EXISTS rnaSeqValidGenes');
#my $validGenesStmt     = $dbh->prepare('CREATE TABLE rnaSeqValidGenes
#                                            (PRIMARY KEY(bgeeGeneId))
#                                            SELECT DISTINCT t1.bgeeGeneId
#                                            FROM rnaSeqResult AS t1
#                                            WHERE t1.reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"');
#
#printf('Identifying set of valid genes for ranking: ');
#$dropValidGenesStmt->execute()  or die $dropValidGenesStmt->errstr;
#$validGenesStmt->execute()  or die $validGenesStmt->errstr;
#printf("Done\n");


##############################################
# COMPUTE RANKS PER LIBRARY                  #
##############################################
sub compute_update_rank_lib_batch {
    my ($libBatchRef, $pid) = @_;
    my $batchLength = scalar @{ $libBatchRef };

    # Connection to database in the parent process must have been closed
    # before calling this sub, otherwise it will generate errors
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289).
    # Get a new database connection for each thread.
    my $dbh_thread = Utils::connect_bgee_db($bgee_connector);
    Utils::start_transaction($dbh_thread);

    for my $k ( 0..$batchLength-1 ) {

        # ======= Compute ranks ========
        my $rnaSeqLibraryId = ${$libBatchRef}[$k];
        my $rnaSeqResultsStmt = $dbh_thread->prepare(
            'SELECT DISTINCT t1.bgeeGeneId, t1.tpm
             FROM rnaSeqResult AS t1 '.
             # join to table rnaSeqValidGenes to force
             # the selection of valid genes
             'INNER JOIN rnaSeqValidGenes AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId
             WHERE t1.rnaSeqLibraryId = ?
             AND t1.reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"
             ORDER BY t1.tpm DESC');

        $rnaSeqResultsStmt->execute($rnaSeqLibraryId) or die $rnaSeqResultsStmt->errstr;

        my @results = map { {'id' => $_->[0], 'val' => $_->[1]} } @{$rnaSeqResultsStmt->fetchall_arrayref};
        my %sorted = Utils::fractionnal_ranking(@results);
        # we get ranks as keys, with reference to an array of gene IDs with that rank as value
        my %reverseHash = Utils::revhash(%sorted);


        # ======= Update ranks ========
        # if several genes at a same rank, we'll update them at once
        # with a 'bgeeGeneId IN (?,?, ...)' clause.
        # If only one gene at a given rank, updated with the prepared statement below.
        my $rankUpdateStart = 'UPDATE rnaSeqResult SET rank = ? WHERE rnaSeqLibraryId = ? and bgeeGeneId ';
        my $rnaSeqResultUpdateStmt = $dbh_thread->prepare($rankUpdateStart.'= ?');

        for my $rank ( keys %reverseHash ){
            my $geneIds_arrRef = $reverseHash{$rank};
            my @geneIds_arr = @$geneIds_arrRef;
            my $geneCount = scalar @geneIds_arr;
            # with the multi-update query, sometimes there are thousands of genes
            # with equal ranks and the query is extremely slow. Apparently it's better
            # to alway use the single-gene update query
            for ( my $i = 0; $i < $geneCount; $i++ ){
                $rnaSeqResultUpdateStmt->execute($rank, $rnaSeqLibraryId, $geneIds_arr[$i])
                    or die $rnaSeqResultUpdateStmt->errstr;
            }
#            if ( $geneCount == 1 ){
#                my $geneId = $geneIds_arr[0];
#                $rnaSeqResultUpdateStmt->execute($rank, $rnaSeqLibraryId, $geneId)
#                    or die $rnaSeqResultUpdateStmt->errstr;
#            } else {
#                my $query = $rankUpdateStart.'IN (';
#                for ( my $i = 0; $i < $geneCount; $i++ ){
#                    if ( $i > 0 ){
#                        $query .= ', ';
#                    }
#                    $query .= '?';
#                }
#                $query .= ')';
#                my $rnaSeqRankMultiUpdateStmt = $dbh_thread->prepare($query);
#                $rnaSeqRankMultiUpdateStmt->execute($rank, $rnaSeqLibraryId, @geneIds_arr)
#                    or die $rnaSeqRankMultiUpdateStmt->errstr;
#            }
        }

        #print status
        printf("Lib: %s - PID: %s - %d/%d\n", $rnaSeqLibraryId, $pid, $k+1, $batchLength);
    }

    $dbh_thread->commit() or die('Failed commit');
    $dbh_thread->disconnect();
}



# Clean potentially already computed ranks
#my $cleanRNASeq = $dbh->prepare("UPDATE rnaSeqResult SET rank = NULL");
#my $cleanLib    = $dbh->prepare("UPDATE rnaSeqLibrary SET libraryMaxRank = NULL,
#                                                              libraryDistinctRankCount = NULL");
#printf("Cleaning existing data: ");
#$cleanRNASeq->execute() or die $cleanRNASeq->errstr;
#$cleanLib->execute() or die $cleanLib->errstr;
#printf("Done\n");



# Queries to compute gene ranks per library.
# We rank all genes that have received at least one read in any condition.
# So we always rank the same set of gene in a given species over all libraries.
# We assume that each library maps to only one species through its contained genes.
my $libSql = 'SELECT t1.rnaSeqLibraryId FROM rnaSeqLibrary AS t1
              WHERE EXISTS (SELECT 1 FROM rnaSeqResult AS t2
                  WHERE t1.rnaSeqLibraryId = t2.rnaSeqLibraryId
                  AND t2.readsCount > 0
              )';
if ($sample_count > 0) {
    $libSql .= ' ORDER BY t1.rnaSeqLibraryId
                 LIMIT '.$sample_offset.', '.$sample_count;
}
my $rnaSeqLibStmt = $dbh->prepare($libSql);

my $t0 = time();
# Get the list of all rna-seq libraries
$rnaSeqLibStmt->execute()  or die $rnaSeqLibStmt->errstr;


my @libs = map { $_->[0] } @{$rnaSeqLibStmt->fetchall_arrayref};
my $l = @libs;
printf("Found %d libraries\n", $l);
# We are going to compute/store ranks using parallelization,
# each child process will be responsible to compute/update ranks for a batch of libraries
my $libBatchSize = 100;
my $iterationCount = ceil($l/$libBatchSize);
my $parallel = $parallel_jobs;
if ($iterationCount < $parallel_jobs) {
    $parallel = $iterationCount;
}

print("Rank computation per library...\n");
# Disconnect the DBI connection open in parent process, otherwise it will generate errors
# (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
$dbh->disconnect();
my $pm = new Parallel::ForkManager($parallel);
while ( my @next_libs = splice(@libs, 0, $libBatchSize) ) {
    # Forks and returns the pid for the child
    my $pid = $pm->start and next;
    print("\nStart batch of $libBatchSize libraries, process ID $pid...\n");
    compute_update_rank_lib_batch(\@next_libs, $pid);
    print("\nDone batch of $libBatchSize libraries, process ID $pid.\n");
    $pm->finish;
}
$pm->wait_all_children;

print("Rank computation per library done\n");


# ##############
# Store max rank and number of distinct ranks per library.
# NOTE: as of Bgee 15, since we can run this script on the cluster for parallelization,
# this table should be produced afterwards (for intance, by the makefile)
# after launching all the jobs
#$dbh = Utils::connect_bgee_db($bgee_connector);
#$dbh->{'AutoCommit'} = 1;
#
#my $sql =
#"UPDATE rnaSeqLibrary AS t0
# INNER JOIN (
#     SELECT t1.rnaSeqLibraryId, MAX(t1.rank) AS maxRank, COUNT(DISTINCT t1.rank) AS distinctRankCount
#     FROM rnaSeqResult AS t1
#     WHERE t1.reasonForExclusion NOT IN ('$Utils::EXCLUDED_FOR_PRE_FILTERED',
#         '$Utils::EXCLUDED_FOR_UNDEFINED', '$Utils::EXCLUDED_FOR_ABSENT_CALLS')
#     GROUP BY t1.rnaSeqLibraryId
# ) AS ranks ON t0.rnaSeqLibraryId = ranks.rnaSeqLibraryId
# SET t0.libraryMaxRank = ranks.maxRank, t0.libraryDistinctRankCount = ranks.distinctRankCount
# WHERE EXISTS (
#     SELECT 1 FROM rnaSeqResult AS t2
#     WHERE t2.expressionId IS NOT NULL AND t2.rnaSeqLibraryId = t0.rnaSeqLibraryId
# )";
#
#$t0 = time();
#printf('Inserting max ranks and distinct rank counts in rnaSeqLibrary table...');
#my $maxRankLibStmt = $dbh->prepare($sql);
#$maxRankLibStmt->execute()  or die $maxRankLibStmt->errstr;
#printf("Done in %.2fs\n", (time() - $t0));
#$dbh->disconnect();

exit 0;