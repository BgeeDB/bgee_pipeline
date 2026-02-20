#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw( time );
use diagnostics;

# Updates the ranks in rnaSeqResult table.
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
my ($number_threads) = (8); # default 8 parallel threads used to compute ranks
my ($libs_per_job)  = (100); # default 100 libraries per thread
my (@lib_ids) = ();
my ($sample_offset) = (0);
my ($sample_count) = (0);
my ($all_datatypes) = (0);
my ($bulk) = (0);
my ($full_length) = (0);
my ($droplet) = (0);
my %opts = ('bgee=s'           => \$bgee_connector, # Bgee connector string
            'number_threads=i' => \$number_threads,
            'libs_per_job=i'   => \$libs_per_job,
            'lib_ids=s'        => \@lib_ids,
            'sample_offset=i'  => \$sample_offset,
            'sample_count=i'   => \$sample_count,
            'all_datatypes' => \$all_datatypes,
            'bulk' => \$bulk,
            'full_length' => \$full_length,
            'droplet' => \$droplet
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\t-number_threads   Number of threads used to compute ranks (optional)
\t-libs_per_job     Number of libraries per thread (optional)
\t-lib_ids          a comma-separated list of library IDs to treat, instead of providing sample_offset and sample_count (optional)
\t-sample_offset    The offset parameter to retrieve libraries to compute ranks for
\t-sample_count     The row_count parameter to retrieve libraries to compute ranks for
\t-all_datatypes    (optional) compute ranks for all datatypes (bulk RNA-Seq, full-length single cell RNA-Seq and droplet-based single cell RNA-Seq libraries)
\t-bulk             (optional) compute ranks for bulk RNA-Seq libraries
\t-full_length       (optional) compute ranks for full-length single cell RNA-Seq libraries
\t-droplet           (optional) compute ranks for droplet-based single cell RNA-Seq libraries
\n";
    exit 1;
}

# Validate number_threads
if ( $number_threads < 0 || $number_threads =~ /\D/ ){
    print "\tError: number_threads must be a non-negative integer\n";
    exit 1;
}
$number_threads = 1 if $number_threads == 0;

if ($sample_offset < 0 || $sample_count < 0) {
    die('sample_offset and sample_count cannot be negative');
}
if ($sample_offset > 0 && $sample_count == 0) {
    die('sample_count must be provided if sample_offset is provided');
}
if ($number_threads <= 0 || $libs_per_job <= 0) {
    die("Invalid argument number_threads/libs_per_job\n");
}

my $datatype_options_provided = 0;
$datatype_options_provided = 3  if ( $all_datatypes );
$datatype_options_provided++  if ( $bulk );
$datatype_options_provided++  if ( $full_length );
$datatype_options_provided++  if ( $droplet );
if ( $datatype_options_provided == 0 ){
    print "\tError: at least one datatype option must be provided among -all_datatypes, -bulk, -full_length, -droplet\n";
    exit 1;
}
if ( $datatype_options_provided > 3 ){
    print "\tError: if -all_datatypes is selected then  -bulk, -full_length, -droplet should not be selected\n";
    exit 1;
}

@lib_ids = split(/,/, join(',', @lib_ids));
if ($sample_offset > 0 && @lib_ids) {
    die("Not possible to provide library IDs and offset parameters at the same time\n");
}

if ($datatype_options_provided > 0 && @lib_ids) {
    die("not possible to provide datatype options and library IDs at the same time");
}
# Reasoning of the computations (as of Bgee 15.0 it's exactly the same reasoning for
# bulk RNA-Seq and single-cell RNA-Seq data):
# * compute gene fractional ranks in table rnaSeqLibraryAnnotatedSampleGeneResult, for each RNA-Seq library,
#   based on TPM values.
# * Store also number of distinct ranks per library



##############################################
# COMPUTE RANKS PER LIBRARY                  #
##############################################
sub compute_update_rank_lib_batch {
    my ($libBatchRef) = @_;
    my $batchLength = scalar @{ $libBatchRef };

    for my $k ( 0..$batchLength-1 ) {
        my $rnaSeqLibraryId = ${$libBatchRef}[$k];
        # Connection to database in the parent process must have been closed
        # before calling this sub, otherwise it will generate errors
        # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289).
        # Get a new database connection for each thread.
        my $dbh_thread = Utils::connect_bgee_db($bgee_connector);
        # ======= Retrieve annotated samples for this library ========
        my $rnaSeqAnnotSamplesStmt = $dbh_thread->prepare(
            'SELECT DISTINCT rnaSeqLibraryAnnotatedSampleId
             FROM rnaSeqLibraryAnnotatedSample
             WHERE rnaSeqLibraryId = ?');
        $rnaSeqAnnotSamplesStmt->execute($rnaSeqLibraryId) or die $rnaSeqAnnotSamplesStmt->errstr;

        my @results = map { $_->[0] } @{$rnaSeqAnnotSamplesStmt->fetchall_arrayref};
        # Disconnect the DBI connection open in this thread, otherwise it will generate errors
        # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
        $dbh_thread->disconnect();

        compute_update_rank_annotated_sample_batch(\@results);
    }
}
sub compute_update_rank_annotated_sample_batch {
    my ($annotatedSampleBatchRef) = @_;
    my $batchLength = scalar @{ $annotatedSampleBatchRef };

    for my $k ( 0..$batchLength-1 ) {
        my $rnaSeqLibraryAnnotatedSampleId = ${$annotatedSampleBatchRef}[$k];

        # Connection to database in the parent process must have been closed
        # before calling this sub, otherwise it will generate errors
        # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289).
        # Get a new database connection for each thread.
        my $dbh_thread = Utils::connect_bgee_db($bgee_connector);
        Utils::start_transaction($dbh_thread);

        # ======= Compute ranks ========
        my $rnaSeqResultsStmt = $dbh_thread->prepare(
            'SELECT DISTINCT t1.bgeeGeneId, t1.abundance
             FROM rnaSeqLibraryAnnotatedSampleGeneResult AS t1
             WHERE t1.rnaSeqLibraryAnnotatedSampleId = ?
             AND t1.expressionId IS NOT NULL
             AND t1.reasonOfExclusion = ?
             ORDER BY t1.abundance DESC');
        #reason for exclusion should always be Utils::$CALL_NOT_EXCLUDED for genes with an expressionId, but we add this condition to be more explicit
        $rnaSeqResultsStmt->execute($rnaSeqLibraryAnnotatedSampleId, Utils::$CALL_NOT_EXCLUDED) or die $rnaSeqResultsStmt->errstr;

        my @results = map { {'id' => $_->[0], 'val' => $_->[1]} } @{$rnaSeqResultsStmt->fetchall_arrayref};
        # some annotatedSamples have no geneResults, check if @results is non-empty
        if (@results) {
            my %sorted = Utils::fractionnal_ranking(@results);
            # we get ranks as keys, with reference to an array of gene IDs with that rank as value
            my %reverseHash = Utils::revhash(%sorted);

            # ======= Update ranks ========
            my $rnaSeqTmpTableStmt = $dbh_thread->prepare(
                    'CREATE TEMPORARY TABLE rnaSeqAnnotSampleRanking
                     SELECT bgeeGeneId, rnaSeqLibraryAnnotatedSampleId, rawRank
                     FROM rnaSeqLibraryAnnotatedSampleGeneResult
                     LIMIT 0');
            $rnaSeqTmpTableStmt->execute() or die $rnaSeqTmpTableStmt->errstr;
            my $sqlInsertTmpRanks = 'INSERT INTO rnaSeqAnnotSampleRanking (bgeeGeneId, rnaSeqLibraryAnnotatedSampleId, rawRank)
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
                    push (@insertTmpRanksValues, ($geneIds_arr[$i], $rnaSeqLibraryAnnotatedSampleId, $rank));
                    $valueCount++;
                }
            }
            my $insertTmpRanksStmt = $dbh_thread->prepare($sqlInsertTmpRanks);
            $insertTmpRanksStmt->execute(@insertTmpRanksValues) or die $insertTmpRanksStmt->errstr;

            my $updateRnaSeqLibsStmt = $dbh_thread->prepare('UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1
                    INNER JOIN rnaSeqAnnotSampleRanking AS t2
                        ON t1.bgeeGeneId = t2.bgeeGeneId AND t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId '.
                        # To be able to use the index on rnaSeqLibraryAnnotatedSampleGeneResult, plus this constraint is correct
                       'AND t1.expressionId IS NOT NULL
                    SET t1.rawRank = t2.rawRank');
            $updateRnaSeqLibsStmt->execute() or die $updateRnaSeqLibsStmt->errstr;

            my $dropRnaSeqTmpTableStmt = $dbh_thread->prepare('DROP TABLE rnaSeqAnnotSampleRanking');
            $dropRnaSeqTmpTableStmt->execute() or die $dropRnaSeqTmpTableStmt->errstr;
        }
        # ======= Commit ========
        $dbh_thread->commit() or die('Failed commit');
        $dbh_thread->disconnect();

        #print status
        printf("AnnotatedSampleId: %s - PID: %s - %d/%d\n", $rnaSeqLibraryAnnotatedSampleId, $$, $k+1, $batchLength);
    }
}

sub get_libraries_to_process {
    my ($single_cell, $sample_multiplexing) = @_;
    my $dbh = Utils::connect_bgee_db($bgee_connector);

    # Queries to compute gene ranks per library.
    # We rank all genes that have received at least one read in any condition.
    # So we always rank the same set of gene in a given species over all libraries.

    #XXX the comment above is WRONG. We do not rank genes that have received at least one read in any condition
    #    because genes with reasonOfExclusion = "biotype not targeted" does not have any expressionId. These genes are not ranked but are however
    #    the same for all libraries of a given species and a given population captured.
    #TODO: Investigate the impact of not always using the same set of genes considering that the score is then normalized by the max rank for the population captured, species, isSingleCell and isSampleMultiplexing (i.e a different max rank for bulk, droplet and full length depending on species and population captured by the library protocol).

    # We assume that each library maps to only one species through its contained genes.
    my $libSql = 'SELECT t1.rnaSeqLibraryId FROM rnaSeqLibrary AS t1 '.
              'WHERE rnaSeqTechnologyIsSingleCell =  ? '.
              'AND rnaSeqTechnologyIsSampleMultiplexing = ? '.
              'AND EXISTS (SELECT 1 FROM rnaSeqLibraryAnnotatedSample AS t2 '.
                  'INNER JOIN rnaSeqLibraryAnnotatedSampleGeneResult AS t3 '.
                  'ON t3.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId '.
                  'WHERE t1.rnaSeqLibraryId = t2.rnaSeqLibraryId '.
                  'AND t3.expressionId IS NOT NULL '.
              ') AND NOT EXISTS (SELECT 1 FROM rnaSeqLibraryAnnotatedSample AS t2 '.
                  'INNER JOIN rnaSeqLibraryAnnotatedSampleGeneResult AS t3 '.
                  'ON t3.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId '.
                  'WHERE t1.rnaSeqLibraryId = t2.rnaSeqLibraryId '.
                  'AND t3.rawRank IS NOT NULL)';
    if ($sample_count > 0) {
        $libSql .= ' ORDER BY t1.rnaSeqLibraryId '.
                    'LIMIT '.$sample_offset.', '.$sample_count;
    }
    my $rnaSeqLibStmt = $dbh->prepare($libSql);

    # Get the list of all rna-seq libraries
    $rnaSeqLibStmt->execute($single_cell, $sample_multiplexing)  or die $rnaSeqLibStmt->errstr;

    my @libs = map { $_->[0] } @{$rnaSeqLibStmt->fetchall_arrayref};
    # Disconnect the DBI connection open in parent process, otherwise it will generate errors
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
    $dbh->disconnect();
    return @libs;
}

my @libs = @lib_ids;
if (!@lib_ids) {
    if ($all_datatypes || $bulk) {
        ## add libraries to @libs
        push(@libs, get_libraries_to_process(0, 0));
    }
    if ($all_datatypes || $full_length) {
        ## add libraries to @libs
        push(@libs, get_libraries_to_process(1, 0));
    }
    if ($all_datatypes || $droplet) {
        ## add libraries to @libs
        push(@libs, get_libraries_to_process(1, 1));
    }
}

my $number_libs = @libs;
if (!$number_libs) {
    print "Nothing to be done, exiting.\n";
    exit 0;
}
if ($number_libs < $libs_per_job) {
    $libs_per_job = $number_libs;
}
printf("Found %d libraries\n", $number_libs);
# We are going to compute/store ranks using parallelization,
# each child process will be responsible to compute/update ranks for a batch of libraries
my $iterationCount = ceil($number_libs/$libs_per_job);
my $parallel = $number_threads;
if ($iterationCount < $number_threads) {
    $parallel = $iterationCount;
}

# To save cores needed, we run the computation in the main thread when only one job is requested.
# And we do that in different loops to not need to use ForManager if not needed.
if ($parallel == 1) {
    while ( my @next_libs = splice(@libs, 0, $libs_per_job) ) {
        print("\nStart batch of $libs_per_job libraries, process ID $$...\n");
        compute_update_rank_lib_batch(\@next_libs);
        print("\nDone batch of $libs_per_job libraries, process ID $$.\n");
    }
} else {
    print("Rank computations per library with $parallel threads and $libs_per_job libraries per thread...\n");
    my $pm = new Parallel::ForkManager($parallel);
    while ( my @next_libs = splice(@libs, 0, $libs_per_job) ) {
        # Forks and returns the pid for the child
        my $pid = $pm->start and next;
        # See https://stackoverflow.com/a/1673011/1768736 about PID
        print("\nStart batch of $libs_per_job libraries, process ID $$...\n");
        compute_update_rank_lib_batch(\@next_libs);
        print("\nDone batch of $libs_per_job libraries, process ID $$.\n");
        $pm->finish;
    }
    $pm->wait_all_children;
}

print("Rank computation per library done\n");

exit 0;
