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
my ($samples_per_job)  = (100); # default 100 libraries per thread
my (@lib_ids) = ();
my ($sample_offset) = (0);
my ($sample_count) = (0);
my ($all_datatypes) = (0);
my ($bulk) = (0);
my ($full_length) = (0);
my ($droplet) = (0);
my %opts = ('bgee=s'           => \$bgee_connector, # Bgee connector string
            'number_threads=i' => \$number_threads,
            'samples_per_job=i'   => \$samples_per_job,
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
\t-samples_per_job  Number of annotated samples per thread (optional)
\t-lib_ids          a comma-separated list of library IDs to treat, instead of providing sample_offset and sample_count (optional)
\t-sample_offset    The offset parameter to retrieve libraries to compute ranks for
\t-sample_count     The row_count parameter to retrieve libraries to compute ranks for
\t-all_datatypes    (optional) compute ranks for all datatypes (bulk RNA-Seq, full-length single cell RNA-Seq and droplet-based single cell RNA-Seq libraries)
\t-bulk             (optional) compute ranks for bulk RNA-Seq libraries
\t-full_length      (optional) compute ranks for full-length single cell RNA-Seq libraries
\t-droplet          (optional) compute ranks for droplet-based single cell RNA-Seq libraries
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
if ($number_threads <= 0 || $samples_per_job <= 0) {
    die("Invalid argument number_threads/samples_per_job\n");
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
# COMPUTE RANKS PER ANNOTATED SAMPLE         #
##############################################

# Batch size for INSERT into temp ranking table, to avoid huge SQL strings
# and stay within max_allowed_packet limits (~50K genes/sample average)
my $INSERT_BATCH_SIZE = 1000;

sub compute_update_rank_annotated_sample_batch {
    my ($annotatedSampleBatchRef) = @_;
    my $batchLength = scalar @{ $annotatedSampleBatchRef };
    return unless $batchLength;

    # Open ONE connection for the entire batch — not per sample.
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
    my $dbh_thread = Utils::connect_bgee_db($bgee_connector);

    # Create the MEMORY temp table once per connection; TRUNCATE before each sample.
    # MEMORY engine avoids disk I/O for the ranking staging table.
    # The PRIMARY KEY on (rnaSeqLibraryAnnotatedSampleId, bgeeGeneId) matches the PK of
    # rnaSeqLibraryAnnotatedSampleGeneResult, so the UPDATE JOIN uses index lookups.
    $dbh_thread->do('CREATE TEMPORARY TABLE IF NOT EXISTS rnaSeqAnnotSampleRanking (
        bgeeGeneId                    mediumint unsigned NOT NULL,
        rnaSeqLibraryAnnotatedSampleId mediumint unsigned NOT NULL,
        rawRank                        decimal(9,2) unsigned NOT NULL,
        PRIMARY KEY (rnaSeqLibraryAnnotatedSampleId, bgeeGeneId)
    ) ENGINE=MEMORY') or die $dbh_thread->errstr;

    # Prepare statements once, reuse across all samples in this batch.
    my $rnaSeqResultsStmt = $dbh_thread->prepare(
        'SELECT DISTINCT t1.bgeeGeneId, t1.abundance
         FROM rnaSeqLibraryAnnotatedSampleGeneResult AS t1
         WHERE t1.rnaSeqLibraryAnnotatedSampleId = ?
         AND t1.expressionId IS NOT NULL
         AND t1.reasonForExclusion = ?
         ORDER BY t1.abundance DESC');

    my $updateRnaSeqLibsStmt = $dbh_thread->prepare(
        'UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1
         INNER JOIN rnaSeqAnnotSampleRanking AS t2
             ON t1.bgeeGeneId = t2.bgeeGeneId
             AND t1.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId
             AND t1.expressionId IS NOT NULL
         SET t1.rawRank = t2.rawRank');

    for my $k ( 0..$batchLength-1 ) {
        my $rnaSeqLibraryAnnotatedSampleId = ${$annotatedSampleBatchRef}[$k];

        # TRUNCATE resets the MEMORY table without DDL overhead (no CREATE/DROP per sample)
        $dbh_thread->do('TRUNCATE TABLE rnaSeqAnnotSampleRanking') or die $dbh_thread->errstr;

        Utils::start_transaction($dbh_thread);

        # ======= Compute ranks ========
        #reason for exclusion should always be $Utils::CALL_NOT_EXCLUDED for genes with an expressionId, but we add this condition to be more explicit
        $rnaSeqResultsStmt->execute($rnaSeqLibraryAnnotatedSampleId, $Utils::CALL_NOT_EXCLUDED)
            or die $rnaSeqResultsStmt->errstr;

        my @results = map { {'id' => $_->[0], 'val' => $_->[1]} } @{$rnaSeqResultsStmt->fetchall_arrayref};
        # some annotatedSamples have no geneResults, check if @results is non-empty
        if (@results) {
            my %sorted = Utils::fractionnal_ranking(@results);
            # we get ranks as keys, with reference to an array of gene IDs with that rank as value
            my %reverseHash = Utils::revhash(%sorted);

            # ======= Insert ranks into temp table in batches to stay within max_allowed_packet ========
            my @all_rows;
            for my $rank ( keys %reverseHash ) {
                for my $geneId ( @{$reverseHash{$rank}} ) {
                    push @all_rows, [$geneId, $rnaSeqLibraryAnnotatedSampleId, $rank];
                }
            }
            while (my @batch = splice(@all_rows, 0, $INSERT_BATCH_SIZE)) {
                my $placeholders = join(', ', ('(?, ?, ?)') x scalar @batch);
                my $insertStmt = $dbh_thread->prepare(
                    "INSERT INTO rnaSeqAnnotSampleRanking ".
                    "(bgeeGeneId, rnaSeqLibraryAnnotatedSampleId, rawRank) VALUES $placeholders");
                $insertStmt->execute(map { @$_ } @batch) or die $insertStmt->errstr;
            }

            # ======= Update ranks via PK-scoped JOIN ========
            $updateRnaSeqLibsStmt->execute() or die $updateRnaSeqLibsStmt->errstr;
        }
        # ======= Commit ========
        $dbh_thread->commit() or die('Failed commit');

        #print status
        printf("AnnotatedSampleId: %s - PID: %s - %d/%d\n", $rnaSeqLibraryAnnotatedSampleId, $$, $k+1, $batchLength);
    }
    $dbh_thread->disconnect();
}

# Retrieve annotated sample IDs that have expressionId set but no rawRank yet,
# for the given datatype. Operates directly at the annotated sample level —
# this is the actual unit of work for rank computation.
#
# The EXISTS/NOT EXISTS subqueries use the PK (rnaSeqLibraryAnnotatedSampleId, bgeeGeneId)
# of rnaSeqLibraryAnnotatedSampleGeneResult: MySQL stops at the first matching row per sample,
# so each check is a fast PK range scan rather than a full-table scan.
sub get_annotated_samples_to_process {
    my ($single_cell, $sample_multiplexing) = @_;

    #XXX We do not rank genes with reasonForExclusion = "biotype not targeted" because they have no expressionId.
    #    These genes are consistent across libraries of the same species and population captured.
    #TODO: Investigate the impact of not always using the same set of genes considering that the score is then normalized
    #      by the max rank for the population captured, species, isSingleCell and isSampleMultiplexing.

    # Step 1: quickly retrieve all candidate annotated sample IDs from the small
    # rnaSeqLibraryAnnotatedSample table (~119K rows for Bgee 16.0), without touching the huge
    # rnaSeqLibraryAnnotatedSampleGeneResult table.
    my $dbh = Utils::connect_bgee_db($bgee_connector);
    my $step1_sql = 'SELECT las.rnaSeqLibraryAnnotatedSampleId'.
                   ' FROM rnaSeqLibraryAnnotatedSample AS las'.
                   ' INNER JOIN rnaSeqLibrary AS lib ON lib.rnaSeqLibraryId = las.rnaSeqLibraryId'.
                   ' WHERE lib.rnaSeqTechnologyIsSingleCell = ?'.
                   ' AND lib.sampleMultiplexing = ?';
    if ($sample_count > 0) {
        $step1_sql .= ' ORDER BY las.rnaSeqLibraryAnnotatedSampleId'.
                     ' LIMIT '.$sample_offset.', '.$sample_count;
    }
    my $stmt = $dbh->prepare($step1_sql);
    $stmt->execute($single_cell, $sample_multiplexing) or die $stmt->errstr;
    my @all_candidates = map { $_->[0] } @{$stmt->fetchall_arrayref};
    # Disconnect the DBI connection open in parent process, otherwise it will generate errors
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
    $dbh->disconnect();
    return () unless @all_candidates;

    printf("Found %d candidate annotated samples, filtering in parallel with %d threads...\n",
        scalar @all_candidates, $number_threads);

    # Step 2: in parallel, for each candidate run two fast PK-prefix range scans on
    # rnaSeqLibraryAnnotatedSampleGeneResult (LIMIT 1 stops at the first matching row):
    #   - EXISTS: at least one row with expressionId set
    #   - NOT EXISTS: no row with rawRank set yet
    # Each check uses the PK (rnaSeqLibraryAnnotatedSampleId, bgeeGeneId) so it
    # reads at most a handful of pages regardless of the 6B row count.
    my $exists_sql   = 'SELECT 1 FROM rnaSeqLibraryAnnotatedSampleGeneResult'.
                       ' WHERE rnaSeqLibraryAnnotatedSampleId = ?'.
                       ' AND expressionId IS NOT NULL LIMIT 1';
    my $norank_sql   = 'SELECT 1 FROM rnaSeqLibraryAnnotatedSampleGeneResult'.
                       ' WHERE rnaSeqLibraryAnnotatedSampleId = ?'.
                       ' AND rawRank IS NOT NULL LIMIT 1';

    my @filtered_samples;
    my $pm_filter = new Parallel::ForkManager($number_threads);
    $pm_filter->run_on_finish(sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
        push @filtered_samples, @{$data_ref} if defined $data_ref;
    });

    my $batch_size = ceil(scalar(@all_candidates) / $number_threads);
    $batch_size = 1 if $batch_size < 1;
    while (my @batch = splice(@all_candidates, 0, $batch_size)) {
        my $pid = $pm_filter->start and next;
        my $dbh_thread = Utils::connect_bgee_db($bgee_connector);
        my $existsStmt = $dbh_thread->prepare($exists_sql);
        my $norankStmt = $dbh_thread->prepare($norank_sql);
        my @valid;
        for my $sample_id (@batch) {
            $existsStmt->execute($sample_id) or die $existsStmt->errstr;
            next unless defined $existsStmt->fetchrow_arrayref;  # no expressionId → skip
            $norankStmt->execute($sample_id) or die $norankStmt->errstr;
            next if defined $norankStmt->fetchrow_arrayref;      # already ranked → skip
            push @valid, $sample_id;
        }
        $existsStmt->finish;
        $norankStmt->finish;
        $dbh_thread->disconnect();
        $pm_filter->finish(0, \@valid);
    }
    $pm_filter->wait_all_children;

    return @filtered_samples;
}

# Collect all annotated samples to process, across requested datatypes.
# lib_ids is kept for backward compatibility: if provided, treat them as annotated sample IDs.
my @samples = ();
if (!@lib_ids) {
    if ($all_datatypes || $bulk) {
        push(@samples, get_annotated_samples_to_process(0, 0));
    }
    if ($all_datatypes || $full_length) {
        push(@samples, get_annotated_samples_to_process(1, 0));
    }
    if ($all_datatypes || $droplet) {
        push(@samples, get_annotated_samples_to_process(1, 1));
    }
} else {
    print "Using provided library IDs to retrieve annotated sample IDs to process: ".join(', ', @samples)."\n";
    my $dbh = Utils::connect_bgee_db($bgee_connector);
    my $sql = 'SELECT las.rnaSeqLibraryAnnotatedSampleId'.
                ' FROM rnaSeqLibraryAnnotatedSample AS las'.
                ' WHERE las.rnaSeqLibraryAnnotatedSampleId IN ('.join(', ', ('?') x scalar @lib_ids).')';
    my $stmt = $dbh->prepare($sql);
    $stmt->execute(@lib_ids) or die $stmt->errstr;
    @samples = map { $_->[0] } @{$stmt->fetchall_arrayref};
    $dbh->disconnect();
}

my $number_samples = scalar @samples;
if (!$number_samples) {
    print "Nothing to be done, exiting.\n";
    exit 0;
}
if ($number_samples < $samples_per_job) {
    $samples_per_job = $number_samples;
}
printf("Found %d annotated samples to process\n", $number_samples);

# Distribute samples across threads; each child handles a batch of $samples_per_job samples
# using a single DB connection and a reusable MEMORY temp table.
my $iterationCount = ceil($number_samples/$samples_per_job);
my $parallel = $number_threads;
if ($iterationCount < $number_threads) {
    $parallel = $iterationCount;
}

# To save cores needed, we run the computation in the main thread when only one job is requested.
# And we do that in different loops to not need to use ForManager if not needed.
if ($parallel == 1) {
    while ( my @next_samples = splice(@samples, 0, $samples_per_job) ) {
        print("\nStart batch of $samples_per_job annotated samples, process ID $$...\n");
        compute_update_rank_annotated_sample_batch(\@next_samples);
        print("\nDone batch of $samples_per_job annotated samples, process ID $$.\n");
    }
} else {
    print("Rank computations with $parallel threads and $samples_per_job annotated samples per thread...\n");
    my $pm = new Parallel::ForkManager($parallel);
    while ( my @next_samples = splice(@samples, 0, $samples_per_job) ) {
        my $pid = $pm->start and next;
        # See https://stackoverflow.com/a/1673011/1768736 about PID
        print("\nStart batch of $samples_per_job annotated samples, process ID $$...\n");
        compute_update_rank_annotated_sample_batch(\@next_samples);
        print("\nDone batch of $samples_per_job annotated samples, process ID $$.\n");
        $pm->finish;
    }
    $pm->wait_all_children;
}

print("Rank computation per annotated sample done\n");

exit 0;
