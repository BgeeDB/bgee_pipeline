#!/usr/bin/env perl

## Julien W. Feb. 2026
## Script used to parallelize delete of gene results for a given RNA-Seq experiment, to speed up the process.
## It is not part of the pipeline but a helper script to be used when needed, e.g. when an insert of gene results for an experiment fails and needs to be re-run.

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use Parallel::ForkManager;
use Time::HiRes qw(time);

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

############################
# CONFIGURATION
############################

# Define arguments & their default value
my ($bgee_connector)        = ('');
my ($experiment_id)         = ('');
my ($parallel_jobs)         = (4);
my ($delete_libraries)      = (0);
my ($is_single_cell)        = (-1);
my ($is_sample_multiplexing) = (-1);
my %opts = (
    'bgee=s'                 => \$bgee_connector,          # Bgee connector string
    'experiment_id=s'        => \$experiment_id,           # RNA-Seq experiment ID
    'parallel_jobs=i'        => \$parallel_jobs,           # Number of parallel jobs (default: 4)
    'delete_libraries'       => \$delete_libraries,        # Also delete libraries
    'isSingleCell=i'         => \$is_single_cell,          # delete only single-cell RNA-Seq data
    'isSampleMultiplexing=i' => \$is_sample_multiplexing,  # delete sample multiplexing RNA-Seq data
);

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $experiment_id eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -experiment_id=SRP467631 -parallel_jobs=12
\t-bgee             Bgee connector string
\t-experiment_id    RNA-Seq experiment ID to delete
\t-parallel_jobs    Number of parallel jobs (default: 4)
\t-delete_libraries Also delete libraries (default: false)
\t-isSingleCell     1: Delete single-cell RNA-Seq data, 0: Delete bulk RNA-Seq data
\t-isSampleMultiplexing 1: Delete sample multiplexing RNA-Seq data, 0: Delete non-sample multiplexing RNA-Seq data
\n";
    exit 1;
}

if ($is_single_cell != -1 && $is_single_cell != 0 && $is_single_cell != 1) {
    print "\n\tError: Invalid value for -isSingleCell. Must be 0 (for bulk RNA-Seq) or 1 (for single-cell RNA-Seq).\n";
    exit 1;
}
if ($is_sample_multiplexing != -1 && $is_sample_multiplexing != 0 && $is_sample_multiplexing != 1) {
    print "\n\tError: Invalid value for -isSampleMultiplexing. Must be 0 (for non-sample multiplexing RNA-Seq) or 1 (for sample multiplexing RNA-Seq).\n";
    exit 1;
}

############################
# Database Connection
############################

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

############################
# Récupération des IDs
############################

print "Fetching IDs...\n";

my $query = "
    SELECT DISTINCT t1.rnaSeqLibraryId, t2.rnaSeqLibraryAnnotatedSampleId
    FROM rnaSeqLibrary AS t1
    INNER JOIN rnaSeqLibraryAnnotatedSample AS t2
        ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId
    WHERE t1.rnaSeqExperimentId = ?";

my @params = ($experiment_id);

if ($is_single_cell != -1) {
    $query .= " AND t1.rnaSeqTechnologyIsSingleCell = ?";
    push @params, $is_single_cell;
}

if ($is_sample_multiplexing != -1) {
    $query .= " AND t1.sampleMultiplexing = ?";
    push @params, $is_sample_multiplexing;
}

print "$query\n";
print "Parameters: " . join(", ", @params) . "\n";

my $sth_ids = $dbh->prepare($query);
$sth_ids->execute(@params);

my %libsToAnnotSampleIds;
while (my ($lib_id, $annot_sample_id) = $sth_ids->fetchrow_array) {
    push @{$libsToAnnotSampleIds{$lib_id}}, $annot_sample_id;
}

$sth_ids->finish;
$dbh->disconnect;

my $total_samples = 0;
for my $lib_id (keys %libsToAnnotSampleIds) {
    $total_samples += scalar(@{$libsToAnnotSampleIds{$lib_id}});
}
my $total_libs = scalar(keys %libsToAnnotSampleIds);
print "Total libraries: $total_libs, Total annotated samples: $total_samples\n";

############################
# Parallel deletion
############################

print "Starting delete with $parallel_jobs threads...\n";

my $start = time();
my $pm = Parallel::ForkManager->new($parallel_jobs);

for my $libId (keys %libsToAnnotSampleIds) {
    for my $annotatedSampleId (@{$libsToAnnotSampleIds{$libId}}) {
        # Fork a process for each annotated sample
        my $pid = $pm->start and next;

        # Create a new database connection for this child process
        my $dbh_child = Utils::connect_bgee_db($bgee_connector);

        # Delete gene results
        my $sth_del = $dbh_child->prepare("
            DELETE FROM rnaSeqLibraryAnnotatedSampleGeneResult
            WHERE rnaSeqLibraryAnnotatedSampleId = ?
        ");
        $sth_del->execute($annotatedSampleId);
        $sth_del->finish;

        # Delete annotated sample if requested
        if ($delete_libraries) {
            my $sth_del_sample = $dbh_child->prepare("
                DELETE FROM rnaSeqLibraryAnnotatedSample
                WHERE rnaSeqLibraryAnnotatedSampleId = ?
            ");
            $sth_del_sample->execute($annotatedSampleId);
            $sth_del_sample->finish;
        }

        $dbh_child->disconnect;

        # End this child process
        $pm->finish;
    }
}

# Wait for all child processes to complete
$pm->wait_all_children;

# Delete libraries after all samples are deleted
if ($delete_libraries) {
    print "Deleting libraries...\n";
    my $dbh_lib = Utils::connect_bgee_db($bgee_connector);
    for my $libId (keys %libsToAnnotSampleIds) {
        my $sth_del_lib = $dbh_lib->prepare("
            DELETE FROM rnaSeqLibrary
            WHERE rnaSeqLibraryId = ?
        ");
        $sth_del_lib->execute($libId);
        $sth_del_lib->finish;
    }
    $dbh_lib->disconnect;
}

my $end = time();
my $duration = $end - $start;

print "Finished in ".sprintf("%.2f", $duration/3600)." hours\n";