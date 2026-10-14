#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Parallel::ForkManager;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/..";
use Utils;


# -----------------------------------------------------------------------------
# Script: fix_degenerate_rna_seq_pvalues.pl
# Author: Julien Wollbrett
# Created: Oct. 2026
#
# Description:
#   Written to fix pvalues before the release of Bgee 16.0. This script can be removed for future
#   releases as BgeeCall and the droplet pipeline will be updated to solve the issue. The file
#   listing number of reference intergenic sequences per species available at 
#   generated_files/post_processing/count_intergenic_per_species_fix_16.tsv should be deleted at
#   the same time.
#   
#   This script replaces the p-values of 0 of the 'rnaSeqLibraryAnnotatedSampleGeneResult' table
#   computed in the annotated samples where all the reference intergenic regions with reads have
#   the same abundance v: the standard deviation of their distribution is 0, so that no z-score
#   can be computed. These samples are identified by the standard deviation stored in
#   sdAbundanceReferenceIntergenicDistribution:
#   - 0, mostly full-length single-cell samples where no reference intergenic region has reads
#     (v = 0);
#   - 1 for droplet-based single-cell samples: the droplet calls store 2^(standard deviation of
#     the log2 abundances), so that a standard deviation of 0 is stored as 1. This happens when
#     each reference intergenic region with reads has a single UMI.
#
#   The genes more expressed than v were given a p-value of 0 with a NULL z-score, which is not
#   a valid p-value. It is replaced by the empirical p-value of the gene against the N reference
#   intergenic regions of its species, (r + 1) / (N + 1), r being the number of regions at least
#   as expressed as the gene: r = 0 for a gene more expressed than v, so that its p-value is
#   1 / (N + 1).
#
#   Only the p-values of 0 with a NULL z-score of these samples are updated:
#   - the other genes keep their p-value: 1 below v, and, for droplet samples, the proportion
#     of regions with reads K / N at v, close to their empirical p-value (K + 1) / (N + 1),
#     K being the number of regions at v, which is not stored in the database;
#   - the p-values of 0 with a finite z-score are lower than the precision of the column (1e-30)
#     and are kept.
#
#   The samples are processed in parallel, each process with its own connection: the primary key
#   of rnaSeqLibraryAnnotatedSampleGeneResult starts with the sample ID, so that a process only
#   reads and updates the rows of its own samples.
#
#   The script can be relaunched, after an interruption for instance: the p-values already
#   replaced are not 0 any more. The p-values of the 'expression' table are averaged from these
#   p-values: the conditions listed at the end of the run have to be recomputed with
#   update_rna_seq_expression.pl.
#
# Usage Example:
#   perl fix_degenerate_rna_seq_pvalues.pl -bgee=<BGEECMD> -intergenicCounts=<file.tsv> -numberThreads=25
#
#   Options:
#     -bgee              Bgee connector string (required)
#     -intergenicCounts  TSV file with 2 columns: species ID, number of reference intergenic
#                        regions of the species (required). The first line can be a header.
#     -numberThreads     Number of processes updating the samples in parallel (default: 1)
#     -speciesIds        Comma-separated IDs of the species to process (default: all the species
#                        with such samples)
#     -debug             Count the p-values to replace, do not update them
# -----------------------------------------------------------------------------


###############################################################################
# Argument parsing
###############################################################################

my $bgee_connector    = '';
my $intergenic_counts = '';
my $number_threads    = 1;
my $species_ids       = '';
my $debug             = 0;

my %opts = (
    'bgee=s'             => \$bgee_connector,
    'intergenicCounts=s' => \$intergenic_counts,
    'numberThreads=i'    => \$number_threads,
    'speciesIds=s'       => \$species_ids,
    'debug'              => \$debug,
);

my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $intergenic_counts eq '' ){
    print "\n\tInvalid or missing argument:
	e.g. $0  -bgee=\$(BGEECMD) -intergenicCounts=<file.tsv> -numberThreads=25
	-bgee              Bgee connector string (required)
	-intergenicCounts  TSV file with 2 columns: species ID, number of reference intergenic regions
	                   of the species (required). The first line can be a header.
	-numberThreads     Number of processes updating the samples in parallel (default 1)
	-speciesIds        Comma-separated IDs of the species to process (default: all the species
	                   with such samples)
	-debug             Count the p-values to replace, do not update them
\n";
    exit 1;
}
if ( $number_threads < 1 ){
    warn "numberThreads < 1; defaulting to 1\n";
    $number_threads = 1;
}
my %requestedSpecies;
for my $speciesId ( split(/,/, $species_ids) ){
    die "Invalid species ID in -speciesIds: $speciesId\n"  if ( $speciesId !~ /^\d+$/ );
    $requestedSpecies{$speciesId} = 1;
}

# Number of samples processed by a same process, for a connection not to be opened per sample
my $SAMPLES_PER_JOB = 50;

# -----------------------------------------------------------------------------
# Number of reference intergenic regions per species
# -----------------------------------------------------------------------------
my %intergenicCount;
open(my $countFile, '<', $intergenic_counts)  or die "Cannot open $intergenic_counts: $!\n";
while ( defined(my $line = <$countFile>) ){
    chomp $line;
    $line =~ s/\r$//;
    next  if ( $line =~ /^\s*$/ );
    my ($speciesId, $count) = split(/\t/, $line);
    # Header line
    next  if ( $. == 1 && (!defined $speciesId || $speciesId !~ /^\d+$/) );
    if ( !defined $count || $speciesId !~ /^\d+$/ || $count !~ /^\d+$/ || $count == 0 ){
        die "Invalid line $. in $intergenic_counts, a species ID and a number of regions are expected: $line\n";
    }
    die "Species $speciesId defined twice in $intergenic_counts\n"  if ( exists $intergenicCount{$speciesId} );
    $intergenicCount{$speciesId} = $count;
}
close $countFile;
die "No species in $intergenic_counts\n"  if ( !%intergenicCount );

# -----------------------------------------------------------------------------
# Samples to fix, and the p-value to give in each species
# -----------------------------------------------------------------------------
my $bgee = Utils::connect_bgee_db($bgee_connector);
$bgee->{RaiseError} = 1;

$| = 1;

print "Retrieving the samples whose reference intergenic regions with reads all have the same abundance...\n";
my $sqlSamples =
    'SELECT t1.rnaSeqLibraryAnnotatedSampleId, t3.speciesId, t3.exprMappedConditionId, '.
    't2.rnaSeqTechnologyIsSingleCell, t2.sampleMultiplexing '.
    'FROM rnaSeqLibraryAnnotatedSample AS t1 '.
    'INNER JOIN rnaSeqLibrary AS t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId '.
    'INNER JOIN cond AS t3 ON t1.conditionId = t3.conditionId '.
    'WHERE (t1.sdAbundanceReferenceIntergenicDistribution = 0 '.
    # Droplet calls store 2^(standard deviation of the log2 abundances)
    'OR (t2.rnaSeqTechnologyIsSingleCell = 1 AND t2.sampleMultiplexing = 1 '.
    'AND t1.sdAbundanceReferenceIntergenicDistribution = 1))';
if ( %requestedSpecies ){
    $sqlSamples .= ' AND t3.speciesId IN ('.join(',', sort { $a <=> $b } keys %requestedSpecies).')';
}
$sqlSamples .= ' ORDER BY t3.speciesId, t1.rnaSeqLibraryAnnotatedSampleId';
my $querySamples = $bgee->prepare($sqlSamples);
$querySamples->execute() or die $querySamples->errstr;
# speciesId -> number of samples
my %samplesPerSpecies;
# [sampleId, speciesId, exprMappedConditionId, datatype]
my @samples;
while ( my @data = $querySamples->fetchrow_array ){
    # Same datatype names as the options of update_rna_seq_expression.pl
    my $datatype = !$data[3] ? 'bulk' : ($data[4] ? 'droplet' : 'fullLength');
    push @samples, [$data[0], $data[1], $data[2], $datatype];
    $samplesPerSpecies{$data[1]}++;
}
$querySamples->finish;
print 'Done, ', scalar(@samples), " samples retrieved.\n";

my @missingSpecies = grep { !exists $intergenicCount{$_} } sort { $a <=> $b } keys %samplesPerSpecies;
if ( @missingSpecies ){
    die "No number of reference intergenic regions in $intergenic_counts for the species of these samples: ",
        join(', ', @missingSpecies), "\n";
}

my %pValuePerSpecies;
for my $speciesId ( sort { $a <=> $b } keys %samplesPerSpecies ){
    # Decimal division, for the value to keep the 30 decimals of the column: a division of
    # integers only keeps 4 decimals (div_precision_increment), 1/(N+1) would then be 0.
    ($pValuePerSpecies{$speciesId}) =
        $bgee->selectrow_array('SELECT CAST(1 AS DECIMAL(31, 30)) / '.($intergenicCount{$speciesId} + 1));
    print "Species $speciesId: $samplesPerSpecies{$speciesId} samples, N = $intergenicCount{$speciesId}, ",
          "p-value = $pValuePerSpecies{$speciesId}\n";
}
# The processes open their own connection: the one of the parent is not shared with them
$bgee->disconnect;

# -----------------------------------------------------------------------------
# Parallel replacement of the p-values, by jobs of $SAMPLES_PER_JOB samples
# -----------------------------------------------------------------------------
my %replacedPerSpecies = map { $_ => 0 } keys %samplesPerSpecies;
# datatype -> exprMappedConditionIds holding at least one replaced p-value
my %conditionsToRecompute;
my $processedSamples = 0;
my $failedJobs = 0;

my $pm = Parallel::ForkManager->new($number_threads);
$pm->run_on_finish(sub {
    my ($pid, $exitCode, $jobIndex, $exitSignal, $coreDump, $result) = @_;
    if ( $exitCode != 0 || $exitSignal != 0 || !defined $result ){
        warn "Job $jobIndex failed (exit code $exitCode, signal $exitSignal)\n";
        $failedJobs++;
        return;
    }
    for my $speciesId ( keys %{$result->{'replaced'}} ){
        $replacedPerSpecies{$speciesId} += $result->{'replaced'}->{$speciesId};
    }
    for my $datatype ( keys %{$result->{'conditions'}} ){
        $conditionsToRecompute{$datatype}{$_} = 1  for keys %{$result->{'conditions'}->{$datatype}};
    }
    my $previous = $processedSamples;
    $processedSamples += $result->{'samples'};
    if ( int($processedSamples / 1000) > int($previous / 1000) || $processedSamples == scalar(@samples) ){
        print "\t$processedSamples/", scalar(@samples), " samples processed\n";
    }
});

my $jobIndex = 0;
for ( my $start = 0; $start < scalar(@samples); $start += $SAMPLES_PER_JOB ){
    my $end = $start + $SAMPLES_PER_JOB - 1;
    $end = $#samples  if ( $end > $#samples );
    my @jobSamples = @samples[$start..$end];
    $jobIndex++;
    $pm->start($jobIndex) and next;

    my $bgee_thread = Utils::connect_bgee_db($bgee_connector);
    $bgee_thread->{RaiseError} = 1;
    my $result = fix_samples($bgee_thread, \@jobSamples, \%pValuePerSpecies, $debug);
    $bgee_thread->disconnect;
    $pm->finish(0, $result);
}
$pm->wait_all_children;

# -----------------------------------------------------------------------------
# Summary, and conditions of the 'expression' table to recompute
# -----------------------------------------------------------------------------
my $totalReplaced = 0;
for my $speciesId ( sort { $a <=> $b } keys %replacedPerSpecies ){
    print "Species $speciesId: $replacedPerSpecies{$speciesId} p-values ", ($debug ? 'to replace' : 'replaced'), "\n";
    $totalReplaced += $replacedPerSpecies{$speciesId};
}
print "Done, $totalReplaced p-values ", ($debug ? 'to replace' : 'replaced'), " in total.\n";

for my $datatype ( sort keys %conditionsToRecompute ){
    my @conditionIds = sort { $a <=> $b } keys %{$conditionsToRecompute{$datatype}};
    print "Conditions of the expression table to recompute for $datatype (", scalar(@conditionIds), "): ",
          join(',', @conditionIds), "\n";
}
if ( %conditionsToRecompute ){
    print "Recompute them with: perl update_rna_seq_expression.pl -bgee=<BGEECMD> ",
          join(' ', map { "-$_" } sort keys %conditionsToRecompute), " -processAll\n";
}
if ( $failedJobs > 0 ){
    print "$failedJobs jobs failed, relaunch the script to process their samples.\n";
    exit 1;
}

exit 0;

# -----------------------------------------------------------------------------
# Subroutine: fix_samples
#   Replaces (or counts, in debug mode) the p-values of 0 with a NULL z-score of the provided
#   samples, and returns the number of p-values replaced per species, the conditions holding them
#   per datatype, and the number of samples processed.
# -----------------------------------------------------------------------------
sub fix_samples {
    my ($bgee_thread, $jobSamples, $pValuePerSpecies, $debugMode) = @_;

    my $updatePValues = $bgee_thread->prepare(
        'UPDATE rnaSeqLibraryAnnotatedSampleGeneResult SET pValue = ? '.
        'WHERE rnaSeqLibraryAnnotatedSampleId = ? AND pValue = 0 AND zScore IS NULL'
    );
    my $countPValues = $bgee_thread->prepare(
        'SELECT COUNT(*) FROM rnaSeqLibraryAnnotatedSampleGeneResult '.
        'WHERE rnaSeqLibraryAnnotatedSampleId = ? AND pValue = 0 AND zScore IS NULL'
    );
    my %result = ('replaced' => {}, 'conditions' => {}, 'samples' => 0);
    for my $sample ( @$jobSamples ){
        my ($sampleId, $speciesId, $exprMappedConditionId, $datatype) = @$sample;
        my $replaced;
        if ( $debugMode ){
            $countPValues->execute($sampleId) or die $countPValues->errstr;
            ($replaced) = $countPValues->fetchrow_array;
            $countPValues->finish;
        } else {
            $replaced = $updatePValues->execute($pValuePerSpecies->{$speciesId}, $sampleId)
                or die $updatePValues->errstr;
            $replaced = 0  if ( $replaced eq '0E0' );
        }
        $result{'replaced'}->{$speciesId} += $replaced;
        $result{'conditions'}->{$datatype}->{$exprMappedConditionId} = 1  if ( $replaced > 0 );
        $result{'samples'}++;
    }
    $updatePValues->finish;
    return \%result;
}
