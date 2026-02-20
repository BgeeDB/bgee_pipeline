#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Parallel::ForkManager;

# Julien W. Feb. 2026
#this script removes pvalues for calls not to consider based on the biotype of the genes and the population captured by the RNA-Seq library e.g. if the library is polyA, we will not consider pvalues for non-polyA genes, and then we will set their pvalue to NA. if a pvalue is set to NA it will not be considered in the downstream analyses (insertion in expression table and propagation) and the reason for exclusion will be "biotype not targeted".

# USAGE: perl rna_seq_pvalues_cleaning.pl -bgee=connection_string -number_threads=number_of_threads -debug

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my ($bgee_connector) = ('');
my $number_threads   = (0);
my ($debug)          = (0);
my ($remove_pvalue)  = (0);
my %opts = ('bgee=s'            => \$bgee_connector,   # Bgee connector string
            'numberThreads=s'  => \$number_threads,   # Number of threads to run in parallel
            'debug'             => \$debug,
            'removePvalue'     => \$remove_pvalue,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $number_threads eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\t-number_threads   number of update threads to run in parallel
\t-debug            printing the update SQL queries, not executing them
\t-remove_pvalue    if set, pvalues will be set to NA for calls that are not to be considered based on the biotype of the genes and the population captured by the RNA-Seq library, otherwise only the reason for exclusion will be updated but pvalues will not be changed.
\n";
    exit 1;
}

# Validate numberThreads
if ( $number_threads < 0 || $number_threads =~ /\D/ ){
    print "\tError: numberThreads must be a non-negative integer\n";
    exit 1;
}
$number_threads = 1 if $number_threads == 0;

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

$| = 1;

print "Start filtering genes with no positive read/UMI counts...\n";

# first update all previously pre-filtered genes to "not excluded". It allows to identify genes that were excluded in a previous release because of zero-count filtering but that are not anymore in the current release because of new calls. This is only useful for minor releases.
my $update_previously_excluded_query = "UPDATE rnaSeqLibraryAnnotatedSampleGeneResult SET reasonForExclusion = ? WHERE reasonForExclusion = ?";
if ( $debug ){
    print "DEBUG: would execute query: $update_previously_excluded_query with reasonForExclusion: ".Utils::$CALL_NOT_EXCLUDED." and previous reasonForExclusion: ".Utils::$EXCLUDED_FOR_PRE_FILTERED."\n";
}
else {
    my $sth_update_previously_excluded = $bgee->prepare($update_previously_excluded_query);
    $sth_update_previously_excluded->execute(Utils::$CALL_NOT_EXCLUDED, Utils::$EXCLUDED_FOR_PRE_FILTERED) or die $sth_update_previously_excluded->errstr;
    $sth_update_previously_excluded->finish();
}

# Get all Bgee species
my $select_species = "SELECT DISTINCT speciesId FROM species";
my $sth_species = $bgee->prepare($select_species);
$sth_species->execute();
my @species_to_process = ();
while ( my @row = $sth_species->fetchrow_array() ){
    push @species_to_process, $row[0];
}
$sth_species->finish();

print "Start to process ".scalar(@species_to_process)." species to exclude genes with only zero-count filtering.\n";

# Parallelize by species for better load balancing
my $pm2 = Parallel::ForkManager->new($number_threads);
for my $speciesId ( @species_to_process ){
    $pm2->start and next;
    # start a new connection for each thread
    my $bgee_thread = Utils::connect_bgee_db($bgee_connector);
    
    print "Processing species $speciesId for zero-count filtering...\n";
    
    # Find genes for this species that never have positive readsCount or UMIsCount
    my $find_genes_query = "SELECT DISTINCT t1.bgeeGeneId 
                            FROM rnaSeqLibraryAnnotatedSampleGeneResult AS t1 
                            INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId 
                            WHERE t2.speciesId = ? 
                            AND NOT EXISTS (
                                SELECT 1
                                FROM rnaSeqLibraryAnnotatedSampleGeneResult AS t3
                                WHERE t3.bgeeGeneId = t1.bgeeGeneId
                                AND (t3.readsCount > 0 OR t3.UMIsCount > 0)
                            )";
    
    my $sth_genes = $bgee_thread->prepare($find_genes_query);
    $sth_genes->execute($speciesId, $speciesId);
    my @genes_to_exclude = ();
    while ( my @row = $sth_genes->fetchrow_array() ){
        push @genes_to_exclude, $row[0];
    }
    $sth_genes->finish();
    
    if ( scalar(@genes_to_exclude) > 0 ){
        print "  Found ".scalar(@genes_to_exclude)." genes with zero counts for species $speciesId\n";
        
        # Update in batches of 1000 genes to avoid very long queries
        my $batch_size = 1000;
        for (my $i = 0; $i < scalar(@genes_to_exclude); $i += $batch_size) {
            my $end = $i + $batch_size - 1;
            $end = $#genes_to_exclude if $end > $#genes_to_exclude;
            my @batch = @genes_to_exclude[$i..$end];
            
            my $placeholders = join(',', ('?') x scalar(@batch));
            my $update_query = "UPDATE rnaSeqLibraryAnnotatedSampleGeneResult 
                               SET reasonForExclusion = ?";
            if ( $remove_pvalue ){
                $update_query .= ", pValue = NULL";
            }
            $update_query .= " WHERE bgeeGeneId IN ($placeholders) 
                               AND reasonForExclusion = ?";
            
            if ( $debug ){
                print "DEBUG: would execute query for batch ".($i/$batch_size + 1)." with ".scalar(@batch)." genes\n";
            }
            else {
                my $sth_update = $bgee_thread->prepare($update_query);
                $sth_update->execute(Utils::$EXCLUDED_FOR_PRE_FILTERED, @batch, Utils::$CALL_NOT_EXCLUDED) or die $sth_update->errstr;
                $sth_update->finish();
            }
        }
    } else {
        print "  No genes with zero counts found for species $speciesId\n";
    }
    
    $bgee_thread->disconnect();
    $pm2->finish; # end the child process
}
$pm2->wait_all_children();

print "finished pre-filtering genes with always zero-count.\n";

print "Start to filter genes based on the biotype and population captured by the RNA-Seq library...\n";
# retrieve rnaSeqLibraryAnnotatedSampleId and corresponding rnaSeqPopulationCaptureId
my $select_libs = "SELECT distinct t2.rnaSeqLibraryAnnotatedSampleId, t1.rnaSeqPopulationCaptureId FROM rnaSeqLibrary AS t1 INNER JOIN rnaSeqLibraryAnnotatedSample t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId where not exists (select 1 from rnaSeqLibraryAnnotatedSampleGeneResult where rnaSeqLibraryAnnotatedSampleGeneResult.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId and exclusionReason = '".Utils::$BIOTYPE_NOT_TARGETED."')";
my $sth_libs = $bgee->prepare($select_libs);
$sth_libs->execute();
my %rna_seq_libraries_to_process = ();
while ( my @row = $sth_libs->fetchrow_array() ){
    my ($annotatedSampleId, $popCapturedId) = @row;
    $rna_seq_libraries_to_process{$annotatedSampleId} = $popCapturedId;
}
$sth_libs->finish();

# default exclusion reason is Utils::$CALL_NOT_EXCLUDED, we will update it to Utils::$BIOTYPE_NOT_TARGETED for genes that are not targeted by the population captured by the library
my $update_query = "UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1 INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId LEFT JOIN rnaSeqPopulationCaptureToBiotype AS t3 ON t3.geneBiotypeId = t2.geneBioTypeId AND t3.rnaSeqPopulationCaptureId = ? SET t1.reasonForExclusion = ?";
if ( $remove_pvalue ){
    $update_query .= ", t1.pValue = NULL";
}
$update_query .= " WHERE t3.geneBiotypeId IS NULL AND t1.rnaSeqLibraryAnnotatedSampleId = ? AND t1.reasonForExclusion = ?";

# parallelize the updates by annotated sample for better load balancing
my $pm = Parallel::ForkManager->new($number_threads);
for my $annotatedSampleId ( keys %rna_seq_libraries_to_process ){
    $pm->start and next; # fork and continue to the next sample
    # start a new connection for each thread
    my $bgee_thread = Utils::connect_bgee_db($bgee_connector);
    
    my $popCapturedId = $rna_seq_libraries_to_process{$annotatedSampleId};
    

    if ( $debug ){
        print "DEBUG: would execute query: $update_query with popCapturedId: $popCapturedId and annotatedSampleId: $annotatedSampleId\n";
    }
    else {
        my $sth_update = $bgee_thread->prepare($update_query);
        $sth_update->execute($popCapturedId, Utils::$BIOTYPE_NOT_TARGETED, $annotatedSampleId, Utils::$CALL_NOT_EXCLUDED) or die $sth_update->errstr;
        $sth_update->finish();
    }
    
    $bgee_thread->disconnect();
    $pm->finish; # end the child process
}
$pm->wait_all_children();

print "biotype filtering done.\n";

