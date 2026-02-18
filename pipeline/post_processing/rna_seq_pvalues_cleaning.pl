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
my $number_threads   = '';
my ($debug)          = (0);
my ($remove_pvalue)  = (0);
my %opts = ('bgee=s'            => \$bgee_connector,   # Bgee connector string
            'number_threads=s'  => \$number_threads,   # Number of threads to run in parallel
            'debug'             => \$debug,
            'remove_pvalue'     => \$remove_pvalue,
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

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

$| = 1;

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

