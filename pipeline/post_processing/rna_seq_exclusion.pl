#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Parallel::ForkManager;

# Julien W. Feb. 2026
#this script removes pvalues for calls not to consider based on the biotype of the genes and the population captured by the RNA-Seq library e.g. if the library is polyA, we will not consider pvalues for non-polyA genes, and then we will set their pvalue to NA. if a pvalue is set to NA it will not be considered in the downstream analyses (insertion in expression table and propagation) and the reason for exclusion will be "biotype not targeted".

#Julien W. March 2026
# Add possibility to process gene exclusion for a subset of species. Useful if few data have to be added and impact only a subset of species, to speed up the processing.
# WARNING: This option also make the script delete all data for the specified species from the expression table if the expressionId does not also have insitu data. I also update expressionId = NULL in the rnaSeqSampleGeneResult table for the specified species.

# USAGE: perl rna_seq_pvalues_cleaning.pl -bgee=connection_string -number_threads=number_of_threads -debug

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my ($bgee_connector) = ('');
my $number_threads   = (0);
my ($debug)          = (0);
my ($remove_pvalue)  = (0);
my ($pre_filtering)  = (0);
my ($biotype_filtering)  = (0);
my ($init_post_processing)  = (-1);
my $species_to_process = '';
my %opts = ('bgee=s'            => \$bgee_connector,   # Bgee connector string
            'numberThreads=s'  => \$number_threads,   # Number of threads to run in parallel
            'debug'             => \$debug,
            'pre_filtering'     => \$pre_filtering,
            'biotype_filtering' => \$biotype_filtering,
            'removePvalue'     => \$remove_pvalue,
            'initPostProcessing=i' => \$init_post_processing,
            'speciesToProcess=s' => \$species_to_process,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $number_threads eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t-bgee              Bgee connector string
\t-number_threads    Number of update threads to run in parallel
\t-debug             printing the update SQL queries, not executing them
\t-pre_filtering     If set, genes with zero counts in all libraries will be excluded with reason \"pre-filtering\"
\t-biotype_filtering If set, genes that are not targeted by the population captured by the library will be excluded with reason \"biotype not targeted\". This filtering is applied after the pre_filtering, so it can update calls that were previously excluded for pre-filtering but that are not anymore in the current release because of new calls. This is only useful for minor releases.
\t-remove_pvalue     If set, pvalues will be set to NA for calls that are not to be considered based on the biotype of the genes and the population captured by the RNA-Seq library, otherwise only the reason for exclusion will be updated but pvalues will not be changed.
\t-init_post_processing If 1: the script will first update all previously excluded genes to \"not excluded\". It will then delete all expression data from the expression table if the expressionId does not also have insitu data. Finally, it will set expressionId to NULL in the rnaSeqLibraryAnnotatedSampleGeneResult. Combined with the species_to_process option, this will reinitialize post_processing for the selected species. This option should always be selected if new RNA-Seq data are added to the database and some postprocessing steps have already been run.
                     If 0: the script will not update previously pre-filtered genes. Only use that option for a brand new database where no postprocessing steps of RNA-Seq data have been run yet, or if you are sure that no new RNA-Seq calls have been added since the last release.
\t-species_to_process Comma-separated list of speciesId to process. If not set, all species will be processed for exclusion.
\n";
    exit 1;
}

# Validate numberThreads
if ( $number_threads < 0 || $number_threads =~ /\D/ ){
    print "\tError: numberThreads must be a non-negative integer\n";
    exit 1;
}
$number_threads = 1 if $number_threads == 0;
if ($init_post_processing != 0 && $init_post_processing != 1) {
    print "\tError: initPostProcessing must be either 0 (not update previously pre-filtered genes) or 1 (update previously pre-filtered genes to not excluded)\n";
    exit 1;
}

# initialize database connection variable
my $bgee = '';

$| = 1;

sub species_list_to_sql {
    my ($list) = @_;
    my @items = split /,/, $list;
    return join(',', map { "'$_'" } @items);
}

sub execute_with_deadlock_retry {
    my ($sth, @params) = @_;
    my $max_retries = 5;
    
    for (my $retry = 0; $retry < $max_retries; $retry++) {
        eval {
            $sth->execute(@params) or die $sth->errstr;
        };
        if ($@) {
            if (($@ =~ /Deadlock found/i || $@ =~ /Lost connection/i || $@ =~ /MySQL server has gone away/i) && $retry < $max_retries - 1) {
                my $error_type = $@ =~ /Deadlock/i ? 'Deadlock' : 'Connection loss';
                warn "$error_type detected, retrying (attempt ".($retry + 2)."/$max_retries)...\n";
                sleep(2 + $retry);  # Exponential backoff
                next;
            }
            die $@;
        }
        return;  # Success
    }
}

# first update all previously excluded genes to "not excluded". It allows to identify genes that wereexcluded in a previous release because of zero-count filtering but that are not anymore in the current releasebecause of new calls. This is only useful for minor releases.
if ( $init_post_processing == 1 ){
    print "undo postprocessing step already done for specified species or all of them if non is provided...\n";
    
    $bgee = Utils::connect_bgee_db($bgee_connector);
    # number of annotated samples to process per batch
    my $BATCH_SIZE = 100;
    
    # Get list of species to process
    my @species_list = ();
    if ($species_to_process) {
        @species_list = split /,/, $species_to_process;
    } else {
        # Get all species
        my $all_species_query = "SELECT DISTINCT speciesId FROM species";
        my $sth_all_species = $bgee->prepare($all_species_query);
        execute_with_deadlock_retry($sth_all_species);
        while (my @row = $sth_all_species->fetchrow_array()) {
            push @species_list, $row[0];
        }
        $sth_all_species->finish();
    }
    $bgee->disconnect();
    print "Processing ".scalar(@species_list)." species to reinitialize post-processing\n";
    
    # Process each species with batch approach
    for my $speciesId (@species_list) {
        print "Processing species $speciesId\n";
        $bgee = Utils::connect_bgee_db($bgee_connector);
        
        # Get list of sample IDs
        my $sample_query = "SELECT DISTINCT t1.rnaSeqLibraryAnnotatedSampleId 
                            FROM rnaSeqLibraryAnnotatedSample AS t1
                            INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId
                            WHERE t2.speciesId = ?
                            ORDER BY t1.rnaSeqLibraryAnnotatedSampleId";
        my $sample_sth = $bgee->prepare($sample_query);
        execute_with_deadlock_retry($sample_sth, $speciesId);
        
        my @sample_ids;
        while (my ($sample_id) = $sample_sth->fetchrow_array()) {
            push @sample_ids, $sample_id;
        }
        $sample_sth->finish();
        $bgee->disconnect(); # disconnect before forking to avoid connection issues in child processes
        
        my $sample_count = scalar(@sample_ids);
        print "Found $sample_count samples\n";
        
        if ($sample_count == 0) {
            next;
        }
        
        if ( $debug ){
            print "DEBUG: would update $sample_count samples in batches of $BATCH_SIZE using $number_threads threads\n";
        } else {
            print "  Processing $sample_count samples using $number_threads threads\n";
            
            # Parallelize batch processing using ForkManager
            my $pm = Parallel::ForkManager->new($number_threads);
            
            # Track errors from child processes
            $pm->run_on_finish(
                sub {
                    my ($pid, $exit_code) = @_;
                    if ($exit_code != 0) {
                        die "Child process $pid exited with error code $exit_code\n";
                    }
                }
            );
            
            for (my $i = 0; $i < $sample_count; $i += $BATCH_SIZE) {
                $pm->start and next; # fork and continue to next batch
                
                # Start a new connection for each thread
                my $bgee_thread = Utils::connect_bgee_db($bgee_connector);
                
                # Get batch of sample IDs
                my $end_idx = ($i + $BATCH_SIZE - 1 < $sample_count) ? $i + $BATCH_SIZE - 1 : $sample_count - 1;
                my @batch_samples = @sample_ids[$i..$end_idx];
                
                # Build and execute update query
                my $placeholders = join(',', ('?') x scalar(@batch_samples));
                my $update_query = "UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1
                                   INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId
                                   SET t1.reasonForExclusion = ?, t1.expressionId = NULL
                                   WHERE t1.rnaSeqLibraryAnnotatedSampleId IN ($placeholders)
                                   AND t2.speciesId = ?";
                
                my $update_sth = $bgee_thread->prepare($update_query);
                execute_with_deadlock_retry($update_sth, $Utils::CALL_NOT_EXCLUDED, @batch_samples, $speciesId);
                $update_sth->finish();
                
                $bgee_thread->disconnect();
                $pm->finish; # end the child process
            }
            $pm->wait_all_children();
        }
    }
    
    # Reconnect to database after parallel processing (forking can corrupt parent connection)
    $bgee = Utils::connect_bgee_db($bgee_connector);
    
    # Delete expression data (not batched, should be relatively fast)
    my $delete_expression_query = "delete t1 from expression as t1 inner join cond as t2 on t1.conditionId = t2.conditionId where t1.inSituNumberObs = 0";
    if ($species_to_process) {
        $delete_expression_query .= " and t2.speciesId IN (".species_list_to_sql($species_to_process).")";
    }
    if ( $debug ){
        print "DEBUG: would execute query: $delete_expression_query with species IDs: ".$species_to_process."\n";
    }
    else {
        my $sth_delete_expression = $bgee->prepare($delete_expression_query);
        execute_with_deadlock_retry($sth_delete_expression);
        $sth_delete_expression->finish();
    }
    $bgee->disconnect();
    
    print "Init post-processing completed\n";
}

if ($pre_filtering) {

    print "Start filtering genes with no positive read/UMI counts...\n";

    $bgee = Utils::connect_bgee_db($bgee_connector);

    # Get all Bgee species
    my $select_species = "SELECT DISTINCT speciesId FROM species";
    if ($species_to_process) {
        $select_species .= " WHERE speciesId IN (".species_list_to_sql($species_to_process).")";
    }
    my $sth_species = $bgee->prepare($select_species);
    execute_with_deadlock_retry($sth_species);
    my @species_to_process = ();
    while ( my @row = $sth_species->fetchrow_array() ){
        push @species_to_process, $row[0];
    }
    $sth_species->finish();

    print "Start to process ".scalar(@species_to_process)." species to exclude genes with only zero-count filtering.\n";

    $bgee->disconnect(); # disconnect before forking to avoid connection issues in child processes

    # Parallelize by species for better load balancing
    my $pm2 = Parallel::ForkManager->new($number_threads);
    for my $speciesId ( @species_to_process ){
        $pm2->start and next;
        # start a new connection for each thread
        my $bgee_thread = Utils::connect_bgee_db($bgee_connector);
        
        print "Processing species $speciesId for zero-count filtering...\n";
        
        # Find genes for this species that never have positive readsCount or UMIsCount
        my $find_genes_query = "SELECT t1.bgeeGeneId
                                FROM rnaSeqLibraryAnnotatedSampleGeneResult t1
                                JOIN gene t2 ON t1.bgeeGeneId = t2.bgeeGeneId
                                WHERE t2.speciesId = ?
                                GROUP BY t1.bgeeGeneId
                                HAVING SUM(CASE WHEN readsCount > 0 THEN 1 ELSE 0 END) = 0
                                AND SUM(CASE WHEN UMIsCount > 0 THEN 1 ELSE 0 END) = 0;";
        
        my $sth_genes = $bgee_thread->prepare($find_genes_query);
        execute_with_deadlock_retry($sth_genes, $speciesId);
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
                    execute_with_deadlock_retry($sth_update, $Utils::EXCLUDED_FOR_PRE_FILTERED, @batch, $Utils::CALL_NOT_EXCLUDED);
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

}

if ( $biotype_filtering ) {
    print "Start to filter genes based on the biotype and population captured by the RNA-Seq library...\n";

    $bgee = Utils::connect_bgee_db($bgee_connector);
    # retrieve all rnaSeqLibraryAnnotatedSampleId and corresponding rnaSeqPopulationCaptureId. Do not check the reason for exclusion yet as the query detecting already filtered ones would be really slow. We will check the reason for exclusion in the update query to only update calls that are not already excluded for another reason. This way we will also update calls that were previously excluded for another reason but that are not anymore in the current release because of new calls. This is only useful for minor releases.
    my $select_libs = "SELECT distinct t2.rnaSeqLibraryAnnotatedSampleId, t1.rnaSeqPopulationCaptureId FROM rnaSeqLibrary AS t1 INNER JOIN rnaSeqLibraryAnnotatedSample t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId";
    if ($species_to_process) {
        $select_libs .= " INNER JOIN cond as t3 on t2.conditionId = t3.conditionId WHERE t3.speciesId IN (".species_list_to_sql($species_to_process).")";
    }
    my $sth_libs = $bgee->prepare($select_libs);
    execute_with_deadlock_retry($sth_libs);
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

    $bgee->disconnect(); # disconnect before forking to avoid connection issues in child processes

    # parallelize the updates by annotated sample for better load balancing
    my $pm = Parallel::ForkManager->new($number_threads);
    for my $annotatedSampleId ( keys %rna_seq_libraries_to_process ){
        $pm->start and next; # fork and continue to the next sample
        # start a new connection for each thread
        my $bgee_thread = Utils::connect_bgee_db($bgee_connector);
        
        my $popCapturedId = $rna_seq_libraries_to_process{$annotatedSampleId};
        

        if ( $debug ){
            print "DEBUG: would execute query: $update_query with popCapturedId: $popCapturedId and annotatedSampleId: $annotatedSampleId. Would set reasonForExclusion to $Utils::BIOTYPE_NOT_TARGETED when reasonForExclusion is $Utils::CALL_NOT_EXCLUDED\n";
        }
        else {
            my $sth_update = $bgee_thread->prepare($update_query);
            execute_with_deadlock_retry($sth_update, $popCapturedId, $Utils::BIOTYPE_NOT_TARGETED, $annotatedSampleId, $Utils::CALL_NOT_EXCLUDED);
            $sth_update->finish();
        }
        
        $bgee_thread->disconnect();
        $pm->finish; # end the child process
    }
    $pm->wait_all_children();

    print "biotype filtering done.\n";
}

