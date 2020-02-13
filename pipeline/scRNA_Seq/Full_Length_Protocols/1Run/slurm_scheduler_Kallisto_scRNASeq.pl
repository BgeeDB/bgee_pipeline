#!/usr/bin/env perl

## Script created by SMoretti and adapted by SFonsecaCosta to run Kallisto for scRNA-Seq data
## This script launches the slurm jobs for Kallisto using R script

# Perl core modules
use strict;
use warnings;
use diagnostics;

use File::Path qw(make_path);
use File::Slurp;
use FindBin qw( $RealBin ); # directory where the script is lying
use Getopt::Long;

# Define arguments & their default value
my ($scrna_seq_sample_info, $raw_cells_folder, $infoFolder, $output_folder, $cluster_kallisto_cmd, $cluster_R_cmd) = ('', '', '', '', '', '');
my %opts = ('scrna_seq_sample_info=s'     => \$scrna_seq_sample_info,
            'raw_cells_folder=s'          => \$raw_cells_folder,
            'infoFolder=s'                => \$infoFolder, # same as GTF folder
            'output_folder'               => \output_folder,
            'cluster_kallisto_cmd=s'      => \$cluster_kallisto_cmd,
            'cluster_R_cmd=s'             => \$cluster_R_cmd,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $scrna_seq_sample_info eq '' || $raw_cells_folder eq '' ||  $infoFolder eq '' || $output_folder eq '' || $cluster_kallisto_cmd eq '' || $cluster_R_cmd eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -scrna_seq_sample_info=\$(SC_RNASEQ_SAMPINFO_FILEPATH) -raw_cells_folder=\$(SC_RNASEQ_FASTQ_FULL_LENGTH) -infoFolder=\$(RNASEQ_CLUSTER_GTF) -output_folder=\$(SC_RNASEQ_CLUSTER_KALLISTO) -cluster_kallisto_cmd=\$(CLUSTER_KALLISTO_CMD) -cluster_R_cmd=\$(CLUSTER_R_CMD)
\t-scrna_seq_sample_info    scrna_seq_sample_info
\t-raw_cells_folder=s       Folder with are all libraries with Fastq files
\t-infoFolder=s             Folder with Kallisto indexes (same as GTF folder)
\t-output_folder=s          Folder where the kallisto output should be saved
\t-cluster_kallisto_cmd=s   Command to load kallisto module on cluster
\t-cluster_R_cmd=s          Command to load R module on cluster
\n";
    exit 1;
}

# Tests
die "Invalid or missing [$scrna_seq_sample_info]: $?\n"  if ( !-e $scrna_seq_sample_info || !-s $scrna_seq_sample_info );

my $main_script = $RealBin.'/kallisto.R';

# kallisto is no multithreaded unless bootstraps are used
my $nbr_processors = 1;
my $memory_usage   = 10;      # in GB
my $user_email     = 'sara.fonsecacosta@unil.ch'; # for email notification
my $account        = 'mrobinso_bgee';
my $queue          = 'normal';

my $job_limit      = 120; # Number of simultaneous jobs running

# reading library infos
my $count = 0;
JOB:
for my $line ( read_file("$scrna_seq_sample_info", chomp=>1) ){
    next JOB  if ( $line =~ /^#/); # header

    my @fields     = split ("\t", $line);
    my $library_id = $fields[0];

    # Test to not re-run already finished jobs
    if ( -s "$output_folder/$library_id/DONE.txt" ){
        print "\n$library_id not launched because already analyzed\n";
        next JOB;
    }
    # Check running jobs to not resubmit them while running
    if ( "`squeue --user=\$USER --account=$account --long | grep ' $library_id '` =~ / (RUNN|PEND)ING /"){
        print "\n$library_id not launched because it is currently being analyzed (see squeue/sacct)\n";
        next JOB;
    }

    # Let's launch this library!
    print "\nLaunching $library_id ...\n";
    $count++;

    # library-specific arguments
    my $output_file = $output_folder.'/'.$library_id.'/'.$library_id.'.out';
    my $rm_output_command = 'rm -f '.$output_file.";\n";

    my $error_file = $output_folder.'/'.$library_id.'/'.$library_id.'.err';
    my $rm_error_command = 'rm -f '.$error_file.";\n";

    my $sbatch_file = $output_folder.'/'.$library_id.'/'.$library_id.'.sbatch';

    my $script_plus_args = "/usr/bin/time -v R CMD BATCH --no-save --no-restore '--args library_id=\"$library_id\" scrna_seq_sample_info=\"$scrna_seq_sample_info\" raw_cells_folder=\"$raw_cells_folder\" infoFolder=\"$infoFolder\" output_folder=\"$output_folder\"' $main_script $output_folder/$library_id/kallisto.Rout";

    # Wait for free places in job queue
    my $running_jobs = check_running_jobs();
    WAIT_FREE_JOB_IN_QUEUE:
    while ( $running_jobs >= $job_limit ){
        print "No more possible slot for the job, waiting and resubmitting\n";
        sleep 30;
        $running_jobs = check_running_jobs();
    }

    # Script can be launched! Construct SLURM sbatch command:
    # First, remove previous .out and .err files
    my $sbatch_command = $rm_output_command.$rm_error_command;
    $sbatch_command .= 'module add Bioinformatics/Software/vital-it'."\n";
    $sbatch_command .= "$cluster_kallisto_cmd\n";
    $sbatch_command .= "$cluster_R_cmd\n\n";
    $sbatch_command .= $script_plus_args;
    print "Command submitted to cluster:\n$sbatch_command\n";

    # Create the SBATCH script
    open (my $OUT, '>', "$sbatch_file")  or die "Cannot write [$sbatch_file]\n";
    print {$OUT} sbatch_template($queue, $account, $nbr_processors, $memory_usage, $output_file, $error_file, $library_id);
    print {$OUT} "$sbatch_command\n";
    close $OUT;

    # Then, run the job
    system("sbatch $sbatch_file")==0  or print "Failed to submit job [$library_id]\n";
}

print "\n######################################################\nAll done. $count jobs submitted.\n######################################################\n";
exit 0;

sub check_running_jobs {
    my $running_jobs = `squeue --user=\$USER --account=$account | grep -v 'JOBID' | wc -l` || 0;
    chomp($running_jobs);
    return $running_jobs;
}

# Add main sbatch command and options
sub sbatch_template {
    my ($queue, $account, $nbr_processors, $memory_usage, $output_file, $error_file, $library_id) = @_;

    my $template="#!/bin/bash

#SBATCH --partition=$queue
#SBATCH --account=$account

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$nbr_processors
#SBATCH --mem=${memory_usage}G
##SBATCH --time=...

#SBATCH --output=$output_file
#SBATCH --error=$error_file
#SBATCH --export=NONE
#SBATCH --job-name=$library_id

";

    return $template;
}
