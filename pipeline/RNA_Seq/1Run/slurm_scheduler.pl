#!/usr/bin/env perl

## Julien Roux, Jan 11, 2016
# This script launches the slurm jobs
# The script rna_seq_mapping_and_analysis.pl is launched for each library
# It is inspired from https://svn.vital-it.ch/svn/Selectome/trunk/scripts/pipeline/bsub_scheduler.pl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use File::Path qw(make_path);
use File::Slurp;
use FindBin qw( $RealBin ); # directory where the script is lying
use Getopt::Long;

# Define arguments & their default value
my ($sample_info_file, $exclude_sample_file, $output_log_folder, $index_folder, $fastq_folder, $kallisto_out_folder, $ens_release, $ens_metazoa_release, $enc_passwd_file, $cluster_kallisto_cmd, $cluster_R_cmd) = ('', '', '', '', '', '', '', '', '', '', '', '', '');
my %opts = ('sample_info_file=s'     => \$sample_info_file,
            'exclude_sample_file=s'  => \$exclude_sample_file,
            'output_log_folder=s'    => \$output_log_folder,
            'index_folder=s'         => \$index_folder, # same as GTF folder
            'fastq_folder=s'         => \$fastq_folder,
            'kallisto_out_folder=s'  => \$kallisto_out_folder,
            'ens_release=s'          => \$ens_release,
            'ens_metazoa_release=s'  => \$ens_metazoa_release,
            'enc_passwd_file=s'      => \$enc_passwd_file,
            'cluster_kallisto_cmd=s' => \$cluster_kallisto_cmd,
            'cluster_R_cmd=s'        => \$cluster_R_cmd,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $sample_info_file eq '' || $output_log_folder eq '' || $index_folder eq '' || $fastq_folder eq '' || $kallisto_out_folder eq '' || $ens_release eq '' || $ens_metazoa_release eq '' || $enc_passwd_file eq '' || $cluster_kallisto_cmd eq '' || $cluster_R_cmd eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -sample_info_file=\$(RNASEQ_SAMPINFO_FILEPATH) -exclude_sample_file=\$(RNASEQ_SAMPEXCLUDED_FILEPATH) -output_log_folder=\$(RNASEQ_CLUSTER_LOG) -index_folder=\$(RNASEQ_CLUSTER_GTF)  -fastq_folder=\$(RNASEQ_SENSITIVE_FASTQ) -kallisto_out_folder=\$(RNASEQ_CLUSTER_ALL_RES) -ens_release=\$(ENSRELEASE) -ens_metazoa_release=\$(ENSMETAZOARELEASE) -enc_passwd_file=\$(ENCRYPT_PASSWD_FILE) -cluster_kallisto_cmd=\$(CLUSTER_KALLISTO_CMD) -cluster_R_cmd=\$(CLUSTER_R_CMD)
\t-sample_info_file       rna_seq_sample_info.txt
\t-exclude_sample_file    rna_seq_sample_excluded.txt
\t-output_log_folder      folder for .out and .err files (produced by queuing system), and .Rout files produced by R
\t-index_folder=s         Folder with Kallisto indexes (same as GTF folder)
\t-fastq_folder=s         Folder with Fastq files on big bgee
\t-kallisto_out_folder=s  Folder with Kallisto output and results
\t-ens_release=s          Ensembl release
\t-ens_metazoa_release=s  Ensembl Metazoa release
\t-enc_passwd_file=s      File with password necessary to decrypt the GTEx data
\t-cluster_kallisto_cmd=s Command to load kallisto module on cluster
\t-cluster_R_cmd=s        Command to load R module on cluster
\n";
    exit 1;
}

# Tests
die "Invalid or missing [$sample_info_file]: $?\n"  if ( !-e $sample_info_file || !-s $sample_info_file );
# Create output folder if not present
make_path "$output_log_folder",   {verbose=>0, mode=>0775};

# Setting up SLURM parameters #################################
my $main_script = $RealBin.'/rna_seq_mapping_and_analysis.pl';
## TODO launch slurm_scheduler.pl from /data/ul/dee/bgee/GIT/pipeline/RNA_Seq/
##        Beware that git pull command should be executed before
##        kallisto_out_folder should be on /scratch/temporary. If too slow, consider using /scratch/local/ + cp of results file to /scratch/temporary/ or /home/bbgee, or /data/ (read-only, but should be fine via scp)

# kallisto is no multithreaded unless bootstraps are used
my $nbr_processors = 1;
# RAM needed: 10GB should be enough
my $memory_usage   = 10;      # in GB
my $user_email     = 'bgee@sib.swiss'; # for email notification
my $account        = 'mrobinso_bgee';
my $queue          = 'ax-long';

my $job_limit      = 120; # Number of simultaneous jobs running

# Sample to exclude if any
my %manually_excluded;
EXCLUSION:
for my $line ( read_file("$exclude_sample_file", chomp=>1) ){
    #libraryId    excluded    comment    annotatorId    lastModificationDate
    next EXCLUSION  if ( $line =~ /^#/); # header or comment

    my ($sampleId, $to_exclude) = split(/\t/, $line);
    next EXCLUSION  if ( $to_exclude ne 'TRUE' );

    $manually_excluded{$sampleId} = 1;
}

# reading library infos
my $count = 0;
JOB:
for my $line ( read_file("$sample_info_file", chomp=>1) ){
    next JOB  if ( $line =~ /^#/); # header

    my @fields     = split ("\t", $line);
    my $library_id = $fields[0];
    # Excluded library
    next JOB  if ( exists $manually_excluded{$library_id} );

    # Test to not re-run already finished jobs
    if ( -s "$output_log_folder/$library_id/DONE.txt" ){
        print "\n$library_id not launched because already analyzed\n";
        next JOB;
    }
    # Check running jobs to not resubmit them while running
    if ( `squeue --user=\$USER --account=$account --long | grep ' $library_id '` =~ / (RUNN|PEND)ING / ){
        print "\n$library_id not launched because it is currently being analyzed (see squeue/sacct)\n";
        next JOB;
    }

    # Let's launch this library!
    print "\nLaunching $library_id ...\n";
    $count++;

    # Create output folder for library
    make_path "$output_log_folder/$library_id", {verbose=>0, mode=>0775};

    # library-specific arguments
    my $output_file = $output_log_folder.'/'.$library_id.'/'.$library_id.'.out';
    my $rm_output_command = 'rm -f '.$output_file.";\n";

    my $error_file = $output_log_folder.'/'.$library_id.'/'.$library_id.'.err';
    my $rm_error_command = 'rm -f '.$error_file.";\n";

    my $sbatch_file = $output_log_folder.'/'.$library_id.'/'.$library_id.'.sbatch';


    #TODO Update for paths in Jura +simplify
    my $script_plus_args = "/usr/bin/time -v perl $main_script -library_id=$library_id -sample_info_file=$sample_info_file -exclude_sample_file=$exclude_sample_file -index_folder=$index_folder -fastq_folder=$fastq_folder -kallisto_out_folder=$kallisto_out_folder -output_log_folder=$output_log_folder -ens_release=$ens_release -ens_metazoa_release=$ens_metazoa_release -enc_passwd_file=$enc_passwd_file";


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
    # Potential other options:
    # #SBATCH --mail-user=$user_email
    # #SBATCH --mail-type=ALL

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
