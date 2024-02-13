#!/usr/bin/env perl

## Julien Wollbrett, Jan 28, 2020
# This script launches the slurm jobs that generate GTF files containing intergenic sequences
# The script 0Before/prepare_GTF.R is launched for each GTF file previously downloaded in Ensembl
# It is inspired from https://svn.vital-it.ch/svn/Selectome/trunk/scripts/pipeline/bsub_scheduler.pl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use File::Path qw(make_path);
use File::Slurp;
use FindBin qw($RealBin);
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;
use Getopt::Long;
use Time::localtime;


# Define arguments & their default value
my ($gtf_dir, $block_size_N, $proportion_N, $account, $partition, $output_gtf_path, $output_log_folder, $cluster_R_cmd) = ( '', '', '', '', '', '', '', '');
my %opts = ('gtf_dir=s'             => \$gtf_dir,
            'block_size_N=s'        => \$block_size_N,
            'proportion_N=s'        => \$proportion_N,
            'account=s'              => \$account,
            'partition=s'            => \$partition,
            'output_gtf_path=s'     => \$output_gtf_path,
            'output_log_folder=s'   => \$output_log_folder,
            'cluster_R_cmd=s'       => \$cluster_R_cmd
            );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $gtf_dir eq '' || $block_size_N eq '' || $proportion_N eq '' || $account eq '' || $partition eq '' || $output_gtf_path eq '' || $output_log_folder eq '' || $cluster_R_cmd eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -gtf_dir=\$(RNASEQ_CLUSTER_GTF) -block_size_N=31 -proportion_N=0.05 -output_gtf_path=\$(RNASEQ_CLUSTER_GTF) -output_log_folder=\$(RNASEQ_CLUSTER_GTF) -cluster_R_cmd=\$(CLUSTER_R_CMD)
\t-gtf_dir=s            Folder with GTF files downloaded from Ensembl
\t-block_size_N=s       max number of consecutive Ns present in intergenic sequences
\t-proportion_N=s       max proportion of Ns in intergenic sequences
\t-account=s            name of the account used to run the script on the cluster
\t-partition=s          slurm partition to use to run the script on the cluster
\t-output_gtf_path=s    Folder where GTF files with intergenic regions are saved
\t-output_log_folder=s  Folder for .out and .err files (produced by queuing system), and .Rout files produced by R
\t-cluster_R_cmd=s      Command used to load R on the cluster
\n";
    exit 1;
}

# Create output folder if not present
make_path "$output_log_folder",   {verbose=>0, mode=>0775};
print "module add : $cluster_R_cmd\n";

# Setting up SLURM parameters #################################
my $main_script = $RealBin.'/prepare_GTF.R';

my $nbr_processors = 1;
# RAM needed: 30GB should be enough
my $memory_usage   = 30;      # in GB

my $jobs_during_day   = 100; # Number of simultaneous jobs during working days
my $jobs_during_night = 120; # Number of simultaneous jobs during week-end & night

# create one unique name for all the jobs
my $job_name = "prepare_GTF";

# reading library infos
my $count = 0;

print "will launch jobs with name : $job_name\n";

opendir (DIR, $gtf_dir) or die $!;

JOB:
while (my $file = readdir(DIR) ){

    next JOB if (! ($file =~ /\.gtf\.gz$/i)); # header

    $count++;

    # Create output folder for library
    my $genome_file = $file=~ s/gtf\.gz$/genome\.fa/r;
    my $output_log_file = $file=~ s/gtf\.gz$/out/r;
    my $err_log_file = $file=~ s/gtf\.gz$/err/r;
    my $sbatch_file_name = $file=~ s/gtf\.gz$/sbatch/r;

    $err_log_file = $output_log_folder.'/'.$err_log_file;
    $output_log_file = $output_log_folder.'/'.$output_log_file;
    my $sbatch_file = "${output_log_folder}/${sbatch_file_name}";

    my $script_plus_args = "R CMD BATCH --vanilla --slave '--args gene_gtf_path=\"$gtf_dir/$file\" output_gtf_path=\"$output_gtf_path\" genome_fasta_path=\"$gtf_dir/$genome_file\" N_block_size=$block_size_N N_proportion=$proportion_N'  0Before/prepare_GTF.R";

    my $sbatch_command = "$cluster_R_cmd\n";
    $sbatch_command .= $script_plus_args;
    print "Command submitted to cluster:\n$sbatch_command\n";

    # Create the SBATCH script
    open (my $OUT, '>', "$sbatch_file")  or die "Cannot write [$sbatch_file]\n";
    print {$OUT} Utils::sbatch_template($partition, $account, $nbr_processors, $memory_usage, $output_log_file, $err_log_file, $job_name);
    print {$OUT} "$sbatch_command\n";
    close $OUT;

    # Then, run the job
    system("sbatch $sbatch_file")==0  or print "Failed to submit job [$file]\n";
}

close(DIR);

# check that all jobs finished
my $active_jobs = Utils::check_active_jobs_number($job_name);
while ($active_jobs > 0) {
    sleep(30);
    $active_jobs = Utils::check_active_jobs_number($job_name);
}
print "delete sbatch files...\n";
system("rm -f $gtf_dir./*.sbatch");
print "\nGeneration of custom GTF files has been done for $count species. \\_o_/\n";
exit 0;


