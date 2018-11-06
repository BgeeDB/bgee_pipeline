#!/usr/bin/env perl

## Julien Roux, Jan 11, 2016
# This script launches the bsub jobs
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
use Time::localtime;

# Define arguments & their default value
my ($sample_info_file, $exclude_sample_file, $output_log_folder, $index_folder, $fastq_folder, $kallisto_out_folder, $ens_release, $ens_metazoa_release, $data_host, $data_login, $enc_passwd_file, $vit_kallisto_cmd, $vit_R_cmd) = ('', '', '', '', '', '', '', '', '', '', '', '', '', '', '');
my %opts = ('sample_info_file=s'     => \$sample_info_file,
            'exclude_sample_file=s'  => \$exclude_sample_file,
            'output_log_folder=s'    => \$output_log_folder,
            'index_folder=s'         => \$index_folder, # same as GTF folder
            'fastq_folder=s'         => \$fastq_folder,
            'kallisto_out_folder=s'  => \$kallisto_out_folder,
            'ens_release=s'          => \$ens_release,
            'ens_metazoa_release=s'  => \$ens_metazoa_release,
            'data_host=s'            => \$data_host,
            'data_login=s'           => \$data_login,
            'enc_passwd_file=s'      => \$enc_passwd_file,
            'vit_kallisto_cmd=s'     => \$vit_kallisto_cmd,
            'vit_R_cmd=s'            => \$vit_R_cmd,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $sample_info_file eq '' || $output_log_folder eq '' || $index_folder eq '' || $fastq_folder eq '' || $kallisto_out_folder eq '' || $ens_release eq '' || $ens_metazoa_release eq '' || $data_host eq '' || $data_login eq '' || $enc_passwd_file eq '' || $vit_kallisto_cmd eq '' || $vit_R_cmd eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -sample_info_file=\$(RNASEQ_SAMPINFO_FILEPATH) -exclude_sample_file=\$(RNASEQ_SAMPEXCLUDED_FILEPATH) -output_log_folder=\$(RNASEQ_VITALIT_LOG) -index_folder=\$(RNASEQ_VITALIT_GTF)  -fastq_folder=\$(RNASEQ_BIGBGEE_FASTQ) -kallisto_out_folder=\$(RNASEQ_VITALIT_ALL_RES) -ens_release=\$(ENSRELEASE) -ens_metazoa_release=\$(ENSMETAZOARELEASE) -data_host=\$(DATAHOST) -data_login=\$(DATA_LOGIN) -enc_passwd_file=\$(ENCRYPT_PASSWD_FILE) -vit_kallisto_cmd=\$(VIT_KALLISTO_CMD) $vit_R_cmd=\$(VIT_R_CMD)
\t-sample_info_file       rna_seq_sample_info.txt
\t-exclude_sample_file    rna_seq_sample_excluded.txt
\t-output_log_folder      folder for .out and .err files (produced by LSF system), and .Rout files produced by R
\t-index_folder=s         Folder with Kallisto indexes (same as GTF folder)
\t-fastq_folder=s         Folder with Fastq files on big bgee
\t-kallisto_out_folder=s  Folder with Kallisto output and results
\t-ens_release=s          Ensembl release
\t-ens_metazoa_release=s  Ensembl Metazoa release
\t-data_host=s            Bigbgee machine with RNA-seq fastq files
\t-data_login=s           Login for bigbgee
\t-enc_passwd_file=s      File with password necessary to decrypt the GTEx data
\t-vit_kallisto_cmd=s     Command to load kallisto module on vital-it
\t-vit_R_cmd=s            Command to load R module on vital-it
\n";
    exit 1;
}

# Tests
die "Invalid or missing [$sample_info_file]: $?\n"  if ( !-e $sample_info_file || !-s $sample_info_file );
# Create output folder if not present
make_path "$output_log_folder",   {verbose=>0, mode=>0775};

# Setting up bsub parameters #################################
my $main_script = $RealBin.'/rna_seq_mapping_and_analysis.pl';
## TODO launch bsub_scheduler.pl from /data/ul/dee/bgee/GIT/pipeline/RNA_Seq/
##        Beware that git pull command should be executed before
##        kallisto_out_folder should be on /scratch/cluster/monthly. if too slow, consider using /scratch/local/ + cp of results file to /scratch/cluster/monthly/ or /home/bbgee, or /data/ (read-only, but should be fine via scp)

# kallisto is no multithreaded unless bootstraps are used
my $nr_processors = 1;
# memory limit: should match the limit of the queue
my $memory_limit  = 160_000_000;
# RAM needed: 4GB should be enough
my $memory_usage  = 4000;
my $queue         = 'bgee';
my $user_email    = 'bgee@sib.swiss'; # for email notification

# my $jobs_during_day   = 250; # Number of simultaneous jobs during working days
# my $jobs_during_night = 300; # Number of simultaneous jobs during week-end & night
# But this is a lot too much for the IO on bigbgee by ssh. Let's limit for now
my $jobs_during_day   = 10;
my $jobs_during_night = 10;

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
    if ( `bjobs -w | grep ' $library_id '` =~ / (RUN|PEND) / ){
        print "\n$library_id not launched because it is currently being analyzed (see bjobs)\n";
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

    my $report_file = $output_log_folder.'/'.$library_id.'/'.$library_id.'.report';

    # remove ending ";" at end of module loading commands, which can mess up job submission command line
    $vit_kallisto_cmd =~ s/\;$//;
    $vit_R_cmd =~ s/\;$//;
    # espace spaces so that the arguments are not cut after first space by GetOpt
    ##$vit_kallisto_cmd =~ s/\s/\\ /g;
    ##$vit_R_cmd =~ s/\s/\\ /g;


    my $script_plus_args = "time perl $main_script -library_id=$library_id -sample_info_file=$sample_info_file -exclude_sample_file=$exclude_sample_file -index_folder=$index_folder -fastq_folder=$fastq_folder -kallisto_out_folder=$kallisto_out_folder -output_log_folder=$output_log_folder -ens_release=$ens_release -ens_metazoa_release=$ens_metazoa_release -data_host=$data_host -data_login=$data_login -enc_passwd_file=$enc_passwd_file -vit_kallisto_cmd=\\\"$vit_kallisto_cmd\\\" -vit_R_cmd=\\\"$vit_R_cmd\\\"";

    # Adjust number of jobs to time and day
    my $job_limit = $jobs_during_day;
    # Get current date and time to set $job_limit upper if during nights or week-end days.
    # Get Week Day
    if ( localtime()->[6]==0 || localtime()->[6]==6 ){
        $job_limit = $jobs_during_night;
        # More jobs at the same time if current day is Sunday (0) or Saturday (6)
    }
    # Get Hour
    else {
        my $hour = localtime()->[2];
        $job_limit = $jobs_during_night  if ( $hour<=6 || $hour>=19 );
        # More jobs at the same time if during night, from 19h00 to 06h59
    }

    # Wait for free places in job queue
    my $running_jobs = check_running_jobs();
    WAIT_FREE_JOB_IN_QUEUE:
    while ( $running_jobs >= $job_limit ){
        print "No more possible slot for the job, waiting and resubmitting\n";
        sleep 30;
        $running_jobs = check_running_jobs();
    }


    # Script can be launched! Construct bsub command:

    # First, remove previous .out and .err files
    my $bsub_command = $rm_output_command.$rm_error_command;
    # Add main bsub command and options
    # Potential other options:
    # -R span[ptile=$nr_processors]
    # -u $user_email

    $bsub_command .= "echo $script_plus_args | bsub -n $nr_processors -M $memory_limit -R rusage[mem=$memory_usage] -o $output_file -e $error_file -q $queue -J \"$library_id\";";
    print "Command submitted to Vital-IT:\n$bsub_command\n";
    # also print the command on .out file
    open (my $OUT, '>', "$report_file")  or die "Cannot write [$report_file]\n";
    print {$OUT} "Vital-IT / Bsub command submitted:\n$bsub_command\n";
    close $OUT;

    # Then, run the job
    system(". /mnt/common/lsf/conf/profile.lsf; $bsub_command")==0
        or print "Failed to submit job [$library_id]\n";
}

print "\n######################################################\nAll done. $count jobs submitted.\n######################################################\n";
exit 0;


sub check_running_jobs {
    my $running_jobs = `. /mnt/common/lsf/conf/profile.lsf; bjobs | grep -v 'JOBID' | wc -l` || 0;
    chomp($running_jobs);
    return $running_jobs;
}

