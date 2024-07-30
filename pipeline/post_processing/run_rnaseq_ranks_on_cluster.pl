#!/usr/bin/env perl

## Julien Wollbrett, April 12, 2021
## Update Frederic Bastian, April 10 2023: add parameter to specify single-cell or not
# This script generates sbatch files to run on cluster with slurm queuing system
# A bash script called "generate_rnaseq_ranks_jobs.sh" is created at the same location than sbatch scripts.
# It is possible to run this bash script to run all jobs

# TODO: integrate in the pipeline by creating a rule in the Makefile
# TODO: create a variable in the Makefile.common to store all slurm information

# Perl core modules
use strict;
use warnings;
use diagnostics;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

use Getopt::Long;

# Define arguments & their default value
my ($script_path, $output_cluster_dir, $bgee_connector) = ('', '', '',  '', '');
my ($queue, $account) = ('', '');
my $is_single_cell = -1;
my $parallelJobs = 1;
my %opts = ('bgee=s'                 => \$bgee_connector,
            'output_cluster_dir=s'   => \$output_cluster_dir,,
            'script_path=s'          => \$script_path,
            'is_single_cell=i'       => \$is_single_cell,
            'queue=s'                => \$queue,
            'account=s'              => \$account,
            'parallelJobs=i'         => \$parallelJobs
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options ||  $script_path eq '' || $output_cluster_dir eq ''
    || $bgee_connector eq '' || $queue eq '' || $parallelJobs == 1 || $account eq '' || ($is_single_cell != 0 && $is_single_cell != 1)){
    print "\n\tInvalid or missing argument:
\t-bgee                 connector allowing to connect to Bgee database
\t-script_path          Absolute path to the directory containing the ranks generation script
\t-output_cluster_dir   path to the directory where sbatch and log files will be written on the cluster
\t-is_single_cell       0: target bulk RNA-Seq libraries; 1: target single-cell RNA-Seq libraries
\t-queue                partition to use to run jobs on the cluster
\t-account              account to use to run jobs on the cluster
\t-parallelJobs         max number of jobs to run in parallel
\n";
    exit 1;
}


####################### FUNCTIONS #######################

sub create_perl_command {
    my ($nbr_processors, $script_path, $bgee_connector,
        $libs_per_thread, $lib_ids) = @_;

    # the ranks_rnaseq script has one master thread to process data except if
    # only 1 thread is asked.
    my $nbr_threads = $nbr_processors;
    if ($nbr_processors > 1) {
         $nbr_threads--;
    }
    my $template = "perl $script_path -bgee=${bgee_connector} -parallel_jobs=${nbr_threads} -libs_per_job=$libs_per_thread -lib_ids=${lib_ids}";
    return $template;
}

####################### Setting up SLURM parameters #######################
my $nbr_processors  = 1;
my $libs_per_thread = 2;
my $memory_usage    = 5;      # in GB

my $script     = 'ranks_rnaseq.pl';
my $log_prefix = 'generateRnaSeqRanks_';
if ($is_single_cell == 1) {
    $script     = 'ranks_scrnaseq.pl';
    $log_prefix = 'generateScRnaSeqRanks_';
}

#concatenate script_path and the name of the script to use
$script_path = "${script_path}/$script";

# Connect to Bgee DB to retrieve all libraries for which no ranks have
# been processed for now.
my @remainingLibraries;
my $dbh = Utils::connect_bgee_db($bgee_connector);
my $libSql = 'SELECT t1.rnaSeqLibraryId FROM rnaSeqLibrary AS t1
              WHERE rnaSeqTechnologyIsSingleCell = ?
              AND EXISTS (SELECT 1 FROM rnaSeqLibraryAnnotatedSample AS t2
                  INNER JOIN rnaSeqLibraryAnnotatedSampleGeneResult AS t3
                  ON t3.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId
                  WHERE t1.rnaSeqLibraryId = t2.rnaSeqLibraryId
                  AND t3.expressionId IS NOT NULL
              ) AND NOT EXISTS (SELECT 1 FROM rnaSeqLibraryAnnotatedSample AS t2
                  INNER JOIN rnaSeqLibraryAnnotatedSampleGeneResult AS t3
                  ON t3.rnaSeqLibraryAnnotatedSampleId = t2.rnaSeqLibraryAnnotatedSampleId
                  WHERE t1.rnaSeqLibraryId = t2.rnaSeqLibraryId
                  AND t3.rawRank IS NOT NULL)';
my $rnaSeqLibStmt = $dbh->prepare($libSql);
$rnaSeqLibStmt->execute($is_single_cell == 1 ? '1' : '0')  or die $rnaSeqLibStmt->errstr;
while ( my @data = $rnaSeqLibStmt->fetchrow_array ){
    push(@remainingLibraries, $data[0]);
}
$rnaSeqLibStmt->finish;

my $number_libraries_per_job = 10;

# list of sbatch files to run
my @sbatchCommands = ();

my $libraries_offset = 0;

my $job_number = 0;
#create a new sbatch file
while (@remainingLibraries != 0) {
    my $jobName = "${log_prefix}${libraries_offset}";
    my $file_name = "${output_cluster_dir}/sbatch/$jobName.sbatch";
    open(my $file_handler, '>', $file_name) or die $!;
    # create template of the sbatch file
    my $output_file = "${output_cluster_dir}/logs/${log_prefix}${libraries_offset}.out";
    my $error_file = "${output_cluster_dir}/logs/${log_prefix}${libraries_offset}.err";
    my $template = Utils::sbatch_template($queue, $account, 1,
          50, "${clusterOutput}${jobName}.out", "${clusterOutput}/${jobName}.err",
          $jobName);
    $template .= "export PATH=/software/bin:\$PATH;\n";
    my $libraries_parameter = '';
    for my $i (1..$number_libraries_per_job){
        if(@remainingLibraries != 0) {
            if($i == 1) {
                $libraries_parameter .= shift(@remainingLibraries);
            } else {
                $libraries_parameter .= ",".shift(@remainingLibraries);
            }
            $libraries_offset++;
        } else {
            $i = $number_libraries_per_job;
        }
    }
    $template .= create_perl_command($nbr_processors, $script_path,
        $bgee_connector, $libs_per_thread, $libraries_parameter);

    print $file_handler $template;
    close($file_handler);
    push @sbatchCommands "sbatch $file_name";
    $job_number++;
}

print "created $job_number sbatch files.\n";

my $numberJobRun = 0;
my $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $log_prefix);
my $startTime = localtime->strftime('%Y-%m-%dT%H:%M:%S');
foreach my $sbatchCommand (@sbatchCommands) {
	while ($jobsRunning >= $parallelJobs) {
        sleep(15);
        $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $log_prefix);
    }
    #system("$sbatchCommand >/dev/null");
    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $log_prefix);
	$numberJobRun++;
    }
}

while ($jobsRunning > 0) {
    sleep(15);
    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $log_prefix);
}

print "all jobs finished. Run $numberJobRun jobs\n";

my %jobs_status = Utils::count_status_finished_jobs($log_prefix, $startTime);
print $jobs_status{"completed"}." jobs completed, ".$jobs_status{"failed"}." jobs failed, ".$jobs_status{"out_of_memory"}." jobs failed with an out of memory issue and ".$jobs_status{"cancelled"}." jobs have been cancelled.\n";

print "all jobs finished\n";

exit 0;

