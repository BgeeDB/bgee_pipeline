#!/usr/bin/env perl

## Julien Wollbrett, April 12, 2021
# This script generates sbatch files to run on cluster with slurm queuing system
# It is possible to directly run all jobs with the parameter "-run_jobs".
# A bash script called "generate_rnaseq_ranks_jobs.sh" is created at the same location than sbatch scripts.
# It is possible to run this bash script to run all jobs

# Perl core modules
use strict;
use warnings;
use diagnostics;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

use Getopt::Long;

# Define arguments & their default value
my ($pipeline_cluster_dir, $script_relative_path, $output_dir, $output_cluster_dir, $bgee_pwd) = ('', '', '', '', '');
my %opts = ('output_dir=s'           => \$output_dir,
            'output_cluster_dir=s'   => \$output_cluster_dir,,
            'pipeline_cluster_dir=s' => \$pipeline_cluster_dir,
            'script_relative_path=s' => \$script_relative_path
            'bgee_pwd=s'             => \$bgee_pwd,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $pipeline_cluster_dir eq '' || $script_relative_path eq '' || $output_dir eq ''|| $output_cluster_dir eq ''
    || $bgee_pwd eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -jar_path=\$(PATH_TO_JAR) -output_dir=\$(PATH_TO_OUTPUT) -output_cluster_dir=\$(OUTPUT_CLUSTER_DIR)
\t-pipeline_cluster_dir path to Bgee pipeline directory
\t-script_relative_path relative path of the script ranks_rnaseq.pl from the root of bgee pipeline
\t-output_dir           path to the directory (somewhere in our home directory of the cluster) where all
\t                      sbatch files and the bash file allowing to run all jobs will be created
\t-bgee_pwd             password to connect to bgee database
\t-output_cluster_dir   path to the directory where log files should be written on the cluster
\t                      (!! Be sure this path exists !!)
\n";
    exit 1;
}

# Setting up SLURM parameters #################################
my $partition       = 'cpu';
my $account         = 'mrobinso_bgee';
my $bgee_db_version = 'bgee_v15_0';
my $bgee_user       = 'root';
my $bgee_pwd        = 'BdD.pw';
my $bgee_port       = 3306;
my $nbr_processors  = 1;
my $libs_per_thread = 10;
my $memory_usage    = 5;      # in GB
my $time_limit      = '12:00:00';
my $log_prefix      = 'generateRnaSeqRanks_';
my $serveur_url     = 'rbioinfo.unil.ch';

my $bgee_connector= get_bgee_connector($bgee_user, $serveur_url, $bgee_pwd, $bgee_port, $bgee_db_version);

# Connect to Bgee DB to retrieve all libraries for which no ranks have
# been processed for now.
my @remainingLibraries;
my $dbh = Utils::connect_bgee_db($bgee_connector);
my $libSql = 'SELECT t1.rnaSeqLibraryId FROM rnaSeqLibrary AS t1
              WHERE EXISTS (SELECT 1 FROM rnaSeqResult AS t2
                  WHERE t1.rnaSeqLibraryId = t2.rnaSeqLibraryId
                  AND t2.expressionId IS NOT NULL
              ) AND NOT EXISTS (SELECT 1 FROM rnaSeqResult AS t2
                  WHERE t1.rnaSeqLibraryId = t2.rnaSeqLibraryId
                  AND t2.rank IS NOT NULL)';
my $rnaSeqLibStmt = $dbh->prepare($libSql);
$rnaSeqLibStmt->execute()  or die $rnaSeqLibStmt->errstr;
while ( my @data = $rnaSeqLibStmt->fetchrow_array ){
    push(@remainingLibraries, $data[0]);
}
$rnaSeqLibStmt->finish;

my $number_rnaSeq_libraries = @remainingLibraries;
my $number_libraries_per_job = 10;




# bash file containing all sbatch to run
my $bash_file = "${output_dir}/run_all_jobs.sh";
open(my $bash_file_handler, '>', $bash_file) or die $!;
print $bash_file_handler "#!/usr/bin/env bash\n";

my $libraries_offset = 0;

my $job_number = 0;
#create a new sbatch file
while (@remainingLibraries != 0) {
    my $file_name = "${output_dir}/${log_prefix}${libraries_offset}.sbatch";
    my $cluster_file_name = "${output_cluster_dir}/${log_prefix}${libraries_offset}.sbatch";
    open(my $file_handler, '>', $file_name) or die $!;
    my $job_name = "ranks_${job_number}";
    # create template of the sbatch file
    my $output_file = "${output_cluster_dir}${log_prefix}${libraries_offset}.out";
    my $error_file = "${output_cluster_dir}${log_prefix}${libraries_offset}.err";
    my $template = sbatch_template($partition, $account, $time_limit, $nbr_processors,
        $memory_usage, $output_file, $error_file, $pipeline_cluster_dir, $job_name);

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
    $template .= create_perl_command($nbr_processors, $pipeline_cluster_dir, $script_relative_path,
        $bgee_connector, $libs_per_thread, $libraries_parameter);

    print $file_handler $template;
    close($file_handler);
    print $bash_file_handler "sbatch $cluster_file_name\n";
    $job_number++;
}
close($bash_file_handler);
exit 0;

sub create_perl_command {
    my ($nbr_processors, $pipeline_cluster_dir, $script_relative_path, $bgee_connector,
        $libs_per_thread, $lib_ids) = @_;

    # the ranks_rnaseq script has one master thread to process data except if
    # only 1 thread is asked.
    my $nbr_threads = $nbr_processors;
    if ($nbr_processors > 1) {
         $nbr_threads--;
    }
    my $template = "perl ${pipeline_cluster_dir}${script_relative_path}ranks_rnaseq.pl -bgee=${bgee_connector} -parallel_jobs=${nbr_threads} -libs_per_job=$libs_per_thread -lib_ids=${lib_ids}";
    return $template;
}

# Add main sbatch command and options
sub sbatch_template {
    my ($partition, $account, $time_limit, $nbr_processors, $memory_usage, $output_file, $error_file, $pipeline_path, $job_name) = @_;

    my $template="#!/bin/bash
#SBATCH --partition=$partition
#SBATCH --account=$account
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$nbr_processors
#SBATCH --mem=${memory_usage}G
#SBATCH --time=$time_limit
#SBATCH --output=$output_file
#SBATCH --error=$error_file
#SBATCH --export=NONE
#SBATCH --job-name=$job_name
module add Bioinformatics/Software/vital-it;
module add Development/Ensembl_API/97;

export PIPELINE_PATH=$pipeline_path

";

    return $template;
}

sub get_bgee_connector {
    my ($bgee_user, $serveur_url, $bgee_pwd, $bgee_port, $bgee_db_version) = @_;
    my $bgee_cmd = "user=${bgee_user}__pass=${bgee_pwd}__host=${serveur_url}__port=${bgee_port}__name=${bgee_db_version}";
    return $bgee_cmd;
}

sub retrieve_all_librari

