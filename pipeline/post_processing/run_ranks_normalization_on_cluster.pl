#!/usr/bin/env perl

## Julien Wollbrett, April 12, 2021
# This script generates sbatch files to run ranks normalization on cluster with slurm queuing system
# It is possible to directly run all jobs with the parameter "-run_jobs".
# A bash script called "run_all_normalization.sh" is created at the same location than sbatch scripts.
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
my ($pipeline_cluster_dir, $script_relative_path, $output_dir, $output_cluster_dir, $bgee_connector, $bgee_species) = ('', '', '', '', '', '');
my %opts = ('output_dir=s'           => \$output_dir,
            'output_cluster_dir=s'   => \$output_cluster_dir,,
            'pipeline_cluster_dir=s' => \$pipeline_cluster_dir,
            'bgee=s'                 => \$bgee_connector,
            'script_relative_path=s' => \$script_relative_path,
            'bgee_species=s'         => \$bgee_species
           );

# Check arguments
my $emptyArg = '-';

my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $pipeline_cluster_dir eq '' || $script_relative_path eq '' || $output_dir eq ''|| $output_cluster_dir eq ''
    || $bgee_connector eq '' || $bgee_species eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -output_dir=\$(PATH_TO_OUTPUT) -output_cluster_dir=\$(OUTPUT_CLUSTER_DIR) -pipeline_cluster_dir=\"/path_to_pipeline/\" -script_relative_path=\"pipeline/post_processing/\" -bgee=\"user=\${bgee_user}__pass=\${bgee_pwd}__host=\${serveur_url}__port=\${bgee_port}__name=\${database_name}\" -bgee_species=\"-\"
\t-pipeline_cluster_dir path to Bgee pipeline directory
\t-script_relative_path relative path of the script ranks_global.pl from the root of bgee pipeline
\t-output_dir           path to the directory (somewhere in our home directory of the cluster) where all
\t                      sbatch files and the bash file allowing to run all jobs will be created
\t-bgee                 Bgee connector
\t-bgee_species         species for which ranks normalization has to be run
\t-output_cluster_dir   path to the directory where log files should be written on the cluster
\t                      (!! Be sure this path exists !!)
\n";
    exit 1;
}

# Setting up SLURM parameters #################################
my $partition       = 'cpu';
my $account         = 'mrobinso_bgee';
my $nbr_processors  = 1;
my $memory_usage    = 5;      # in GB
my $time_limit      = '3-00:00:00';
my $log_prefix      = 'nr_';

my @speciesList = ();
if ($bgee_species ne $emptyArg) {
    @speciesList = split(',', $bgee_species);
} else {
    my $dbh = Utils::connect_bgee_db($bgee_connector);
    my $speciesSql = 'select speciesId from species';
    my $speciesStmt = $dbh->prepare($speciesSql);
    $speciesStmt -> execute()  or die $speciesStmt->errstr;
    while ( my @data = $speciesStmt -> fetchrow_array ){
        push(@speciesList, $data[0]);
    }
}

# bash file containing all sbatch to run
my $bash_file = "${output_dir}/run_all_normalization.sh";
open(my $bash_file_handler, '>', $bash_file)  or die $!;
print $bash_file_handler "#!/usr/bin/env bash\n";
foreach my $species_id (@speciesList){
    my $file_name = "${output_dir}/${log_prefix}${species_id}.sbatch";
    open(my $file_handler, '>', $file_name)  or die $!;
    my $job_name = "${log_prefix}${species_id}";
    # create template of the sbatch file
    my $output_file = "${output_cluster_dir}${log_prefix}${species_id}.out";
    my $error_file = "${output_cluster_dir}${log_prefix}${species_id}.err";
    my $template = sbatch_template($partition, $account, $time_limit, $nbr_processors,
            $memory_usage, $output_file, $error_file, $job_name);
    $template .= create_perl_command($pipeline_cluster_dir, $script_relative_path,
        $bgee_connector, $species_id);

    print $file_handler $template;
    close($file_handler);
    print $bash_file_handler "sbatch $file_name\n";
}
close($bash_file_handler);
exit 0;

sub create_perl_command {
    my ($pipeline_cluster_dir, $script_relative_path, $bgee_connector,
        $species_id) = @_;

    # the ranks_rnaseq script has one master thread to process data except if
    # only 1 thread is asked.
    my $nbr_threads = $nbr_processors;
    if ($nbr_processors > 1) {
         $nbr_threads--;
    }
    my $template = "perl ${pipeline_cluster_dir}${script_relative_path}normalize_ranks.pl -bgee=${bgee_connector} -species=${species_id}";
    return $template;
}

# Add main sbatch command and options
sub sbatch_template {
    my ($partition, $account, $time_limit, $nbr_processors, $memory_usage, $output_file, $error_file, $job_name) = @_;

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
module use /software/module/;
module add Development/Ensembl_API/97;

export PATH=/software/bin:\$PATH
";

    return $template;
}

