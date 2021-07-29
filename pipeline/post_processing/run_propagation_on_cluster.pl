#!/usr/bin/env perl

## Julien Wollbrett, April 6, 2021
# This script generates sbatch files to run on cluster with slurm queuing system
# It is possible to directly run all jobs with the parameter "-run_jobs".
# A bash script called "run_all_jobs.sh" is created at the same location than sbatch scripts.
# It is possible to run this bash script to run all jobs

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Tie::IxHash;
use Getopt::Long;

# Define arguments & their default value
my ($jar_path, $output_dir, $output_cluster_dir, $bgee_pwd) = ('', '', '', '','');
my $run_jobs;
my %opts = ('jar_path=s'            => \$jar_path,
            'output_dir=s'          => \$output_dir,
            'output_cluster_dir=s'  => \$output_cluster_dir,
            'bgee_pwd=s'            => \$bgee_pwd,        
            'run_jobs'              => \$run_jobs
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);

if ( !$test_options || $jar_path eq '' || $output_dir eq ''|| $output_cluster_dir eq '' || $bgee_pwd eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -jar_path=\$(PATH_TO_JAR) -output_dir=\$(PATH_TO_OUTPUT) -output_cluster_dir=\$(OUTPUT_CLUSTER_DIR)
\t-jar_path           path to jar with all dependancies (somewhere in your home directory of the cluster)
\t-output_dir         path to the directory (somewhere in our home directory of the cluster) where all
\t                    sbatch files and the bash file allowing to run all jobs will be created
\t-output_cluster_dir path to the directory where log files should be written on the cluster
\t                    (!! Be sure this path exists !!)
\t-bgee_pwd           password to connect to bgee database
\t-run_jobs           boolean allowing to directly run all jobs if present. If not present
\t                    a bash file allowing to run all jobs at the same time is created
\n";
    exit 1;
}

# Setting up SLURM parameters #################################
my $partition      = 'cpu';
my $account        = 'mrobinso_bgee';
my $nbr_processors = 6;
# RAM needed: 10GB should be enough
my $memory_usage   = 29;      # in GB
my $time_limit     = '3-00:00:00';
my $log_prefix     = 'generatePropagatedCalls_';

# In the future these info should be retrieved from the database (with a query retrieving number
# of genes per species IDs provided as input) the query will also retrieve the number of globalCond
# per species and adapt the number of genes per job accordindly

my $serveur_url = "rbioinfo.unil.ch";
tie my %species_to_gene_number, 'Tie::IxHash';
%species_to_gene_number = (
6239 => 26529,
7237 => 15820,
7240 => 15047,
7740 => 31630,
7897 => 14614,
7918 => 21861,
7936 => 29144,
7994 => 26360,
8010 => 26351,
8030 => 47473,
8049 => 20691,
8081 => 21371,
8090 => 24044,
8154 => 25404,
8355 => 34459,
9031 => 21655,
9103 => 14715,
9258 => 18699,
9483 => 25428,
9531 => 24567,
9541 => 26482,
9544 => 34471,
9545 => 25864,
9555 => 28226,
9593 => 24883,
9597 => 22923,
9598 => 32286,
9615 => 29676,
9685 => 27238,
9796 => 29504,
9823 => 31055,
9925 => 23583,
9940 => 25959,
9974 => 21412,
9986 => 26392,
10116 => 29580,
10141 => 24342,
10181 => 28022,
13616 => 33507,
28377 => 24793,
30608 => 24190,
32507 => 22221,
52904 => 20656,
60711 => 27427,
69293 => 21148,
105023 => 22021);

# if we retrieve data from the database, $gene_row_count should be calculated
# depending on the number of global conditions
my $gene_row_count = 1000;
#my $species_order = 1;

# bash file containing all sbatch to run
my $bash_file = "${output_dir}/run_all_jobs.sh";
open(my $bash_file_handler, '>', $bash_file) or die $!;
print $bash_file_handler "#!/usr/bin/env bash\n";
foreach my $species_id (keys %species_to_gene_number){
    my $gene_number = $species_to_gene_number{$species_id};
    my $gene_offset = 0;

    #create a new sbatch file
    while ($gene_offset < $gene_number) {
        my $file_name = "${output_dir}/${log_prefix}${species_id}_${gene_offset}.sbatch";
        my $cluster_file_name = "${output_cluster_dir}/${log_prefix}${species_id}_${gene_offset}.sbatch";
        open(my $file_handler, '>', $file_name) or die $!;
        my $job_name = "propCalls_${species_id}_${gene_offset}";
        # create template of the sbatch file
        my $output_file = "${output_cluster_dir}${log_prefix}${species_id}_${gene_offset}.out";
        my $error_file = "${output_cluster_dir}${log_prefix}${species_id}_${gene_offset}.err";
        my $template = sbatch_template($partition, $account, $time_limit, $nbr_processors,
            $memory_usage, $output_file, $error_file, $jar_path, $job_name);
        if($gene_offset + $gene_row_count <= $gene_number) {
            $template .= create_java_command($jar_path, $output_dir, $species_id, $gene_offset,
                $gene_row_count, $memory_usage, $bgee_pwd, $serveur_url);
            $gene_offset = $gene_offset + $gene_row_count;
        }elsif ($gene_offset + $gene_row_count > $gene_number) {
            my $genes_remaining = $gene_number - $gene_offset;
            $template .= create_java_command($jar_path, $output_dir, $species_id, $gene_offset,
                $genes_remaining, $memory_usage, $bgee_pwd, $serveur_url);
            $gene_offset = $gene_number;
        }
        print $file_handler $template;
        close($file_handler);
        if($run_jobs) {
            # Then, run the job
            system("sbatch $file_name")==0  or print "Failed to submit job for species".
            "[$species_id] and gene offset [$gene_offset]\n";
        }
        print $bash_file_handler "sbatch $cluster_file_name\n";
    }
    #$species_order = $species_order + 1;
}
close($bash_file_handler);
exit 0;

sub create_java_command {
    my ($jar_path, $output_dir, $species_id, $gene_offset, $gene_row_count, $memory_usage, $bgee_pwd, $serveur_url) = @_;
    my $template = "java -Xmx${memory_usage}g -Dbgee.dao.jdbc.username=root -Dbgee.dao.jdbc.password=${bgee_pwd} -Dbgee.dao.jdbc.driver.names=com.mysql.jdbc.Driver,net.sf.log4jdbc.sql.jdbcapi.DriverSpy -Dbgee.dao.jdbc.url='jdbc:log4jdbc:mysql://${serveur_url}:3306/bgee_v15_0?useSSL=false&enableQueryTimeouts=false&sessionVariables=net_write_timeout=260000,net_read_timeout=260000,wait_timeout=260000' -jar \${JAR_PATH}bgee-pipeline-15.0-with-dependencies.jar InsertPropagatedCalls insertCalls $species_id - $gene_offset $gene_row_count 0";
}

# Add main sbatch command and options
sub sbatch_template {
    my ($partition, $account, $time_limit, $nbr_processors, $memory_usage, $output_file, $error_file, $jar_path, $job_name) = @_;

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
module use /software/module/
module add Development/java/1.8.0_242;

export JAR_PATH=$jar_path

";

    return $template;
}

