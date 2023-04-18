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
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my ($jar_path, $output_dir, $sbatch_cluster_dir, $log_cluster_dir, $database_name, $bgee_pwd, $bgee_species) = ('', '', '', '','', '', '');
my $run_jobs;
my %opts = ('jar_path=s'            => \$jar_path,
            'output_dir=s'          => \$output_dir,
            'log_cluster_dir=s'     => \$log_cluster_dir,
            'sbatch_cluster_dir=s'  => \$sbatch_cluster_dir,
            'database_name=s'       => \$database_name,
            'bgee_pwd=s'            => \$bgee_pwd,
            'bgee_species=s'        => \$bgee_species,
            'run_jobs'              => \$run_jobs,
           );

# Check arguments
my $emptyArg = '-';

my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $jar_path eq '' || $bgee_species eq '' || $output_dir eq ''|| $sbatch_cluster_dir eq '' || $log_cluster_dir eq '' || $database_name eq '' || $bgee_pwd eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -jar_path=\$(PATH_TO_JAR) -bgee_species='-' -output_dir=\$(PATH_TO_OUTPUT) -output_dir=\$(PATH_TO_OUTPUT) -sbatch_cluster_dir=\$(SBATCH_CLUSTER_DIR) -log_cluster_dir=\$LOG_CLUSTER_DIR)
\t-jar_path           path to jar with all dependancies (somewhere in your home directory of the cluster)
\t-output_dir         path to the directory where all sbatch files and the bash file will be created. This parameter allow to
                      generate sbatch files not on the cluster and move them after. If run on the cluster, then output_dir and
                      sbatch_cluster_dir should have the same value
\t-sbatch_cluster_dir path to the directory where sbatch files should be written on the cluster
\t                    (!! Be sure this path exists !!)
\t-log_cluster_dir    path to the directory where log files should be written on the cluster
\t                    (!! Be sure this path exists !!)
\t-database_name      name of the database (e.g bgee_v15_0)
\t-bgee_pwd           password to connect to bgee database
\t-bgee_species       list of species propagation has to be run on or '-' if all species
\t-run_jobs           boolean allowing to directly run all jobs if present. If not present
\t                    a bash file allowing to run all jobs at the same time is created
\n";
    exit 1;
}

# Setting up SLURM parameters #################################
my $partition      = 'cpu';
my $account        = 'mrobinso_bgee';
my $nbr_processors = 6;
my $memory_usage   = 80;      # in GB
my $time_limit     = '1-00:00:00'; #in days
my $log_prefix     = 'generatePropagatedCalls_';

# In the future these info should be retrieved from the database (with a query retrieving number
# of genes per species IDs provided as input) the query will also retrieve the number of globalCond
# per species and adapt the number of genes per job accordindly
# TODO : REALLY NEED TO RETRIEVE THESE INFO FROM THE DATABASE WITH THE QUERY :
# select speciesId, count(distinct bgeeGeneId) as geneNumber from cond as t1 inner join expression as t2 on t1.conditionId = t2.conditionIdg roup by speciesId order by speciesId;

my $serveur_url = 'localhost';
my $bgee_user   = 'root';
my $bgee_port   = 3306;


my $bgee_connector= get_bgee_connector($bgee_user, $serveur_url, $bgee_pwd, $bgee_port, $database_name);
my $dbh = Utils::connect_bgee_db($bgee_connector);

# if no species provided then retrieve all species from the database_name
# species list
my @speciesList = ();
if ($bgee_species ne $emptyArg) {
    @speciesList = split(',', $bgee_species);
}

tie my %species_to_gene_number, 'Tie::IxHash';

my $genesSql = 'select speciesId, count(distinct bgeeGeneId) as geneNumber
              from cond as t1 inner join expression as t2 on t1.conditionId = t2.conditionId';
if (@speciesList) {
    $genesSql .= ' WHERE speciesId IN (';
    for my $i (0 .. $#speciesList) {
        if ($i > 0) {
            $genesSql .= ', ';
        }
        $genesSql .= $speciesList[$i];
    }
    $genesSql .= ')';
}

$genesSql .= ' group by speciesId order by speciesId';

my $genesStmt = $dbh->prepare($genesSql);
$genesStmt->execute()  or die $genesStmt->errstr;
while ( my @data = $genesStmt->fetchrow_array ){
    $species_to_gene_number{$data[0]} = $data[1];
}
$genesStmt->finish;

# $gene_row_count could be calculated depending on the number of global conditions
my $gene_row_count = 50;

# bash file containing all sbatch to run
my $bash_file = "${output_dir}/run_all_jobs.sh";
open(my $bash_file_handler, '>', $bash_file)  or die $!;
print $bash_file_handler "#!/usr/bin/env bash\n";
foreach my $species_id (keys %species_to_gene_number){
    my $gene_number = $species_to_gene_number{$species_id};
    my $gene_offset = 0;

    #create a new sbatch file
    while ($gene_offset < $gene_number) {
        my $file_name = "${output_dir}/${log_prefix}${species_id}_${gene_offset}.sbatch";
        my $cluster_file_name = "${sbatch_cluster_dir}/${log_prefix}${species_id}_${gene_offset}.sbatch";
        open(my $file_handler, '>', $file_name)  or die $!;
        my $job_name = "propCalls_${species_id}_${gene_offset}";
        # create template of the sbatch file
        my $output_file = "${log_cluster_dir}${log_prefix}${species_id}_${gene_offset}.out";
        my $error_file = "${log_cluster_dir}${log_prefix}${species_id}_${gene_offset}.err";
        my $template = sbatch_template($partition, $account, $time_limit, $nbr_processors,
            $memory_usage, $output_file, $error_file, $jar_path, $job_name);
        if($gene_offset + $gene_row_count <= $gene_number) {
            $template .= create_java_command($jar_path, $output_dir, $species_id, $gene_offset,
                $gene_row_count, $memory_usage, $bgee_pwd, $serveur_url, $bgee_user, $bgee_port, $nbr_processors);
            $gene_offset = $gene_offset + $gene_row_count;
        }elsif ($gene_offset + $gene_row_count > $gene_number) {
            my $genes_remaining = $gene_number - $gene_offset;
            $template .= create_java_command($jar_path, $output_dir, $species_id, $gene_offset,
                $genes_remaining, $memory_usage, $bgee_pwd, $serveur_url, $bgee_user, $bgee_port, $nbr_processors);
            $gene_offset = $gene_number;
        }
        print $file_handler $template;
        close($file_handler);
        if($run_jobs) {
            # Then, run the job
            system("sbatch $file_name")==0  or print 'Failed to submit job for species'.
            "[$species_id] and gene offset [$gene_offset]\n";
        }
        print $bash_file_handler "sbatch ${cluster_file_name}\n";
    }
    #$species_order = $species_order + 1;
}
close($bash_file_handler);
exit 0;

sub create_java_command {
    my ($jar_path, $output_dir, $species_id, $gene_offset, $gene_row_count, $memory_usage, 
        $bgee_pwd, $serveur_url, $bgee_user, $bgee_port, $nbr_processors) = @_;
    my $template = "java -Xmx${memory_usage}g -Djava.util.concurrent.ForkJoinPool.common.parallelism=${nbr_processors} -Dbgee.dao.jdbc.username=${bgee_user} -Dbgee.dao.jdbc.password=${bgee_pwd} -Dbgee.dao.jdbc.driver.names=com.mysql.jdbc.Driver,net.sf.log4jdbc.sql.jdbcapi.DriverSpy -Dbgee.dao.jdbc.url='jdbc:log4jdbc:mysql://${serveur_url}:${bgee_port}/${database_name}?useSSL=false&enableQueryTimeouts=false&sessionVariables=net_write_timeout=520000,net_read_timeout=520000,wait_timeout=520000' -jar \${JAR_PATH} InsertPropagatedCalls insertCalls $species_id - $gene_offset $gene_row_count 0";
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
module load openjdk/11.0.15_10;

export JAR_PATH=$jar_path

";

    return $template;
}

sub get_bgee_connector {
    my ($bgee_user, $serveur_url, $bgee_pwd, $bgee_port, $database_name) = @_;
    my $bgee_cmd = "user=${bgee_user}__pass=${bgee_pwd}__host=${serveur_url}__port=${bgee_port}__name=${database_name}";
    return $bgee_cmd;
}
