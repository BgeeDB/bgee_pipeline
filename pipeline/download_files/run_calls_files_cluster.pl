#!/usr/bin/env perl

# Julien Wollbrett, January 2022
# Allows to generate calls donwload files on cluster. One job is created for each species

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
my ($jar_path, $output_dir, $output_cluster_dir, $database_name, $bgee_pwd, $bgee_species) = ('', '', '', '','', '', '');
my $run_jobs;
my %opts = ('jar_path=s'            => \$jar_path,
            'output_dir=s'          => \$output_dir,
            'output_cluster_dir=s'  => \$output_cluster_dir,
            'database_name=s'       => \$database_name,
            'bgee_pwd=s'            => \$bgee_pwd,
            'bgee_species=s'        => \$bgee_species,
            'run_jobs'              => \$run_jobs,
           );

# Check arguments
my $emptyArg = '-';

my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $jar_path eq '' || $bgee_species eq '' || $output_dir eq ''|| $output_cluster_dir eq '' || $database_name eq '' || $bgee_pwd eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -jar_path=\$(PATH_TO_JAR) -bgee_species='-' -output_dir=\$(PATH_TO_OUTPUT) -output_cluster_dir=\$(OUTPUT_CLUSTER_DIR)
\t-jar_path           path to jar with all dependancies (somewhere in your home directory of the cluster)
\t-output_dir         path to the directory (somewhere in our home directory of the cluster) where all
\t                    sbatch files and the bash file allowing to run all jobs will be created
\t-output_cluster_dir path to the directory where log files should be written on the cluster
\t                    (!! Be sure this path exists !!)
\t-database_name      name of the database (e.g bgee_v15_0)
\t-bgee_pwd           password to connect to bgee database
\t-bgee_species       list of species or '-' for all species
\t-run_jobs           boolean allowing to directly run all jobs if present. If not present
\t                    a bash file allowing to run all jobs at the same time is created
\n";
    exit 1;
}

# Setting up SLURM parameters #################################
my $partition      = 'cpu';
my $account        = 'mrobinso_bgee';
my $nbr_processors = 1;
my $memory_usage   = 100;      # in GB
my $time_limit     = '2-00:00:00'; #in days
my $log_prefix     = 'cf_';
my $serveur_url    = 'rbioinfo.unil.ch';
my $bgee_user      = 'bgee';
my $bgee_port      = 3306;
my $max_jobs       = 20;

# hash containing all combination of condition parameter to process
my $organ_param  = "ANAT_ENTITY_ID";
my $stage_param  = "DEV_STAGE_ID";
my $sex_param    = "SEX_ID";
my $strain_param = "STRAIN_ID";

# for now we only generate files for anat. entities and all cond. parameters.
# If files containing other combinations have to be generated please update the
# following hash accordingly
my %cond_param_comb = (
    AE_    => $organ_param , 
    ALL_   => "$organ_param,$stage_param,$sex_param,$strain_param"
);


my $bgee_connector= get_bgee_connector($bgee_user, $serveur_url, $bgee_pwd, $bgee_port, $database_name);

# if no species provided then retrieve all species from the database_name
# species list
my @species_ids = ();
if ($bgee_species ne $emptyArg) {
    @species_ids = split(',', $bgee_species);
}

if (@species_ids == 0) {
    my $dbh = Utils::connect_bgee_db($bgee_connector);
    my $species_sql = 'select speciesId from species order by speciesId';
    my $species_stmt = $dbh->prepare($species_sql);
    $species_stmt->execute()  or die $species_stmt->errstr;
    while ( my @data = $species_stmt->fetchrow_array ){
        push(@species_ids, $data[0]);
    }
    $species_stmt->finish;
}

# bash file containing all sbatch to run
my $bash_file = "${output_dir}/run_all_jobs.sh";
open(my $bash_file_handler, '>', $bash_file)  or die $!;
print $bash_file_handler "#!/usr/bin/env bash\n";
foreach my $species_id (@species_ids){
    foreach my $cond_key (keys %cond_param_comb){
        my $file_name = "${output_dir}/${log_prefix}${species_id}_${cond_key}.sbatch";
        open(my $file_handler, '>', $file_name)  or die $!;
        my $job_name = "${log_prefix}${species_id}_${cond_key}";
        # create template of the sbatch file
        my $output_file = "${output_cluster_dir}${log_prefix}${species_id}_${cond_key}.out";
        my $error_file = "${output_cluster_dir}${log_prefix}${species_id}_${cond_key}.err";
        my $template = sbatch_template($partition, $account, $time_limit, $nbr_processors,
            $memory_usage, $output_file, $error_file, $jar_path, $job_name);
        $template .= create_java_command($jar_path, $output_cluster_dir, $species_id, $cond_param_comb{$cond_key},
            $memory_usage, $bgee_pwd, $serveur_url, $bgee_user, $bgee_port);

        print $file_handler $template;
        close($file_handler);
        if($run_jobs) {
            # Then, run the job
            system("sbatch $file_name")==0  or print 'Failed to submit job for species'.
            "[$species_id]\n";
        }
        print $bash_file_handler "sbatch $file_name\n";
    }
    #$species_order = $species_order + 1;
}
close($bash_file_handler);

my $manage_number_jobs_file = "${output_dir}/manage_number_jobs.pl";
open(my $jobs_file_handler, '>', $manage_number_jobs_file)  or die $!;
print $jobs_file_handler Utils::limit_number_jobs_cluster($max_jobs, $bash_file, $log_prefix, $account);
close($jobs_file_handler);

exit 0;

sub create_java_command {
    my ($jar_path, $output_cluster_dir, $species_id, $cond_params, $memory_usage, 
        $bgee_pwd, $serveur_url, $bgee_user, $bgee_port) = @_;

    my $expr_files = "expr_simple,expr_advanced";
    my $template = "java -Xmx${memory_usage}g -Dbgee.dao.jdbc.username=${bgee_user} -Dbgee.dao.jdbc.password=${bgee_pwd} 
                   -Dbgee.dao.jdbc.driver.names=com.mysql.jdbc.Driver,net.sf.log4jdbc.sql.jdbcapi.DriverSpy 
                   -Dbgee.dao.jdbc.url='jdbc:log4jdbc:mysql://${serveur_url}:${bgee_port}/${database_name}?useSSL=false&enableQueryTimeouts=false&sessionVariables=net_write_timeout=520000,net_read_timeout=520000,wait_timeout=520000' 
                   -jar \${JAR_PATH} GenerateBasicExprFile $species_id $expr_files $output_cluster_dir $cond_params";
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

sub get_bgee_connector {
    my ($bgee_user, $serveur_url, $bgee_pwd, $bgee_port, $database_name) = @_;
    my $bgee_cmd = "user=${bgee_user}__pass=${bgee_pwd}__host=${serveur_url}__port=${bgee_port}__name=${database_name}";
    return $bgee_cmd;
}
