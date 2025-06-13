#!/usr/bin/env perl

## Julien Wollbrett, Nov 1, 2019
# This script launches the slurm jobs creating kallisto indexes for all species.


# Perl core modules
use strict;
use warnings;
use diagnostics;
use FindBin;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;
use File::Path qw(make_path);
use Getopt::Long;
use Time::localtime;

# Define arguments & their default value
my ($transcriptome_folder, $metadata_file, $output_log_folder, $ref_intergenic_folder, $account, $partition, $short_index_length, $container_cmd) = ('', '', '', '', '', '', '', '', '', '');
my %opts = ('transcriptome_folder=s' => \$transcriptome_folder, # same as GTF folder
            'metadata_file=s'        => \$metadata_file, # same as output folder
            'output_log_folder=s'    => \$output_log_folder,
            'ref_intergenic_folder=s'=> \$ref_intergenic_folder,
            'account=s'              => \$account,
            'partition=s'            => \$partition,
            'container_cmd=s'        => \$container_cmd,
           );

my $test_options = Getopt::Long::GetOptions(%opts);

# TO IMPLEMENT
# kallisto index generation is no multithreaded
my $nbr_processors = 1;
# RAM needed: 10GB should be enough
my $memory_usage   = 40;  # in GB
my $time_limit     = '12:00:00';
my $job_prefix = 'ref_intergenic_index';

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# Info of processed libraries coming from the pipeline
my %processedLibraries = get_processed_libraries_info($metadata_file, 1);
my %speciesId_to_name;
foreach my $experimentId (keys %processedLibraries){
    foreach my $libraryId (keys %{$processedLibraries{$experimentId}}){
        $speciesId_to_name{$processedLibraries{$experimentId}{$libraryId}{'speciesId'}} = $processedLibraries{$experimentId}{$libraryId}{'speciesName'};
    }
}

my $sbatch_folder = $output_log_folder."/sbatch_ref_intergenic_index/";
make_path($sbatch_folder);
my $clusterOutput_folder = $output_log_folder."/clusterOutput_ref_intergenic_index/";
make_path($clusterOutput_folder);
foreach my $speciesId (keys %speciesId_to_name) {
    print "Species ID: $speciesId, Species Name: $speciesId_to_name{$speciesId}\n";
    # retrieve files from transcriptome_folder that start with the species name and end with .gtf_all
    my $species_name = $speciesId_to_name{$speciesId};
    # retrieve the species name replacing space by underscore
    $species_name =~ s/ /_/g;
    my @files = glob("$transcriptome_folder/${species_name}*.gtf_all");
    foreach my $file (@files) {
        if($file =~ ('.gtf_all$')) {
            my $transcriptome_ref_intergenic_index_path = $file =~ s/gtf_all/transcriptome_ref_intergenic.idx/r;
            if (-e $transcriptome_ref_intergenic_index_path || -e $transcriptome_ref_intergenic_index_path.'.xz') {
                print("index $transcriptome_ref_intergenic_index_path already exists for species $species_name");
		next;
            }

            my $jobName = "${job_prefix}_${speciesId}";
            # name of the file to create that contain both the transcriptome and the intergenic sequences
            my $transcriptome_ref_intergenic_file_path = $file =~ s/gtf_all/transcriptome_ref_intergenic.fa/r;
            if (-e $transcriptome_ref_intergenic_file_path.'.gz') {
                system("gunzip -c $transcriptome_ref_intergenic_file_path.gz > $transcriptome_ref_intergenic_file_path") == 0
                    or die "Failed to gunzip $transcriptome_ref_intergenic_file_path.gz";
            } elsif (-e $transcriptome_ref_intergenic_file_path.'xz') {
                system("unxz $transcriptome_ref_intergenic_file_path.xz") == 0
                    or die "Failed to unxz $transcriptome_ref_intergenic_file_path.xz";
            }
            my $sbatch_commands = "";
            if (! -e $transcriptome_ref_intergenic_file_path) {
                # ref intergenic sequences are always gzipped
                my $ref_intergenic_path = $ref_intergenic_folder.'/'.$speciesId.'_intergenic.fa.gz';
                # transcriptome file without any intergenic sequences
                my $transcriptome_wo_intergenic_path = $file =~ s/gtf_all/transcriptome_wo_intergenic.fa/r;
                if (! -e $ref_intergenic_path) {
                    warn "Reference intergenic file $ref_intergenic_path does not exist. Please provide it. No transcriptome with reference intergenic file will be generated for $species_name.\n";
                    next;
                }
                if (-e $transcriptome_wo_intergenic_path.'.gz') {
                    system("gunzip -c $transcriptome_wo_intergenic_path.gz > $transcriptome_wo_intergenic_path") == 0
                        or die "Failed to gunzip $transcriptome_wo_intergenic_path.gz";
                } elsif (-e $transcriptome_wo_intergenic_path.'xz') {
                    system("unxz $transcriptome_wo_intergenic_path.xz") == 0
                        or die "Failed to unxz $transcriptome_wo_intergenic_path.xz";
                }
                if (! -e $transcriptome_wo_intergenic_path) {
                    warn "Transcriptome file without intergenic sequences $transcriptome_wo_intergenic_path does not exist. Please provide it. No transcriptome with reference intergenic file will be generated for $species_name.\n";
                    next;
                }
                system("cat $transcriptome_wo_intergenic_path <(zcat $ref_intergenic_path) > $transcriptome_ref_intergenic_file_path") == 0
                    or die "Failed to concatenate $transcriptome_wo_intergenic_path and $ref_intergenic_path into $transcriptome_ref_intergenic_file_path";
                next if (! -e $transcriptome_ref_intergenic_file_path);
            }
            $sbatch_commands .= "# generate index with default kmer size\n";
            $sbatch_commands .= "$container_cmd kallisto index -i $transcriptome_ref_intergenic_index_path $transcriptome_ref_intergenic_file_path &&\n";
            $sbatch_commands .= "$container_cmd gzip $transcriptome_ref_intergenic_file_path\n";

            # create the sbatch file
            my $sbatch_file_path = "$sbatch_folder/${species_name}.sbatch";
            my $output_file_path = "$clusterOutput_folder/${species_name}.out";
            my $error_file_path = "$clusterOutput_folder/${species_name}.err";
            open (my $OUT, '>', "$sbatch_file_path")  or die "Cannot write [$sbatch_file_path]\n";
            print {$OUT} Utils::sbatch_template($partition, $account, 1, $memory_usage, $output_file_path, $error_file_path,
                $jobName);
            print {$OUT} $sbatch_commands;
            close $OUT;
            system("sbatch $sbatch_file_path")==0  or print "Failed to submit job [$file]\n";
        }
    }
}
# count number of jobs running
my $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $job_prefix);
while ($jobsRunning > 0) {
    sleep(15);
    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $job_prefix);
}


