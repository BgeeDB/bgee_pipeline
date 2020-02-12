#!/usr/bin/env perl

## Julien Wollbrett, Nov 1, 2019
# This script launches the slurm jobs creating kallisto indexes for all species.


# Perl core modules
use strict;
use warnings;
use diagnostics;

use File::Path qw(make_path);
use FindBin qw( $RealBin ); # directory where the script is lying
use Getopt::Long;
use Time::localtime;

# Define arguments & their default value
my ($transcriptome_folder, $output_log_folder, $ens_release, $ens_metazoa_release, $short_index_length, $cluster_kallisto_cmd, $cluster_tophat_cmd) = ('', '', '', '', '', '', '');
my %opts = ('transcriptome_folder=s' => \$transcriptome_folder, # same as GTF folder
			'output_log_folder=s'    => \$output_log_folder,
            'ens_release=s'          => \$ens_release,
            'ens_metazoa_release=s'  => \$ens_metazoa_release,
			'short_index_length=s'   => \$short_index_length,
            'cluster_kallisto_cmd=s' => \$cluster_kallisto_cmd,
            'cluster_tophat_cmd=s'   => \$cluster_tophat_cmd
           );

my $test_options = Getopt::Long::GetOptions(%opts);

# TO IMPLEMENT
# kallisto index generation is no multithreaded
my $nbr_processors = 1;
# RAM needed: 10GB should be enough
my $memory_usage   = 90;      # in GB
my $user_email     = 'julien.wollbrett@unil.ch'; # for email notification
my $account        = 'mrobinso_bgee';
my $queue          = 'ax-normal';
my $time_limit	   = '12:00:00';

# retrieve path to all transcriptome
# for each
#	create string for sbatch file with header, command to generate index,
#	check number of jobs running
#   run jobs for one species ( k31 and k15)

opendir (DIR, $transcriptome_folder) or die "cannot open directory [$transcriptome_folder]";
while (my $file = readdir(DIR)) {
	if($file =~ ($ens_release.'.gtf_all$') || $file =~ ($ens_metazoa_release.'.gtf_all$')) {

		# initialise path to all files
		my $genome_file_path = $transcriptome_folder.'/'.($file =~ s/gtf_all/genome.fa/r);
		my $transcriptome_file_path = $transcriptome_folder.'/'.($file =~ s/gtf_all/transcriptome.fa/r);
		my $transcriptome_index_path = $transcriptome_folder.'/'.($file =~ s/gtf_all/transcriptome.idx/r);
		my $short_transcriptome_index_path = $transcriptome_folder.'/'.($file =~ s/gtf_all/transcriptome_k$short_index_length.idx/r);
		my $sbatch_file_path = $transcriptome_folder.'/'.($file =~ s/gtf_all/index.sbatch/r);
		my $output_file_path = $output_log_folder.'/'.($file =~ s/gtf_all/index.out/r);
		my $error_file_path = $output_log_folder.'/'.($file =~ s/gtf_all/index.err/r);

		# generate sbatch commands

		next if (-e $short_transcriptome_index_path && -e $transcriptome_index_path);

		# load vital-it softwares
		my $sbatch_commands = "module add Bioinformatics/Software/vital-it\n";

		if (!-e $transcriptome_file_path) {
			if (-e "$transcriptome_file_path.xz") {
				$sbatch_commands .= "# unxz already existing transcriptome file\n";
				$sbatch_commands .= "unxz $transcriptome_file_path.xz\n";
			} else {
				if (!-e $genome_file_path) {
					die "can not acces to genome file $genome_file_path";
				}
				# generate transcriptome with tophat gtf_to_fasta
				# load tophat module
				$sbatch_commands .= "$cluster_tophat_cmd\n";
				# generate transcriptome
				$sbatch_commands .= "gtf_to_fasta $transcriptome_folder/$file $genome_file_path $transcriptome_file_path\n";
				$sbatch_commands .= "# update transcriptome file to correct the header of each sequence\n";
				$sbatch_commands .= 'perl -i -pe \'s/^>\\d+ +/>/\' '.$transcriptome_file_path."\n";
			}
		}

		# load kallisto module
		$sbatch_commands .= "$cluster_kallisto_cmd\n";

		# generate index with default kmer size
		if (!-e $transcriptome_index_path) {
			$sbatch_commands .= "generate index with default kmer size\n";
			$sbatch_commands .= "kallisto index -i $transcriptome_index_path $transcriptome_file_path\n";
		}

		# generate index with short kmer size
		if (!-e $short_transcriptome_index_path) {
			$sbatch_commands .= "#generate index with a kmer size of $short_index_length\n";
			$sbatch_commands .= "kallisto index -k $short_index_length -i $short_transcriptome_index_path $transcriptome_file_path\n";
		}

		# delete genome file not useful for next steps of the pipeline
		$sbatch_commands .= "# rm -rf $genome_file_path\n";

		# compress transcriptome file
		$sbatch_commands .= "xz --threads=2 -9 $transcriptome_file_path\n";

		# compress tlst transcriptome file
		$sbatch_commands .= "if [ -f \"$transcriptome_file_path.tlst\" ]; then\n";
		$sbatch_commands .= "	xz --threads=2 -9 $transcriptome_file_path.tlst\n";
		$sbatch_commands .= "fi\n";

		# generate sbatch file
		open (my $OUT, '>', "$sbatch_file_path")  or die "Cannot write [$sbatch_file_path]\n";
		print {$OUT} sbatch_header($queue, $account, $nbr_processors, $memory_usage, $output_file_path, $error_file_path, $file, $time_limit);
		print {$OUT} $sbatch_commands;
		close $OUT;

		# Then, run the job
		print "run sbatch script $sbatch_file_path => ";
    	system("sbatch $sbatch_file_path")==0  or print "Failed to submit job [$file]\n";
    	# TODO: should check number of jobs running in order to remove .sbatch files once all indexes have been generate
	}

}

############ Functions #############

# Add main sbatch command and options
sub sbatch_header {
    my ($queue, $account, $nbr_processors, $memory_usage, $output_file, $error_file, $job_name, $time_limit) = @_;
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
#SBATCH --time=$time_limit

#SBATCH --output=$output_file
#SBATCH --error=$error_file
#SBATCH --export=NONE
#SBATCH --job-name=$job_name

";

    return $template;
}

