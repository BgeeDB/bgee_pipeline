#!/usr/bin/env perl

## Highly inspired from the script slurm_index_creation.pl used in the bulk RNASeq pipeline.
##TODO: refactor the code to avoid code redundancy
## Julien Wollbrett, May 9 2023
# This script launches the slurm jobs creating both single nucleus transcriptome and kallisto index.


# Perl core modules
use strict;
use warnings;
use diagnostics;

use File::Path qw(make_path);
use FindBin qw( $RealBin ); # directory where the script is lying
use Getopt::Long;
use Time::localtime;

# Define arguments & their default value
my ($transcriptome_folder, $output_log_folder, $account, $partition, $cluster_kallisto_cmd, $cluster_tophat_cmd) = ('', '', '', '', '', '', '', '', '');
my %opts = ('transcriptome_folder=s' => \$transcriptome_folder, # same as GTF folder
            'output_log_folder=s'    => \$output_log_folder,
            'account=s'              => \$account,
            'partition=s'            => \$partition,
            'cluster_kallisto_cmd=s' => \$cluster_kallisto_cmd,
            'cluster_tophat_cmd=s'   => \$cluster_tophat_cmd
           );

my $test_options = Getopt::Long::GetOptions(%opts);

# TO IMPLEMENT
# kallisto index generation is no multithreaded
my $nbr_processors = 1;
# RAM needed: 10GB should be enough
my $memory_usage   = 90;      # in GB
my $time_limit     = '12:00:00';

# retrieve path to all transcriptome
# for each
#   create string for sbatch file with header, command to generate index,
#   check number of jobs running
#   run jobs for one species ( k31 and k15)

opendir (DIR, $transcriptome_folder) or die "cannot open directory [$transcriptome_folder]";
while (my $file = readdir(DIR)) {
    if($file =~ ('.nascent.gtf$')) {

        # initialise path to all files
        my $genome_file_path = $transcriptome_folder.'/'.($file =~ s/.nascent.gtf/.genome.fa/r);
        my $transcriptome_file_path = $transcriptome_folder.'/'.($file =~ s/.nascent.gtf/.transcriptome.fa/r);
        my $transcriptome_nascent_file_path = $transcriptome_folder.'/'.($file =~ s/.nascent.gtf/.nascent_transcriptome.fa/r);
        my $transcriptome_single_nucleus_file_path = $transcriptome_folder.'/'.($file =~ s/.nascent.gtf/.single_nucleus_transcriptome.fa/r);
        my $transcriptome_single_nucleus_index_path = $transcriptome_folder.'/'.($file =~ s/.nascent.gtf/.single_nucleus_transcriptome.idx/r);
        my $sbatch_file_path = $transcriptome_folder.'/'.($file =~ s/.nascent.gtf/.index_single_nucleus.sbatch/r);
        my $output_file_path = $output_log_folder.'/'.($file =~ s/.nascent.gtf/.index_single_nucleus.out/r);
        my $error_file_path = $output_log_folder.'/'.($file =~ s/.nascent.gtf/.index_single_nucleus.err/r);

        # generate sbatch commands

        next if (-e $transcriptome_single_nucleus_index_path);

        # load vital-it softwares
        my $sbatch_commands = "module use /software/module/\n";

        # generate transcriptome with all intergenic sequences
        if (!-e $transcriptome_single_nucleus_file_path) {
            if (-e "$transcriptome_single_nucleus_file_path.xz") {
                $sbatch_commands .= "# unxz already existing transcriptome file\n";
                $sbatch_commands .= "unxz $transcriptome_single_nucleus_file_path.xz\n";
            } else {
                if (!-e $genome_file_path) {
                    if (-e "$genome_file_path.xz") {
                        $sbatch_commands .= "# unxz already existing genome file\n";
                        $sbatch_commands .= "unxz $genome_file_path.xz\n";
                    } else {
                        die "can not acces to genome file $genome_file_path";
                    }
                }
                if (!-e $transcriptome_file_path) {
                    if (-e "$transcriptome_file_path.xz") {
                        $sbatch_commands .= "# unxz already existing transcriptome file\n";
                        $sbatch_commands .= "unxz $transcriptome_file_path.xz\n";
                    } else {
                        die "can not acces to transcriptome file $transcriptome_file_path";
                    }
                }
                # generate transcriptome with tophat gtf_to_fasta
                # load tophat module
                $sbatch_commands .= "$cluster_tophat_cmd\n";
                # generate transcriptome
                $sbatch_commands .= "gtf_to_fasta $transcriptome_folder/$file $genome_file_path $transcriptome_nascent_file_path\n";
                $sbatch_commands .= "# update transcriptome file to correct the header of each sequence\n";
                $sbatch_commands .= 'perl -i -pe \'s/^>\\d+ +/>/\' '.$transcriptome_nascent_file_path."\n";
                # concatenate transcriptome containing both matured and intergenic transcripts and the transcriptome contining ascent transcripts
                $sbatch_commands .= "cat $transcriptome_file_path $transcriptome_nascent_file_path > $transcriptome_single_nucleus_file_path\n";
                # remove nascent transcripts file
                $sbatch_commands .= "rm $transcriptome_nascent_file_path\n";
                # compress transcriptome file
                $sbatch_commands .= "xz --threads=2 -9 $transcriptome_file_path\n";
                # compress genome file
                $sbatch_commands .= "xz --threads=2 -9 $genome_file_path\n";
            }
        }

        # load kallisto module
        $sbatch_commands .= "$cluster_kallisto_cmd\n";

        # generate index with default kmer size
        if (!-e $transcriptome_single_nucleus_index_path) {
            $sbatch_commands .= "# generate index with default kmer size\n";
            $sbatch_commands .= "kallisto index -i $transcriptome_single_nucleus_index_path $transcriptome_single_nucleus_file_path\n";
        }

        # compress single nucleus transcriptome file
        $sbatch_commands .= "xz --threads=2 -9 $transcriptome_single_nucleus_file_path\n";

        # compress tlst transcriptome file
        $sbatch_commands .= "if [ -f \"$transcriptome_file_path.tlst\" ]; then\n";
        $sbatch_commands .= "   xz --threads=2 -9 $transcriptome_file_path.tlst\n";
        $sbatch_commands .= "fi\n";

        # generate sbatch file
        open (my $OUT, '>', "$sbatch_file_path")  or die "Cannot write [$sbatch_file_path]\n";
        print {$OUT} sbatch_header($partition, $account, $nbr_processors, $memory_usage, $output_file_path, $error_file_path, $file, $time_limit);
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
    my ($partition, $account, $nbr_processors, $memory_usage, $output_file, $error_file, $job_name, $time_limit) = @_;

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

";

    return $template;
}

