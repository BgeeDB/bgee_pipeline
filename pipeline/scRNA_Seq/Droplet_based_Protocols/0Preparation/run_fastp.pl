#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use File::Slurp;
use FindBin;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;

# Julien Wollbrett, created March 2023
# Script used to run fastp and generate R.stat file for all downloaded target-based run

# Commands to run on the front node of the cluster
# 
# module use /software/module/
# module load Development/Ensembl_API/97;
# module load UHTS/Quality_control/fastp/0.19.5;
# export PATH=/software/bin:$PATH
#
#####################################################################

# Define arguments & their default value
my ($metadata_file, $fastq_dir, $fastp_path)  = ('', '', '');
my %opts = ('metadata_file=s'    => \$metadata_file,   # File containing metadata of all run to process
            'fastq_dir=s'        => \$fastq_dir,       # Directory where all target-based FASTQ files are donwloaded
            'fastp_path=s'       => \$fastp_path       # Path to the fastp tool
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadata_file || $fastq_dir eq ''|| $fastp_path eq ''){
    print "\n\tInvalid or missing argument:
\te.g., $0  -metadata_file=\$(METADATA_FILE) -fastq_dir=R\$(FASTQ_DIR) > $@.tmp 2>warnings.$@
\t-metadata_file          File containing metadata of all run to process
\t-fastq_dir              Directory where all target-based FASTQ files are donwloaded
\t-fastp_path             Path to the fastp tool\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../target_base_utils.pl");

my %metadata = get_processed_libraries_info($metadata_file);

for my $experiment_id (sort keys %metadata) {
    for my $library_id (sort keys %{$metadata{$experiment_id}}) {
        for my $run_id (sort keys %{$metadata{$experiment_id}{$library_id}}) {
            ## speciesId is a key used to retrieve the species ID of the library
            next if ($run_id eq "speciesId");
	    my $run_path = "$fastq_dir/EXPERIMENTS/$experiment_id/$library_id/$run_id/";
            my $fastq_path = "${run_path}FASTQ/";
            my $fastq_fastp = '';
            my $fastq_R     = '';
            if ( -s "${fastq_path}${run_id}_R1.fastq.gz" 
                && -s "${fastq_path}${run_id}_R2.fastq.gz" ){
                $fastq_fastp = "${fastq_path}${run_id}_R1.fastq.gz -I ${fastq_path}${run_id}_R2.fastq.gz";
                $fastq_R     = "${fastq_path}${run_id}_R1.fastq.gz    ${fastq_path}${run_id}_R2.fastq.gz";
            } else {
                warn "FASTQ files were not properly downloaded for experiment : ", $experiment_id,
                " library : ", $library_id, " run : ", $run_id, "\n";
                next;
            }
            # Run FastP (A quality control tool for high throughput sequence data) for ALL SRR (runs)
            # as well as basic read length statistics with R
            #NOTE Would be nice to have all basic stats from FastP (currently some are done in R)
            if ( !-e "${run_path}${run_id}.fastp.html.xz" || !-e "${run_path}${run_id}.fastp.json.xz" ){
                system("$fastp_path -i $fastq_fastp --json ${run_path}${run_id}.fastp.json --html ${run_path}${run_id}.fastp.html --thread 2  > ${run_path}${run_id}.fastp.log 2>&1")==0
                    or do { warn "\tfastp failed for [${run_path}${run_id}]\n" };
                system("xz -9 ${run_path}${run_id}.fastp.html ${run_path}${run_id}.fastp.json");
            }
            if ( !-e "${run_path}${run_id}.R.stat" ){
                system("/bin/echo \"#min\tmax\tmedian\tmean\" > ${run_path}${run_id}.R.stat");
                #NOTE for cases like SRX1372530 with paired-end files coming with a single-end file in the same run, use ${prefix}*.fastq.gz ???
                system("zcat $fastq_R | sed -n '2~4p' | awk '{print length(\$0)}' | Rscript -e 'd<-scan(\"stdin\", quiet=TRUE);cat(min(d), max(d), median(d), mean(d), sep=\"\\t\");cat(\"\\n\")' >> ${run_path}${run_id}.R.stat");
            }
        }
    }
}


