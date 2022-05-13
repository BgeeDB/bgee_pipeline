#!usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;

# Julien Wollbrett, created November 2020

# This script transform the scrna_seq_sample_info.txt file into a file used as input to run BgeeCall.
#TODO use a tsv reader module to filter not on column position but on column names

my ($sc_sample_info_file)    = ('');
my ($transcriptome_dir) = ('');
my ($annotation_dir) = ('');
my $fastq_dir = '';
my $output_dir = '';
my $bgeecall_file = '';
my $ref_intergenic_dir = '';

my %opts = ('sc_sample_info_file=s' => \$sc_sample_info_file,           # path to scrna_seq_sample_info file
            'transcriptome_dir=s'   => \$transcriptome_dir,             # path to directory containing all transcriptomes
            'annotation_dir=s'      => \$annotation_dir,                # path to directory containing all annotations
            'fastq_dir=s'           => \$fastq_dir,                     # path to directory containing all fastq files
            'output_dir=s'            => \$output_dir,                    # path to the directory where BgeeCall results should be stored
            'bgeecall_file=s'       => \$bgeecall_file,                 # path to the output file compatible with BgeeCall
            'ref_intergenic_dir=s'  => \$ref_intergenic_dir             # path to directory containing all reference intergenic sequences
           );

# test arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $sc_sample_info_file eq '' || $transcriptome_dir eq '' || $annotation_dir eq '' || $fastq_dir eq '' || $output_dir eq '' || $bgeecall_file eq '' || $ref_intergenic_dir eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0  -sc_sample_info_file=\$(RNASEQ_SAMPINFO_FILEPATH) -transcriptome_dir=\$(RNASEQ_CLUSTER_GTF) -annotation_dir=\$(RNASEQ_CLUSTER_GTF) -fastq_dir=\$(RNASEQ_SENSITIVE_FASTQ) -output_dir=\$(OUTPUT_DIR) -bgeecall_file=\$(RNASEQ_BGEECALL_FILE) -ref_intergenic_dir=\$(CLUSTER_REF_INTERGENIC_FOLDER)
\t-sc_sample_info_file     Path to rna_seq_sample_info file
\t-transcriptome_dir    Path to directory containing all transcriptomes
\t-annotation_dir       Path to directory containing all annotations
\t-fastq_dir            Path to directory containing all FASTQ
\t-output_dir           Path to directory containing BgeeCall results
\t-bgeecall_file        Path to the output file compatible with BgeeCall
\t-ref_intergenic_dir   Path to directory containing all reference intergenic sequences
\n";
    exit 1;
}

open(FH, '>', $bgeecall_file) or die $!;
# write header
print FH "species_id\trun_ids\treads_size\trnaseq_lib_path\ttranscriptome_path\tannotation_path\toutput_directory\tcustom_intergenic_path\n";


open(my $sample_info, $sc_sample_info_file) || die "failed to read sample info file: $!";
while (my $line = <$sample_info>) {
    chomp $line;
     ## skip comment lines
    next  if ( ($line =~ m/^libraryId/) or ($line =~ m/^\"libraryId/) );
    my @line = split(/\t/, $line);
    my $number_columns = 13;
    if (! scalar @line eq $number_columns) {
        die "all lines of full length single cell sample info file should have $number_columns columns";
    }
    #column organism
    my $species_name = $line[15];
    opendir(my $dir, $transcriptome_dir);
    my @transcriptome_file = grep(/$species_name.*\.transcriptome_wo_intergenic.fa/, readdir($dir));
    closedir($dir);
    my $transcriptome_path = "$transcriptome_dir$transcriptome_file[1]";
    my $annotation_path = $transcriptome_path =~ s/transcriptome_wo_intergenic\.fa/transcriptome\.gtf/r;
    my $fastq_path = "$fastq_dir$line[0]";
    my $library_output_dir = "$output_dir$line[0]";
    my $intergenic_file = "$ref_intergenic_dir$line[4]_intergenic.fa.gz";
    my $output_line = "$line[4]\t\t$line[14]\t$fastq_path\t$transcriptome_path\t$annotation_path\t$library_output_dir\t$intergenic_file\n";

    print FH $output_line;

}
close(FH);
