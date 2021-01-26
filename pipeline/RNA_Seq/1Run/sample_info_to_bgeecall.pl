#!usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;

# Julien Wollbrett, created March 2020

# This script transform the rna_seq_sample_info.txt file into a file used as input to run BgeeCall.

my ($sample_info_file)    = ('');
my ($transcriptome_dir) = ('');
my ($annotation_dir) = ('');
my $fastq_dir = '';
my $bgeecall_file = '';
my $ref_intergenic_dir = '';
my $sample_excluded = '';
my $output_dir = '';
my %opts = ('sample_info_file=s'    => \$sample_info_file,      # path to rna_seq_sample_info file
            'sample_excluded=s'     => \$sample_excluded,       # path to rna_seq_sqmple_excluded file
            'transcriptome_dir=s'   => \$transcriptome_dir,     # path to directory containing all transcriptomes
            'annotation_dir=s'      => \$annotation_dir,        # path to directory containing all annotations
            'fastq_dir=s'           => \$fastq_dir,             # path to directory containing all fastq files
            'bgeecall_file=s'       => \$bgeecall_file,         # path to the output file compatible with BgeeCall
            'ref_intergenic_dir=s'  => \$ref_intergenic_dir,    # path to directory containing all reference intergenic sequences
            'output_dir=s'          => \$output_dir             # path to the directory where all library results will be saved   
	);

# test arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $sample_info_file eq '' || $output_dir eq '' || $sample_excluded eq '' || $transcriptome_dir eq '' || $annotation_dir eq '' || $fastq_dir eq '' || $bgeecall_file eq '' || $ref_intergenic_dir eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0  -sample_info_file=\$(RNASEQ_SAMPINFO_FILEPATH) -sample_excluded=$(RNASEQ_SAMPEXCLUDED_FILEPATH) -transcriptome_dir=\$(RNASEQ_CLUSTER_GTF) -annotation_dir=\$(RNASEQ_CLUSTER_GTF) -fastq_dir=\$(RNASEQ_SENSITIVE_FASTQ) -bgeecall_file=\$(RNASEQ_BGEECALL_FILE) -ref_intergenic_dir=$(CLUSTER_REF_INTERGENIC_FOLDER)
\t-sample_info_file     Path to rna_seq_sample_info file
\t-sample_excluded     Path to rna_seq_sample_excluded file
\t-transcriptome_dir    Path to directory containing all transcriptomes
\t-annotation_dir       Path to directory containing all annotations
\t-fastq_dir            Path to directory containing all FASTQ
\t-output_dir           Path to directory that will contain calls for all libraries
\t-bgeecall_file        Path to the output file compatible with BgeeCall
\t-ref_intergenic_dir   Path to directory containing all reference intergenic sequences
\n";
    exit 1;
}

open(my $FH, '>', $bgeecall_file)  or die $!;
# write header
print {$FH} "species_id\trun_ids\treads_size\trnaseq_lib_path\ttranscriptome_path\tannotation_path\toutput_directory\tcustom_intergenic_path\n";

open(my $excluded, $sample_excluded) || die "failed to read sample excluded file: $!";
my @excluded_libraries;
while (my $line = <$excluded>) {
    chomp $line;
     ## skip comment lines
    next  if ( ($line =~ m/^#/) or ($line =~ m/^\"#/) );
    my @line = split(/\t/, $line);
    if ($line[1] eq "TRUE") {
        push(@excluded_libraries, $line[0])
    }
}

open(my $sample_info, $sample_info_file) || die "failed to read sample info file: $!";
while (my $line = <$sample_info>) {
    chomp $line;
     ## skip comment lines
    next  if ( ($line =~ m/^#/) or ($line =~ m/^\"#/) );
    my @line = split(/\t/, $line);
    # skip libraries present in the sample excluded file
    next if (grep( /^$line[0]/, @excluded_libraries));
    my $number_columns = 11;
    if (! scalar @line eq $number_columns) {
        die "all lines of sample info file should have $number_columns columns";
    }
    my $genomeFilePath = $line[4];
    $line[4] =~ m/.+\/(.+)/;
    my $prefixFilePath = $1;
    my $transcriptome_path = "$transcriptome_dir$prefixFilePath.transcriptome_wo_intergenic.fa";
    my $annotation_path = "$transcriptome_dir$prefixFilePath.transcriptome.gtf";
    my $fastq_path = "$fastq_dir$line[0]";
    my $intergenic_file = "$ref_intergenic_dir$line[2]_intergenic.fa.gz";
    my $lib_output_dir ="$output_dir/all_results/$line[0]";
    my $output_line = "$line[2]\t\t$line[9]\t$fastq_path\t$transcriptome_path\t$annotation_path\t$lib_output_dir\t$intergenic_file\n";

    print {$FH} $output_line;

}
close($FH);

