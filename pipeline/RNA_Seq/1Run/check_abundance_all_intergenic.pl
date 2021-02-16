#!usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Find qw(find);

# Julien Wollbrett, created Feb. 2021

# This script check .err and .txt files created for each library during the run of kallisto with all intergenic

my $sample_info_file = '';
my $sample_excluded = '';
my $result_dir = '';
my $output_file = '';
my %opts = ('sample_info_file=s'     => \$sample_info_file,      # path to rna_seq_sample_info file
            'sample_excluded=s'      => \$sample_excluded,       # path to rna_seq_sqmple_excluded file
            'result_dir=s'          => \$result_dir,            # path to the directory where all library results are saved after running kallisto with all intergenic
            'output_file=s'          => \$output_file            # path to the file where summary of results for all libraries is created
    );

# test arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ($sample_info_file eq '' || $output_file eq '' || $sample_excluded eq '' || $result_dir eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0  -sample_info_file=\$(RNASEQ_SAMPINFO_FILEPATH) -sample_excluded=$(RNASEQ_SAMPEXCLUDED_FILEPATH) -output_file=\"output_file.txt\" -result_dir=\$(RNASEQ_CLUSTER_ABUNDANCE_ALL)
\t-sample_info_file     Path to rna_seq_sample_info file
\t-sample_excluded      Path to rna_seq_sample_excluded file
\t-result_dir          Path to the directory where all library results are saved after running kallisto with all intergenic 
\t-output_file          Path to the file where summary of results for all libraries is created
\n";
    exit 1;
}

open my $fh_output, '>',  $output_file or die "Could not open $output_file: $!";

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
close $excluded;

open(my $sample_info, $sample_info_file) || die "failed to read sample info file: $!";
while (my $line = <$sample_info>) {
    chomp $line;
     ## skip comment lines
    next  if ( ($line =~ m/^#/) or ($line =~ m/^\"#/) );
    my @line = split(/\t/, $line);
    # skip libraries present in the sample excluded file
    next if (grep( /^$line[0]/, @excluded_libraries));
    my $library_path = "$result_dir/$line[0]";
    if (! -e "$library_path/DONE.txt") {
        print {$fh_output} "DONE.txt file is absent for library $line[0]\n";
        next;
    }
    my $library_txt = "$library_path/$line[0].txt";
    my @problem_lines;
    if(-e $library_txt) {
        open (my $fh_txt, $library_txt) || die "Could not open $library_txt: $!";
        while (my $problem_line = <$fh_txt>) {
            if($problem_line =~ /Problem/) {
                push(@problem_lines, $problem_line);
            }
        }
        while (my $warning_line = <$fh_txt>) {
            next if($warning_line =~ /Warning/);
            print {$fh_output} $warning_line;
        }
        close $fh_txt;
    } else {
        print {$fh_output} "file $library_txt does not exist!\n";
    }
    foreach my $problem (@problem_lines) {
        if($problem =~ /Read length in FASTQ file \[([^\]]+)\] is not consistent with SRA record \[([^\]]+)\].*/) {
            my $fq_size = $1;
            my $sra_size = $2;
            next if (abs($fq_size - $sra_size) < 1);
            next if (abs($fq_size - 2*$sra_size) < 2);
            next if (abs($sra_size - 2*$fq_size) < 2);
            print {$fh_output} $problem;
        } else {
            print {$fh_output} $problem;
        }
    }
    my $library_err = "$library_path/$line[0].err";
    if(-e $library_err) {
        open (my $fh_err, $library_err) || die "Could not open $library_err: $!";
        while (my $err_line = <$fh_err>) {
            next if($err_line =~ /Unable to locate/);
            print {$fh_output} $err_line;
        }
        close $fh_err;
    } else {
        print {$fh_output} "file $library_err does not exist!\n";
    } 
}
close $sample_info;
close $fh_output;
