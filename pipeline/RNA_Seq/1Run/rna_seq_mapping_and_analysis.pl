#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use FindBin qw( $RealBin ); # directory where the script is lying
use File::Path qw(make_path);
use File::Slurp;
use List::Util qw(min max);
use Getopt::Long;
#use JSON::XS; # See http://blogs.perl.org/users/e_choroba/2018/03/numbers-and-strings-in-json.html
use Cpanel::JSON::XS;

my $GTEX_exp_id = 'SRP012682';

## Julien Roux 05/01/2016
## This script executes the main RNA-seq analysis, including pseudo-mapping of reads with Kallisto, summing the counts at the transcripts level to the gene level, ...
## It is inspired by Marta's previous script which used Tophat to map reads, but it is greatly simplified: no config file, download of files removed, no generation of BAM files, no HT-seq counting step
## The script should be launched for a given libraryId (e.g. SRX081872), and the analysis will be executed only for this library
## All output files are written in the results folder, as well as log files (e.g. SRX081872.out, SRX081872.err, and SRX081872.Rout)

# Define arguments & their default value
my ($library_id, $sample_info_file, $exclude_sample_file, $index_folder, $fastq_folder, $kallisto_out_folder, $output_log_folder, $enc_passwd_file) = ('', '', '', '', '', '', '', '', '', '', '', '');
my %opts = ('library_id=s'           => \$library_id,
            'sample_info_file=s'     => \$sample_info_file,
            'exclude_sample_file=s'  => \$exclude_sample_file,
            'index_folder=s'         => \$index_folder, # same as GTF folder
            'fastq_folder=s'         => \$fastq_folder,
            'kallisto_out_folder=s'  => \$kallisto_out_folder,
            'output_log_folder=s'    => \$output_log_folder,
            'enc_passwd_file=s'      => \$enc_passwd_file,
           );
# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $library_id eq '' || $sample_info_file eq '' || $index_folder eq '' || $fastq_folder eq '' || $kallisto_out_folder eq '' || $output_log_folder eq '' || $enc_passwd_file eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -library_id=... -sample_info_file=\$(RNASEQ_SAMPINFO_FILEPATH) -exclude_sample_file=\$(RNASEQ_SAMPEXCLUDED_FILEPATH) -index_folder=\$(RNASEQ_CLUSTER_GTF)  -fastq_folder=\$(RNASEQ_SENSITIVE_FASTQ) -kallisto_out_folder=\$(RNASEQ_CLUSTER_ALL_RES) -enc_passwd_file=\$(ENCRYPT_PASSWD_FILE)
\t-library_id=s           Library to process
\t-sample_info_file=s     TSV with information on species and runs for each library
\t-exclude_sample_file=s  rna_seq_sample_excluded.txt
\t-index_folder=s         Folder with Kallisto indexes (same as GTF folder)
\t-fastq_folder=s         Folder with Fastq files on big bgee
\t-kallisto_out_folder=s  Folder with Kallisto output and results
\t-output_log_folder=s    Folder where R output is written
\t-enc_passwd_file=s      File with password necessary to decrypt the GTEx data
\n";
    exit 1;
}

die "Invalid or missing [$sample_info_file]: $?\n"  if ( !-e $sample_info_file || !-s $sample_info_file );

# Sample to exclude if any
my %manually_excluded;
EXCLUSION:
for my $line ( read_file("$exclude_sample_file", chomp=>1) ){
    #libraryId    excluded    comment    annotatorId    lastModificationDate
    next EXCLUSION  if ( $line =~ /^#/); # header or comment

    my ($sampleId, $to_exclude) = split(/\t/, $line);
    next EXCLUSION  if ( $to_exclude ne 'TRUE' );

    $manually_excluded{$sampleId} = 1;
}
die "Excluded library [$library_id]\n"  if ( exists $manually_excluded{$library_id} );


# Reading $sample_info_file
my $libraryExists = 'FALSE';
my $exp_id;
my $organism;
my $genomeFilePath;
my $database;
my $libraryType,
my $libraryInfo,
my $readLength,
my @run_ids;

my %run_paths;
my %run_types;
my %run_orientations;
my $k = 0;
#libraryId  experimentId  speciesId  organism  genomeFilePath  platform  libraryType  libraryInfo  readLength   runIds
for my $line ( read_file($sample_info_file, chomp => 1) ){
    # TODO use the getAllRnaSeqLibrariesInfo subroutine in rna_seq_utils.pl?
    next  if ( ($line =~ m/^#/) or ($line =~ m/^\"#/) );
    my @fields = split("\t", $line);
    if ( $fields[0] eq $library_id ){
        $libraryExists  = 'TRUE';
        $exp_id         = $fields[1];
        $organism       = $fields[3];
        $genomeFilePath = $fields[4];
        $database       = $fields[5];
        $libraryType    = $fields[7];
        $readLength     = $fields[9];
        my @tmp_runs;
        for my $run ( split(',', $fields[10]) ){
            push @tmp_runs, $run  if ( !exists $manually_excluded{$run} );
        }
        if ( scalar @tmp_runs >= 1 ){
            push @run_ids, @tmp_runs;
        }
        else {
            die "All runs are excluded for [$library_id]!\n";
        }
    }
}
die "Library [$library_id] is not present in [$sample_info_file]\n"  unless ( $libraryExists );
print "Treating library: $library_id from $exp_id in $organism (Run ID(s): ", join(', ', @run_ids), ")\n";

# creating output directories if necessary:
# Do a mkdir -p  +  set group & mode at each intermediary folder
$kallisto_out_folder .= '/'.$library_id;
make_path $kallisto_out_folder, {verbose=>0, mode=>0775};
# Log file with commands submitted
my $report_file = $output_log_folder.'/'.$library_id.'/'.$library_id.'.report';

#########################################################################################################
## Checking the presence of fastq.gz files. Download script get_SRA.pl verifies that all files are zipped

# defining the path to folder with fastq.gz files for run
my $fastqSamplePath = $fastq_folder.'/'.$library_id;

print "\tChecking presence of Fastq files for all runs\n";
for my $run ( @run_ids ){
    # check if fastq file exists, and if there is a single one for single-end library / two for paired-end library)
    if ( $libraryType eq 'SINGLE' ){
        # GTEx files are encrypted. Files have .enc
        if ( $exp_id eq $GTEX_exp_id ){
            if ( system('test -f '.$fastqSamplePath.'/'.$run.'.fastq.gz.enc') eq 0 ){
                print "\tFound fastq.gz.enc file for single-end library (run $run)\n";
            }
            else {
                die "Problem: fastq.gz.enc file not found for single-end library (run $run) [$fastqSamplePath/$run.fastq.gz.enc]\n";
            }
        }
        else {
            if ( system('test -f '.$fastqSamplePath.'/'.$run.'.fastq.gz') eq 0 ){
                print "\tFound fastq.gz file for single-end library (run $run)\n";
            }
            else {
                die "Problem: fastq.gz file not found for single-end library (run $run) [$fastqSamplePath/$run.fastq.gz]\n";
            }
        }
    }
    elsif ( $libraryType eq 'PAIRED' ){
        # GTEx files are encrypted. Files have .enc
        if ( $exp_id eq $GTEX_exp_id ){
            if ( system('test -f '.$fastqSamplePath.'/'.$run.'_1.fastq.gz.enc') eq 0 ){
                print "\tFound fastq.gz.enc file for left reads for paired-end library (run $run)\n";
            }
            else {
                die "Problem: fastq.gz.enc file for left reads not found for paired-end library (run $run) [$fastqSamplePath/$run\_1.fastq.gz.enc]\n";
            }
            if ( system('test -f '.$fastqSamplePath.'/'.$run.'_2.fastq.gz.enc') eq 0 ){
                print "\tFound fastq.gz.enc file for right reads for paired-end library (run $run)\n";
            }
            else {
                die "Problem: fastq.gz.enc file for right reads not found for paired-end library (run $run) [$fastqSamplePath/$run\_2.fastq.gz.enc]\n";
            }
        }
        else {
            if ( system('test -f '.$fastqSamplePath.'/'.$run.'_1.fastq.gz') eq 0 ){
                print "\tFound fastq.gz file for left reads for paired-end library (run $run)\n";
            }
            else {
                die "Problem: fastq.gz file for left reads not found for paired-end library (run $run) [$fastqSamplePath/$run\_1.fastq.gz]\n";
            }
            if ( system('test -f '.$fastqSamplePath.'/'.$run.'_2.fastq.gz') eq 0 ){
                print "\tFound fastq.gz file for right reads for paired-end library (run $run)\n";
            }
            else {
                die "Problem: fastq.gz file for right reads not found for paired-end library (run $run) [$fastqSamplePath/$run\_2.fastq.gz]\n";
            }
        }
    }
}
#############################################################################################
# Before launching Kallisto, we need to verify read length
# We extract the first read/read pair from fastq files, and calculate its length
# Note: potential problems if all reads do not have same length in fastq files (e.g., reads were trimmed)

print "\tExtracting read length for all runs\n";
my $shortReads   = 0;
my $lengthCutoff = 50; # Read shorter than $lengthCutoff nt will be pseudo-mapped using shorter k-mer length. This parameter can be changed, but I assumed most recent RNA-seq data produce at least $lengthCutoff nt-long reads.
my @allLengths;
my $maxLength;

for my $run ( @run_ids ){
    my $read = `grep -v '^#' $fastqSamplePath\/$run\.R.stat | cut -f3`;
    chomp($read);
    if ( (!defined $read) or ($read eq '') ){
        die "Problem: Read length could not be extracted for run [$run]\n";
    }
    print "\tMedian read length = $read for run [$run]\n";

    # record min and max lengths across runs
    push @allLengths, $read;
    # verify that extracted read length is consistent with SRA info
    if ( ($read != $readLength) and ($readLength ne '') ){
        warn "\nProblem: Read length in FASTQ file [$read] is not consistent with SRA record [$readLength]. Please check [$run]\n";
    }
    # reads too short for Kallisto index with default k-mer length
    if ( $read < $lengthCutoff ){
        $shortReads = 1;
        warn "\nWarning: Length of reads [$read] too short for pseudo-mapping on index with k-mer length of 31nt. Library will be pseudo-mapped on index with k-mer length of 15nt [$run]\n";
    }
}

# report min and max read length across runs
open (my $REPORT, '>>', "$report_file")  or die "Cannot write [$report_file]\n";
print {$REPORT} "\n";
print {$REPORT} 'Minimum read length across runs: ', min(@allLengths), "\n";
print {$REPORT} 'Maximum read length across runs: ', max(@allLengths), "\n";
close $REPORT;


#TODO reads too short for Kallisto index with default k-mer length, use index built with k-mer size = 15nt
#  See: https://www.biostars.org/p/104321/
#       https://groups.google.com/forum/#!topic/kallisto-sleuth-users/clOeSROnnFI
my $index = $index_folder.'/';
$genomeFilePath =~ m/.+\/(.+)/;
if ( $shortReads == 1 ){
    $index .= $1.'.transcriptome_k15.idx';
}
else {
    $index .= $1.'.transcriptome.idx';
}
print "\tKallisto index $index will be used for pseudo-mapping\n";


#############################################################################################
# Before launching Kallisto, check FastP run. Each fastq file has its own FastP!
# This can help investigate problems (low number of reads pseudo-mapped, etc).
#TODO For future pipeline, if we trim adapters, this can also be useful
#
# Extract the number of reads in FastP files
print "\tExtracting total number of reads from FastP report...\n";
my %number_reads;
for my $run ( @run_ids ){
    open(my $FASTP_REPORT, "xzcat $fastqSamplePath/${run}.fastp.json.xz |")
        or warn "\nWarning: missing FastP report for run $run\n";

    my $json_report = '';
    $json_report .= $_  while <$FASTP_REPORT>;
    close $FASTP_REPORT;

    my $datastructure = decode_json($json_report);
    $number_reads{$run} = $datastructure->{'summary'}->{'before_filtering'}->{'total_reads'};
}

open (my $REPORT2, '>>', "$report_file")  or die "Cannot write [$report_file]\n";
print {$REPORT2} "\nNumber of reads (from FastP reports):\n";
my $total_reads_fastp = 0;
for my $run ( keys %number_reads ){
    print {$REPORT2} "\t", $run, "\t", $number_reads{$run}, "\n";
    $total_reads_fastp += $number_reads{$run};
}
print {$REPORT2} "\tTotal number reads\t", $total_reads_fastp, "\n";
close $REPORT2;

## TODO For each SRR, we should also store a file with the number of lines (wc -l). On this side, we can verify it is consistent with the number of reads in FastP reports (and number of reads processed by Kalisto = sum of all runs $total_reads_fastp too).


#############################################################################################
# Pseudo-mapping to transcriptome with kallisto:
# - No bootstraping (if needed, just add option -b 100 for 100 bootstraps)
# - No multithreading (not used unless boostraps are done, if needed, add -t ... option)
# - Sequence bias correction enabled
# - for single-end libraries we should provide fragment length, but it is not possible to estimate since we do not have BioAnalyzer results. So we give 180 bp by default, which should be close to real value, see https://groups.google.com/forum/#!topic/kallisto-sleuth-users/h5LeAlWS33w
# - We also need to provide sd for fragment size, but I have found very little info on this. Based on this post (https://groups.google.com/forum/#!topic/rsem-users/S31Rx01Xd18) I plugged sd=30. This can be changed later
# - output folder looks something like: all_results_bgee_v15/
# - for next Bgee releases, if for a species the genome version and annotation did not change, it is possible to copy paste the results for the the concerned libraries, so that these are not rerun. But be careful with this! In most cases folders of it is probably faster to rerun everything.

#my $local_not_remote_path = '...';
# First check that this library was not previously successfully processed
if ( ( -s $kallisto_out_folder.'/abundance.tsv' ) && ( -s $kallisto_out_folder.'/run_info.json' ) ){
    print "\tKallisto was already successfully run for this library. Skipping this step.\n";
}
# if Kallisto needs to be launched
else {
    # creating and invoking the full kallisto command
    my $kallisto_command = '';
    if ( $libraryType eq 'SINGLE' ){
        $kallisto_command .= "kallisto quant -i $index -o $kallisto_out_folder -t 1 --single -l 180 -s 30 --bias ";
        for my $run ( @run_ids ){
            if ( $exp_id ne $GTEX_exp_id ){
                $kallisto_command .= $fastqSamplePath.'/'.$run.'.fastq.gz ';
            }
            if ( $exp_id eq $GTEX_exp_id ){
                $kallisto_command .= '<(cat '.$fastqSamplePath.'/'.$run.'.fastq.gz.enc | openssl enc -aes-128-cbc -d -pass file:'.$enc_passwd_file.') ';
            }
        }
    }
    elsif ( $libraryType eq 'PAIRED' ){
        $kallisto_command .= "kallisto quant -i $index -o $kallisto_out_folder -t 1 --bias ";
        for my $run ( @run_ids ){
            if ( $exp_id ne $GTEX_exp_id ){
                $kallisto_command .= $fastqSamplePath.'/'.$run.'_1.fastq.gz '.$fastqSamplePath.'/'.$run.'_2.fastq.gz ';
            }
            if ( $exp_id eq $GTEX_exp_id ){
                $kallisto_command .= '<(cat '.$fastqSamplePath.'/'.$run.'_1.fastq.gz.enc | openssl enc -aes-128-cbc -d -pass file:'.$enc_passwd_file.') <(cat '.$fastqSamplePath.'/'.$run.'_2.fastq.gz.enc | openssl enc -aes-128-cbc -d -pass file:'.$enc_passwd_file.') ';
            }
        }
    }
    print "\tKallisto command was built: \"", $kallisto_command, "\"\n\tNow launching Kallisto...\n";
    # Print command to .report file
    open (my $REPORT3, '>>', "$report_file")  or die "Cannot write [$report_file]\n";
    print {$REPORT3} "\nComputing node / Kallisto command submitted:\n$kallisto_command\n";
    close $REPORT3;
    # Submit command. Bash syntax is needed to be able to have parentheses and pipes. See http://stackoverflow.com/questions/571368/how-can-i-use-bash-syntax-in-perls-system
    my @args = ( "bash", "-c",  $kallisto_command );
    system(@args)==0  or die "Problem: System call to Kallisto failed\n";
    print "\tDone\n";
}

############################################################################################
# Parsing of Kallisto's output

# Check that this library was successfully processed
if ( ( -s $kallisto_out_folder.'/abundance.tsv' ) && ( -s $kallisto_out_folder.'/run_info.json' ) ){
    open(my $IN1, '<', $kallisto_out_folder.'/abundance.tsv')  or die "could not read kallisto abundance file\n";
    my $line = <$IN1>; # header removed
    # number of TPMs abundances that are NAs
    my $nan = 0;
    # number of pseudo-aligned reads = sum of estimated counts
    my $aligned = 0;
    while ( defined ($line = <$IN1>) ){
        my @tmp = split(/\t/, $line);
        $nan++  if ( $tmp[4] eq '-nan' );
        $aligned += $tmp[3];
    }
    close $IN1;

    if ( $nan > 1000 ){
        # rename abundance file
        system('mv ', $kallisto_out_folder, '/abundance.tsv ', $kallisto_out_folder, '/abundance.tsv.problem');
        die "Problem: Kallisto TPM results include numerous \"-nan\", please check for a problem [$library_id]\n";
    }

    my $total_reads_kallisto;
    open(my $IN2, '<', $kallisto_out_folder.'/run_info.json')  or die "could not read kallisto JSON file\n";
    while ( defined ($line = <$IN2>) ){
        if ( $line =~ m/\"n\_processed\"\:\s(\d+),/ ){
            $total_reads_kallisto = $1;
        }
    }
    close $IN2;
    if ( defined $total_reads_kallisto ){
        # Export % pseudoaligned reads to .report file
        open (my $REPORT4, '>>', "$report_file")  or die "Cannot write [$report_file]\n";
        print {$REPORT4} "\nKallisto pseudo-aligned ", int($aligned), " reads out of $total_reads_kallisto (", $aligned / $total_reads_kallisto * 100, "%)\n";
        close $REPORT4;

        if ( $aligned / $total_reads_kallisto < 0.2 ){
            warn "\nProblem: It seems that less than 20% of the reads were pseudo-aligned by Kallisto, please check for a problem [$library_id]\n";
        }
        #NOTE this arbitrary threshold of 20% can be changed

        #NOTE for paired-end if fastp reads both files at the same time, the read number is two times the real value!
        if ( $total_reads_kallisto != $total_reads_fastp && $total_reads_kallisto != ($total_reads_fastp / 2) ){
            warn "\nProblem: The number of reads processed by FastP and Kallisto differs. Please check for a problem [$library_id] [$total_reads_kallisto vs $total_reads_fastp]\n";
        }
        ## TODO add step to check that the number of reads is also consistent with the original fastq files on storage (see TODO above)
    }
    if ( $aligned < 1_000_000 ){
        warn "\nProblem: Less than 1,000,000 reads were pseudo-aligned by Kallisto, please check for a problem [$library_id]\n";
    }
    # Note: this arbitrary threshold of 1,000,000 reads can be changed

    # TODO? Also export other info, such as fragment length distribution:
    # h5dump -d /aux/fld all_results_bgee_v15/SRX567038/abundance.h5 | less
}
else {
    die "Problem: No abundance.tsv or run_info.json file found for this library. Kallisto run was probably not successful [$library_id]\n";
}

############################################################################################
## Treating kallisto's results with rna_seq_analysis.R

# names of output files:
# Gene-level file (counts and TPMs) added to the kallisto output folder
my $count_info_file = $kallisto_out_folder.'/'.'abundance_gene_level.tsv';
print "\tGene-leve count and TPM file: $count_info_file\n";
# defining R_log file for rna_seq_analysis.R script
my $R_log_file = $output_log_folder.'/'.$library_id.'/'.$library_id.'.Rout';
print "\tR log file: $R_log_file\n";
# defining gene to transcript mapping file
my $gene2transcript = $index_folder.'/';
$genomeFilePath =~ m/.+\/(.+)/;
$gene2transcript .= $1.'.tx2gene';
# defining gene to gene biotype file
my $gene2biotype = $index_folder.'/';
$genomeFilePath =~ m/.+\/(.+)/;
$gene2biotype .= $1.'.gene2biotype';

# First check that this library was not previously successfully analyzed
if ( -s $count_info_file && -s $R_log_file ){
    print "\tR analysis script was already successfully run for this library. Skipping this step.\n";
}
else {
    #creating and invoking command that launches rna_seq_analysis.R script
    my $analyze_count_command = "R CMD BATCH --no-save --no-restore \'--args".
                              ' kallisto_count_folder="'.$kallisto_out_folder.'"'.
                              ' tx2gene_file="'.$gene2transcript.'"'.
                              ' gene2biotype_file="'.$gene2biotype.'"'.
                              ' gene_count_file="'.$count_info_file.'"'.
                              ' library_id="'.$library_id.'"\' '.
                              $RealBin.'/rna_seq_analysis.R '.$R_log_file;
    print "\tAnalysis script rna_seq_analysis.R command was built: \"", $analyze_count_command, "\"\n\tNow launching R script...\n";
    # Print command to .report file
    open (my $REPORT5, '>>', "$report_file")  or die "Cannot write [$report_file]\n";
    print {$REPORT5} "\nComputing node / R command submitted:\n$analyze_count_command\n";
    close $REPORT5;
    # Submit command
    system($analyze_count_command)==0  or die "Problem: System call to analyze_count_command failed\n";
    print "\tDone\n";
}

############################################################################################
# If everything was successful, we write a DONE.txt file in results folder
open (my $REPORT6, '>', $kallisto_out_folder.'/DONE.txt')  or die "Cannot write DONE.txt file\n";
print {$REPORT6} "$library_id was succesfully processed!\n";
close $REPORT6;
############################################################################################

exit 0;

