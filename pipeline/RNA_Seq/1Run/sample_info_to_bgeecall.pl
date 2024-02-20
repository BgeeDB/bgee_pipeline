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

my $renameFastqScript = "./rename_fastq.sh";
open(my $FH, '>', $bgeecall_file)  or die $!;
open(my $FH_missing, '>', "./missing_fastq.sh")  or die $!;
open(my $FH_processed, '>', "./already_processed")  or die $!;
open(my $FH_rename, '>', "$renameFastqScript")  or die $!;
# write header
print {$FH} "species_id\trun_ids\treads_size\trnaseq_lib_path\ttranscriptome_path\tannotation_path\toutput_directory\tcustom_intergenic_path\n";
print {$FH_missing} "#!/usr/bin/env bash\n";
print {$FH_rename} "#!/usr/bin/env bash\n";
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
my %missingFastqSpeciesDir = ();
LIBRARY: while (my $line = <$sample_info>) {
    chomp $line;
     ## skip comment lines
    next LIBRARY if ( ($line =~ m/^#/) or ($line =~ m/^\"#/) );
    my @line = split(/\t/, $line);
    # skip libraries present in the sample excluded file
    next LIBRARY if (grep( /^$line[0]/, @excluded_libraries));
    my $number_columns = 11;
    if (! scalar @line eq $number_columns) {
        die "all lines of sample info file should have $number_columns columns";
    }
    my $genomeFilePath = $line[4];
    $line[4] =~ m/.+\/(.+)/;
    my $prefixFilePath = $1;
    my $transcriptome_path = "$transcriptome_dir$prefixFilePath.transcriptome_wo_intergenic.fa";
    my $annotation_path = "$transcriptome_dir$prefixFilePath.transcriptome.gtf";
    my $fastq_path = "$fastq_dir$line[2]/$line[0]";
    #retrieve mean read length from R.stat file. If more than one run then only the first one
    #is considered to define mean read length
    my @runIds = split(/,/,$line[10]);
    my $rstatFilePath = "$fastq_path/$runIds[0].R.stat";
    my $meanReadLength = ();
    if (! -e $rstatFilePath) {
        $meanReadLength = $line[9];
    }else {
        open RSTAT, "< $rstatFilePath";
        chomp (my @linesRstatFile = <RSTAT>);
        if (exists $linesRstatFile[1]) {
            my @readStats = split(/\t/, $linesRstatFile[1]);
            $meanReadLength = $readStats[3];
        } else {
	    $meanReadLength = $line[9];
        }
        close RSTAT;
    }
    my $intergenic_file = "$ref_intergenic_dir$line[2]_intergenic.fa.gz";
    my $lib_output_dir ="$output_dir/all_results/$line[0]";
    if ( -e "$lib_output_dir/gene_level_abundance+calls.tsv") {
            print {$FH_processed} "$line[2] - $line[0]\n";
            next LIBRARY;
    }
    if (! -e "$fastq_path") {
        my $mkdirCmd = "mkdir -p $fastq_dir$line[2]";
	if (! exists($missingFastqSpeciesDir{$mkdirCmd})) {
            #TODO: generate a bash script updating name of RUN_ID.fastq.gz file when both RUN_ID_1.fastq.gz file and RUN_ID.fastq.gz files exists AND the library is paired
	    #the new name should be RUN_ID.fastq.gz.EXTRA in order to allow BgeeCall to detect how to run kallisto otherwise BgeeCall throw an error
	    $missingFastqSpeciesDir{$mkdirCmd} = 1;
	    print {$FH_missing} "$mkdirCmd\n";
        }
        print {$FH_missing} "cp -r /archive/FAC/FBM/DEE/mrobinso/bgee_sensitive/FASTQ/RNAseq/$line[2]/$line[0] $fastq_path\n";
        next LIBRARY;
    }
    my $findFastqIssue = 0;
    for my $runId (@runIds) {
      #check if there is both _1.fastq.gz and fastq.gz file for each run of a library. It causes an error while running BgeeCall so if the library is
      # paired end and both files exist then we create a bash script allowing to rename the fastq.gz file to fastq.gz.EXTRA
      if ($line[7] eq "PAIRED" and -e "$fastq_path/${runId}_1.fastq.gz" and -e "$fastq_path/${runId}.fastq.gz") {
        $findFastqIssue = 1;
	print {$FH_rename} "mv $fastq_path/${runId}.fastq.gz $fastq_path/${runId}.fastq.gz.EXTRA\n";
      }
      #now check that both R.stat file and .fastp.json.xz files exist"
      if (! -e "$fastq_path/${runId}.R.stat" || ! -e "$fastq_path/${runId}.fastp.json.xz") {
        warn "R.stat or fastp.json.xz file does not exist for the run $runId of library $line[0] (taxId: $line[2]). Generate those files to generate calls for that library";
	next LIBRARY;
      }
    }
    if ($findFastqIssue) {
      warn "Found FASTQ files named both as single-end run and paired-end run for library $line[0](taxId: $line[2]). To process this library please first run the script $renameFastqScript and then generate again the generation of the BgeeCall input file";
      next LIBRARY;
    }
    my $output_line = "$line[2]\t\t$meanReadLength\t$fastq_path\t$transcriptome_path\t$annotation_path\t$lib_output_dir\t$intergenic_file\n";
    print {$FH} $output_line;

}
close($FH);
close($FH_missing);
close($FH_processed);

