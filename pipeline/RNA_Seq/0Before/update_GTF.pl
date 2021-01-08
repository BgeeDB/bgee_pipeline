#!/usr/bin/env perl

## Julien Wollbrett, Jan 7, 2021
# This script modify gtf files of species coming from ncbi
# These modifications are mandatory to generate a fasta transcriptome file using the gtf and corresponding genome.
# It is a hack to allow to run gtf_to_fasta from TopHat (and in the future gffread from cufflinks) for species coming from ncbi
# modifications to the gtf file are :
#- delete empty transcript_id "" attribute : in ensembl transcript_id attribut is not present for gene feature. This attribute is present in ncbi gene feature with an empty value. It has to be removed otherwise gffread and gtf_to_fasta detect the line as transcript info and generate the sequence for the full gene.
#- delete line when transcript_id attribute is associated to value "unknown_transcript_1". As for Bgee 15 unknown transcript ids are always tag with unknown_transcript_1. they can be correspond to transcript
#    1) "tRNA" or "rRNA" 
#    2) linked to the exception "rearrangement required for product" (ebi.ac.uk/ena/WebFeat/qualifiers/exception.html)
#    3) created using The Vertebrate Mitochondrial Code (transl_table=2)

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use IO::Compress::Gzip qw(gzip $GzipError) ;

## Define arguments & their default value
my ($path_to_gtf_folder) = ('',);
my %opts = ('path_to_gtf_folder=s'                => \$path_to_gtf_folder
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $path_to_gtf_folder eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0 -path_to_gtf_folder=CLUSTER_GTF_DIR >> $@.tmp 2> $@.warn
\t-path_to_gtf_folder               Path to the directory containing all gtf annotation files
\n";
    exit 1;
}

#XXX Could filter on source of data to only select species coming from ncbi but it requires to query the database. It is maybe better not to force species info to already be in the database at this point (to confirm)
#retrieve all gtf files
my @files = glob "$path_to_gtf_folder/*.gtf.gz";
my $tempFile = "$path_to_gtf_folder/temp.gtf";
foreach my $file (@files) {
    #create temp file where updated GTF will be store
    open (my $OUT, '>', "$tempFile")  or die "Cannot write [$tempFile]\n";
    # start reading orginal GTF
    open(IN, "gunzip -c $file |") || die "can’t open pipe to $file";
    while(my $line = <IN>) {
        # remove lines containing "unknown_transcript_1"
        next if (index($line, "unknown_transcript_1") != -1);
        # remove transcript_id = "" when exist
        $line =~ s/transcript_id ""; //;
        print {$OUT} $line;

    }
    close(IN);
    close $OUT;
    #gzip temp fil
    my $tempFileGz = "$tempFile.gz";
    gzip $tempFile => $tempFileGz or die "gzip failed: $GzipError\n";
    # move temp file to original file
    rename $tempFileGz, $file or die "can not rename file $tempFileGz to $file";
    #delete uncompressed temp file
    unlink $tempFile;
}



