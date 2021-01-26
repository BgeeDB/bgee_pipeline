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
#- add gene_biotype attribute at the end of lines corresponding to exon: in RefSeq gtf gene_biotype attribute is only present in the gene lines. Those lines are not kept in the gtf_all file generated later in the pipeline. It is then mandatory to add this info to exon lines (in order to be able to retrieve the mapping between gene and biotypes)

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;
use File::Basename;
use Data::Dumper;

## Define arguments & their default value
my ($path_to_gtf_folder, $sample_info_file) = ('','');
my %opts = ('path_to_gtf_folder=s'    => \$path_to_gtf_folder,
            'sample_info_file=s'      => \$sample_info_file  
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $path_to_gtf_folder eq '' || $sample_info_file eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0 -path_to_gtf_folder=CLUSTER_GTF_DIR -sample_info_file=$(SAMPLE_INFO_PATH)>> $@.tmp 2> $@.warn
\t-path_to_gtf_folder     Path to the directory containing all gtf annotation files
\t-sample_info_file       Path to the sample_info file
\n";
    exit 1;
}

# detect species coming from RefSeq/NCBI using information in the sample_info file
my %ref_seq_path;
open(my $sample_info, $sample_info_file) || die "failed to read sample_info file: $!";
while (my $line = <$sample_info>) {
    chomp $line;
    next if( ($line =~ m/^#/) or ($line =~ m/^\"#/) );
    my @fields = split "\t" , $line;
    my $number_columns = 11;
    if (! scalar @fields eq $number_columns) {
        die "all lines of sample info file should have $number_columns columns";
    }
    if($fields[5] eq 'RefSeq') {
        $ref_seq_path{$fields[4]}= "";
    }
}

print keys %ref_seq_path;
print "\n";


foreach my $key (keys %ref_seq_path) {
    my $file_prefix = basename($key);
    my @files = glob "$path_to_gtf_folder/$file_prefix*.gtf.gz";

    if(!defined($files[0])) {
        warn "no gtf.gz file map regex $path_to_gtf_folder/$file_prefix*.gtf.gz";
        next;
    } elsif(scalar @files > 1) {
        die "more than one gtf.gz file map regex $path_to_gtf_folder/$file_prefix*.gtf.gz";
    }
    my $file = $files[0];
    print "parse file $file\n";
    #create temp file where updated GTF will be store
    my $tempFile = "$path_to_gtf_folder/temp.gtf";
    
    #preliminary read of gtf file to check that gtf file was not already updated
    open(IN, "gunzip -c $file |") || die "can’t open pipe to $file";
    my $already_modified = 0;
    while (my $line = <IN>) {
        next if ($line =~ /^#/);
        # check if firest exon line already has a gene_biotype attribute
        if ($line =~ /\texon\t/) {
            if ($line =~ /(gene_biotype [^;]+)/) {
                $already_modified = 1;
            }
            last;
        }
    }
    close(IN);

    if ($already_modified) {
        print "file $file already modified\n";
        next;
    } else {
        print "file $file needs to be modified\n";
    }

    #first read of gtf file to retrieve mapping gene<-biotype
    open(IN, "gunzip -c $file |") || die "can’t open pipe to $file";
    my %geneToBiotype;
    while(my $line = <IN>) {
        next if ($line =~ /^#/);
        # retrieve mapping gene_id <- gene_biotype
        my @gene_id = $line =~ /\t(gene_id [^;]+)/;
        if ($line =~ /\tgene\t/) {
            my @gene_biotype = $line =~ /(gene_biotype [^;]+)/;
            $geneToBiotype{$gene_id[0]}{'biotype'} = $gene_biotype[0];
            my @gene_xref = $line =~ /"GeneID:([^"]+)/;
            $geneToBiotype{$gene_id[0]}{'xref'} = $gene_xref[0];
        }
    }
    close(IN);

    print Dumper (\%geneToBiotype);

    #second read of gtf file to update it
    open(IN, "gunzip -c $file |") || die "can’t open pipe to $file";
    open (my $OUT, '>', "$tempFile")  or die "Cannot write [$tempFile]\n";
    while(my $line = <IN>) {
        # remove lines containing "unknown_transcript_1"
        next if (index($line, "unknown_transcript_1") != -1);
        # remove lines starting with #
        next if ($line =~ /^#/);
        # remove transcript_id = "" when exist
        $line =~ s/transcript_id ""; //;
        #add gene_biotype info to lines corresponding to exon
        my @gene_id = $line =~ /\t(gene_id [^;]+)/;
        chomp($line);
        # remove gene name as gene_id and put the real gene_id
        my $gene_xref_tag = "gene_id $geneToBiotype{$gene_id[0]}{'xref'}";
        $line =~ s/$gene_id[0]/$gene_xref_tag/g;
        my @gene_name = split "\"", $gene_id[0];
        #split at each " character in order to only keep the gene name itself
        my $gene_name_tag = "gene_name \"".$gene_name[1]."\";";
        $line .= $gene_name_tag;
        if ($line =~ /\texon\t/) {
            $line .= " $geneToBiotype{$gene_id[0]}{'biotype'};";
        }
        $line .= "\n";
        print {$OUT} $line;
    }
    close(IN);
    close $OUT;
    #gzip temp file
    my $tempFileGz = "$tempFile.gz";
    gzip $tempFile => $tempFileGz or die "gzip failed: $GzipError\n";
    # move temp file to original file
    rename $tempFileGz, $file or die "can not rename file $tempFileGz to $file";
    #delete uncompressed temp file
    unlink $tempFile;
    print "finished updating file $file\n";
}



