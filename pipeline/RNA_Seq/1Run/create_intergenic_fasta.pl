#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use File::Path qw(make_path);
use File::Slurp;
use Getopt::Long;
use Bio::SeqIO;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;


## Julien Wollbrett 16/09/2019
## This script allows creation of intergenic fasta files for all species.
## It generates reference and non reference (called other) intergenic fasta files.
## Basically, this script will :
##	1. Uncompress transcriptome files (keep compressed files)
##  2. Use gaussian choice and sum by species files to distinguish reference and other intergenic regions
##	3. create tsv files with coordinates of intergenic regions
##	4. create intergenic fasta files
##	5. remove uncompressed transcriptome fasta files

## Alessandro Brandulas Cammarata 30/10/2024
## Changed step 2 to use the classification column of the sum abundance file instead of the gaussian choice file to filter reference intergenic regions

## Define arguments & their default values
my ($sample_info_path, $transcriptomes_folder, $transcriptome_compression_ext, $sum_abundance_file_path, $ref_intergenic_dir, $other_intergenic_dir) = ('', '', '', '', '', '');
my %opts = (
    'sample_info_path=s'              => \$sample_info_path,
    'transcriptomes_folder=s'         => \$transcriptomes_folder,
    'transcriptome_compression_ext=s' => \$transcriptome_compression_ext,
    'sum_abundance_file_path=s'       => \$sum_abundance_file_path,
    'ref_intergenic_dir=s'            => \$ref_intergenic_dir,
    'other_intergenic_dir=s'          => \$other_intergenic_dir
);

######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $sample_info_path eq '' || $transcriptomes_folder eq '' || $transcriptome_compression_ext eq '' || $sum_abundance_file_path eq '' || $ref_intergenic_dir eq '' || $other_intergenic_dir eq '') {
    print "\n\tInvalid or missing argument:
    \te.g. $0 -sample_info_path=RNASEQ_SAMPLE_INFO -transcriptomes_folder=RNASEQ_VITALIT_GTF -transcriptome_compression_ext=TRANSCRIPTOME_FILE_EXT -sum_abundance_file_path=SUM_ABUNDANCE_FILE_PATH -ref_intergenic_dir=VITALIT_REF_INTERGENIC_FOLDER -other_intergenic_dir=VITALIT_OTHER_INTERGENIC_FOLDER >> \$@.tmp 2> \$@.warn
    \t-sample_info_path                Path to the rna_seq_sample_info.txt file
    \t-transcriptomes_folder           Folder where transcriptome FASTA files are stored.
    \t-transcriptome_compression_ext   Extension of compressed transcriptome files (e.g., .xz)
    \t-sum_abundance_file_path         Path to the abundance files created during the sum_by_species step of the RNA-Seq pipeline. SPECIES_ID pattern in the file name will be automatically changed to the corresponding species ID
    \t-ref_intergenic_dir              Path to the directory where reference intergenic FASTA files are created
    \t-other_intergenic_dir            Path to the directory where non-reference (other) intergenic FASTA files are created
    \n";
    exit 1;
}

######################## Preliminary steps ########################

## Retrieve sample_info information
my %speNameToSpeId;
for my $line ( read_file($sample_info_path, chomp => 1) ){
    # Do not parse the header
    next if ( $line =~ m/^#.*/);
    my @fields = split("\t", $line);
    my $speciesId = $fields[2];
    my $genomeFilePath = $fields[4];
    $genomeFilePath =~ m/.+\/(.+)/;
    $speNameToSpeId{$1} = $speciesId;
}

## Open directory containing all compressed transcriptome FASTA files
die "Invalid or missing [$transcriptomes_folder]: $?\n" if ( !-e $transcriptomes_folder || !-d $transcriptomes_folder );
opendir (DIR, $transcriptomes_folder) or die $!;

########################          Functions            ########################

## Function that creates coordinate TSV file of intergenic regions
## This file will contain 3 columns: chr, start, end
sub create_coordinate_file {
    my ($coordinate_file, @intergenic_ids) = @_;
    print "Start creation of: $coordinate_file\n";
    open(my $df, '>', $coordinate_file);
    print $df "chr\tstart\tend\n";
    for my $intergenic_id (0 .. $#intergenic_ids) {
        my @splitted = split("_", $intergenic_ids[$intergenic_id]);
        # If the array contains less than 3 values, it means it is not an intergenic region
        next if (scalar @splitted < 3);
        my $end = $splitted[$#splitted];
        my $start = $splitted[$#splitted - 1];
        # Remove 2 last elements of the array
        splice @splitted, -2;
        my $chr_name = join('_', @splitted);
        print $df "$chr_name\t$start\t$end\n";
    }
    close $df;
}

######################## Generate intergenic sequences ########################

# For each compressed transcriptome FASTA file
foreach my $species_name (keys %speNameToSpeId){
    print "Start generation of intergenic sequences for $species_name\n";
    my $file_path = $transcriptomes_folder . "/" . $species_name . ".transcriptome.fa" . $transcriptome_compression_ext;

    # Remove archive extension to get uncompressed file name
    my $uncompressed_file = substr($file_path, 0, length($file_path) - length($transcriptome_compression_ext));
    # If file already exists, remove it
    if (-e $uncompressed_file) {
        unlink $uncompressed_file;
    }

    ## Uncompress xz archive
    if ($transcriptome_compression_ext =~ /\.xz/) {
        print("Uncompress archive: $file_path\n");
        # Uncompress
        if (system("unxz -k $file_path") != 0) {
            die "Cannot uncompress xz file: $file_path\n";
        }
    } else {
        die "The script only uncompresses .xz archives\n";
    }

    my $speciesId = $speNameToSpeId{$species_name};

    ## Retrieve IDs of reference and other intergenic sequences from the sum_abundance_file.
    (my $sum_abundance_file_species = $sum_abundance_file_path) =~ s/SPECIES_ID/$speciesId/g;
    my @ref_intergenic_ids;
    my @other_intergenic_ids;

    # Parse sum abundance file to detect ref and other intergenic IDs
    for my $line ( read_file($sum_abundance_file_species, chomp => 1) ) {
        next if ( ($line =~ m/^speciesId/) || ($line =~ m/^gene_id/));  # Skip headers
        my @fields = split("\t", $line);
        # Check that the line has at least 7 columns
        if (scalar @fields < 7) {
            die "Each line of file $sum_abundance_file_species should have at least 7 tab-separated columns: $?\n";
        }
        my $gene_id = $fields[0];
        my $type = $fields[4];
        my $classification = $fields[6];
        if ($type eq "intergenic") {
            if ($classification eq "Reference_intergenic") {
                push @ref_intergenic_ids, $gene_id;
            } elsif ($classification eq "Outlier_intergenic") {
                push @other_intergenic_ids, $gene_id;
            }
        }
    }

    print "Number of reference intergenic sequences: " . scalar @ref_intergenic_ids . "\n";
    print "Number of other intergenic sequences: " . scalar @other_intergenic_ids . "\n";

    ## Create coordinate TSV files for reference and other intergenic sequences.
    my $ref_intergenic_coordinate_file = $ref_intergenic_dir . $speciesId . "_coordinates.tsv";
    my $other_intergenic_coordinate_file = $other_intergenic_dir . $speciesId . "_coordinates.tsv";
    create_coordinate_file($ref_intergenic_coordinate_file, @ref_intergenic_ids);
    create_coordinate_file($other_intergenic_coordinate_file, @other_intergenic_ids);

    ## Create intergenic FASTA files
    my $ref_intergenic_file = $ref_intergenic_dir . $speciesId . "_intergenic.fa";
    my $other_intergenic_file = $other_intergenic_dir . $speciesId . "_other_intergenic.fa";

    print "Start creation of intergenic sequence files for species $speciesId\n";
    my $seq_in  = Bio::SeqIO->new(-file => "$uncompressed_file", -format => "fasta");
    my $ref_intergenic_out = Bio::SeqIO->new(-file => ">$ref_intergenic_file", -format => "fasta");
    my $other_intergenic_out = Bio::SeqIO->new(-file => ">$other_intergenic_file", -format => "fasta");

    my %ref_intergenic_hash = map { $_ => 1 } @ref_intergenic_ids;
    my %other_intergenic_hash = map { $_ => 1 } @other_intergenic_ids;

    while(my $seq = $seq_in->next_seq) {
        my $transcript_id = $seq->primary_id;
        # Keep only intergenic regions
        if ( exists $ref_intergenic_hash{$transcript_id} ) {
            $ref_intergenic_out->write_seq($seq);
        } elsif ( exists $other_intergenic_hash{$transcript_id} ) {
            $other_intergenic_out->write_seq($seq);
        }
    }
    print "Intergenic files have been created successfully for species $speciesId\n";
    # Remove uncompressed transcriptome file
    unlink $uncompressed_file;
}
closedir DIR;