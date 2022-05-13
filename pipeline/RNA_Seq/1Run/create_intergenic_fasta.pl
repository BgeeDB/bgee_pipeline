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


## Define arguments & their default value
my ($sample_info_path, $transcriptomes_folder, $transcriptome_compression_ext, $sum_abundance_file_path, $gaussian_file_path, $ref_intergenic_dir, $other_intergenic_dir) = ('', '', '', '', '', '',  '');
my %opts = ('sample_info_path=s'       			=> \$sample_info_path,
            'transcriptomes_folder=s'  			=> \$transcriptomes_folder,
            'transcriptome_compression_ext=s'	=> \$transcriptome_compression_ext,
            'sum_abundance_file_path=s'			=> \$sum_abundance_file_path,
			'gaussian_file_path=s'				=> \$gaussian_file_path,
			'ref_intergenic_dir=s'				=> \$ref_intergenic_dir,
			'other_intergenic_dir=s'			=> \$other_intergenic_dir
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $sample_info_path eq '' || $transcriptomes_folder eq '' || $transcriptome_compression_ext eq '' || $sum_abundance_file_path eq '' || $gaussian_file_path eq '' || $ref_intergenic_dir eq '' || $other_intergenic_dir eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0 -sample_info_path=RNASEQ_SAMPLE_INFO-transcriptomes_folder=RNASEQ_VITALIT_GTF -transcriptome_compression_ext=TRANSCRIPTOME_FILE_EXT -sum_abundance_file_path=SUM_ABUNDANCE_FILE_PATH -gaussian_file_path=RNASEQ_VITALIT_GAUSSIAN_CHOICE -ref_intergenic_dir=VITALIT_REF_INTERGENIC_FOLDER -other_intergenic_dir=VITALIT_OTHER_INTERGENIC_FOLDER >> $@.tmp 2> $@.warn
\t-sample_info_path      			Path to the rna_seq_sample_info.txt file
\t-transcriptomes_folder=s     		Folder where transcriptome fasta files are stored.
\t-transcriptome_compression_ext=s  Extension of compressed transcriptome files. (e.g .xz)
\t-sum_abundance_file_path=s     	Path to the abundance files created during the sum_by_species step of the RNA-Seq Bgee pipeline. SPECIES_ID pattern in the file name will be automatically changed to the corresponding species id
\t-gaussian_file_path=s				Path to the gaussian file
\t-ref_intergenic_dir=s'			Path to the directory where reference intergenic fasta files are created
\t-other_intergenic_dir=s'			Path to the directory where non reference (other) intergenic fasta files are created
\n";
    exit 1;
}

######################## Preliminary steps ########################


## Retrieve sample_info information
my %speNameToSpeId;
for my $line ( read_file($sample_info_path, chomp => 1) ){
	# do not parse the header
	next  if ( $line =~ m/^#.*/);
	my @fields = split("\t", $line);
	my $speciesId = $fields[2];
	my $genomeFilePath = $fields[4];
	$genomeFilePath =~ m/.+\/(.+)/;
	$speNameToSpeId{$1} = $speciesId;
}

## Retrieve gaussian choice information
my %gaussianChoice;
for my $line ( read_file($gaussian_file_path, chomp => 1) ){
	# do not parse the header
	next  if ( $line =~ m/^speciesId.*/);
	my @fields = split("\t", $line);
	# if not all fields are present in the gaussian choice file
	if(scalar @fields != 9) {
		die "Each line of the gaussian file should have 9 tabular seprarated columns : $?\n";
	}
	my $speciesId = $fields[0];
	$gaussianChoice{$speciesId}->{'numberGaussiansIntergenic'} = $fields[2];
	$gaussianChoice{$speciesId}->{'selectedGaussiansIntergenic'} = $fields[5];
	$gaussianChoice{$speciesId}->{'selectionSideIntergenic'} = $fields[6];

}

## open directory containing all compressed transcriptome fasta files
die "Invalid or missing [$transcriptomes_folder]: $?\n"  if ( !-e $transcriptomes_folder || !-d $transcriptomes_folder );
opendir (DIR, $transcriptomes_folder) or die $!;


########################          Functions            ########################

## function that create coordinate tsv file of intergenic regions
## This file will contain 3 columns : chr, start, end
sub create_coordinate_file {
	my ($coordinate_file, @intergenic_ids)= @_;
	print "Start creation of : $coordinate_file\n";
	open(my $df, '>', $coordinate_file);
	print $df "chr\tstart\tend\n";
	for my $intergenic_id (0 .. $#intergenic_ids) {
    	my @splitted = split("_", $intergenic_ids[$intergenic_id]);
	# if the array contains less than 3 values it means it is not an intergenic region
    	next if(scalar @splitted < 3);
	my $end = $splitted[$#splitted];
    	my $start = $splitted[$#splitted - 1];
    	#remove 2 last elements of the array
    	splice @splitted, -2;
    	my $chr_name = join('_', @splitted);
    	print $df "$chr_name\t$start\t$end\n";
	}
	close $df;
}

######################## Generate intergenic sequences ########################

# for each compressed transcriptome fasta file
foreach my $species_name (keys %speNameToSpeId){
	print "start generation of intergenic sequences for $species_name";
	my $file_path = $transcriptomes_folder."\/".$species_name."\.transcriptome\.fa".$transcriptome_compression_ext;

	# remove archive extension to file name
	my $uncompressed_file = substr($file_path, 0, length($file_path) - length($transcriptome_compression_ext));
	# if file already exist we remove it
	if (-e $uncompressed_file) {
		unlink $uncompressed_file;
	}


	## uncompress xz archive
	if ($transcriptome_compression_ext =~ /\.xz/) {
		print("Uncompress archive : $file_path\n");
		# uncompress
		if (system("unxz -k $file_path") != 0) {
			exit "Can not uncompress xz file : $file_path";
		}
	} else {
		exit("The script only uncompress .xz archives");
	}

	my $speciesId = $speNameToSpeId{$species_name};

	## retrieve ids of reference and other intergenic sequences from gaussian
	## choice file and sum by species file.
	(my $sum_abundance_file_species = $sum_abundance_file_path) =~ s/SPECIES_ID/$speciesId/g;
	my @ref_intergenic_ids;
	my @other_intergenic_ids;
	my $selectedRefIntergenic = $gaussianChoice{$speciesId}{'selectedGaussiansIntergenic'};
	my $selectionSideIntergenic = $gaussianChoice{$speciesId}{'selectionSideIntergenic'};
	## init with value allowing to detect tpmThreshold
	my $tpmCutoff = 9**9**9;
	if ($selectionSideIntergenic eq "Left") {
		$tpmCutoff = -$tpmCutoff;
	}
	# parse sum abundance file to calculate tpm cutoff
	for my $line ( read_file($sum_abundance_file_species, chomp => 1) ) {
		next  if ( ($line =~ m/^speciesId/));
		my @fields = split("\t", $line);
		# if not all fields are present in the gaussian choice file
		if(scalar @fields != 7) {
			die "Each line of file $sum_abundance_file_species should have 7 tabular seprarated columns : $?\n";
		}
		if ($fields[6] eq "intergenic_$selectedRefIntergenic") {
			if ($selectionSideIntergenic eq "Left" && $tpmCutoff < $fields[2]) {
				$tpmCutoff = $fields[2];
			} elsif ($selectionSideIntergenic eq "Right" && $tpmCutoff > $fields[2]) {
				$tpmCutoff = $fields[2];
			}
		}
	}
	# parse sum abundance file to detect ref and other intergenic ids
	for my $line ( read_file($sum_abundance_file_species, chomp => 1) ) {
		next  if ( ($line =~ m/^speciesId/));
		my @fields = split("\t", $line);
		if ($fields[4] eq "intergenic") {
			if ($selectionSideIntergenic eq "Right") {
				if ($fields[2] < $tpmCutoff) {
					push @ref_intergenic_ids, $fields[0];
				} else {
					push @other_intergenic_ids, $fields[0];
				}
			} elsif ($selectionSideIntergenic eq "Left") {
				if ($fields[2] <= $tpmCutoff) {
					push @ref_intergenic_ids, $fields[0];
				} else {
					push @other_intergenic_ids, $fields[0];
				}
			} else {
				exit "should select Right or Left side of the gaussian but selected $selectionSideIntergenic";
			}
		}
	}
	print "number ref intergenic sequences : ".scalar @ref_intergenic_ids."\n";
	print "number other intergenic sequences : ".scalar @other_intergenic_ids."\n";

	## Create coordinate tsv file for ref. intergenic and other intergenic sequences.
	my $ref_intergenic_coordinate_file = $ref_intergenic_dir.$speciesId."_coordinates.tsv";
	my $other_intergenic_coordinate_file = $other_intergenic_dir.$speciesId."_coordinates.tsv";
	create_coordinate_file($ref_intergenic_coordinate_file, @ref_intergenic_ids);
	create_coordinate_file($other_intergenic_coordinate_file, @other_intergenic_ids);

	## Create intergenic fasta files
	my $ref_intergenic_file = $ref_intergenic_dir.$speciesId."_intergenic.fa";
	my $other_intergenic_file = $other_intergenic_dir.$speciesId."_other_intergenic.fa";

	print "Start creation of intergenic sequence files for species $speciesId\n";
	my $seq_in  = Bio::SeqIO->new(-file => "$uncompressed_file", -format => "fasta");
	my $ref_intergenic_out = Bio::SeqIO->new(-file => ">$ref_intergenic_file", -format => "fasta");
	my $other_intergenic_out = Bio::SeqIO->new(-file => ">$other_intergenic_file", -format => "fasta");
	while(my $seq = $seq_in->next_seq) {
		my $transcript_id = $seq->primary_id;
		# keep only intergenic regions
		if ( grep {$_ eq $transcript_id} @ref_intergenic_ids) {
			$ref_intergenic_out->write_seq($seq);
  		} elsif ( grep {$_ eq $transcript_id} @other_intergenic_ids) {
			$other_intergenic_out->write_seq($seq);
  		}
 	}
  	print "intergenic files have been created successfully for species $speciesId\n";
  	# remove uncompressed transcriptomic file
  	unlink $uncompressed_file;
}
closedir DIR;
