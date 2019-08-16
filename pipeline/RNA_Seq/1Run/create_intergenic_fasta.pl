#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Data::Dumper;
use File::Path qw(make_path);
use Getopt::Long;
use Bio::SeqIO;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Data::Dumper;



## Julien Wollbrett 16/09/2019
## This script allows creation of intergenic fasta files for all species.
## It generates reference and non reference (called other) intergenic fasta files.
## Basically, this script will :
##	1. Uncompress transcriptome files (keep compressed files)
##  2. Use gaussian choice and sum by species files to distinguish reference and other intergenic regions
##	3. create intergenic fasta files
##	4. remove uncompressed transcriptome fasta files


## Define arguments & their default value
my ($bgee_connector, $ens_release, $transcriptomes_folder, $transcriptome_suffix, $transcriptome_compression_ext, $ref_intergenic_dir) = ('', '', '', '', '');
my %opts = ('bgee=s'       						=> \$bgee_connector,     # Bgee connector string
			'ens_release=s'      				=> \$ens_release,
            'transcriptomes_folder=s'  			=> \$transcriptomes_folder, 
            'transcriptome_compression_ext=s'	=> \$transcriptome_compression_ext, 
            'sum_abundance_file_path=s'			=> \$sum_abundance_file_path
			'gaussian_file_path=s'				=> \$gaussian_file_path
			'ref_intergenic_dir=s'				=> \$ref_intergenic_dir
			'other_intergenic_dir=s'			=> \$other_intergenic_dir
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $ens_release eq '' || $transcriptomes_folder eq '' || $transcriptome_compression_ext eq '' || $sum_abundance_file_path eq '' ||
	$gaussian_file_path eq '' || $ref_intergenic_dir eq '' || $other_intergenic_dir eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0 -ens_release=... -transcriptomes_folder=\$(OUTPUT_DIR) -transcriptome_suffix=\$(TRANSCRIPTOME_FILE_SUFFIX)  -transcriptome_compression_ext=\$(TRANSCRIPTOME_FILE_EXT) -intergenic_pattern=\$(INTERGENIC_PATTERN)
\t-bgee      						Bgee connector string
\t-ens_release=s          			Release version of ensembl (e.g 84).
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

## Retrieve species information from Bgee database

my $dbh = Utils::connect_bgee_db($bgee_connector);
# retrieve a map where key correspond to species and genus of species (e.g Homo_sapiens) 
# and value correspond to ensembl id of the species
my $bgeeSpecies = $dbh->prepare('SELECT genus, species, speciesId FROM species');
$bgeeSpecies->execute()  or die $bgeeSpecies->errstr;
my %speNameToSpeId;
while(my @row = $bgeeSpecies->fetchrow_array()){
	# Some species names in bgee db contain more than one word (e.g lupus familiaris). 
	# For bgee transcriptome files, the species name always correspond to the last word of the name (e.g familiaris)
	# That's why we only keep the last word of species name
	my @speciesNames = split(/ /,$row[1]);
	my $species = $speciesNames[(scalar @speciesNames)-1];
	my $speciesFileName = "$row[0]_$species";
	$speNameToSpeId {$speciesFileName} = $row[2];
}

## Retrieve gaussian choice information
my %gaussiansInfo;
my $header = 0;
for my $line ( read_file($gaussian_file_path, chomp => 1) ){
	# do not parse the header
	next  if ( ($line =~ m/^speciesId/);
	my @fields = split("\t", $line);
	# if not all fields are present in the gaussian choice file
	if(scalar @fields != 10) {
		die "Each line of the gaussian file should have 10 tabular seprarated columns : $?\n";
	}
	my $speciesId = $fields[0];
	$rnaSeqLibrary{$speciesId}->{'numberGaussiansIntergenic'} = $fields[4];
	$rnaSeqLibrary{$speciesId}->{'selectedGaussiansIntergenic'} = $fields[6];
	$rnaSeqLibrary{$speciesId}->{'selectionSideIntergenic'} = $fields[7];
	
}

## open directory containing all compressed transcriptome fasta files
die "Invalid or missing [$transcriptomes_folder]: $?\n"  if ( !-e $transcriptomes_folder || !-d $transcriptomes_folder );
opendir (DIR, $transcriptomes_folder) or die $!;

######################## Generate intergenic sequences ########################

# for each compressed transcriptome fasta file
while (my $file = readdir(DIR)) {
	if ($file =~ /.+$ens_release\.transcriptome\.fa$transcriptome_compression_ext$/) {
		my $file_path = "$transcriptomes_folder$file";
		
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
		
		## retrieve speciesId from file name thanks to species info from the DB
		$file =~ /^([^_]+_[a-zA-Z]+).*/;
		my $speciesId = $speNameToSpeId{$1};
		
		
		## retrieve ids of reference and other intergenic sequences from gaussian
		## choice file and sum by species file.
		$sum_abundance_file_path =~ s/SPECIES_ID/$speciesId/g;
		my @ref_intergenic_ids;
		my @other_intergenic_ids;
		my $selectedRefIntergenic = $rnaSeqLibrary{$speciesId}{'selectedGaussiansIntergenic'}
		my $selectionSideIntergenic = $rnaSeqLibrary{$speciesId}{'selectionSideIntergenic'}
		## init with value allowing to detect tpmThreshold
		my $tpmCutoff = 9**9**9;
		if ($selectionSideIntergenic eq "Left") {
			$tpmCutoff = -$tpmCutoff;
		}
		# parse sum abundance file to calculate tpm cutoff
		for my $line ( read_file($sum_abundance_file_path, chomp => 1) ) {
			next  if ( ($line =~ m/^speciesId/));
			my @fields = split("\t", $line);
			# if not all fields are present in the gaussian choice file
			if(scalar @fields != 7) {
				die "Each line of file $sum_abundance_file_path should have 7 tabular seprarated columns : $?\n";
			}
			if ($fields[6] eq "intergenic_$selectedRefIntergenic") {
				if ($selectionSideIntergenic eq "Left" && $tpmCutoff < $fields[2]) {
					$tpmCutoff = $fields[2];
				} else if ($selectionSideIntergenic eq "Right" && $tpmCutoff > $fields[2]) {
					$tpmCutoff = $fields[2];
				}
			}
		}
		# parse sum abundance file to detect ref and other intergenic ids
		for my $line ( read_file($sum_abundance_file_path, chomp => 1) ) {
			next  if ( ($line =~ m/^speciesId/));
			my @fields = split("\t", $line);
			if ($fields[5] eq "intergenic") {
				if ($selectionSideIntergenic eq "Right") {
					if ($fields[2] < $tpmCutoff) {
						push @ref_intergenic_ids, fields[0];
					} else {
						push @other_intergenic_ids, fields[0];
					}
				} else if ($selectionSideIntergenic eq "Left") {
					if ($fields[2] <= $tpmCutoff) {
						push @ref_intergenic_ids, fields[0];
					} else {
						push @other_intergenic_ids, fields[0];
					}
				}
			}
		}
		
		## Create intergenic fasta files
		if(!-e $ref_intergenic_dir) {mkdir $ref_intergenic_dir;}
		if(!-e $other_intergenic_dir) {mkdir $other_intergenic_dir;}
		my $ref_intergenic_file = $ref_intergenic_dir.$speciesId."_intergenic.fa";
		my $other_intergenic_file = $other_intergenic_dir.$speciesId."_other_intergenic.fa";

		print "Start creation of : $intergenic_file\n";
		my $seq_in  = Bio::SeqIO->new(-file => "$uncompressed_file", -format => "fasta");
		my $ref_intergenic_out = Bio::SeqIO->new(-file => ">$ref_intergenic_file", -format => "fasta");
		my $other_intergenic_out = Bio::SeqIO->new(-file => ">$other_intergenic_file", -format => "fasta");
		while(my $seq = $seq_in->next_seq) {
			my $transcript_id = $seq->primary_id;
			# keep only intergenic regions
			if ( grep {$_ eq $transcript_id} @ref_intergenic_ids) {
				$ref_intergenic_out->write_seq($seq);
  			} else if ( grep {$_ eq $transcript_id} @other_intergenic_ids) {
				$other_intergenic_out->write_seq($seq);
  			}   			
 		}
  		print "$intergenic files have been created successfully for species $speciesId\n";
  		# remove uncompressed transcriptomic file
  		unlink $uncompressed_file;
  		
	}
}
closedir DIR;