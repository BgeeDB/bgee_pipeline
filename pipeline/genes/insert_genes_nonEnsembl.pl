#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use List::MoreUtils qw{uniq};
use List::Compare;
use LWP::Simple;
use File::Slurp;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Utils_insert_genes;

use Data::Dumper qw(Dumper);

# Patrick Tran Van, updated May 2019
# Insertion of non Ensembl genes
#############################################################

### Function : Formatting the GO annotation

sub go_formatting {

	# Sort the GO annotation:
	
	#$VAR1 = [
    #      'GO:0006412',
    #      'GO:0003735',
    #      'GO:0022627'
    #    ];
        
    # become:
    
    #$VAR1 = [
    #      'GO:0003735',
    #      'GO:0006412',
    #      'GO:0022627'
    #    ];
    		
	my $go_description = shift;
	my $altid_go = shift;
	my $obs_go = shift;
	
	# Transform the GO description into array
	my @go = split /=/, $go_description;
	my @go_list = split /,/, $go[1];
	
	#print join(", ", @go_list) . "\n\n";
	
	my @final_go_list;	# Array that take account of obsolete and alternative GO.
	my $go_id;

	foreach $go_id (@go_list) {
 	
	 	my $final_go_id = "";
		my $go_lower_id = $go_id;
		$go_lower_id =~ s/GO/go/;	# Replace GO by go
		
		if (not defined ($obs_go->{$go_lower_id})) {	# Skip obsolete GO, not inserted in Bgee previously
			
			if (defined ($altid_go->{$go_lower_id}))	# Replace alt_id GO by main GO if any
			{
				
				$final_go_id = $altid_go->{$go_lower_id};
				$final_go_id =~ s/go/GO/;
			}
			
			else {
				$final_go_id = $go_id;
			}
			
			#print $final_go_id . "\n";
			push @final_go_list, $final_go_id;
		}

	}


	# Sort the GO array
	my @sorted_go_list = sort @final_go_list;	
	#print join(", ", @sorted_go_list) . "\n\n";	
	return (@sorted_go_list);
	
}

#################### MAIN

# Define arguments & their default value
my ($species, $bgee_connector) = ('', '');
my ($obsGO, $annot_dir) = ('', '');
my ($debug) = (0);
my %opts = ('species=s'    => \$species,            # speciesCommonName from TSV for or Bgee db
            'bgee=s'       => \$bgee_connector,     # Bgee connector string
            'obsGO=s'      => \$obsGO,              # go.obsolete file
            'annot_dir=s'  => \$annot_dir,          # path to the directory containing fasta and GFF files of all non ensembl species
            'debug'        => \$debug,              # debug mode, do not insert/update in database
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $species eq '' || $bgee_connector eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0  -species=61478__0__Timema_douglasi__non_Ensembl  -bgee=\$(BGEECMD)  -obsGO=go.obsolete 
\t-species         speciesId from Bgee db with the genomeSpeciesId concatenated
\t-bgee            Bgee connector string
\t-obsGO           go.obsolete file
\tannot_dir'  	   Directory containing fasta and GFF files of all non ensembl species
\t-debug           Debug mode, do not insert/update in database
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

# Factor 1: Need to map to another genomeSpeciesId?

my ($species_tmp, $speciesBgee, $newSpecies, $scientific_name, $ensSource) = Utils_insert_genes::map_species($dbh, $species);

# 7227__0__Timema__Maker
# become 7227 0 Timema Maker

$species = $species_tmp;

## Factor 2: Get alternatives GO

my %altid_go = Utils_insert_genes::get_alt_go($dbh);

## Factor 3: Get Obsolete GO

my %obs_go = Utils_insert_genes::get_obs_go($dbh, $obsGO);

## Factor 4: Get used dataSources

my %InsertedDataSources = Utils_insert_genes::get_source_id($dbh);

#print Dumper \%InsertedDataSources;

## Gene info (id, description)
# Get individual gene info
my $geneDB       = $dbh->prepare('INSERT INTO gene (geneId, geneName, geneDescription, geneBioTypeId, speciesId)
                                  VALUES (?, ?, ?, (SELECT geneBioTypeId FROM geneBioType WHERE geneBioTypeName=?), ?)');
my $synonymDB    = $dbh->prepare('INSERT INTO geneNameSynonym (bgeeGeneId, geneNameSynonym)
                                  VALUES (?, ?)');
my $xrefDB       = $dbh->prepare('INSERT INTO geneXRef (bgeeGeneId, XRefId, XRefName, dataSourceId)
                                  VALUES (?, ?, ?, (SELECT dataSourceId FROM dataSource WHERE dataSourceName=?))');
my $goDB         = $dbh->prepare('INSERT INTO geneToGeneOntologyTerm (bgeeGeneId, goId, goEvidenceCode)
                                  VALUES (?, ?, ?)');
my $geneToTermDB = $dbh->prepare('INSERT INTO geneToTerm (bgeeGeneId, term)
                                  VALUES (?, ?)');
print "Inserting gene info...\n";
#GENE:

# Create an empty hash (dictionary)

my %gene_hash;

# Get the annotation file

my $get_annotation_path  = $dbh->prepare('SELECT genomeFilePath FROM species WHERE speciesId = (?)'); 
$get_annotation_path->execute($speciesBgee);
my $ann_file = "";
$ann_file = $get_annotation_path->fetchrow_array(); 
$get_annotation_path->finish();

$ann_file =~ s/.fasta/_vbgee.gff/; 
$ann_file = "$annot_dir$ann_file";

print "Parsing the annotation file: " . $ann_file . " ...\n";

# Go through the annotation file

open (FILE, $ann_file);

while (<FILE>) {
	chomp;
	
	my @field = split("\t+");
	my $type = $field[2];	# gene, mRNA, exon, CDS, five_prime_UTR, three_prime_UTR
	
	# Look only for genes
	if ($type eq "gene") {
		
		my @gene_description = split /;/, $field[8];

		my @id = split /=/, $gene_description[0];
		my $gene_id = $id[1];
		
		# Parsing ...
		
		my $target;
		my $idx;
		my $idx_taxon;
		my $id_species;
		my $gene_description_ref  = \@gene_description;
						
		# ... the GO annotation (1)

		$target = "ontology_term";
		$idx = undef;
		
		# Get the index of "ontology_term" in the array @gene_description
		
		foreach ( 0 .. $#gene_description ) {
			if ( index ( $gene_description_ref->[ $_ ], $target ) >= 0 ) {
				$idx = $_;
				last;
			}
		}
		
		if (defined $idx){
			
			my @annotation = split /=/, $gene_description[$idx];
			# ontology_term=GO:0005634,GO:0003677
			
			# become
			
			# $VAR1 = [
            # 'ontology_term',
            # 'GO:0005634,GO:0003677'
            # ];
            
			if ( scalar @annotation == 2 ) {			
				my @go = go_formatting($gene_description[$idx], \%altid_go, \%obs_go);
				$gene_hash{$gene_id}{GO} = \@go;	# Create a reference to the array by using backslash
			}
			
			else {
				$gene_hash{$gene_id}{GO} = "undef";
			}
			
		}			 
		
		else {	# Gene has no GO or "ontology_term" not present
			$gene_hash{$gene_id}{GO} = "undef";
		}	
		
		# ... the functional annotation 
		# Gene name (2)
	
		$target = "topblasthit_gene";
		$idx = undef;
		
		# Get the index of "topblasthit_gene" in the array @gene_description
		
		foreach ( 0 .. $#gene_description ) {
			if ( index ( $gene_description_ref->[ $_ ], $target ) >= 0 ) {
				$idx = $_;
				last;
			}
		}
		
		if (defined $idx){
			
			my @annotation = split /=/, $gene_description[$idx];

			if ( scalar @annotation == 2 ) {			
				$gene_hash{$gene_id}{Name} = $annotation[1];
			}
			
			else {
				$gene_hash{$gene_id}{Name} = "undef";
			}
			
		}			 
		
		else {	# Gene has no gene name or "topblasthit_gene" not present
			$gene_hash{$gene_id}{Name} = "undef";
		}	

		# Gene accession (3)
	
		$target = "topblasthit_xref";
		$idx = undef;
		
		# Get the index of "topblasthit_xref" in the array @gene_description
		
		foreach ( 0 .. $#gene_description ) {
			if ( index ( $gene_description_ref->[ $_ ], $target ) >= 0 ) {
				$idx = $_;
				last;
			}
		}
		
		if (defined $idx){
			
			my @annotation = split /=/, $gene_description[$idx];

			if ( scalar @annotation == 2 ) {			
				$gene_hash{$gene_id}{Accession} = $annotation[1];
			}
			
			else {
				$gene_hash{$gene_id}{Accession} = "undef";
			}
			
		}			 
		
		else {	# Gene has no gene name or "topblasthit_gene" not present
			$gene_hash{$gene_id}{Accession} = "undef";
		}	


		# Gene description (4) + Taxon best hit
	
		$target = "topblasthit_description";
		$idx = undef;
		
		# Get the index of "topblasthit_description" in the array @gene_description
		
		foreach ( 0 .. $#gene_description ) {
			if ( index ( $gene_description_ref->[ $_ ], $target ) >= 0 ) {
				$idx = $_;
				last;
			}
		}
		
		if (defined $idx){
			
			$target = "topblasthit_taxon";
			$idx_taxon = undef;
		
			# Get the index of "topblasthit_description" in the array @gene_description
		
			foreach ( 0 .. $#gene_description ) {
				if ( index ( $gene_description_ref->[ $_ ], $target ) >= 0 ) {
					$idx_taxon = $_;
					last;
				}
			}
			
			if (defined $idx_taxon){
				my @taxon_annotation = split /=/, $gene_description[$idx_taxon];
				
				if ( scalar @taxon_annotation == 2 ) {	
					$id_species = $taxon_annotation[1];
				}
				
				else {
					$id_species = "undef";
				}
			
			}
			
			else {
				$id_species = "undef";
			}
			
			my @annotation = split /=/, $gene_description[$idx];

			if ( scalar @annotation == 2 ) {			
				$gene_hash{$gene_id}{Description} = $annotation[1] . " [Hit tax id: " . $id_species . ", Acc: " . $gene_hash{$gene_id}{Accession} ."]";
				# RNA-binding protein fusilli isoform X1 [Hit tax id: Zootermopsis nevadensis, Acc: gi|1227956332|ref|XP_021916099.1|]
				# Include functional description, taxon hit and ncbi accessions.
			}
			
			else {
				$gene_hash{$gene_id}{Description} = "undef";
			}
			
		}			 
		
		else {	# Gene has no gene name or "topblasthit_description" not present
			$gene_hash{$gene_id}{Description} = "undef";
		}	
	 
	} 

}
close (FILE);


my $gene_name;
my $description;
my $accession;
my $biotype = "protein_coding";

#print Dumper \%gene_hash;


#=for comment

foreach my $id2insert (sort {lc $a cmp lc $b} keys %gene_hash) {	# Sort to always get the same order

	my $bgeeGeneId;

	if (ref($gene_hash{$id2insert}) eq 'HASH') {	# Gene has a function

		$gene_name = $gene_hash{$id2insert}{Name};
    	$description = $gene_hash{$id2insert}{Description};
    	$accession = $gene_hash{$id2insert}{Accession};
    	
#=for comment
		## Insert gene info

		if ( ! $debug ){
			$geneDB->execute($id2insert, $gene_name, $description, $biotype, $speciesBgee)  or die $geneDB->errstr;
			$bgeeGeneId = $dbh->{mysql_insertid};
			die "Cannot get bgeeGeneId [$bgeeGeneId]\n"  if ( $bgeeGeneId !~ /^\d+$/ );
		}
		
		else {
			print "\n[$id2insert] [] [$description]   [$biotype] [$speciesBgee]\n";
		}

		
		# Get UniProt ID and gene name, missing in B2G annotation
		
		#print Dumper \$gene_hash{$id2insert};
		
		my $uniprot_gene = "";
		my $content = "";
		my $uniprot_result = "";
		
		$content = get("https://www.uniprot.org/uniprot/?query=" . $accession . "&format=tab&columns=entry%20name");

		if ( defined $content ){
			(undef, $uniprot_result) = split("\n", $content, -1);
			
			if ( defined $uniprot_result ){ $uniprot_gene = $uniprot_result;}
			else { $uniprot_gene = "undef"; }
		}	
		else {
			warn "\tCannot get UniProt ID for [$accession]\n";
		}

		# Fill geneXRef table with uniprot information
		
		my $dbname = "UniProtKB/Swiss-Prot";
		
		#print $uniprot_gene;
		
		if ( ! $debug ){
	    	$xrefDB->execute($bgeeGeneId, $id2insert, $uniprot_gene, $dbname)  or die $xrefDB->errstr;
		}
		else {
		    print "xref: [$id2insert] [$id2insert] [$uniprot_gene] [$dbname]\n";
		}
 
        	
    	if (ref($gene_hash{$id2insert}{GO}) eq 'ARRAY') {	# Gene has GO terms
    		
			foreach my $go_id (@{$gene_hash{$id2insert}{GO}}){
								
					if ( ! $debug ){
					    $goDB->execute($bgeeGeneId, $go_id, "IEA")  or do {warn "[$id2insert] [$go_id] [Evidence code]\n"; die $goDB->errstr};
					}
					else {
					    print "go:   [$id2insert] [$go_id] [Evidence code]\n";
					}
			}
		}				    			
	}
	
	else {
		
		if ( ! $debug ){
			$geneDB->execute($id2insert, "", "", $biotype, $speciesBgee)  or die $geneDB->errstr;
			$bgeeGeneId = $dbh->{mysql_insertid};
			die "Cannot get bgeeGeneId [$bgeeGeneId]\n"  if ( $bgeeGeneId !~ /^\d+$/ );
		}
		
		else {
			print "\n[$id2insert] [] [$description]   [$biotype] [$speciesBgee]\n";
		}

	}
	
}
#=cut


# Update ensemblGene to 0
my $updt_geneDB  = $dbh->prepare('UPDATE gene SET ensemblGene = 0 WHERE speciesId = (?)');                             
$updt_geneDB->execute($speciesBgee);

$updt_geneDB->finish();
$geneDB->finish;
$synonymDB->finish;
$xrefDB->finish;
$goDB->finish;
$geneToTermDB->finish;
print "Gene nbr for $scientific_name: ", scalar keys %gene_hash, "\n\n";



# Close db connections
$dbh->disconnect;

exit 0;

