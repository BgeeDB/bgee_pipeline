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

###############################################################################
# Patrick Tran Van, created April 2019
#
# Insertion of non Ensembl genes.
###############################################################################

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

### Function : Formatting the functional annotation

sub func_formatting {
	
	my ($func_description) = @_; 
	
	my @annotation = split /=/, $func_description;
	# 'arth_topblasthit=gi|1330931529|gb|PNF42263.1|cGMP-dependent protein kinase, isozyme 1 [Cryptotermes secundus]'
	
	# become
	
	# $VAR1 = [
    #      'arth_topblasthit',
    #      'gi|1330931529|gb|PNF42263.1|cGMP-dependent protein kinase, isozyme 1 [Cryptotermes secundus]'
    #    ];

	my $original_annotation = $annotation[1];
	# 'gi|1330931529|gb|PNF42263.1|cGMP-dependent protein kinase, isozyme 1 [Cryptotermes secundus]'
	
	my $index_sep_acc_annot = rindex($original_annotation,"|");
	
	# rindex STR,SUBSTR :
	# Returns the position of the last occurrence of SUBSTR in STR. If POSITION is specified, returns the last occurrence beginning at or before that position.

	my $all_acc = substr($original_annotation, 0, $index_sep_acc_annot+1);
	my $annot = substr($original_annotation, $index_sep_acc_annot+1);
	
	# substr EXPR,INDEX,LENGTH:
	# Extracts a substring out of EXPR and returns it.

	my $intermed_annotation = $annot =~ s/\Q[/[Taxon hit id: /r;	# Replace '[' by '[Taxon hit id:'
	
	# RNA-binding protein fusilli isoform X1 [Taxon hit id: Zootermopsis nevadensis]
	
	my $final_annotation = $intermed_annotation =~ s/\Q]/, Acc: $all_acc]/r;	# Replace ']' by ', Acc: ...]'
	
	# RNA-binding protein fusilli isoform X1 [Hit tax id: Zootermopsis nevadensis, Acc: gi|1227956332|ref|XP_021916099.1|]
	# Include functional description, taxon hit and ncbi accessions.
	
	# Get the last accession (searcheable on uniprot)
	my @all_acc_list = split /\|/, $all_acc;
	my $acc = $all_acc_list[3];
	
	# Return description and accession
	return ($final_annotation, $acc);
}

#################### MAIN

# Define arguments & their default value
my ($species, $bgee_connector) = ('', '');
my ($obsGO) = ('');
my ($debug) = (0);
my %opts = ('species=s'    => \$species,            # speciesCommonName from TSV for or Bgee db
            'bgee=s'       => \$bgee_connector,     # Bgee connector string
            'obsGO=s'      => \$obsGO,              # go.obsolete file
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
\t-debug           Debug mode, do not insert/update in database
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

# Factor 1: Need to map to another genomeSpeciesId?

my ($prefix, $species_tmp, $speciesBgee, $newSpecies, $scientific_name, $ensSource) = Utils_insert_genes::map_species($dbh, $species);

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

print "Parsing the annotation file: " . $ann_file . " ...\n";

# Go through the annotation file

open (FILE, $ann_file);

while (<FILE>) {
	chomp;
	
	
	my @field = split("\t+");
	my $type = $field[2];	# gene, mRNA, exon, CDS, five_prime_UTR, three_prime_UTR
	
	# Look only for gene
	if ($type eq "gene") {
		
		my @gene_description = split /;/, $field[8];

		my @id = split /=/, $gene_description[0];
		my $gene_id = $id[1];
		
		# Parsing ...
		
		my $target;
		my $idx;
		my $gene_description_ref  = \@gene_description;
			
		# ... the functional annotation 
	
		$target = "topblasthit";
		$idx = undef;
		
		# Get the index of "xxx_topblasthit" in the array @gene_description
		
		foreach ( 0 .. $#gene_description ) {
		    if ( index ( $gene_description_ref->[ $_ ], $target ) >= 0 ) {
		        $idx = $_;
		        last;
		    }
		}
		
		if (defined $idx){
			
			my @result = func_formatting($gene_description[$idx]);
			$gene_hash{$gene_id}{Function} = $result[0];
			$gene_hash{$gene_id}{Accession} = $result[1];
			
			# ... the GO annotation 
	
			$target = "ontology_term";
			$idx = undef;
			
			# Get the index of "xxx_ontology_term" in the array @gene_description
			
			foreach ( 0 .. $#gene_description ) {
			    if ( index ( $gene_description_ref->[ $_ ], $target ) >= 0 ) {
			        $idx = $_;
			        last;
			    }
			}
			
			if (defined $idx){
				
				#print Dumper \%obs_go;
				my @go = go_formatting($gene_description[$idx], \%altid_go, \%obs_go);
				$gene_hash{$gene_id}{GO} = \@go;	# Create a reference to the array by using backslash
				
			}
			
		} 
		
		else {	# Gene has no function 
			$gene_hash{$gene_id} = "undef";
		}		
	 
	} 

}
close (FILE);



my $description;
my $accession;
my $biotype = "protein_coding";

#=for comment

foreach my $id2insert (sort {lc $a cmp lc $b} keys %gene_hash) {	# Sort to always get the same order

	my $bgeeGeneId;
		       
	if (ref($gene_hash{$id2insert}) eq 'HASH') {	# Gene has a function

    	$description = $gene_hash{$id2insert}{Function};
    	$accession = $gene_hash{$id2insert}{Accession};
    	
#=for comment
		## Insert gene info

		if ( ! $debug ){
			$geneDB->execute($id2insert, "", $description, $biotype, $speciesBgee)  or die $geneDB->errstr;
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

