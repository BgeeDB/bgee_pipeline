#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Text::CSV;
use Data::Dumper;

my ($bgee_connector, $paralogs_dir_path, $orthologs_dir_path) = ('', '', '');
my $debug = 0;
my %opts = ('bgee=s'                    => \$bgee_connector,            # Bgee connector string
            'paralogs_dir_path=s'       => \$paralogs_dir_path,              # path to directory containing all paralogy files
            'orthologs_dir_path=s'      => \$orthologs_dir_path,             # path to directory containing all orthology files
            'debug'                     => \$debug,                     # debug mode, do not insert/update in database
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( $bgee_connector eq '' || $paralogs_dir_path eq '' || $orthologs_dir_path eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. perl $0 -bgee=\$(BGEECMD) -paralogs_dir_path=path/to/fpara/dir -orthologs_dir_path=path/to/ortho/dir
\t-bgee                     Bgee connector string
\t-paralogs_dir_path        path to directory containing all paralogy files
\t-orthologs_dir_path       path to directory containing all orthology files
\t-debug                    Debug mode, do not insert/update in database
\n";
    exit 1;
}


# TODO: insert symetric homology relations only once and modify API
# LOGIC OF INSERTION:
# * load all Bgee species
# * for each species
#       * load all genes
#       * detect files containing speciesId in their name
#       * for each file
#           * load genes of second species
#           * insert ortholog/paralog info with taxonId (each homology relation is symmetric, always store twice for the moment)
#           * tag file as already visited

# load species in Bgee database
my $load_species_query = "select speciesId from species";
# load all taxon that are Least Common Ancestor in Bgee. Mandatory to filter LCA taxon and to map OMA taxon name to taxonId.
my $load_taxon_query = "select taxonId, taxonScientificName from taxon where bgeeSpeciesLCA = 1";
# load all genes of one species
my $load_genes_query = "select bgeeGeneId, geneId from gene where speciesId = ?";
# insert homology for one species (and its symetric relations)
my $insert_orthology = "insert into geneOrthologs (bgeeGeneId, targetGeneId, taxonId) VALUES(?, ?, ?)";
my $insert_paralogy = "insert into geneParalogs (bgeeGeneId, targetGeneId, taxonId) VALUES(?, ?, ?)";



# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

# load bgee species
my @bgee_species_array;
my $sth = $dbh->prepare($load_species_query);
$sth->execute()  or die $sth->errstr;
my $speciesId;
$sth->bind_columns(\$speciesId);
while($sth->fetch()) {
   push @bgee_species_array, $speciesId;
}

# load bgee taxon
my %taxonName_to_taxonId;
$sth = $dbh->prepare($load_taxon_query);
$sth->execute()  or die $sth->errstr;
my ($taxonId, $taxonName);
$sth->bind_columns(\$taxonId, \$taxonName);
while($sth->fetch()) {
   $taxonName_to_taxonId{$taxonName} = $taxonId;
}

# retrieve all ortholog files
opendir my $ortho_dir, $orthologs_dir_path or die "Cannot open directory: $!";
my @ortho_files = readdir $ortho_dir;

# array used to keep note of already parsed files
my %already_parsed_files;

# for each bgee species
foreach my $species_id (@bgee_species_array) {
    #load genes of this species
    my %ensemblId_to_bgeeId = load_genes_map($species_id);
    print "$species_id\n";

    #parse all files where current species in first position in file name
    my @first_species_files = grep /_$species_id-/, @ortho_files;
    foreach my $first_species_file (@first_species_files) {
        next if (exists($already_parsed_files{$first_species_file}));
        $first_species_file =~ /orthologs_bgee_\w*-(\w*)\.csv/;
        my $second_species_id = $1;
        my %ensemblId_to_bgeeId_homologous_species = load_genes_map($second_species_id);
        my @values_to_insert;
        open my $homo_file_handler, "$orthologs_dir_path/$first_species_file" or die "failed to read input file: $!";
        while (my $line = <$homo_file_handler>) {
            
            chomp $line;
            #skip header
            next if $line =~ /^gene1,gene2,tax_level$/;
            #split each line into array
            my @line = split(/,/, $line);
            #print "$line[0] $line[1] $line[2]\n";
            my $current_species_bgee_id = $ensemblId_to_bgeeId{$line[0]};
            my $homolog_species_bgee_id = $ensemblId_to_bgeeId_homologous_species{$line[1]};
            my $lca_taxon_id = $taxonName_to_taxonId{$line[2]};
            #print "$current_species_bgee_id - $homolog_species_bgee_id - $lca_taxon_id\n";
            if ($current_species_bgee_id ='' || $homolog_species_bgee_id = '' || $lca_taxon_id='') {
                warn "missing values : $current_species_bgee_id, $homolog_species_bgee_id, $lca_taxon_id";
            } else {
                push @values_to_insert, "($current_species_bgee_id, $homolog_species_bgee_id, $lca_taxon_id)";
                push @values_to_insert, "($homolog_species_bgee_id, $current_species_bgee_id, $lca_taxon_id)";
            }
            

        }
        $dbh->{AutoCommit} = 0;
        $sth = $dbh->prepare_cached( $insert_orthology );
        $sth->execute( @$_ ) for @values_to_insert;
        $sth->finish;
        $dbh->{AutoCommit} = 1;
        #print "number of rows in $first_species_file : $#values_to_insert\n";
        $already_parsed_files{$second_species_file} = 1;
        close $homo_file_handler;
    }

    #parse all files where current species in second position in file name
    my @second_species_files = grep /-$species_id\./, @ortho_files; 
    foreach my $second_species_file (@second_species_files) {
        next if (exists($already_parsed_files{$second_species_file}));
        $second_species_file =~ /orthologs_bgee_(\w*)-\w*\.csv/;
        my $first_species_id = $1;
        my %ensemblId_to_bgeeId_homologous_species = load_genes_map($first_species_id);
        my @values_to_insert;
        open my $homo_file_handler, "$orthologs_dir_path/$second_species_file" or die "failed to read input file: $!";
        while (my $line = <$homo_file_handler>) {
            
            chomp $line;
            #skip header
            next if $line =~ /^gene1,gene2,tax_level$/;
            #split each line into array
            my @line = split(/,/, $line);
            #print "$line[0] $line[1] $line[2]\n";
            my $current_species_bgee_id = $ensemblId_to_bgeeId{$line[1]};
            my $homolog_species_bgee_id = $ensemblId_to_bgeeId_homologous_species{$line[0]};
            my $lca_taxon_id = $taxonName_to_taxonId{$line[2]};
            #print "$current_species_bgee_id - $homolog_species_bgee_id - $lca_taxon_id\n";
            if ($current_species_bgee_id ='' || $homolog_species_bgee_id = '' || $lca_taxon_id='') {
                warn "missing values : $current_species_bgee_id, $homolog_species_bgee_id, $lca_taxon_id";
            } else {
                push @values_to_insert, "($current_species_bgee_id, $homolog_species_bgee_id, $lca_taxon_id)";
                push @values_to_insert, "($homolog_species_bgee_id, $current_species_bgee_id, $lca_taxon_id)";
            }
            

        }
        $dbh->{AutoCommit} = 0;
        $sth = $dbh->prepare_cached( $insert_orthology );
        $sth->execute( @$_ ) for @values_to_insert;
        $sth->finish;
        $dbh->{AutoCommit} = 1;
        #print "number of rows in $first_species_file : $#values_to_insert\n";
        $already_parsed_files{$second_species_file} = 1;
        close $homo_file_handler;
    }

}


# FUNCTIONS

sub load_genes_map {
    my ($species_id) = @_;
    my %ensemblId_to_bgeeId;
    my $sth2 = $dbh->prepare($load_genes_query);
    $sth2->execute($species_id)  or die $sth2->errstr;
    my ($bgeeGeneId, $geneId);
    $sth2->bind_columns(\$bgeeGeneId, \$geneId);
    while($sth2->fetch()) {
        $ensemblId_to_bgeeId{$geneId} = $bgeeGeneId;
    }
    return %ensemblId_to_bgeeId;
}