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
            'paralogs_dir_path=s'       => \$paralogs_dir_path,         # path to directory containing all paralogy files
            'orthologs_dir_path=s'      => \$orthologs_dir_path,        # path to directory containing all orthology files
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


# LOGIC OF INSERTION:
# * load all Bgee species
# * for each species
#       * load all genes
#       * detect files containing speciesId in their name
#       * for each file
#           * load genes of second species
#           * insert ortholog/paralog info with taxonId (each homology relation is symmetric)
#           * tag file as already visited

# For now we do not insert paralogy relations at taxon level not present in Bgee.

# load species in Bgee database
my $load_species_query = "select speciesId from species";
# load all taxon that are descendant of the oldest Bgee LCA taxon (Bilateria) to map OMA taxon name to taxonId
# only LCA bgee taxon are required to insert orthologs but... paralogy can be at any taxon level.
my $load_taxon_query = "select taxonId, taxonScientificName FROM taxon where taxonLeftBound >= (select MIN(taxonLeftBound) from taxon where bgeeSpeciesLCA = 1) AND taxonRightBound <= (select MAX(taxonRightBound) from taxon where bgeeSpeciesLCA = 1);";

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

#insert_homology_type($orthologs_dir_path, \@bgee_species_array,  \%taxonName_to_taxonId, $dbh, "ortholog");
insert_homology_type($paralogs_dir_path, \@bgee_species_array,  \%taxonName_to_taxonId, $dbh, "paralog");

############# FUNCTIONS #############

sub load_genes_map {
    my ($species_id) = @_;
    # load all genes of one species
    my $load_genes_query = "select bgeeGeneId, geneId from gene where speciesId = ?";
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

# insert orthologs OR paralogs in Bgee DB depending on the value of the homology type provided as 
# argument. Homology type can be : "orthologs" or "paralogs"
sub insert_homology_type {

    my $homology_dir_path = $_[0];
    my %taxonName_to_taxonId = %{$_[2]};
    my $dbh = $_[3];
    my $homology_type = $_[4];
    my @species_array = @{$_[1]};
    
    #init variables mandatory for insertion
    my ($homo_file_prefix, $insert_homo_query);
    if ($homology_type eq "ortholog") {
        $homo_file_prefix = "orthologs";
        $insert_homo_query = "insert into geneOrthologs (bgeeGeneId, targetGeneId, taxonId) VALUES(?, ?, ?)";
    } elsif ($homology_type eq "paralog") {
        $homo_file_prefix = "paralogs";
        $insert_homo_query = "insert into geneParalogs (bgeeGeneId, targetGeneId, taxonId) VALUES(?, ?, ?)";
    } else {
        die "homology type should be \"ortholog\" or \"paralog\" but was $homology_type.";
    }

    # retrieve all ortholog files

    opendir my $homo_dir, $homology_dir_path or die "Cannot open directory: $!";
    my @homo_files = readdir $homo_dir;

    print "start insertion of $homology_type\n";
    # array used to keep note of already parsed files
    my %already_parsed_files;

    # for each bgee species
    foreach my $species_id (@species_array) {
        #load genes of this species
        my %ensemblId_to_bgeeId = load_genes_map($species_id);
        print "start insertion of species $species_id\n";

        #parse all files where current species in first position in file name
        my @first_species_files = grep /_$species_id-/, @homo_files;
        %already_parsed_files = insert_from_file_name(\@first_species_files, $insert_homo_query, $homo_file_prefix, 
            $homology_dir_path, \%taxonName_to_taxonId, $dbh, \%ensemblId_to_bgeeId, \%already_parsed_files, 1);

        my @second_species_files = grep /-$species_id\./, @homo_files;
        %already_parsed_files = insert_from_file_name(\@second_species_files, $insert_homo_query, $homo_file_prefix, 
            $homology_dir_path, \%taxonName_to_taxonId, $dbh, \%ensemblId_to_bgeeId, \%already_parsed_files, 2);
    }
}

sub insert_from_file_name {

    #retrieve arguments
    my @species_files = @{$_[0]};
    my $query = $_[1];
    my $file_prefix = $_[2];
    my $homology_dir_path = $_[3];
    my %taxonName_to_taxonId = %{$_[4]};
    my $dbh = $_[5];
    my %ensemblId_to_bgeeId = %{$_[6]};
    my %already_parsed_files = %{$_[7]};
    my $species_position = $_[8];


    foreach my $species_file (@species_files) {
        # disable auto commit in order to insert all tuples of one species at the same time
        $dbh->{AutoCommit} = 0;
        $sth = $dbh->prepare_cached($query);
        next if (exists($already_parsed_files{$species_file}));
        my $homo_species_id;
        $species_file =~ /[$file_prefix]_bgee_([0-9]*)-([0-9]*)\.csv/;
        my %ensemblId_to_bgeeId_homologous_species;
        # if paralogs of same species do not need to reload genes
        if ($1 == $2) {
            %ensemblId_to_bgeeId_homologous_species = %ensemblId_to_bgeeId;
        } else {
            if ($species_position == 1) {
                %ensemblId_to_bgeeId_homologous_species = load_genes_map(hack_gorilla_tax_id($2));
            } elsif ($species_position == 2) {
                %ensemblId_to_bgeeId_homologous_species = load_genes_map(hack_gorilla_tax_id($1));
            } else {
                die "species_position should be 1 or 2 but was $species_position";
            }
        }

        my @values_to_insert;
        open my $homo_file_handler, "$homology_dir_path/$species_file" or die "failed to read input file: $!";
        while (my $line = <$homo_file_handler>) {
            
            chomp $line;
            #skip header
            next if $line =~ /^gene1,gene2,tax_level$/;
            #split each line into array
            my @line = split(/,/, $line);
            my ($current_species_bgee_id, $homolog_species_bgee_id) = ('', '');
            if($species_position == 1) {
                $current_species_bgee_id = $ensemblId_to_bgeeId{$line[0]};
                $homolog_species_bgee_id = $ensemblId_to_bgeeId_homologous_species{$line[1]};
            } else { # no need to check values again, we already tested before that values where only 1 or 2
                $current_species_bgee_id = $ensemblId_to_bgeeId{$line[1]};
                $homolog_species_bgee_id = $ensemblId_to_bgeeId_homologous_species{$line[0]};
            }
            
            my $lca_taxon_id = $taxonName_to_taxonId{$line[2]};
            do {
                no warnings 'uninitialized'; # for the do block only
                    
                if ($current_species_bgee_id eq '' || $homolog_species_bgee_id eq '' || $lca_taxon_id eq '') {
                    warn "Can not map OMA data to Bgee : [$line[0], $line[1], $line[2]] became : [$current_species_bgee_id, $homolog_species_bgee_id, $lca_taxon_id]. From file : $species_file";
                } else {
                    $sth->execute($current_species_bgee_id, $homolog_species_bgee_id, $lca_taxon_id);
                }
            }
            

        }

        $sth->finish;
        $dbh->{AutoCommit} = 1;
        $already_parsed_files{$species_file} = 1;
        close $homo_file_handler;
    }

    return %already_parsed_files;

}

# in bgee we do not have the proper taxonId for Gorilla
sub hack_gorilla_tax_id {
    my $tax_id = $_[0];
    if($tax_id == 9595) {
        $tax_id = 9593;
    }
    return $tax_id;
}