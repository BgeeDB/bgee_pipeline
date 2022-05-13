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


# load species in Bgee database
my $load_species_query = "select speciesId, taxonId from species";
# load all Bgee taxon that
my $load_taxon_query = "select taxonId, taxonScientificName FROM taxon;";

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

# load bgee species
my %bgee_species_hash;
my $sth = $dbh->prepare($load_species_query);
$sth->execute()  or die $sth->errstr;
my ($speciesId, $speciesTaxonId) = ('','');
$sth->bind_columns(\$speciesId, \$speciesTaxonId);
while($sth->fetch()) {
    $bgee_species_hash{$speciesId} = $speciesTaxonId;
}

# load bgee taxon
my %taxonName_to_taxonId;
my %taxonId_to_taxonName;
$sth = $dbh->prepare($load_taxon_query);
$sth->execute()  or die $sth->errstr;
my ($taxonId, $taxonName) = ('','');
$sth->bind_columns(\$taxonId, \$taxonName);
while($sth->fetch()) {
   $taxonId_to_taxonName{$taxonId} = $taxonName;
   $taxonName_to_taxonId{$taxonName} = $taxonId;
}

insert_homology_type($orthologs_dir_path, \%bgee_species_hash, $dbh, "ortholog", \%taxonName_to_taxonId);
insert_homology_type($paralogs_dir_path, \%bgee_species_hash, $dbh, "paralog", \%taxonName_to_taxonId);

$dbh->disconnect or warn "Disconnection failed: $DBI::errstr\n";
exit;

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
    my $dbh = $_[2];
    my $homology_type = $_[3];
    my %bgee_species_hash = %{$_[1]};
    my %taxonName_to_taxonId = %{$_[4]};
    
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
    foreach my $species_id (keys %bgee_species_hash) {
        #load genes of this species
        my %ensemblId_to_bgeeId = load_genes_map($species_id);
        my $modified_species_id = hack_tax_id_bgee_to_oma($species_id);
        print "start insertion of species $modified_species_id\n";

        #parse all files where current species in first position in file name
        my @first_species_files = grep /_${modified_species_id}-/, @homo_files;
        %already_parsed_files = insert_from_file_name(\@first_species_files, $insert_homo_query, $homo_file_prefix, 
            $homology_dir_path, $dbh, \%ensemblId_to_bgeeId, \%already_parsed_files, \%bgee_species_hash, 1, \%taxonName_to_taxonId);

        my @second_species_files = grep /-${modified_species_id}\.csv/, @homo_files;
        %already_parsed_files = insert_from_file_name(\@second_species_files, $insert_homo_query, $homo_file_prefix, 
            $homology_dir_path, $dbh, \%ensemblId_to_bgeeId, \%already_parsed_files, \%bgee_species_hash, 2, \%taxonName_to_taxonId);
        if (!@first_species_files && !@second_species_files) {
            warn "No homology files for species $species_id.";
        }
    }
}

sub insert_from_file_name {

    #retrieve arguments
    my @species_files = @{$_[0]};
    my $query = $_[1];
    my $file_prefix = $_[2];
    my $homology_dir_path = $_[3];
    my $dbh = $_[4];
    my %ensemblId_to_bgeeId = %{$_[5]};
    my %already_parsed_files = %{$_[6]};
    my %bgee_species_hash = %{$_[7]};
    my $species_position = $_[8];
    my %taxonName_to_taxonId = %{$_[9]};

    # some orthologs are present twice in OMA files. We do not want to insert them twice
    # in order to save disk space. the symmetry is managed in the Java API. In order not
    # to insert symmetric homology relations gene names and taxonId are stored in a hash
    # allowing to check if geneA_geneB or geneB_geneA has already been inserted
    my %already_inserted;


    foreach my $species_file (@species_files) {


        # init variables used to log homology inserted or not inserted
        my $inserted_homology = 0;
        my $not_mapped_homology = 0;
        my $duplicated_homology = 0;

        # disable auto commit in order to insert all tuples of one species at the same time
        $dbh->{AutoCommit} = 0;
        $sth = $dbh->prepare_cached($query);
        next if (exists($already_parsed_files{$species_file}));

        $species_file =~ /[$file_prefix]_([0-9]*)-([0-9]*)\.csv/;
        my $species1 = $1;
        my $species2 = $2;
        # verify the species for which homology has to be insert is present in the Bgee database
        if($species_position == 1 && !exists($bgee_species_hash{hack_tax_id_oma_to_bgee($species2)})) {
            warn "species 2 does not exist. Do not consider file [$species_file]";
            next;
        }
        if($species_position == 2 && !exists($bgee_species_hash{hack_tax_id_oma_to_bgee($species1)})) {
            warn "species 1 does not exist. Do not consider file [$species_file]";
            next;
        }
        # if paralogs of same species do not need to reload genes
        my %ensemblId_to_bgeeId_homologous_species;
        if ($species1 == $species2) {
            %ensemblId_to_bgeeId_homologous_species = %ensemblId_to_bgeeId;
        } else {
            if ($species_position == 1) {
                %ensemblId_to_bgeeId_homologous_species = load_genes_map(hack_tax_id_oma_to_bgee($species2));
            } elsif ($species_position == 2) {
                %ensemblId_to_bgeeId_homologous_species = load_genes_map(hack_tax_id_oma_to_bgee($species1));
            } else {
                die "species_position should be 1 or 2 but was $species_position";
            }
        }
        open my $homo_file_handler, "$homology_dir_path/$species_file" or die "failed to read input file ".
            "$homology_dir_path/$species_file: $!";
        while (my $line = <$homo_file_handler>) {
            
            chomp $line;
            #skip header
            next if $line =~ /^gene1,gene2,tax_level/;
            next if $line =~ /^gene2,gene1,tax_level/;
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
            my $lca_taxon_name = $line[2];
            my $lca_taxon_id = $line[3];
            
            if (!$current_species_bgee_id || !$homolog_species_bgee_id || !$lca_taxon_id) {
                $not_mapped_homology++;
                #warn "Can not map OMA data to Bgee : [$line[0], $line[1], $line[2]] became : ".
                #"[$current_species_bgee_id, $homolog_species_bgee_id, $lca_taxon_id]. From file : $species_file";
            } else {
                # for Cetartiodactyla taxon name is correct but the taxonId is equal to -3
                # in OMA database. This taxon is not present in the Bgee database. We map it 
                # to Artiodactyla (https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=91561).
                if($lca_taxon_id == -3 && $lca_taxon_name eq "Cetartiodactyla") {
                    $lca_taxon_name = "Artiodactyla";
                    $lca_taxon_id = 91561;
                }
                #before query execution we check that this couple of genes has not already been inserted
                if (!exists($already_inserted{"${current_species_bgee_id}_${homolog_species_bgee_id}"})
                    && !exists($already_inserted{"${homolog_species_bgee_id}_${current_species_bgee_id}"})) {

                    # In OMA the taxonId of a duplication is defined as the closest descendant taxon. Then it happens 
                    # that the taxonId of an in-species duplication is the species itself. In Bgee species and taxon 
                    # are 2 different concepts stored in 2 different tables .The RDB table that stores paralogy
                    # has a foreign key between its column taxonId and the column taxon.taxonId.
                    # In order to link such paralogy relations to a taxonId we map it to the direct parent of the
                    # speciesId (e.g Chlorocebus sabaeus => Chlorocebus). This parent information comes from the 
                    # column species.taxonId of the RDB.
                    if($species1 == $species2 && $lca_taxon_id == $species1) {
                        $lca_taxon_id = $bgee_species_hash{hack_tax_id_oma_to_bgee($lca_taxon_id)};
                    }
                    $sth->execute($current_species_bgee_id, $homolog_species_bgee_id, $lca_taxon_id) 
                        or warn "can not insert $file_prefix [$line[0], $line[1], ${lca_taxon_id}] ".
                        "from file [$species_file]";
                    $already_inserted{"${current_species_bgee_id}_${homolog_species_bgee_id}"} = ();
                    $inserted_homology++;
                } else {
                    $duplicated_homology++;
                }
            }
        }

        $sth->finish;
        $dbh->{AutoCommit} = 1;
        $already_parsed_files{$species_file} = 1;
        close $homo_file_handler;
        print "finished to parse ${species_file}. $inserted_homology homologies inserted. $duplicated_homology ".
        "homologies not inserted because duplicated in the file. $not_mapped_homology homologies not inserted because ".
        "no mapping between OMA gene ID and Bgee gene ID\n";
    }

    return %already_parsed_files;

}

# difference of taxonId used for some species between OMA and Bgee
sub hack_tax_id_oma_to_bgee {
    my $tax_id = $_[0];
    # Gorilla gorilla gorilla => Gorilla gorilla
    if($tax_id == 9595) {
        $tax_id = 9593;
    }
    # Drosophila pseudoobscura pseudoobscura => Drosophila pseudoobscura
    if($tax_id == 46245) {
        $tax_id = 7237;
    }
    return $tax_id;
}

sub hack_tax_id_bgee_to_oma {
    my $tax_id = $_[0];
    # Gorilla gorilla => Gorilla gorilla gorilla
    if($tax_id == 9593) {
        $tax_id = 9595;
    }
    # Drosophila pseudoobscura => Drosophila pseudoobscura pseudoobscura
    if($tax_id == 7237) {
        $tax_id = 46245;
    }
    return $tax_id;
}