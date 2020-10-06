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
#           * insert ortholog/paralog info with taxonId (each homology relation is symmetric, always store twice for the moment)
#           * tag file as already visited


# load species in Bgee database
my $load_species_query = "select speciesId from species";
# load all taxon that are descendant of the oldest Bgee LCA taxon (Bilateria) to map OMA taxon name to taxonId
# only LCA bgee taxon are required to insert orthologs but... paralogy can be at any taxon level.
my $load_taxon_query = "select taxonId, taxonScientificName, taxonLeftBound FROM taxon where taxonLeftBound >= (select MIN(taxonLeftBound) from taxon where bgeeSpeciesLCA = 1) AND taxonRightBound <= (select MAX(taxonRightBound) from taxon where bgeeSpeciesLCA = 1);";

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
my ($taxonId, $taxonName, $taxonLeftBound);
$sth->bind_columns(\$taxonId, \$taxonName, \$taxonLeftBound);
while($sth->fetch()) {
    # if 2 taxon have the same name, keep the one closest to the root
    # TODO remove once OMA provide taxonIds
    if(exists($taxonName_to_taxonId{$taxonName}) && $taxonLeftBound > $taxonName_to_taxonId{$taxonName}{"taxonLeftBound"}) {
        next;
    }
    $taxonName_to_taxonId{$taxonName}{"taxonId"} = $taxonId;
    $taxonName_to_taxonId{$taxonName}{"taxonLeftBound"} = $taxonLeftBound;
    
}
# TODO should check if homologGroup table is empty. If not die and ask to delete table first
my $total_number_group = 0;
$total_number_group = insert_homology_type($orthologs_dir_path, \@bgee_species_array,  \%taxonName_to_taxonId, $dbh, $total_number_group, "ortholog");
$total_number_group = insert_homology_type($paralogs_dir_path, \@bgee_species_array,  \%taxonName_to_taxonId, $dbh, $total_number_group, "paralog");

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
    my $homology_type = $_[5];
    my $total_group_number = $_[4]
    my @species_array = @{$_[1]};
    
    #init variables mandatory for insertion
    my ($homo_file_prefix, $insert_homo_query);
    if ($homology_type eq "ortholog") {
        $homo_file_prefix = "orthologs";
        $insert_homo_query = "insert into homologGroup (homologGroupId, orthologGroup) VALUES(?, 0)";
    } elsif ($homology_type eq "paralog") {
        $homo_file_prefix = "paralogs";
        $insert_homo_query = "insert into homologGroup (homologGroupId, orthologGroup) VALUES(?, 1)";
    } else {
        die "homology type should be \"ortholog\" or \"paralog\" but was $homology_type.";
    }

    # retrieve all ortholog files

    opendir my $homo_dir, $homology_dir_path or die "Cannot open directory: $!";
    my @homo_files = readdir $homo_dir;

    print "start insertion of $homology_type\n";
    # hash used to keep note of already parsed files
    my %already_parsed_files;
    # hash used to know in each homolog group is a gene
    my %gene_to_group;
    # for each bgee species
    foreach my $species_id (@species_array) {
        #load genes of this species
        my %ensemblId_to_bgeeId = load_genes_map($species_id);
        print "start insertion of species $species_id\n";

        #parse all files where current species in first position in file name
        my @first_species_files = grep /_$species_id-/, @homo_files;
        (%already_parsed_files, %gene_to_group, $total_group_number) = insert_from_file_name(\@first_species_files, $insert_homo_query, $homo_file_prefix, 
            $homology_dir_path, \%taxonName_to_taxonId, $dbh, \%ensemblId_to_bgeeId, \%already_parsed_files, \%gene_to_group, $total_group_number, 1);

        my @second_species_files = grep /-$species_id\./, @homo_files;
        (%already_parsed_files, %gene_to_group, $total_group_number) = insert_from_file_name(\@second_species_files, $insert_homo_query, $homo_file_prefix, 
            $homology_dir_path, \%taxonName_to_taxonId, $dbh, \%ensemblId_to_bgeeId, \%already_parsed_files, \%gene_to_group, $total_group_number, 2);

        # Remove bgee gene Ids of current species from %gene_to_group
        foreach my $bgee_gene_id (values %ensemblId_to_bgeeId) {
            delete $gene_to_group{$bgee_gene_id};
        }
    }
    return $total_group_number;
}

sub insert_from_file_name {

    #retrieve arguments
    my @species_files = @{$_[0]};
    my $group_query = $_[1];
    my $file_prefix = $_[2];
    my $homology_dir_path = $_[3];
    my %taxonName_to_taxonId = %{$_[4]};
    my $dbh = $_[5];
    my %ensemblId_to_bgeeId = %{$_[6]};
    my %already_parsed_files = %{$_[7]};
    my %gene_to_group = %{$_[8]};
    my $total_group_number = $_[9];
    my $species_position = $_[10];


    foreach my $species_file (@species_files) {
        # disable auto commit in order to insert all tuples of one species at the same time
        $dbh->{AutoCommit} = 0;
        my $sth_group = $dbh->prepare_cached($group_query);
        my $sth_genes = $dbh->prepare_cached("insert into homologGroupToGene (homologGroupId, bgeeGeneId, taxonId) VALUES(?, ?, ?)");
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

            # move to next line if taxon does not exist in bgee
            my $lca_taxon_id = '';
            if(exists($taxonName_to_taxonId{$line[2]})) {
                $lca_taxon_id = $taxonName_to_taxonId{$line[2]};
            } else {
                warn "can not map OMA taxon name $line[2] to a bgee taxon Id. The line is skipped.";
                next;
            }

            # retrieve genes that have a corresponding bgee gene ID
            my @bgee_genes_to_add = ();
            if($species_position == 1) {
                if(exists($ensemblId_to_bgeeId{$line[0]})) {
                    push(@bgee_genes_to_add, $ensemblId_to_bgeeId{$line[0]});
                } else {
                    warn "can not map OMA gene ID $line[0] to a bgee gene Id. The gene is skipped.";
                }
                if (exists($ensemblId_to_bgeeId_homologous_species{$line[0]})) {
                    push(@bgee_genes_to_add, $ensemblId_to_bgeeId_homologous_species{$line[1]});
                } else {
                    warn "can not map OMA gene ID $line[0] to a bgee gene Id. The gene is skipped.";
                }
            } else { # no need to check values again, we already tested before that values where only 1 or 2
                if(exists($ensemblId_to_bgeeId{$line[1]})) {
                    push(@bgee_genes_to_add, $ensemblId_to_bgeeId{$line[1]});
                } else {
                    warn "can not map OMA gene ID $line[1] to a bgee gene Id. The gene is skipped.";
                }
                if (exists($ensemblId_to_bgeeId_homologous_species{$line[0]})) {
                    push(@bgee_genes_to_add, $ensemblId_to_bgeeId_homologous_species{$line[0]});
                } else {
                    warn "can not map OMA gene ID $line[0] to a bgee gene Id. The gene is skipped.";
                }
            }

            # move to next line if none of the gene Ids are present in bgee
            if (scalar @bgee_genes_to_add == 0) {
                warn "can not map OMA gene ids $line[0] and $line[1] to bgee gene Ids. The line is skipped.";
                next;
            }
            


            my $gene_group = -1;

            # retrieve corresponding homolog group if already exist
            foreach my $bgee_gene_to_add (@bgee_genes_to_add) {
                if (exists($gene_to_group{$bgee_gene_to_add})) {
                    my $current_gene_group = $gene_to_group{$bgee_gene_to_add};
                    if ($gene_group != -1 && $gene_group != $current_gene_group) {
                        die "gene present in 2 different homolog gene groups!!!";
                    }
                    $gene_group = $current_gene_group;
                }
            }
            #create new homolog group if necessary
            if($gene_group == -1) {
                $total_group_number = $total_group_number + 1;
                $gene_group = $total_group_number;
                #create new group in the RDB
                $sth_group->execute($gene_group);
            }
            
            # add data to homologGroupToGene table
            foreach my $bgee_gene_to_add (@bgee_genes_to_add) {
                if (!exists($gene_to_group{$bgee_gene_to_add})) {
                    $sth_genes->execute($gene_group, $bgee_gene_to_add, $lca_taxon_id);
                }
            }
        }

        $sth_genes->finish;
        $sth_group->finish;
        $dbh->{AutoCommit} = 1;
        $already_parsed_files{$species_file} = 1;
        close $homo_file_handler;
    }

    return (%already_parsed_files, %gene_to_group, $total_group_number);

}

# in bgee we do not have the proper taxonId for Gorilla
sub hack_gorilla_tax_id {
    my $tax_id = $_[0];
    if($tax_id == 9595) {
        $tax_id = 9593;
    }
    return $tax_id;
}