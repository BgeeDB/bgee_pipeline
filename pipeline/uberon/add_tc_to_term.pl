#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Wollbrett, created Mar. 2026
# Add taxon constraints to an existing term.
# The script takes as input a TSV file with one term per line,
# in the format:
# OBO term ID\tspecies ID for which TC should be added
# e.g.:
# anatEntityId    speciesId
# UBERON:0001000	9103
# UBERON:0001301	9103,9031
#
# The script will then
# - retrieve the original taxon constraints of the term as they are stored in the database
# - retrieve parent terms and filter them to only include those that have at least the same TC as the term to update
# - retrieve parent relations and filter them to only include those that have at least the same TC as the term to update
# - add the new TC to the term
# - add the new TC to the previously filtered parent terms
# - add the new TC to the previously filtered parent relations
# 
# example of command to run the script:
# perl add_tc_to_term.pl -bgee=user=usr__pass=pwd__host=host__port=port__name=db_name -file_path=../../generated_files/uberon/add_tc.tsv > ../../generated_files/uberon/log_add_tc_to_term.txt
#############################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
# Define arguments & their default value
my ($bgee_connector) = ('');
my $file_path        = '';
my %opts = ('file_path=s'=> \$file_path,
            'bgee=s'     => \$bgee_connector   # Bgee connector string
           );
# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $file_path eq '' || $bgee_connector eq '' ) {
    print "\n\tInvalid or missing argument:
\te.g. $0 -file_path=/path/to/term/file -bgee=connection_string
\t-file_path        Path to file listing anatEntity IDs and species IDs for which TC should be added
\t-bgee             Bgee connector string
\n";
    exit 1;
}

# Read input file
my $fh;
my %terms_to_update = ();
open($fh, '<', $file_path) || die "failed to read input file: $!";
while (my $line = <$fh>) {
    chomp $line;
    # skip the title line if it is present
    next if $line =~ /^anatEntityId/;
    ## skip comments and blank lines and optional repeat of title line
    next if $line =~ /^\#/ || $line =~ /^\s*$/ || $line =~ /^\+/;
    #split each line into array
    my @line = split(/\s+/, $line);
    # the species ID column could contain multiple species IDs separated by commas, so we split it into an array
    my @species_ids = split(/,/, $line[1]);
    # we store the species IDs in the hash with the term ID as key
    $terms_to_update{$line[0]} = \@species_ids;
}
close($fh);

print " connector: $bgee_connector\n";
# Connect to the database
my $dbh = Utils::connect_bgee_db($bgee_connector);

my $IS_A_PART_OF = 'is_a part_of';

# retrieve number of species in Bgee to write NULL as speciesId when all speciesId are present as TC of one term or relation
my $sql_count_species = "SELECT COUNT(speciesId) FROM species";
my $sth_count_species = $dbh->prepare($sql_count_species);
$sth_count_species->execute() or die "Failed to execute query: " . $sth_count_species->errstr();
my ($nb_species) = $sth_count_species->fetchrow_array();
print "Number of species in Bgee: $nb_species\n";
$sth_count_species->finish();

# for each term to update
my $sql_current_TC = "SELECT speciesId FROM anatEntityTaxonConstraint WHERE anatEntityId = ?";
my $sth_current_TC = $dbh->prepare($sql_current_TC);

my $sql_parent_rel_with_TC = "select anatEntityRelationId, GROUP_CONCAT(speciesId) from anatEntityRelationTaxonConstraint where anatEntityRelationId IN (select anatEntityRelationId from anatEntityRelation where anatEntitySourceId IN (select anatEntityTargetId from anatEntityRelation where anatEntitySourceId = ? and relationType = ?) and relationType = ?) group by anatEntityRelationId;";
my $sth_parent_rel_with_TC = $dbh->prepare($sql_parent_rel_with_TC);
my $sql_parent_terms_with_TC = "select anatEntityId, GROUP_CONCAT(speciesId order by speciesId) from anatEntityTaxonConstraint where anatEntityId IN (select anatEntityTargetId from anatEntityRelation where anatEntitySourceId = ? and relationType = ?) group by anatEntityId;";
my $sth_parent_terms_with_TC = $dbh->prepare($sql_parent_terms_with_TC);

# sql queres to add taxon constraints
my $sql_add_anatEntityTaxonConstraint = "INSERT INTO anatEntityTaxonConstraint (anatEntityId, speciesId) VALUES (?, ?)";
my $sth_add_anatEntityTaxonConstraint = $dbh->prepare($sql_add_anatEntityTaxonConstraint);
my $sql_add_anatEntityRelationTaxonConstraint = "INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES (?, ?)";
my $sth_add_anatEntityRelationTaxonConstraint = $dbh->prepare($sql_add_anatEntityRelationTaxonConstraint);
# if after adding the new TC to the term to update, all species IDs are present as TC of the term, we write NULL
# as speciesId in the database and we don't add the new TC to parent terms and relations. It means we first need
# to delete existing TC of the term or the relation to update. The queries below will allow this deletion.
my $sql_delete_anatEntityTaxonConstraint = "DELETE FROM anatEntityTaxonConstraint WHERE anatEntityId = ?";
my $sth_delete_anatEntityTaxonConstraint = $dbh->prepare($sql_delete_anatEntityTaxonConstraint);
my $sql_delete_anatEntityRelationTaxonConstraint = "DELETE FROM anatEntityRelationTaxonConstraint WHERE anatEntityRelationId = ?";
my $sth_delete_anatEntityRelationTaxonConstraint = $dbh->prepare($sql_delete_anatEntityRelationTaxonConstraint);

foreach my $term_id (keys %terms_to_update) {
    my @species_ids = @{$terms_to_update{$term_id}};
    print "\nUpdating term $term_id with new TC for species IDs: " . join(", ", @species_ids) . "\n";
    # first we retrieve the TC of the term to update
    $sth_current_TC->execute($term_id) or die "Failed to execute query: " . $sth_current_TC->errstr();
    my @current_TC = ();
    while (my ($species_id) = $sth_current_TC->fetchrow_array()) {
        push @current_TC, $species_id;
    }
    print "Current TC of term $term_id: " . join(", ", @current_TC) . "\n";

    # then we retrieve parent relations with their TC and filter them to only keep those that have at least the same TC as the term to update
    $sth_parent_rel_with_TC->execute($term_id, $IS_A_PART_OF, $IS_A_PART_OF) or die "Failed to execute query: " . $sth_parent_rel_with_TC->errstr();
    my %parent_rel_with_TC = ();
    while (my ($rel_id, $species_ids_str) = $sth_parent_rel_with_TC->fetchrow_array()) {
        next if !defined($species_ids_str); # if the relation has no TC, we skip it
        #print "Parent relation $rel_id has TC for species IDs: $species_ids_str\n";
        my @rel_species_ids = split(/,/, $species_ids_str);
        # we check if the relation has at least the same TC as the term to update, i.e. if all species IDs of the term to update are present in the relation TC
        my $has_at_least_same_TC = 1;
        foreach my $species_id (@rel_species_ids) {
            if (!grep { $_ == $species_id } @current_TC) {
                $has_at_least_same_TC = 0;
                last;
            }
        }
        if ($has_at_least_same_TC) {
            $parent_rel_with_TC{$rel_id} = \@rel_species_ids;
        }
    }
    print "Parent relations with at least the same TC as term $term_id: " . join(", ", keys %parent_rel_with_TC) . "\n";

    # then we retrieve parent terms with their TC and filter them to only keep those that have at least the same TC as the term to update
    $sth_parent_terms_with_TC->execute($term_id, $IS_A_PART_OF) or die "Failed to execute query: " . $sth_parent_terms_with_TC->errstr();
    my %parent_terms_with_TC = ();
    while (my ($parent_term_id, $species_ids_str) = $sth_parent_terms_with_TC->fetchrow_array()) {
        next if !defined($species_ids_str); # if the parent term has no TC, we skip it
        my @parent_term_species_ids = split(/,/, $species_ids_str);
        # we check if the parent term has at least the same TC as the term to update, i.e. if all species IDs of the term to update are present in the parent term TC
        my $has_at_least_same_TC = 1;
        # check that @parent_term_species_ids contains all @current_TC species IDs

        foreach my $species_id (@parent_term_species_ids) {
            if (!grep { $_ == $species_id } @current_TC) {
                $has_at_least_same_TC = 0;
                last;
            }
        }
        if ($has_at_least_same_TC) {
            $parent_terms_with_TC{$parent_term_id} = \@parent_term_species_ids;
        }
    }
    print "Parent terms with at least the same TC as term $term_id: " . join(", ", keys %parent_terms_with_TC) . "\n";

    # start by updating the term to update and its parent terms.
    #XXX the term to update is part of %parent_terms_with_TC as each term has a reflexive relation to itself
    foreach my $parent_term_id (keys %parent_terms_with_TC) {
        my @parent_term_species_ids = @{$parent_terms_with_TC{$parent_term_id}};
        # if uniq speciesIds in @parent_term_species_ids and @species_ids are the same as the number of species in Bgee, we write NULL as speciesId in the database and we don't add the new TC to parent terms and relations
        my %uniq_species_ids = map { $_ => 1 } (@parent_term_species_ids, @species_ids);
        if (scalar(keys %uniq_species_ids) == $nb_species) {
            print "All species IDs are present as TC of parent term $parent_term_id after adding new TC, we will write NULL as speciesId in the database and we won't add the new TC to parent terms and relations\n";
            # we delete existing TC of the parent term to update
            $sth_delete_anatEntityTaxonConstraint->execute($parent_term_id) or die "Failed to execute query: " . $sth_delete_anatEntityTaxonConstraint->errstr();
            print "Deleted existing TC of parent term $parent_term_id\n";
            $sth_add_anatEntityTaxonConstraint->execute($parent_term_id, undef) or die "Failed to execute query: " . $sth_add_anatEntityTaxonConstraint->errstr();
            print "Added TC with NULL speciesId to parent term $parent_term_id\n";
        } else {
            # we add the new TC to the parent term
            foreach my $species_id (@species_ids) {
                # next if the TC is already present in the parent term TC

                if (!grep { $_ == $species_id } @parent_term_species_ids) {
                    $sth_add_anatEntityTaxonConstraint->execute($parent_term_id, $species_id) or die "Failed to execute query: " . $sth_add_anatEntityTaxonConstraint->errstr();
                    print "Added TC for species ID $species_id to parent term $parent_term_id\n";
                }
            }
        }
    }

    # then we update parent relations
    foreach my $rel_id (keys %parent_rel_with_TC) {
        my @rel_species_ids = @{$parent_rel_with_TC{$rel_id}};
        # if uniq speciesIds in @rel_species_ids and @species_ids are the same as the number of species in Bgee, we write NULL as speciesId in the database and we don't add the new TC to parent relations
        my %uniq_species_ids = map { $_ => 1 } (@rel_species_ids, @species_ids);
        if (scalar(keys %uniq_species_ids) == $nb_species) {
            print "All species IDs are present as TC of parent relation $rel_id after adding new TC, we will write NULL as speciesId in the database and we won't add the new TC to parent relations\n";
            # we delete existing TC of the parent relation to update
            $sth_delete_anatEntityRelationTaxonConstraint->execute($rel_id);
            print "Deleted existing TC of parent relation $rel_id\n";
            $sth_add_anatEntityRelationTaxonConstraint->execute($rel_id, undef);
            print "Added TC with NULL speciesId to parent relation $rel_id\n";
        } else {
            # we add the new TC to the parent relation
            foreach my $species_id (@species_ids) {
                if (!grep { $_ == $species_id } @rel_species_ids) {
                    $sth_add_anatEntityRelationTaxonConstraint->execute($rel_id, $species_id);
                    print "Added TC for species ID $species_id to parent relation $rel_id\n";
                }
            }
        }
    }
}
# stop sth
$sth_current_TC->finish();
$sth_parent_rel_with_TC->finish();
$sth_parent_terms_with_TC->finish();
$sth_add_anatEntityTaxonConstraint->finish();
$sth_add_anatEntityRelationTaxonConstraint->finish();
$sth_delete_anatEntityTaxonConstraint->finish();
$sth_delete_anatEntityRelationTaxonConstraint->finish();
# close database connection
$dbh->disconnect();