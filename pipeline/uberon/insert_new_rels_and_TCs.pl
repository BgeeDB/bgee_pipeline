#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Data::Dump qw(dump);
use Tie::IxHash;

# Frederic Bastian, created May 2024
# Insert relations for new anatomical terms added to the database,
# and infer their taxon constraints from their parents.
# The script takes as input a TSV file with one relation per line,
# in the format:
# OBO child ID\tOBO parent ID
# No header to the TSV file.

# Julien Wollbrett, modified Apr. 2026
# A term keeps the taxon constraints it already has in the database; otherwise it gets the intersection
# of the taxon constraints of its parents. The new terms are processed parents first, so that a new term
# whose parent is also new inherits from that parent.
# The closure of the relations links each term to itself (reflexive relation), to its parents (direct
# relations) and to its other ancestors (indirect relations). Bgee stores it in anatEntityRelation
# to get all the ancestors of a term at once, without following the relations from parent to parent.
# With -indirect_relations, the script also inserts the indirect relations of the closure: from each
# new term to the ancestors of its parents, and from the descendants of the new term already in the
# database to its new ancestors. A relation is valid in the species where every relation of a chain
# producing it is valid; when several chains produce it, it is valid in the union of their species.
# A term taxon constraint, relation or relation taxon constraint that already exists is kept and completed.
#############################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my ($bgee_connector) = ('');
my $file_path        = '';
my $indirect_relations = 0; # whether to insert the indirect relations of the closure
my %opts = ('file_path=s'        => \$file_path,
            'bgee=s'             => \$bgee_connector,   # Bgee connector string
            'indirect_relations' => \$indirect_relations
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $file_path eq '' || $bgee_connector eq '' ) {
    print "\n\tInvalid or missing argument:
\te.g. $0 -file_path=/path/to/relation/file -bgee=connection_string
\t-file_path            Path to relation file
\t-bgee                 Bgee connector string
\t-indirect_relations   Whether to insert the indirect relations of the closure (1) or not (0, default)
\n";
    exit 1;
}

if ($indirect_relations != 0 && $indirect_relations != 1) {
    print "\n\tInvalid value for -indirect_relations argument: $indirect_relations
              \t-indirect_relations should be either 0 or 1\n";
}
if ($indirect_relations) {
    print "Indirect relations of the closure will be inserted\n";
} else {
    print "Indirect relations of the closure will NOT be inserted: they must be inserted before the propagation ",
          "of the conditions, the step conditions_not_existing of post_processing fails without them\n";
}
$| = 1;

##########################################
# RETRIEVE RELATIONS
##########################################
my $fh;
# want to keep the order in %relations with Tie::IxHash
my %relations;
tie %relations, 'Tie::IxHash';
open($fh, '<', $file_path) || die "failed to read input file: $!";
while (my $line = <$fh>) {
    chomp $line;
    ## skip comments and blank lines and optional repeat of title line
    next if $line =~ /^\#/ || $line =~ /^\s*$/ || $line =~ /^\+/;
    #split each line into array
    my @line = split(/\s+/, $line);

    if (!$relations{$line[0]}) {
        my @parents = ();
        $relations{$line[0]} = \@parents;
    }
    push @{$relations{$line[0]}}, $line[1];
}

dump(%relations);

##########################################
# COMPUTE TAXON CONSTRAINTS AND INSERT REL
##########################################
# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

# First we need to get all species in Bgee for when a term exists in all species
my $querySpecies = $bgee->prepare('SELECT speciesId FROM species');
$querySpecies->execute()  or die $querySpecies->errstr;
my @allSpeciesIds = ();
while ( my @data = $querySpecies->fetchrow_array ){
    push(@allSpeciesIds, $data[0]);
}
$querySpecies->finish();
# important to sort for further comparisons
@allSpeciesIds = sort { $a <=> $b } @allSpeciesIds;

my $queryTermTCs = $bgee->prepare('SELECT speciesId FROM anatEntityTaxonConstraint WHERE anatEntityId = ?');
my $queryRel = $bgee->prepare('SELECT anatEntityRelationId FROM anatEntityRelation
                               WHERE anatEntitySourceId = ? AND anatEntityTargetId = ?
                               AND relationType = "is_a part_of" AND relationStatus = ?');
my $queryRelTCs = $bgee->prepare('SELECT speciesId FROM anatEntityRelationTaxonConstraint WHERE anatEntityRelationId = ?');
my $queryAncestors = $bgee->prepare('SELECT anatEntityRelationId, anatEntityTargetId FROM anatEntityRelation
                                     WHERE anatEntitySourceId = ? AND anatEntityTargetId != anatEntitySourceId
                                     AND relationType = "is_a part_of" AND relationStatus IN ("direct", "indirect")');
my $queryDescendants = $bgee->prepare('SELECT anatEntityRelationId, anatEntitySourceId FROM anatEntityRelation
                                       WHERE anatEntityTargetId = ? AND anatEntitySourceId != anatEntityTargetId
                                       AND relationType = "is_a part_of" AND relationStatus IN ("direct", "indirect")');
my $insertAnatTC = $bgee->prepare('INSERT INTO anatEntityTaxonConstraint (anatEntityId, speciesId) VALUES (?, ?)');
my $insertRel = $bgee->prepare('INSERT INTO anatEntityRelation (
                                    anatEntityRelationId,
                                    anatEntitySourceId,
                                    anatEntityTargetId,
                                    relationType,
                                    relationStatus)
                                    VALUES (?, ?, ?, "is_a part_of", ?)');
my $insertRelTC = $bgee->prepare('INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES (?, ?)');

my $maxRelIdStmt = $bgee->prepare('SELECT MAX(anatEntityRelationId) FROM anatEntityRelation');
$maxRelIdStmt->execute()  or die $maxRelIdStmt->errstr;
my ($maxRelId) = $maxRelIdStmt->fetchrow_array;
$maxRelIdStmt->finish();
print "Max relationId: $maxRelId\n";

# Taxon constraints are handled as hash references with species IDs as keys,
# a NULL species ID in the database (valid in all species) giving all the species
sub to_species_set {
    my %set = ();
    for my $speciesId (@_) {
        if (!defined $speciesId) {
            %set = map { $_ => 1 } @allSpeciesIds;
            last;
        }
        $set{$speciesId} = 1;
    }
    return \%set;
}
sub intersect_species {
    my ($set1, $set2) = @_;
    return { map { $_ => 1 } grep { $set2->{$_} } keys %$set1 };
}
sub union_species {
    my ($set1, $set2) = @_;
    return { %$set1, %$set2 };
}
sub is_all_species {
    my ($set) = @_;
    return scalar(grep { $set->{$_} } @allSpeciesIds) == scalar(@allSpeciesIds);
}
sub species_to_string {
    my ($set) = @_;
    return is_all_species($set)? 'all species': join(', ', sort { $a <=> $b } keys %$set);
}

sub term_species {
    my ($termId) = @_;
    $queryTermTCs->execute($termId)  or die $queryTermTCs->errstr;
    return to_species_set(map { $_->[0] } @{$queryTermTCs->fetchall_arrayref()});
}
sub relation_species {
    my ($relId) = @_;
    $queryRelTCs->execute($relId)  or die $queryRelTCs->errstr;
    return to_species_set(map { $_->[0] } @{$queryRelTCs->fetchall_arrayref()});
}
# Returns the [relation ID, related term ID] of the relations retrieved by $query for $termId
sub related_terms {
    my ($query, $termId) = @_;
    $query->execute($termId)  or die $query->errstr;
    return @{$query->fetchall_arrayref()};
}

# Gives $termId the taxon constraints of $species it does not have yet
sub add_term_species {
    my ($termId, $species) = @_;
    my $existing = term_species($termId);
    return if is_all_species($existing);
    if (!%$existing && is_all_species($species)) {
        print "INSERT INTO anatEntityTaxonConstraint (anatEntityId, speciesId) VALUES ($termId, null)\n";
        $insertAnatTC->execute($termId, undef)  or die $insertAnatTC->errstr;
        return;
    }
    for my $speciesId (sort { $a <=> $b } keys %$species) {
        next if $existing->{$speciesId};
        print "INSERT INTO anatEntityTaxonConstraint (anatEntityId, speciesId) VALUES ($termId, $speciesId)\n";
        $insertAnatTC->execute($termId, $speciesId)  or die $insertAnatTC->errstr;
    }
}
# Makes the relation valid in the species of $species, inserting it if it does not exist
sub add_relation {
    my ($sourceId, $targetId, $status, $species) = @_;
    return if !%$species;
    $queryRel->execute($sourceId, $targetId, $status)  or die $queryRel->errstr;
    my ($relId) = $queryRel->fetchrow_array;
    $queryRel->finish();
    my $existing = {};
    if (defined $relId) {
        $existing = relation_species($relId);
    } else {
        $maxRelId++;
        $relId = $maxRelId;
        print "INSERT INTO anatEntityRelation $relId, $sourceId, $targetId, '$status'\n";
        $insertRel->execute($relId, $sourceId, $targetId, $status)  or die $insertRel->errstr;
    }
    return if is_all_species($existing);
    if (!%$existing && is_all_species($species)) {
        print "INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES ($relId, null)\n";
        $insertRelTC->execute($relId, undef)  or die $insertRelTC->errstr;
        return;
    }
    for my $speciesId (sort { $a <=> $b } keys %$species) {
        next if $existing->{$speciesId};
        print "INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES ($relId, $speciesId)\n";
        $insertRelTC->execute($relId, $speciesId)  or die $insertRelTC->errstr;
    }
}

# The new terms ordered parents first: a new term is only processed once its new parents
# have their taxon constraints and relations
my @newTerms = ();
my %visitState = ();
sub visit_new_term {
    my ($termId) = @_;
    my $state = $visitState{$termId} // 0;
    return if $state == 2;
    die "Cycle between the new terms, through $termId\n" if $state == 1;
    $visitState{$termId} = 1;
    for my $parent (@{$relations{$termId}}) {
        visit_new_term($parent) if exists $relations{$parent};
    }
    $visitState{$termId} = 2;
    push @newTerms, $termId;
}
visit_new_term($_) for keys %relations;

print "Retrieve taxon constraints\n";
# The taxon constraints of a new term: those it already has, otherwise the intersection
# of the taxon constraints of its parents
my %newTermTCs = ();
for my $newTerm (@newTerms) {
    print "New term $newTerm\n";
    my $species = term_species($newTerm);
    if (!%$species) {
        $species = undef;
        for my $parent (@{$relations{$newTerm}}) {
            my $parentSpecies = exists $newTermTCs{$parent}? $newTermTCs{$parent}: term_species($parent);
            $species = defined $species? intersect_species($species, $parentSpecies): $parentSpecies;
        }
    }
    if (!%$species) {
        print "Could not compute taxon constraints for term $newTerm\n";
        exit 1;
    }
    print "Resulting TCs: ", species_to_string($species), "\n";
    $newTermTCs{$newTerm} = $species;
}
print "------------------\n";

for my $newTerm (@newTerms) {
    my $species = $newTermTCs{$newTerm};
    add_term_species($newTerm, $species);
    add_relation($newTerm, $newTerm, 'reflexive', $species);

    # The species of the relations from the new term to each of its new ancestors,
    # its parents included: the union, over the chains reaching the ancestor,
    # of the species where every relation of the chain is valid
    my %ancestorSpecies = ();
    for my $parent (@{$relations{$newTerm}}) {
        my $parentSpecies = exists $newTermTCs{$parent}? $newTermTCs{$parent}: term_species($parent);
        my $directSpecies = intersect_species($species, $parentSpecies);
        add_relation($newTerm, $parent, 'direct', $directSpecies);
        $ancestorSpecies{$parent} = union_species($ancestorSpecies{$parent} // {}, $directSpecies);
        next unless $indirect_relations;

        for my $relation (related_terms($queryAncestors, $parent)) {
            my ($relId, $ancestor) = @$relation;
            die "Cycle: $newTerm is an ancestor of its parent $parent\n" if $ancestor eq $newTerm;
            my $chainSpecies = intersect_species($directSpecies, relation_species($relId));
            next unless %$chainSpecies;
            add_relation($newTerm, $ancestor, 'indirect', $chainSpecies);
            $ancestorSpecies{$ancestor} = union_species($ancestorSpecies{$ancestor} // {}, $chainSpecies);
        }
    }
    next unless $indirect_relations;

    # The descendants of the new term already in the database reach its new ancestors through it
    for my $relation (related_terms($queryDescendants, $newTerm)) {
        my ($relId, $descendant) = @$relation;
        my $descendantSpecies = relation_species($relId);
        for my $ancestor (keys %ancestorSpecies) {
            die "Cycle: $descendant is an ancestor of $newTerm\n" if $ancestor eq $descendant;
            add_relation($descendant, $ancestor, 'indirect',
                intersect_species($descendantSpecies, $ancestorSpecies{$ancestor}));
        }
    }
}

$queryTermTCs->finish();
$queryRel->finish();
$queryRelTCs->finish();
$queryAncestors->finish();
$queryDescendants->finish();
$insertAnatTC->finish();
$insertRel->finish();
$insertRelTC->finish();

$bgee->disconnect();
exit 0;
