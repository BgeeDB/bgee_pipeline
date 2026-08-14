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
# Add the option to also insert indirect relations to all ancestors of the parent term, with the same TCs as the relation to the parent term. This option is managed with the -indirect_relations argument, which should be set to 1 to insert indirect relations or 0 (default) to not insert them.
# Added IGNORE in the SQL insert statements to ignore the insertion if a term TC, relation or relation TC already exists for a given term or relation.
#############################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my ($bgee_connector) = ('');
my $file_path        = '';
my $indirect_relations = 0; # whether to insert indirect relations to all ancestors of the parent term (with the same TCs as the relation to the parent term)
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
\t-indirect_relations   Whether to insert indirect relations(1) to all ancestors of the parent term (with the same TCs as the relation to the parent term) or not (0, default)
\n";
    exit 1;
}

if ($indirect_relations != 0 && $indirect_relations != 1) {
    print "\n\tInvalid value for -indirect_relations argument: $indirect_relations
              \t-indirect_relations should be either 0 or 1\n";
}
if ($indirect_relations) {
    print "Indirect relations to ancestors of the parent term will be inserted\n";
} else {
    print "Indirect relations to ancestors of the parent term will NOT be inserted\n";
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

print "Retrieve taxon constraints\n";
# Now for each new term, we retrieve the intersection of the taxon constraints of its parents
my %newTermTCs = ();
my $queryTCs = $bgee->prepare('SELECT speciesId FROM anatEntityTaxonConstraint WHERE anatEntityId = ?');
while(my ($newTerm, $parents) = each(%relations)) {
    print "New term $newTerm\n";
    my @TCParents = ();
    foreach my $parent (@$parents) {
        $queryTCs->execute($parent)  or die $queryTCs->errstr;
        my @TCs = ();
        while ( my @data = $queryTCs->fetchrow_array ){
            if (!$data[0]) {
                print "Push all species\n";
                push(@TCs, @allSpeciesIds);
            } else {
                push(@TCs, $data[0]);
            }
        }
        push(@TCParents, \@TCs);
    }
    dump(@TCParents);
    my(%count,@res);
    for(map @$_, @TCParents){$count{$_}++}
    for(keys %count){push @res, $_ if $count{$_}==@TCParents}

    @res = sort { $a <=> $b } @res;
    if (scalar(@res) == 0) {
        print "Could not compute taxon constraints for term $newTerm\n";
        exit 1;
    }
    if ( @res ~~ @allSpeciesIds ) {
        # Empty array means "valid in all species".
        @res = ();
    }
    print "Resulting TCs:\n";
    dump(@res);
    $newTermTCs{$newTerm} = \@res;
}
$queryTCs->finish();
    print "------------------\n";
    dump(%newTermTCs);
    


my $maxRelIdStmt = $bgee->prepare('SELECT MAX(anatEntityRelationId) FROM anatEntityRelation');
my $maxRelId;
$maxRelIdStmt->execute()  or die $querySpecies->errstr;
if ( my @data = $maxRelIdStmt->fetchrow_array ){
    $maxRelId = $data[0];
}
$maxRelIdStmt->finish();
print "Max relationId: $maxRelId\n";

# in case a term TC, relation or relation TC already exists we ignore the insertion.
my $insertAnatTCs = $bgee->prepare('INSERT IGNORE INTO anatEntityTaxonConstraint (anatEntityId, speciesId) VALUES (?, ?)');
my $insertRel = $bgee->prepare('INSERT IGNORE INTO anatEntityRelation (
                                    anatEntityRelationId,
                                    anatEntitySourceId,
                                    anatEntityTargetId,
                                    relationType,
                                    relationStatus)
                                    VALUES (?, ?, ?, "is_a part_of", ?)');
my $insertRelTCs = $bgee->prepare('INSERT IGNORE INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES (?, ?)');

my $queryAncestors = $bgee->prepare('SELECT anatEntityRelationId, anatEntityTargetId FROM anatEntityRelation WHERE anatEntitySourceId = ? AND relationType IN ("is_a part_of") AND relationStatus IN ("direct", "indirect")');
my $queryRelTCsAncestors = $bgee->prepare('SELECT speciesId FROM anatEntityRelationTaxonConstraint WHERE anatEntityRelationId = ?');

my $queryAncestors = $bgee->prepare('SELECT anatEntityRelationId, anatEntityTargetId FROM anatEntityRelation WHERE anatEntitySourceId = ? AND relationType IN ("is_a part_of") AND relationStatus IN ("direct", "indirect")');
my $queryRelTCsAncestors = $bgee->prepare('SELECT speciesId FROM anatEntityRelationTaxonConstraint WHERE anatEntityRelationId = ?');

while(my ($newTerm, $tcs) = each(%newTermTCs)) {    
    $maxRelId++;
    print "INSERT INTO anatEntityRelation $maxRelId, $newTerm, $newTerm, 'reflexive'\n";
    $insertRel->execute($maxRelId, $newTerm, $newTerm, 'reflexive')  or die $insertRel->errstr;
	if (scalar(@$tcs) == 0) {
		print "INSERT INTO anatEntityTaxonConstraint (anatEntityId, speciesId) VALUES ($newTerm, null)\n";
		$insertAnatTCs->execute($newTerm, undef)  or die $insertAnatTCs->errstr;
		print "INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES ($maxRelId, null)\n";
		$insertRelTCs->execute($maxRelId, undef)  or die $insertRelTCs->errstr;
	} else {
	    foreach my $tc (@$tcs) {
		    print "INSERT INTO anatEntityTaxonConstraint (anatEntityId, speciesId) VALUES ($newTerm, $tc)\n";
	        $insertAnatTCs->execute($newTerm, $tc)  or die $insertAnatTCs->errstr;
		    print "INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES ($maxRelId, $tc)\n";
		    $insertRelTCs->execute($maxRelId, $tc)  or die $insertRelTCs->errstr;
	    }
    }
    foreach my $parent (@{$relations{$newTerm}}) {
        $maxRelId++;
        print "INSERT INTO anatEntityRelation $maxRelId, $newTerm, $parent, 'direct'\n";
        $insertRel->execute($maxRelId, $newTerm, $parent, 'direct')  or die $insertRel->errstr;
           if (scalar(@$tcs) == 0) {
            print "INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES ($maxRelId, null)\n";
            $insertRelTCs->execute($maxRelId, undef)  or die $insertRelTCs->errstr;
        } else {
            foreach my $tc (@$tcs) {
                print "INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES ($maxRelId, $tc)\n";
                $insertRelTCs->execute($maxRelId, $tc)  or die $insertRelTCs->errstr;
            }
        }

        # we also need to insert indirect relations to all ancestors of the parent term, with the same TCs as the relation to the parent term
        if ($indirect_relations) {
            $queryAncestors->execute($parent)  or die $queryAncestors->errstr;
            while ( my @data = $queryAncestors->fetchrow_array ){
                my $relId = $data[0];
                my $ancestor = $data[1];
                $maxRelId++;
                print "INSERT INTO anatEntityRelation $maxRelId, $newTerm, $ancestor, 'indirect'\n";
                $insertRel->execute($maxRelId, $newTerm, $ancestor, 'indirect')  or die $insertRel->errstr;
                # now we retrieve the taxon constraints of the relation from the parent term to the ancestor to insert the same TCs for the new relation between the new term and the ancestor
                $queryRelTCsAncestors->execute($relId)  or die $queryRelTCsAncestors->errstr;
                my @relTCs = ();
                while ( my @data = $queryRelTCsAncestors->fetchrow_array ){
                    if (!$data[0]) {
                        print "Push all species\n";
                        push(@relTCs, @allSpeciesIds);
                    } else {
                        push(@relTCs, $data[0]);
                    }
                }
                if (scalar(@relTCs) == scalar(@allSpeciesIds)) {
                    print "INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES ($maxRelId, null)\n";
                    $insertRelTCs->execute($maxRelId, undef)  or die $insertRelTCs->errstr;
                } else {
                    foreach my $tc (@relTCs) {
                        print "INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES ($maxRelId, $tc)\n";
                        $insertRelTCs->execute($maxRelId, $tc)  or die $insertRelTCs->errstr;
                    }
                }
            }
        }
       }
}
$insertAnatTCs->finish();
$insertRel->finish();
$insertRelTCs->finish();
$queryAncestors->finish();
$queryRelTCsAncestors->finish();

$bgee->disconnect();
exit 0;