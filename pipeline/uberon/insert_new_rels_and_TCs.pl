#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Data::Dump qw(dump);

# Frederic Bastian, created May 2024
# Insert relations for new anatomical terms added to the database,
# and infer their taxon constraints from their parents.
# The script takes as input a TSV file with one relation per line,
# in the format:
# OBO child ID\tOBO parent ID
# No header to the TSV file.
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
\te.g. $0 -file_path=/path/to/relation/file -bgee=connection_string
\t-file_path        Path to relation file
\t-bgee             Bgee connector string
\n";
    exit 1;
}

$| = 1;

##########################################
# RETRIEVE RELATIONS                     
##########################################
my $fh;
my %relations = ();
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

my $insertAnatTCs = $bgee->prepare('INSERT INTO anatEntityTaxonConstraint (anatEntityId, speciesId) VALUES (?, ?)');
my $insertRel = $bgee->prepare('INSERT INTO anatEntityRelation (
                                    anatEntityRelationId,
                                    anatEntitySourceId,
                                    anatEntityTargetId,
                                    relationType,
                                    relationStatus)
                                    VALUES (?, ?, ?, "is_a part_of", ?)');
my $insertRelTCs = $bgee->prepare('INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId) VALUES (?, ?)');

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
   	}
}
$insertAnatTCs->finish();
$insertRel->finish();
$insertRelTCs->finish();
$bgee->disconnect();
exit 0;