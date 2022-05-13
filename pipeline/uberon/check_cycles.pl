#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Frederic Bastian, created November 2016
# Identifies cycls in Uberon relations inserted in Bgee
#############################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my ($bgee_connector) = ('');
my $debug            = 0;
my %opts = ('bgee=s'     => \$bgee_connector,   # Bgee connector string
            'debug|v'    => \$debug,            # more verbose output
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

$| = 1;

##########################################
# RETRIEVE CYCLES                        #
##########################################
# Search cycles
my $queryCycles = $bgee->prepare('SELECT t1.anatEntitySourceId, t3.anatEntityName, t1.relationType,
                                         t1.anatEntityTargetId, t4.anatEntityName,
                                         GROUP_CONCAT(DISTINCT CONCAT_WS("-",
                                                      IF(t2.speciesId IS NULL, "NULL", t2.speciesId),
                                                      IF(t6.speciesId IS NULL, "NULL", t6.speciesId))
                                                      ORDER BY CONCAT_WS("-",
                                                               IF(t2.speciesId IS NULL, "NULL", t2.speciesId),
                                                               IF(t6.speciesId IS NULL, "NULL", t6.speciesId))
                                                      SEPARATOR ",") AS speciesIds
                                  FROM anatEntityRelation AS t1
                                  INNER JOIN anatEntityRelationTaxonConstraint AS t2
                                      ON t1.anatEntityRelationId = t2.anatEntityRelationId
                                  INNER JOIN anatEntity AS t3
                                      ON t1.anatEntitySourceId = t3.anatEntityId
                                  INNER JOIN anatEntity AS t4
                                      ON t1.anatEntityTargetId = t4.anatEntityId
                                  INNER JOIN anatEntityRelation AS t5
                                      ON t5.anatEntityTargetId = t1.anatEntitySourceId
                                      AND t5.anatEntitySourceId = t1.anatEntityTargetId
                                      AND t1.relationType = t5.relationType
                                  INNER JOIN anatEntityRelationTaxonConstraint AS t6
                                      ON t5.anatEntityRelationId = t6.anatEntityRelationId
                                  WHERE (t2.speciesId IS NULL OR t6.speciesId IS NULL OR t2.speciesId = t6.speciesId)
                                      AND t1.relationStatus != "reflexive" AND t5.relationStatus != "reflexive"
                                  GROUP BY t1.anatEntitySourceId, t1.anatEntityTargetId, t1.relationType;');
$queryCycles->execute()  or die $queryCycles->errstr;
my @cycles = ();
while ( my @data = $queryCycles->fetchrow_array ){
    my %cycle;
    $cycle{'sourceId'}     = $data[0];
    $cycle{'relationType'} = $data[2];
    $cycle{'targetId'}     = $data[3];
    $cycle{'speciesIds'}   = $data[5];
    push(@cycles, \%cycle);
}
$queryCycles->finish();


##########################################
# RETRIEVE ALL DIRECT RELATIONS          #
# TO BE ABLE TO BUILD A PATH IN CYCLES   #
##########################################
my $queryRels = $bgee->prepare('SELECT t1.anatEntitySourceId, t1.relationType, t1.anatEntityTargetId,
                                       t1.relationStatus,
                                       GROUP_CONCAT(DISTINCT t2.speciesId
                                                    ORDER BY t2.speciesId
                                                    SEPARATOR ",")
                                FROM anatEntityRelation AS t1
                                INNER JOIN anatEntityRelationTaxonConstraint AS t2
                                    ON t1.anatEntityRelationId = t2.anatEntityRelationId
                                WHERE t1.relationStatus != "reflexive"
                                GROUP BY t1.anatEntityRelationId
                                ORDER BY t1.anatEntitySourceId, t1.relationType, t1.anatEntityTargetId,
                                    t1.relationStatus');
$queryRels->execute()  or die $queryRels->errstr;
# 'direct/indirect' -> 'sourceId' -> relation
my %rels;
while ( my @data = $queryRels->fetchrow_array ){
    my %rel;
    $rel{'sourceId'}       = $data[0];
    $rel{'relationType'}   = $data[1];
    $rel{'targetId'}       = $data[2];
    $rel{'relationStatus'} = $data[3];
    $rel{'speciesIds'}     = $data[4];

    if (!defined $rels{$data[3]}->{$data[0]}) {
        $rels{$data[3]}->{$data[0]} = ();
    }
    push(@{$rels{$data[3]}->{$data[0]}}, \%rel);
}
$queryRels->finish();

##########################################
# RETRIEVE ANAT ENTITIES                 #
##########################################
my $queryAnat = $bgee->prepare('SELECT anatEntityId, anatEntityName FROM anatEntity');
$queryAnat->execute()  or die $queryAnat->errstr;
# ID = name
my %organs;
while ( my @data = $queryAnat->fetchrow_array ){
    $organs{$data[0]} = $data[1];
}
$queryAnat->finish();

$bgee->disconnect;

#############################################
# NOW TRIES TO IDENTIFY NON-REDUNDANT PATHS #
# THAT INCLUDES ALL CYCLES                  #
#############################################
# Note that, for simplicity, we don't verify the species in which relations are valid here,
# so maybe we'll identify more paths for a given cycle that there actually are.
# We will store all paths producing cycles at the end, to count number of times they are involved.
# key -> 'path'/'count'
my %allPaths;
my %subPaths;
CYCLE: for my $cycle ( @cycles ){
    #For now, consider only cycles over is_a part_of
    if ($cycle->{'relationType'} ne 'is_a part_of') {
        next CYCLE;
    }
    print "Cycle: $cycle->{'sourceId'} $cycle->{'relationType'} $cycle->{'targetId'}\n";
    my $pathFound = 0;
    # Start with the source of the cycle, walk until we reach the target of the cycle,
    # and keep walking until we get back to the source.
    # @pathsToWalk: array storing each path going from the source to the target and back to the source.
    # a path is tored as an array of term IDs, so we have an array of arrays here.
    my @pathsToWalk = ();
    # Initialize the walk with the source of the cycle, in an array that initialize the walk
    my @init = ($cycle->{'sourceId'});
    push(@pathsToWalk, \@init);
    # Will store the valid paths source -> target -> source
    my @validatedPaths = ();

    # Keep iterating all possible paths until we have reached the target back to the source,
    # or reached the top of the ontology. We'll stop when there will be no more paths to walk
    # in @pathsToWalk
    while (@pathsToWalk) {
        # retrieve each path being walked, and remove them from the paths being walk;
        # we'll add them back if there is more to walk for this path (like a Dequeue).
        my $walk = pop @pathsToWalk;

        # retrieve the last element of the walk (negative index returns last element)
        my $lastElement = $walk->[-1];
        # did we already get through the target of the cycle, trying to get back to the source?
        my $targetSeen = 0;
        ELEMENT: for my $element (@{$walk}) {
            if ($element eq $cycle->{'targetId'}) {
                $targetSeen = 1;
                last ELEMENT;
            }
        }

        # retrieve the next elements to walk, create a new path for each one of them
        REL: for my $rel (@{$rels{'direct'}->{$lastElement}}) {
            if ($debug) {
                print "\tIterating relation for $lastElement $rel->{'relationType'} $rel->{'targetId'}: ";
            }
            if (!validateRelationType($cycle, $rel)) {
                print "INVALID for cycle\n";
                next REL;
            }
            print "VALID for cycle\n";

            my $nextElement = $rel->{'targetId'};

            # copy the current path walked to create a new path
            my @newWalk = @{$walk};
            push(@newWalk, $nextElement);
            if ($debug) {
                print "\t\ttest next element: $nextElement (";
                for my $element (@newWalk) {
                    print "$element - ";
                }
                print ")\n";
            }

            # if we get back to the source after having visited the target,
            # then we validate this path.
            if ($nextElement eq $cycle->{'sourceId'} && $targetSeen) {
                # different relations can lead to a same path,
                # verify that we didn't already store this path
                my $alreadyStored = 0;
                STORED_PATH: for my $storedPath (@validatedPaths) {
                    if (scalar @newWalk != scalar @{$storedPath}) {
                        next STORED_PATH;
                    }
                    my $i = 0;
                    for my $newElement (@newWalk) {
                        if ($newElement ne $storedPath->[$i]) {
                            next STORED_PATH;
                        }
                        $i++;
                    }
                    $alreadyStored = 1;
                    last STORED_PATH;
                }
                if (!$alreadyStored) {
                    push(@validatedPaths, \@newWalk);
                    $pathFound = 1;
                    if ($debug) {
                        print "\tPath validated!\n";
                    }
                } elsif ($debug) {
                    print "\tPath validated but already stored!\n";
                }
                next REL;
            }

            # protection against cycles: maybe we are walking a different cycle than the one
            # currently being studied
            my $quotedNextElement = quotemeta($nextElement);
            if (grep( /^$quotedNextElement$/, @{$walk})) {
                if ($debug) {
                    print "\t\tinvalid cycle detected\n";
                }
                next REL;
            }


            # check whether the target of the cycle is on the path of the next element,
            # (or the source of the cycle if we already went through the target of the cycle)
            my $onPath = 0;
            if ($nextElement eq $cycle->{'targetId'}) {
                $onPath = 1;
            } else {
                # consider both the direct and indirect rels of the next element
                my @nextRels = ();
                if (defined $rels{'indirect'}->{$nextElement}) {
                    push(@nextRels, @{$rels{'indirect'}->{$nextElement}});
                }
                if (defined $rels{'direct'}->{$nextElement}) {
                    push(@nextRels, @{$rels{'direct'}->{$nextElement}});
                }
                NEXT_REL: for my $nextRel (@nextRels) {
                    if ($debug) {
                        print "\t\t\ttest is on path: $nextRel->{'relationType'} $nextRel->{'targetId'}: ";
                    }
                    if (validateRelationType($cycle, $nextRel) &&
                        (!$targetSeen && $nextRel->{'targetId'} eq $cycle->{'targetId'} ||
                          $targetSeen && $nextRel->{'targetId'} eq $cycle->{'sourceId'})) {
                        $onPath = 1;
                        if ($debug) {
                            print "yes\n";
                        }
                        last NEXT_REL;
                    } elsif ($debug) {
                        print "no\n";
                    }
                }
            }
            # on path, add this path to the paths to walk
            if ($onPath) {
                # add the new path to walk in first position, like a Dequeue
                unshift(@pathsToWalk, \@newWalk);
                if ($debug) {
                    print "\t\t\t=> $nextElement added to path\n";
                }
            } elsif ($debug) {
                print "\t\tLast element not on path to (targetSeen=$targetSeen? target $cycle->{'targetId'}: source $cycle->{'sourceId'})\n";
            }
        }
    }
    if (!$pathFound) {
        print "\tNo path found for cycle $cycle->{'sourceId'} $cycle->{'relationType'} $cycle->{'targetId'}, maybe you can delete the indirect relation\n";
    } else {
        print "\tAll paths for cycle:\n";
        for my $path (@validatedPaths) {
            print "\t\t";
            my $i = 0;
            my $key = "";
            for my $element (@{$path}) {
                print "$element \"$organs{$element}\" - ";
                $key .= $element.' - ';
                if ($i >= 2) {
                    my $subPathKey = $path->[$i-2].' - '.$path->[$i-1].' - '.$path->[$i];
                    my @subPath = ($path->[$i-2], $path->[$i-1], $path->[$i]);
                    $subPaths{$subPathKey}->{'path'} = \@subPath;
                    $subPaths{$subPathKey}->{'count'}++;
                    if (!defined $subPaths{$subPathKey}->{'containingPaths'}) {
                        $subPaths{$subPathKey}->{'containingPaths'} = ();
                    }

                    my $containingPathAlreadyStored = 0;
                    CONT_PATH: for my $containingPath (@{$subPaths{$subPathKey}->{'containingPaths'}}) {
                        if (scalar @{$containingPath} != scalar @{$path}) {
                            next CONT_PATH;
                        }
                        my $i = 0;
                        for my $newElement (@{$path}) {
                            if ($newElement ne $containingPath->[$i]) {
                                next CONT_PATH;
                            }
                            $i++;
                        }
                        $containingPathAlreadyStored = 1;
                        last CONT_PATH;
                    }
                    if (!$containingPathAlreadyStored) {
                        push(@{$subPaths{$subPathKey}->{'containingPaths'}}, $path);
                    }
                }
                $i++;
            }
            $allPaths{$key}->{'path'} = $path;
            $allPaths{$key}->{'count'}++;
            print "\n";
        }
        print "\n";
    }
}


print "\nNow display most represented subgraphs in cycles\n";
for my $key (sort { $subPaths{$b}->{'count'} <=> $subPaths{$a}->{'count'} } keys %subPaths) {
    print "\t$subPaths{$key}->{'count'}: ";
    for my $element (@{$subPaths{$key}->{'path'}}) {
        print "$element \"$organs{$element}\" - ";
    }
    print "\n";
}

print "\nNow display paths ordered by number of times they were involved in a cycle\n";
for my $key (sort { $allPaths{$b}->{'count'} <=> $allPaths{$a}->{'count'}
                    or scalar $allPaths{$a}->{'path'} <=> scalar $allPaths{$b}->{'path'} } keys %allPaths) {
    print "\t$allPaths{$key}->{'count'}: ";
    for my $element (@{$allPaths{$key}->{'path'}}) {
        print "$element \"$organs{$element}\" - ";
    }
    print "\n";
}

exit 0;

# Determine whether the 'relationType' of a relation is compatible with the 'relationType'
# of a cycle (meaning, the composed object property over which there is a cycle).
# Warning, this method assumes that only 3 types of relations exist in Bgee:
# 'is_a part_of', 'develops_from', 'transformation_of'.
sub validateRelationType {
    my ($cycle, $rel) = @_;
    if ($cycle->{'relationType'} eq 'is_a part_of' &&
            $rel->{'relationType'} ne $cycle->{'relationType'}) {
        return 0;
    }
    # 'transformation_of' is a sub-property of 'develops_from',
    # so if the cycle is over a 'develops_from' property,
    # a 'transformation_of' property is also valid for the relation (but not the other way around).
    #
    # And also, because for storage in Bgee database we do not distinguish is_a/part_of relations,
    # a develops_from/transformation_of relation could be propagated through an is_a relation.
    # So, for cycles over 'transformation_of' or 'develops_from' relationTypes,
    # an 'is_a part_of' relation is also valid.
    #
    # So the only invaid remaining case is when cycle is over 'transformation_of'
    # and relation is 'develops_from'.
    if ($cycle->{'relationType'} eq 'transformation_of' &&
            $rel->{'relationType'} eq 'develops_from') {
        return 0;
    }
    return 1;
}
