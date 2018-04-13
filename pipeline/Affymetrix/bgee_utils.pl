#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Frederic Bastian, created November 2012
# list of subs useful for several parts of the pipeline
#############################################################

$| = 1;

# Check whether a given organ does exist during a given developmental stage.
# return 1 if it does, 0 otherwise
sub organExistsAtStage {
    # $bgeeConnection allows to pass an already opened connection to the Bgee database.
    my ($anatEntityId, $stageId, $bgeeConnection) = @_;

    my $dbh = $bgeeConnection->prepare('SELECT 1 FROM anatEntity AS t1 '.
        'INNER JOIN stage AS t2    ON t1.startStageId = t2.stageId '.
        'INNER JOIN stage AS t2bis ON t2.leftBound    <= t2bis.rightBound AND t2.speciesId  = t2bis.speciesId '.
        'INNER JOIN stage AS t3    ON t1.endStageId   = t3.stageId '.
        'INNER JOIN stage AS t3bis ON t3.rightBound   >= t3bis.leftBound  AND t3.speciesId  = t3bis.speciesId '.
        'INNER JOIN stage AS t0    ON t2bis.stageId   = t0.stageId        AND t3bis.stageId = t0.stageId '.
        'WHERE (t1.anatEntityId = ? AND t0.stageId = ?) '.
        'OR t0.stageName = "unknown" LIMIT 1');
    $dbh->execute($anatEntityId, $stageId)  or die $dbh->errstr;
    if ( my @data = $dbh->fetchrow_array ){
        return 1;
    }
    return 0;
}

# Does it really need a comment?
# prefix the sub name to not collision with perl 6...
sub bgeeTrim {
    my ($stringToTrim) = @_;
    if ( defined $stringToTrim ){
        $stringToTrim =~ s{^\s+|\s+$}{}g;
        $stringToTrim =~ s{\s\s+}{ }g;
        $stringToTrim =~ s{^"|"$}{}g; # For quote-protected cells
    }
    return $stringToTrim;
}

1;

