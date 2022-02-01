#!/usr/bin/env perl

# Philippe Moret, created Nov 2015.
# Frederic Bastian, last updated June 2016.
# Normalize the mean rank in the expression table across the different data types.
# Frederic Bastian, last updated Feb. 2017: adapt to new conditions and new schema in Bgee 14
# Julien Wollbrett, August 2021, allows to normalize ranks for a subset of species. Allows to parallelize ranks generation per species
# Julien Wollbrett, January 2022, allows to select a subset of Bgee condition IDs when exactly one speciesId is provided. Avoid running long transactions updating huge number of rows (e.g 2 billion of rows for mouse as for Bgee 15.0)
use strict;
use warnings;
use Time::HiRes qw( time );
use diagnostics;

use FindBin;
use lib "$FindBin::Bin/..";    # Get lib path for Utils.pm
use Utils;
use Getopt::Long;

$| = 1;

# Define arguments and their default value
my ($bgee_connector, $bgee_species, $condition_ids) = ('', '', '');
my %opts = ( 'bgee=s'       => \$bgee_connector, # Bgee connector string
             'species=s'    => \$bgee_species, # list of species separated by a comma
             'conditions=s' => \$condition_ids); #list of Bgee condition Ids separated by a comma
# Check arguments
my $emptyArg = '-';

my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $bgee_species eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -bgee_species='-'
\t-bgee            Bgee connector string
\t-species         Bgee species for which ranks have to be normalized or '-' for all species
\t-conditions      Subset of Bgee condition IDs for which ranks have to be normalized or '-'
                   for all conditions. Subset of condition IDs can only be provided when exactly one species is selected
\n";
    exit 1;
}

# As of Bgee 15.1, these blacklisted terms are directly remapped to the root of the anatEntities,
# so there's no need to update their ranks to the max rank anymore
#my $blacklisted = "('XAO:0003003', 'ZFA:0001093')";

my @speciesList = ();
if ($bgee_species ne $emptyArg) {
    @speciesList = split(',', $bgee_species);
}
if(scalar @speciesList != 1 && $condition_ids ne $emptyArg) {
    print "Condition IDs can be provided only when exactly one species is selected";
    exit 1;
}

my $dbh = Utils::connect_bgee_db($bgee_connector);

#we get the absolute max rank across all conditions for each species
# If condition IDs are provided they are not used in this query. We really want the max rank across all conditions
my $absoluteMaxSql = "SELECT speciesId, GREATEST(
             IFNULL(max(rnaSeqMaxRank), 0.0),
             IFNULL(max(estMaxRank),0.0),
             IFNULL(max(inSituMaxRank),0.0),
             IFNULL(max(affymetrixMaxRank),0.0),
             IFNULL(max(scRnaSeqFullLengthMaxRank),0.0),
             IFNULL(max(rnaSeqGlobalMaxRank), 0.0),
             IFNULL(max(estGlobalMaxRank),0.0),
             IFNULL(max(inSituGlobalMaxRank),0.0),
             IFNULL(max(affymetrixGlobalMaxRank),0.0),
             IFNULL(max(scRnaSeqFullLengthGlobalMaxRank),0.0)) AS max
             FROM globalCond ";
if (@speciesList) {
    $absoluteMaxSql .= ' WHERE speciesId IN (';
    for my $i (0 .. $#speciesList) {
        if ($i > 0) {
            $absoluteMaxSql .= ', ';
        }
        $absoluteMaxSql .= $speciesList[$i];
    }
    $absoluteMaxSql .= ')';
}
#    WHERE anatEntityId NOT IN $blacklisted
$absoluteMaxSql .=" GROUP BY speciesId";

my $getAbsoluteMax = $dbh->prepare($absoluteMaxSql);
$getAbsoluteMax->execute()  or die $getAbsoluteMax->errstr;
my %maxRanks = ();
while ( my @data = $getAbsoluteMax->fetchrow_array ){
    $maxRanks{$data[0]} = $data[1];
    print("Max rank for species $data[0]: $data[1]\n");
}

# update the expression table with normalized mean ranks
my $updateExpressionQuery = "UPDATE gene
STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
STRAIGHT_JOIN globalCond ON globalCond.globalConditionId = globalExpression.globalConditionId
SET rnaSeqMeanRankNorm     = (rnaSeqMeanRank + (rnaSeqMeanRank * ? / rnaSeqMaxRank))/2,
    estRankNorm            = (estRank + (estRank * ? / estMaxRank))/2,
    inSituRankNorm         = (inSituRank + (inSituRank * ? / inSituMaxRank))/2,
    affymetrixMeanRankNorm = (affymetrixMeanRank + (affymetrixMeanRank * ? / affymetrixMaxRank))/2,
    scRnaSeqFullLengthMeanRankNorm = (scRnaSeqFullLengthMeanRank + (scRnaSeqFullLengthMeanRank * ? / scRnaSeqFullLengthMaxRank))/2,
    rnaSeqGlobalMeanRankNorm     = (rnaSeqGlobalMeanRank + (rnaSeqGlobalMeanRank * ? / rnaSeqGlobalMaxRank))/2,
    estGlobalRankNorm            = (estGlobalRank + (estGlobalRank * ? / estGlobalMaxRank))/2,
    inSituGlobalRankNorm         = (inSituGlobalRank + (inSituGlobalRank * ? / inSituGlobalMaxRank))/2,
    affymetrixGlobalMeanRankNorm = (affymetrixGlobalMeanRank + (affymetrixGlobalMeanRank * ? / affymetrixGlobalMaxRank))/2,
    scRnaSeqFullLengthGlobalMeanRankNorm = (scRnaSeqFullLengthGlobalMeanRank + (scRnaSeqFullLengthGlobalMeanRank * ? / scRnaSeqFullLengthGlobalMaxRank))/2
WHERE gene.speciesId = ?";
if ($condition_ids ne $emptyArg) {
    $updateExpressionQuery .= " AND globalConditionId IN (?)";
}
my $updateExpression = $dbh->prepare($updateExpressionQuery);


    # As of Bgee 15.1, these blacklisted terms are directly remapped to the root of the anatEntities,
    # so there's no need to update their ranks to the max rank anymore
#    my $blacklistUnspecifiedEST = $dbh->prepare( "
#    UPDATE gene
#    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
#    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
#        AND globalCond.globalConditionId = globalExpression.globalConditionId
#    SET estRankNorm = ?
#    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND estRank IS NOT NULL" );
#    my $blacklistGlobalUnspecifiedEST = $dbh->prepare( "
#    UPDATE gene
#    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
#    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
#        AND globalCond.globalConditionId = globalExpression.globalConditionId
#    SET estGlobalRankNorm = ?
#    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND estGlobalRank IS NOT NULL" );
#
#    my $blacklistUnspecifiedInSitu = $dbh->prepare( "
#    UPDATE gene
#    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
#    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
#        AND globalCond.globalConditionId = globalExpression.globalConditionId
#    SET inSituRankNorm = ?
#    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND inSituRank IS NOT NULL" );
#    my $blacklistGlobalUnspecifiedInSitu = $dbh->prepare( "
#    UPDATE gene
#    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
#    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
#        AND globalCond.globalConditionId = globalExpression.globalConditionId
#    SET inSituGlobalRankNorm = ?
#    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND inSituGlobalRank IS NOT NULL" );
#
#    my $blacklistUnspecifiedAffymetrix = $dbh->prepare( "
#    UPDATE gene
#    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
#    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
#        AND globalCond.globalConditionId = globalExpression.globalConditionId
#    SET affymetrixMeanRankNorm = ?
#    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND affymetrixMeanRank IS NOT NULL" );
#    my $blacklistGlobalUnspecifiedAffymetrix = $dbh->prepare( "
#    UPDATE gene
#    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
#    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
#        AND globalCond.globalConditionId = globalExpression.globalConditionId
#    SET affymetrixGlobalMeanRankNorm = ?
#    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND affymetrixGlobalMeanRank IS NOT NULL" );
#
#    my $blacklistUnspecifiedRnaSeq = $dbh->prepare( "
#    UPDATE gene
#    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
#    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
#        AND globalCond.globalConditionId = globalExpression.globalConditionId
#    SET rnaSeqMeanRankNorm = ?
#    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND rnaSeqMeanRank IS NOT NULL" );
#    my $blacklistGlobalUnspecifiedRnaSeq = $dbh->prepare( "
#    UPDATE gene
#    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
#    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
#        AND globalCond.globalConditionId = globalExpression.globalConditionId
#    SET rnaSeqGlobalMeanRankNorm = ?
#    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND rnaSeqGlobalMeanRank IS NOT NULL" );
#
#    my $blacklistUnspecifiedScRnaSeqFL = $dbh->prepare( "
#    UPDATE gene
#    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
#    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
#        AND globalCond.globalConditionId = globalExpression.globalConditionId
#    SET scRnaSeqFullLengthMeanRankNorm = ?
#    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND scRnaSeqFullLengthMeanRank IS NOT NULL" );
#    my $blacklistGlobalUnspecifiedScRnaSeqFL = $dbh->prepare( "
#    UPDATE gene
#    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
#    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
#        AND globalCond.globalConditionId = globalExpression.globalConditionId
#    SET scRnaSeqFullLengthGlobalMeanRankNorm = ?
#    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND scRnaSeqFullLengthGlobalMeanRank IS NOT NULL" );

for my $speciesId ( keys %maxRanks ){
    Utils::start_transaction($dbh);
    my $absMax = $maxRanks{$speciesId};
    print "Updating ranks for species $speciesId with max rank $absMax\n";

    my $t0 = time();
    printf('Update expression table with normalized mean ranks per type:   ');
    if ($condition_ids ne $emptyArg) {
         $updateExpression->execute( $absMax, $absMax, $absMax, $absMax, $absMax, $absMax, $absMax, $absMax,
                    $absMax, $absMax, $speciesId, $condition_ids) or die $updateExpression->errstr;
    } else {
            $updateExpression->execute( $absMax, $absMax, $absMax, $absMax, $absMax, $absMax, $absMax, $absMax,
                    $absMax, $absMax, $speciesId) or die $updateExpression->errstr;

    }

    printf( "OK in %.2fs\n", ( time() - $t0 ) );

    # As of Bgee 15.1, these blacklisted terms are directly remapped to the root of the anatEntities,
    # so there's no need to update their ranks to the max rank anymore
#    $t0 = time();
#    printf('Blacklisting unspecified anat entities:    ');
#    my $blacklistAbsMax = $absMax;
#    $blacklistUnspecifiedEST->execute($blacklistAbsMax, $speciesId)
#      or die $blacklistUnspecifiedEST->errstr;
#    $blacklistGlobalUnspecifiedEST->execute($blacklistAbsMax, $speciesId)
#      or die $blacklistUnspecifiedEST->errstr;
#    $blacklistUnspecifiedInSitu->execute($blacklistAbsMax, $speciesId)
#      or die $blacklistUnspecifiedInSitu->errstr;
#    $blacklistGlobalUnspecifiedInSitu->execute($blacklistAbsMax, $speciesId)
#      or die $blacklistUnspecifiedInSitu->errstr;
#    $blacklistUnspecifiedAffymetrix->execute($blacklistAbsMax, $speciesId)
#      or die $blacklistUnspecifiedAffymetrix->errstr;
#    $blacklistGlobalUnspecifiedAffymetrix->execute($blacklistAbsMax, $speciesId)
#      or die $blacklistUnspecifiedAffymetrix->errstr;
#    $blacklistUnspecifiedRnaSeq->execute($blacklistAbsMax, $speciesId)
#      or die $blacklistUnspecifiedRnaSeq->errstr;
#    $blacklistGlobalUnspecifiedRnaSeq->execute($blacklistAbsMax, $speciesId)
#      or die $blacklistUnspecifiedRnaSeq->errstr;
#    $blacklistUnspecifiedScRnaSeqFL->execute($blacklistAbsMax, $speciesId)
#      or die $blacklistUnspecifiedScRnaSeqFL->errstr;
#    $blacklistGlobalUnspecifiedScRnaSeqFL->execute($blacklistAbsMax, $speciesId)
#      or die $blacklistGlobalUnspecifiedScRnaSeqFL->errstr;
#    printf( "OK in %.2fs\n", ( time() - $t0 ) );

    $dbh->commit() or die('Failed commit');
}
$getAbsoluteMax->finish or die('Failed finish');
$updateExpression->finish or die('Failed finish');
$dbh->disconnect();

exit 0;

