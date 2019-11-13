#!/usr/bin/env perl

# Philippe Moret, created Nov 2015.
# Frederic Bastian, last updated June 2016.
# Normalize the mean rank in the expression table across the different data types.
# Frederic Bastian, last updated Feb. 2017: adapt to new conditions and new schema in Bgee 14

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
my ($bgee_connector) = ('');
my %opts = ( 'bgee=s' => \$bgee_connector );    # Bgee connector string

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD)
\t-bgee      Bgee    connector string
\n";
    exit 1;
}

my $blacklisted = "('XAO:0003003', 'ZFA:0001093')";
my $auto = 1;

    my $dbh = Utils::connect_bgee_db($bgee_connector);
    if ( $auto == 0 ) {
        $dbh->{'AutoCommit'} = 0;
    }

    #we get the absolute max rank across all conditions for each species
    my $getAbsoluteMax = $dbh->prepare(
    "SELECT speciesId, GREATEST(IFNULL(max(rnaSeqMaxRank), 0.0),
                 IFNULL(max(estMaxRank),0.0),
                 IFNULL(max(inSituMaxRank),0.0),
                 IFNULL(max(affymetrixMaxRank),0.0),
                 IFNULL(max(rnaSeqGlobalMaxRank), 0.0),
                 IFNULL(max(estGlobalMaxRank),0.0),
                 IFNULL(max(inSituGlobalMaxRank),0.0),
                 IFNULL(max(affymetrixGlobalMaxRank),0.0)) AS max
    FROM globalCond
    WHERE anatEntityId NOT IN $blacklisted
    GROUP BY speciesId" );

    $getAbsoluteMax->execute()  or die $getAbsoluteMax->errstr;
    my %maxRanks = ();
    while ( my @data = $getAbsoluteMax->fetchrow_array ){
        $maxRanks{$data[0]} = $data[1];
        print("Max rank for species $data[0]: $data[1]\n");
    }

    # update the expression table with normalized mean ranks
    my $updateExpression = $dbh->prepare( "
    UPDATE gene
    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
    STRAIGHT_JOIN globalCond ON globalCond.globalConditionId = globalExpression.globalConditionId
    SET rnaSeqMeanRankNorm     = (rnaSeqMeanRank + (rnaSeqMeanRank * ? / rnaSeqMaxRank))/2,
        estRankNorm            = (estRank + (estRank * ? / estMaxRank))/2,
        inSituRankNorm         = (inSituRank + (inSituRank * ? / inSituMaxRank))/2,
        affymetrixMeanRankNorm = (affymetrixMeanRank + (affymetrixMeanRank * ? / affymetrixMaxRank))/2,
        rnaSeqGlobalMeanRankNorm     = (rnaSeqGlobalMeanRank + (rnaSeqGlobalMeanRank * ? / rnaSeqGlobalMaxRank))/2,
        estGlobalRankNorm            = (estGlobalRank + (estGlobalRank * ? / estGlobalMaxRank))/2,
        inSituGlobalRankNorm         = (inSituGlobalRank + (inSituGlobalRank * ? / inSituGlobalMaxRank))/2,
        affymetrixGlobalMeanRankNorm = (affymetrixGlobalMeanRank + (affymetrixGlobalMeanRank * ? / affymetrixGlobalMaxRank))/2
    WHERE gene.speciesId = ?" );

    # Quick and dirty blacklisting of unspecified anatEntities: set the mean rank to the max value
    my $blacklistUnspecifiedEST = $dbh->prepare( "
    UPDATE gene
    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
        AND globalCond.globalConditionId = globalExpression.globalConditionId
    SET estRankNorm = ?, estGlobalRankNorm = ?
    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND estRank IS NOT NULL" );

    my $blacklistUnspecifiedInSitu = $dbh->prepare( "
    UPDATE gene
    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
        AND globalCond.globalConditionId = globalExpression.globalConditionId
    SET inSituRankNorm = ?, inSituGlobalRankNorm = ?
    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND inSituRank IS NOT NULL" );

    my $blacklistUnspecifiedAffymetrix = $dbh->prepare( "
    UPDATE gene
    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
        AND globalCond.globalConditionId = globalExpression.globalConditionId
    SET affymetrixMeanRankNorm = ?, affymetrixGlobalMeanRankNorm = ?
    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND affymetrixMeanRank IS NOT NULL" );

    my $blacklistUnspecifiedRnaSeq = $dbh->prepare( "
    UPDATE gene
    STRAIGHT_JOIN globalCond ON gene.speciesId = globalCond.speciesId
    STRAIGHT_JOIN globalExpression ON gene.bgeeGeneId = globalExpression.bgeeGeneId
        AND globalCond.globalConditionId = globalExpression.globalConditionId
    SET rnaSeqMeanRankNorm = ?, rnaSeqGlobalMeanRankNorm = ?
    WHERE gene.speciesId = ? AND globalCond.anatEntityId IN $blacklisted AND rnaSeqMeanRank IS NOT NULL" );

for my $speciesId ( keys %maxRanks ){
    my $absMax = $maxRanks{$speciesId};
    print "Updating ranks for species $speciesId with max rank $absMax\n";

    my $t0 = time();
    printf('Update expression table with normalized mean ranks per type:   ');
    $updateExpression->execute( $absMax, $absMax, $absMax, $absMax, $absMax, $absMax, $absMax, $absMax,
                                $speciesId)
      or die $updateExpression->errstr;
    printf( "OK in %.2fs\n", ( time() - $t0 ) );

    $t0 = time();
    printf('Blacklisting unspecified anat entities:    ');
    my $blacklistAbsMax = $absMax + 1;
    $blacklistUnspecifiedEST->execute($blacklistAbsMax, $blacklistAbsMax, $speciesId)
      or die $blacklistUnspecifiedEST->errstr;
    $blacklistUnspecifiedInSitu->execute($blacklistAbsMax, $blacklistAbsMax, $speciesId)
      or die $blacklistUnspecifiedInSitu->errstr;
    $blacklistUnspecifiedAffymetrix->execute($blacklistAbsMax, $blacklistAbsMax, $speciesId)
      or die $blacklistUnspecifiedAffymetrix->errstr;
    $blacklistUnspecifiedRnaSeq->execute($blacklistAbsMax, $blacklistAbsMax, $speciesId)
      or die $blacklistUnspecifiedRnaSeq->errstr;
    printf( "OK in %.2fs\n", ( time() - $t0 ) );

    if ($auto == 0) {
        $dbh->commit()  or die('Failed commit');
    }
}
$dbh->disconnect();

exit 0;

