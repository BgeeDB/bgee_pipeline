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
if ( !$test_options || $bgee_connector eq '' || $bgee_species eq '' || $condition_ids eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -bgee_species='-'
\t-bgee            Bgee connector string
\t-species         Bgee species for which ranks have to be normalized or '-' for all species
\t-conditions      Subset of Bgee condition IDs for which ranks have to be normalized or '-'
                   for all conditions. Subset of condition IDs can only be provided when exactly one species is selected
\n";
    exit 1;
}

#Sanity check that only condition provided as argument of the script AND actually
#corresponding to the expected species will be processed
sub retrieveConditionFromDatabase {
    my ($conditionIdsRef, $speciesId, $dbh) = @_;
    my @filteredCondition = ();
    my $conditionSql = "SELECT globalConditionId FROM globalCond WHERE ";
    $conditionSql .= 'globalConditionId IN (';
    for my $i (0 .. $#$conditionIdsRef) {
        if ($i > 0) {
            $conditionSql .= ', ';
        }
        $conditionSql .= $$conditionIdsRef[$i];
    }
    $conditionSql .= ') AND ';
    $conditionSql .= "speciesId = $speciesId";

    my $getCondition = $dbh->prepare($conditionSql);
    $getCondition->execute()  or die $getCondition->errstr;

    while ( my @data = $getCondition->fetchrow_array ){
        push(@filteredCondition, $data[0]);
    }
    return @filteredCondition;
}

#Do not filter on species in order not to lock the gene table when updating globalExpression
#In order to avoid issues we previously verified that all conditionIds come from the expected species
sub createUdpateRanksQuery {
    my ($conditionListRef, $speciesId) = @_;
    # update the expression table with normalized mean ranks
    my $updateExpressionQuery = "UPDATE globalExpression
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
        scRnaSeqFullLengthGlobalMeanRankNorm = (scRnaSeqFullLengthGlobalMeanRank + (scRnaSeqFullLengthGlobalMeanRank * ? / scRnaSeqFullLengthGlobalMaxRank))/2";

    #add as many question marks as number of globalCondition IDs
    if (@$conditionListRef) {
        $updateExpressionQuery .= " WHERE globalCond.globalConditionId IN (";
        foreach (@$conditionListRef[1 .. $#$conditionListRef]) {
            $updateExpressionQuery .= "?,";
        }
        $updateExpressionQuery .= "?)";
    } else {
        $updateExpressionQuery .= " WHERE globalCond.speciesId = ?";
    }
    return $updateExpressionQuery;
}

my $dbh = Utils::connect_bgee_db($bgee_connector);

# if no speciesId provided then query the database to retrieve all speciesIds
# otherwise normalize ranks only for provided speciesIds
my @speciesList = ();
if ($bgee_species ne $emptyArg) {
    @speciesList = split(',', $bgee_species);
} else {
    my $speciesIdsSql = "select distinct speciesId from species order by speciesId";
    my $getSpeciesIds = $dbh->prepare($speciesIdsSql);
    $getSpeciesIds->execute()  or die $getSpeciesIds->errstr;
    while ( my @data = $getSpeciesIds->fetchrow_array ){
        push(@speciesList, $data[0]);
    }
}

my @conditionList = ();
if ($condition_ids ne $emptyArg) {
    @conditionList = split(',', $condition_ids);
}

if(scalar @speciesList != 1 && $condition_ids ne $emptyArg) {
    print "Condition IDs can be provided only when exactly one species is selected";
    exit 1;
}

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
$getAbsoluteMax->finish or die('Failed finish');

for my $speciesId ( keys %maxRanks ){
    Utils::start_transaction($dbh);
    #if condition are provided we first verify that they correspond to the expected species by
    #querying the database
    my @verifiedConditionList = ();
    if (@conditionList) {
        @verifiedConditionList = retrieveConditionFromDatabase(\@conditionList, $speciesId, $dbh);
    }
    my $absMax = $maxRanks{$speciesId};
    print "Updating ranks for species $speciesId with max rank $absMax\n";

    my $updateExpression = $dbh->prepare(createUdpateRanksQuery(\@verifiedConditionList, $speciesId));
    my $t0 = time();
    printf('Update expression table with normalized mean ranks per type:');
    if (@verifiedConditionList) {
	$updateExpression->execute( $absMax, $absMax, $absMax, $absMax, $absMax, $absMax, $absMax, $absMax,
	        $absMax, $absMax, @verifiedConditionList) or die $updateExpression->errstr;
    } else {
	$updateExpression->execute( $absMax, $absMax, $absMax, $absMax, $absMax, $absMax, $absMax, $absMax,
	        $absMax, $absMax, $speciesId) or die $updateExpression->errstr;
    }
    $updateExpression->finish or die('Failed finish');
    printf( "OK in %.2fs\n", ( time() - $t0 ) );

    $dbh->commit() or die('Failed commit');
}
$dbh->disconnect();

exit 0;
