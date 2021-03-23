#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Roux, created 16/01/09
# USAGE: perl insert_expression_in_situ.pl expr/no_expr/both
# insert in situ data in expression table (ZFIN, Xenbase and MGD)
# 
# modified Jan. 2017, FBB: use of new quality computation and new inSituExperimentExpression table
# modified Mar. 2021, JW:  do not anymore count present/absent low/high calls. As for Bgee 15 pvalues are used to summarize quality
#################################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

$| = 1;

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug) = (0);
my %opts = ('debug'      => \$debug,            # printing the update/insert SQL queries, not executing them
            'bgee=s'     => \$bgee_connector   # Bgee connector string
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD)
\t-bgee      Bgee    connector string
\t-debug     printing the update/insert SQL queries, not executing them
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

##########################################
# GET CONDITIONS STUDIED WITH IN SITU    #
##########################################
print "Retrieving conditions...\n";

my $queryConditions = $bgee->prepare('SELECT DISTINCT t2.exprMappedConditionId FROM inSituSpot AS t1
                                      INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId');
$queryConditions->execute()  or die $queryConditions->errstr;
my @exprMappedConditions = ();
while ( my @data = $queryConditions->fetchrow_array ){
    push(@exprMappedConditions, $data[0]);
}

print "Done, ", scalar(@exprMappedConditions), " conditions retrieved.\n";

############################
# ADD PVALUES              #
############################

my $PVALUE_ABSENT_HIGH = 0.5;
my $PVALUE_ABSENT_LOW = 0.1;
my $PVALUE_PRESENT_LOW = 0.01;
my $PVALUE_PRESENT_HIGH = 0.0004;

if ($debug) {
    print 'UPDATE inSituSpot SET pValue = CASE '.
                          ' WHEN inSituData = "'.$Utils::LOW_QUAL.'" AND detectionFlag = "'.$Utils::ABSENT_CALL.'" THEN '.$PVALUE_ABSENT_LOW.
                          ' WHEN inSituData = "'.$Utils::HIGH_QUAL.'" AND detectionFlag = "'.$Utils::ABSENT_CALL.'" THEN '.$PVALUE_ABSENT_HIGH.
                          ' WHEN inSituData = "'.$Utils::LOW_QUAL.'" AND detectionFlag = "'.$Utils::PRESENT_CALL.'" THEN '.$PVALUE_PRESENT_LOW.
                          ' WHEN inSituData = "'.$Utils::HIGH_QUAL.'" AND detectionFlag = "'.$Utils::PRESENT_CALL.'" THEN '.$PVALUE_PRESENT_HIGH.
                          ' ELSE NULL '.
                        ' END '."\n";
} else {
    my $updPValues   = $bgee->prepare('UPDATE inSituSpot SET pValue = CASE '.
                          ' WHEN inSituData = "'.$Utils::LOW_QUAL.'" AND detectionFlag = "'.$Utils::ABSENT_CALL.'" THEN '.$PVALUE_ABSENT_LOW.
                          ' WHEN inSituData = "'.$Utils::HIGH_QUAL.'" AND detectionFlag = "'.$Utils::ABSENT_CALL.'" THEN '.$PVALUE_ABSENT_HIGH.
                          ' WHEN inSituData = "'.$Utils::LOW_QUAL.'" AND detectionFlag = "'.$Utils::PRESENT_CALL.'" THEN '.$PVALUE_PRESENT_LOW.
                          ' WHEN inSituData = "'.$Utils::HIGH_QUAL.'" AND detectionFlag = "'.$Utils::PRESENT_CALL.'" THEN '.$PVALUE_PRESENT_HIGH.
                          ' ELSE NULL '.
                        ' END ');
    $updPValues->execute()  or die $updPValues->errstr;
    $updPValues->finish();
}



##########################################
# PREPARE QUERIES                        #
##########################################

# insert/update expression
my $insUpExpr   = $bgee->prepare('INSERT INTO expression (bgeeGeneId, conditionId) VALUES (?, ?)
                                  ON DUPLICATE KEY UPDATE expressionId=LAST_INSERT_ID(expressionId)');

# Query to update inSituSpot with the expressionId and reasonForExclusion
my $updResult = $bgee->prepare('UPDATE inSituSpot AS t1
                                INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId
                                SET expressionId = ?, reasonForExclusion = ?
                                WHERE t1.bgeeGeneId = ? AND t2.exprMappedConditionId = ?
                                    AND t1.reasonForExclusion != "'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"');

# query to get all the in situ spots for a condition
# (inSituSpotIds are unique over all experiments so we don't need to also retrieve the evidence IDs)
my $queryResults = $bgee->prepare("SELECT t1.bgeeGeneId, t2.inSituExperimentId, t1.inSituSpotId,
                                          t1.detectionFlag, t1.inSituData
                                   FROM inSituSpot AS t1
                                   INNER JOIN inSituEvidence AS t2 ON t1.inSituEvidenceId = t2.inSituEvidenceId
                                   INNER JOIN cond AS t3 ON t1.conditionId = t3.conditionId
                                   WHERE t3.exprMappedConditionId = ? AND t1.reasonForExclusion = '".$Utils::CALL_NOT_EXCLUDED."'");

##########################################
# SUBROUTINES                            #
##########################################
# subroutine to insert/update in expression table
# return the relevant expressionId
sub add_expression {
    my ($bgeeGeneId, $exprMappedConditionId) = @_;

    $insUpExpr->execute($bgeeGeneId, $exprMappedConditionId)  or die $insUpExpr->errstr;
    return $bgee->{'mysql_insertid'}; # expr_id
}

##########################################
# ITERATING CONDITIONS TO INSERT DATA    #
##########################################

my $reasonForExclusion = $Utils::CALL_NOT_EXCLUDED;

print "Processing conditions...\n";
for my $exprMappedConditionId ( @exprMappedConditions ){
    print "\tconditionId: $exprMappedConditionId\n";

    # retrieve genes results for this condition
    print "\t\tRetrieving genes results...\n";
    $queryResults->execute($exprMappedConditionId)  or die $queryResults->errstr;
    my %results = ();
    while ( my @data = $queryResults->fetchrow_array ){
        #data[0] = bgeeGeneId, data[1] = experimentId, data[2] = spotId, data[3] = detectionFlag, data[4] = inSituData
        # (inSituSpotIds are unique over all experiments so we don't need to also retrieve the evidence IDs)
        $results{$data[0]}->{$data[1]}->{$data[2]}->{'call'}    = $data[3];
        $results{$data[0]}->{$data[1]}->{$data[2]}->{'quality'} = $data[4];
    }
    print "\t\tDone, ", scalar(keys %results), ' genes retrieved. ',
          "Generating expression summary...\n";

    # now iterating the genes to insert expression data
    # (one row for a gene-condition)
    for my $geneId ( keys %results ){ 
        
        my $expressionId   = undef;

        if ( $debug ){
            print "INSERT INTO expression (bgeeGeneId, conditionId) VALUES ($geneId, $exprMappedConditionId)...\n";
        } else {
            $expressionId = add_expression($geneId, $exprMappedConditionId);
        }

        # Now update the related inSituSpots
        if ( $debug ){
            my $printExprId = 'notRetrievedForLog';
            print "UPDATE inSituSpot AS t1
                   INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId
                   SET expressionId = $printExprId, reasonForExclusion = $reasonForExclusion
                   WHERE t1.bgeeGeneId = $geneId AND t2.exprMappedConditionId = $exprMappedConditionId
                   AND t1.reasonForExclusion != '".$Utils::EXCLUDED_FOR_PRE_FILTERED."'\n";
        } else {
            $updResult->execute($expressionId, $reasonForExclusion, $geneId, $exprMappedConditionId)
                or die $updResult->errstr;
        }
    }
}
$insUpExpr->finish();
$updResult->finish();
$queryResults->finish();
$bgee->disconnect;
print "Done\n"  if ( $debug );

exit 0;
