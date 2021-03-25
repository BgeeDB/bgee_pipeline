#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Wollbrett, created March 2021

# USAGE: perl insert_scrna_seq_expression.pl -bgee=connection_string <OPTIONAL: -debug>
# After the insertion of full length single cell RNA-Seq data, this script inserts the data
# into the expression table and update the scRnaSeqFullLengthResult table.
# -debug: if provided, run in verbose mode (print the update/insert SQL queries, not executing them)

#############################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug)          = (0);
my %opts = ('bgee=s'     => \$bgee_connector,   # Bgee connector string
            'debug'      => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\t-debug            printing the update/insert SQL queries, not executing them
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

$| = 1;


##########################################
# GET CONDITIONS STUDIED WITH RNA-SEQ    #
##########################################
print "Retrieving conditions...\n";

my $queryConditions = $bgee->prepare('SELECT DISTINCT t2.exprMappedConditionId FROM scRnaSeqFullLengthLibrary AS t1
                                      INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId');
$queryConditions->execute()  or die $queryConditions->errstr;
my @exprMappedConditions = ();
while ( my @data = $queryConditions->fetchrow_array ){
    push(@exprMappedConditions, $data[0]);
}

print 'Done, ', scalar(@exprMappedConditions), " conditions retrieved.\n";


#################################################
# EXAMINE EXPRESSION CALLS ACROSS ALL SAMPLES   #
# TO KNOW WHICH GENES ARE NEVER SEEN AS PRESENT #
#################################################
print "Examining all gene results to update genes never seen as 'present'...\n";
my $findPresentGenes = $bgee->prepare('CREATE TEMPORARY TABLE tempPresentscRnaSeqFullLength (PRIMARY KEY (bgeeGeneId)) ENGINE=InnoDB
                                        AS (
                                            SELECT DISTINCT bgeeGeneId
                                            FROM scRnaSeqFullLengthResult WHERE detectionFlag = "'.$Utils::PRESENT_CALL.'"
                                        )');
if ( $debug ){
    print 'CREATE TEMPORARY TABLE tempPresentscRnaSeqFullLength (PRIMARY KEY (bgeeGeneId)) ENGINE=InnoDB
           AS (SELECT DISTINCT bgeeGeneId FROM scRnaSeqFullLengthResult WHERE detectionFlag = "'.$Utils::PRESENT_CALL.'")'."\n";
} else {
    $findPresentGenes->execute()  or die $findPresentGenes->errstr;
}

# At this point reasonForExclusion is always $Utils::CALL_NOT_EXCLUDED for $Utils::PRESENT_CALL and can be $Utils::CALL_NOT_EXCLUDED
# or $Utils::EXCLUDED_FOR_ABSENT_CALLS for $Utils::ABSENT_CALL.
my $preFilteringUp = $bgee->prepare('UPDATE scRnaSeqFullLengthResult AS t1
                                     LEFT OUTER JOIN tempPresentscRnaSeqFullLength AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId
                                     SET reasonForExclusion = "'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"
                                     WHERE t2.bgeeGeneId IS NULL AND reasonForExclusion != "'.$Utils::EXCLUDED_FOR_ABSENT_CALLS.'"');
if ( $debug ){
    print 'UPDATE scRnaSeqFullLengthResult AS t1
          LEFT OUTER JOIN tempPresentscRnaSeqFullLength AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId
          SET reasonForExclusion = "'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"
          WHERE t2.bgeeGeneId IS NULL AND reasonForExclusion != "'.$Utils::EXCLUDED_FOR_ABSENT_CALLS.'"'."\n";
} else {
    $preFilteringUp->execute()  or die $preFilteringUp->errstr;
}

my $dropTempPresent = $bgee->prepare('DROP TABLE tempPresentscRnaSeqFullLength');
if ( $debug ){
    print "DROP TABLE tempPresentscRnaSeqFullLength\n";
} else {
    $dropTempPresent->execute()  or die $dropTempPresent->errstr;
}

print "Done\n";

##########################################
# PREPARE QUERIES                        #
##########################################

# Insert/update expression
my $insUpExpr   = $bgee->prepare('INSERT INTO expression (bgeeGeneId, conditionId) VALUES (?, ?)
                                  ON DUPLICATE KEY UPDATE expressionId=LAST_INSERT_ID(expressionId)');


# Query to update scRnaSeqFullLengthResult with the expressionId and reasonForExclusion
my $updResult = $bgee->prepare('UPDATE scRnaSeqFullLengthResult AS t1
                                INNER JOIN scRnaSeqFullLengthLibrary AS t2 ON t1.scRnaSeqFullLengthLibraryId = t2.scRnaSeqFullLengthLibraryId
                                INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                SET expressionId = ?, reasonForExclusion = ?
                                WHERE bgeeGeneId = ? and t3.exprMappedConditionId = ?
                                AND t1.reasonForExclusion != "'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"
                                AND t1.pValue IS NOT NULL');

# query to get all the RNA-Seq results for a condition
my $queryResults = $bgee->prepare("SELECT t1.bgeeGeneId, t2.scRnaSeqFullLengthExperimentId, t1.scRnaSeqFullLengthLibraryId,
                                          t1.detectionFlag, t1.pValue, t1.reasonForExclusion
                                   FROM scRnaSeqFullLengthResult AS t1
                                   INNER JOIN scRnaSeqFullLengthLibrary AS t2 ON t1.scRnaSeqFullLengthLibraryId = t2.scRnaSeqFullLengthLibraryId
                                   INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                   WHERE t3.exprMappedConditionId = ? AND t1.reasonForExclusion != '".$Utils::EXCLUDED_FOR_PRE_FILTERED."'
                                   AND t1.pValue IS NOT NULL");


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

print "Processing conditions...\n";
for my $exprMappedConditionId ( @exprMappedConditions ){
    print "\tconditionId: $exprMappedConditionId\n";

    # retrieve genes results for this condition
    print "\t\tRetrieving genes results...\n";
    $queryResults->execute($exprMappedConditionId)  or die $queryResults->errstr;
    my %results = ();
    while ( my @data = $queryResults->fetchrow_array ){
        #data[0] = bgeeGeneId, data[1] = scRnaSeqFullLengthExperimentId, data[2] = scRnaSeqFullLengthLibraryId, data[3] = detectionFlag, data[4] = pvalue
        #XXX here it could be enough to create a hash like that : $results{$data[0]} = 1. However we kept 
        # information about experiments/libraries in case we have to insert them as a JSON field in the 
        # future. To update if pvalue is useless once full pipeline has been run for Bgee 15
        $results{$data[0]}->{$data[1]}->{$data[2]}->{'call'}                = $data[3];
        $results{$data[0]}->{$data[1]}->{$data[2]}->{'pvalue'}              = $data[4];
        $results{$data[0]}->{$data[1]}->{$data[2]}->{'reasonForExclusion'}  = $data[5];
    }
    print "\t\tDone, ", scalar(keys %results), ' genes retrieved. ',
          "Generating expression summary...\n";

    # now iterating the genes to insert expression data
    # (one row for a gene-condition)
    for my $geneId ( keys %results ){

        #check that at least 

        my $reasonForExclusion =  $Utils::CALL_NOT_EXCLUDED;
        my $expressionId = undef;

        # insert or update the expression table if not only undefined calls
        #TODO remove this condition if default exlusion type is still $Utils::CALL_NOT_EXCLUDED once bgee 15.0 is created
        if ( $debug ){
            print "INSERT INTO expression (bgeeGeneId, conditionId) VALUES ($geneId, $exprMappedConditionId)...\n";
        } else {
            $expressionId = add_expression($geneId, $exprMappedConditionId);
        }

        # Now update the related scRnaSeqFullLengthResults
        if ( $debug ){
            my $printExprId = 'notRetrievedForLog';
            print "UPDATE scRnaSeqFullLengthResult AS t1 ",
                  "INNER JOIN scRnaSeqFullLengthLibrary AS t2 ON t1.scRnaSeqFullLengthLibraryId = t2.scRnaSeqFullLengthLibraryId ",
                  "INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId ",
                  "SET expressionId=$printExprId, reasonForExclusion=$reasonForExclusion ",
                  "WHERE bgeeGeneId=$geneId and t3.exprMappedConditionId = $exprMappedConditionId ",
                  "and t1.reasonForExclusion != '".$Utils::EXCLUDED_FOR_PRE_FILTERED."'\n";
        } else {
            $updResult->execute($expressionId, $reasonForExclusion, $geneId, $exprMappedConditionId)
                or die $updResult->errstr;
        }
    }
}

$bgee->disconnect;
print "Done\n"  if ( $debug );

exit 0;

