#!/usr/bin/env perl
# Julien Roux, created 24/07/09
# modified 15/08/16, JR: revamped whole script for bgee v15 with condition Ids
# modified Nov 2016, Frederic Bastian: split former create_files_affy.pl script between insert_affy.pl and insert_affy_expression.pl.
# last modified March 2021 for bgee 15.0

# Perl core modules
use strict;
use warnings;
use diagnostics;
use File::Slurp;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

$| = 1; # no buffering of output


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


##########################################
# GET CONDITIONS STUDIED WITH AFFYMETRIX #
##########################################
print "Retrieving conditions...\n";

my $queryConditions = $bgee->prepare('SELECT DISTINCT t2.exprMappedConditionId FROM affymetrixChip AS t1
                                      INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId');
$queryConditions->execute()  or die $queryConditions->errstr;
my @exprMappedConditions = ();
while ( my @data = $queryConditions->fetchrow_array ){
    push(@exprMappedConditions, $data[0]);
}

print "Done, ", scalar(@exprMappedConditions), " conditions retrieved.\n";

#####################################################
# EXAMINE EXPRESSION CALLS ACROSS ALL SAMPLES       #
# TO KNOW WHICH PROBESETS ARE NEVER SEEN AS PRESENT #
#####################################################
print "Examining all results to update probesets never seen as 'present'...\n";

if($debug) { 
    print "CREATE TEMPORARY TABLE tempPresentAffy (PRIMARY KEY (affymetrixProbesetId, chipTypeId)) ENGINE=InnoDB ",
          "AS (",
          "SELECT DISTINCT t1.affymetrixProbesetId, t2.chipTypeId ",
          "FROM affymetrixProbeset AS t1 ",
          "INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId ",
          'WHERE t1.rawDetectionFlag IN ("'.$Utils::PRESENT_CALL.'", "'.$Utils::MARGINAL_CALL.'"))'."\n";
} else {
    my $findPresentGenes = $bgee->prepare('CREATE TEMPORARY TABLE tempPresentAffy (PRIMARY KEY (affymetrixProbesetId, chipTypeId)) ENGINE=InnoDB
                    AS (
                        SELECT DISTINCT t1.affymetrixProbesetId, t2.chipTypeId
                        FROM affymetrixProbeset AS t1
                        INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                        WHERE t1.rawDetectionFlag IN ("'.$Utils::PRESENT_CALL.'", "'.$Utils::MARGINAL_CALL.'")
                    )');
    $findPresentGenes->execute()  or die $findPresentGenes->errstr;
}

# The "not excluded" status has to be set first, for the next query to properly take into account "undefined" status
if($debug) { 
    print 'UPDATE affymetrixProbeset AS t1',
          'INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId ',
          'INNER JOIN tempPresentAffy AS t3 ',
          'ON t1.affymetrixProbesetId = t3.affymetrixProbesetId AND t2.chipTypeId = t3.chipTypeId ',
          'SET reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"'."\n";
} else {
    my $presentUp = $bgee->prepare('UPDATE affymetrixProbeset AS t1
                                INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                                INNER JOIN tempPresentAffy AS t3
                                ON t1.affymetrixProbesetId = t3.affymetrixProbesetId AND t2.chipTypeId = t3.chipTypeId
                                SET reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"');
    $presentUp->execute()  or die $presentUp->errstr;
}
if($debug) {
    print 'UPDATE affymetrixProbeset AS t1 ',
          'INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId ',
          'LEFT OUTER JOIN tempPresentAffy AS t3 ',
          'ON t1.affymetrixProbesetId = t3.affymetrixProbesetId AND t2.chipTypeId = t3.chipTypeId ',
          'SET reasonForExclusion = "'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'" ',
          'WHERE t3.affymetrixProbesetId IS NULL ',
          'AND t1.reasonForExclusion != "'.$Utils::EXCLUDED_FOR_UNDEFINED.'"'."\n";
} else {
    my $preFilteringUp = $bgee->prepare('UPDATE affymetrixProbeset AS t1
                                     INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                                     LEFT OUTER JOIN tempPresentAffy AS t3
                                     ON t1.affymetrixProbesetId = t3.affymetrixProbesetId AND t2.chipTypeId = t3.chipTypeId
                                     SET reasonForExclusion = "'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"
                                     WHERE t3.affymetrixProbesetId IS NULL AND t1.reasonForExclusion != "'.$Utils::EXCLUDED_FOR_UNDEFINED.'"');
    $preFilteringUp->execute()  or die $preFilteringUp->errstr;
}
if ($debug) {
    print "DROP TABLE tempPresentAffy\n";
} else {
    my $dropTempPresent = $bgee->prepare('DROP TABLE tempPresentAffy');
    $dropTempPresent->execute()  or die $dropTempPresent->errstr;
}
print "Done\n";

##########################################
# PREPARE QUERIES                        #
##########################################

# insert/update expression
my $insUpExpr   = $bgee->prepare('INSERT INTO expression (bgeeGeneId, conditionId) VALUES (?, ?)
                                  ON DUPLICATE KEY UPDATE expressionId=LAST_INSERT_ID(expressionId)');

# Query to update affymetrixProbeset with the expressionId and reasonForExclusion
my $updResult = $bgee->prepare('UPDATE affymetrixProbeset AS t1
                                INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                                INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                SET expressionId = ?, reasonForExclusion = ?
                                WHERE t1.bgeeGeneId = ? AND t3.exprMappedConditionId = ?
                                    AND t1.reasonForExclusion != "'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"');

# query to get all the Affymetrix probesets for a condition
my $queryResults = $bgee->prepare("SELECT t1.bgeeGeneId, t2.microarrayExperimentId, t1.bgeeAffymetrixChipId,
                                          t1.affymetrixProbesetId, t1.rawDetectionFlag, t1.pValue
                                   FROM affymetrixProbeset AS t1
                                   INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId
                                   INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
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

print "Processing conditions...\n";
for my $exprMappedConditionId ( @exprMappedConditions ){
    print "\tconditionId: $exprMappedConditionId\n";

    # retrieve genes results for this condition
    print "\t\tRetrieving genes results...\n";
    $queryResults->execute($exprMappedConditionId)  or die $queryResults->errstr;
    my %results = ();
    while ( my @data = $queryResults->fetchrow_array ){
        #data[0] = bgeeGeneId, data[1] = experimentId, data[2] = chipId, data[3] = probesetId, data[4] = rawDetectionFlag, data[5] = pValue
        #important to retrieve the probeset IDs for first per-chip reconciliation
        $results{$data[0]}->{$data[1]}->{$data[2]}->{$data[3]}->{'call'}    = $data[4];
        $results{$data[0]}->{$data[1]}->{$data[2]}->{$data[3]}->{'pvalue'} = $data[5];
    }
    print "\t\tDone, ", scalar(keys %results), ' genes retrieved. ',
          "Generating expression summary...\n";

    # now iterating the genes to insert expression data
    # (one row for a gene-condition)
    for my $geneId ( keys %results ){


        my $reasonForExclusion = $Utils::CALL_NOT_EXCLUDED;
        my $expressionId   = undef;

        # insert or update the expression table if not only undefined calls
        if ( $reasonForExclusion eq $Utils::CALL_NOT_EXCLUDED ){
            if ( $debug ){
                print "INSERT INTO expression (bgeeGeneId, conditionId) VALUES ($geneId, $exprMappedConditionId)...\n";
            } else {
                $expressionId = add_expression($geneId, $exprMappedConditionId);
            }
        } else {
            warn "No calls available for geneId $geneId in exprMappedConditionId $exprMappedConditionId\n";
        }

        # Now update the related affymetrixProbesets
        if ( $debug ){
            my $printExprId = 'notRetrievedForLog';
            print "UPDATE affymetrixProbeset AS t1 ",
                  "INNER JOIN affymetrixChip AS t2 ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId ",
                  "INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId ",
                  "SET expressionId = $printExprId, reasonForExclusion = '$reasonForExclusion' ",
                  "WHERE t1.bgeeGeneId = $geneId AND t3.exprMappedConditionId = $exprMappedConditionId ",
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
