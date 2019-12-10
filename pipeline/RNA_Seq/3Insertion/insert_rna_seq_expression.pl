#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Frederic Bastian, created November 2012, last update Dec. 2016

# USAGE: perl insert_rna_seq_expression.pl -bgee=connection_string <OPTIONAL: -debug>
# After the insertion of RNA-Seq data, this script inserts the data
# into the expression table, rnaSeqEvidenceCount table, and update the rnaSeqResult table.
# -debug: if provided, run in verbose mode (print the update/insert SQL queries, not executing them)

#############################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
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

my $queryConditions = $bgee->prepare('SELECT DISTINCT t2.exprMappedConditionId FROM rnaSeqLibrary AS t1
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
my $findPresentGenes = $bgee->prepare('CREATE TEMPORARY TABLE tempPresentRnaSeq (PRIMARY KEY (bgeeGeneId)) ENGINE=InnoDB 
                                        AS (
                                            SELECT DISTINCT bgeeGeneId
                                            FROM rnaSeqResult WHERE detectionFlag = "'.$Utils::PRESENT_CALL.'"
                                        )');
$findPresentGenes->execute()  or die $findPresentGenes->errstr;

# The "not excluded" status has to be set first, for the next query to properly take into account "undefined" status
my $presentUp = $bgee->prepare('UPDATE rnaSeqResult AS t1
                                INNER JOIN tempPresentRnaSeq AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId
                                SET reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"');
$presentUp->execute()  or die $presentUp->errstr;
my $preFilteringUp = $bgee->prepare('UPDATE rnaSeqResult AS t1
                                     LEFT OUTER JOIN tempPresentRnaSeq AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId
                                     SET reasonForExclusion = "'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"
                                     WHERE t2.bgeeGeneId IS NULL AND reasonForExclusion != "'.$Utils::EXCLUDED_FOR_UNDEFINED.'"');
$preFilteringUp->execute()  or die $preFilteringUp->errstr;

my $dropTempPresent = $bgee->prepare('DROP TABLE tempPresentRnaSeq');
$dropTempPresent->execute()  or die $dropTempPresent->errstr;

print "Done\n";

##########################################
# PREPARE QUERIES                        #
##########################################

# Insert/update expression
my $insUpExpr   = $bgee->prepare('INSERT INTO expression (bgeeGeneId, conditionId) VALUES (?, ?)
                                  ON DUPLICATE KEY UPDATE expressionId=LAST_INSERT_ID(expressionId)');


# Query to update rnaSeqResult with the expressionId and reasonForExclusion
my $updResult = $bgee->prepare('UPDATE rnaSeqResult AS t1
                                INNER JOIN rnaSeqLibrary AS t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId
                                INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                SET expressionId = ?, reasonForExclusion = ?
                                WHERE bgeeGeneId = ? and t3.exprMappedConditionId = ?
                                    AND t1.reasonForExclusion != "'.$Utils::EXCLUDED_FOR_PRE_FILTERED.'"');

# Insertion into rnaSeqExperimentExpression
my $insExpSummary = $bgee->prepare('INSERT INTO rnaSeqExperimentExpression
                                    (expressionId, rnaSeqExperimentId,
                                        presentHighRNASeqLibraryCount, presentLowRNASeqLibraryCount,
                                        absentHighRNASeqLibraryCount, absentLowRNASeqLibraryCount,
                                        rnaSeqExperimentCallDirection, rnaSeqExperimentCallQuality)
                                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)');

# query to get all the RNA-Seq results for a condition
my $queryResults = $bgee->prepare("SELECT t1.bgeeGeneId, t2.rnaSeqExperimentId, t1.rnaSeqLibraryId,
                                          t1.detectionFlag, t1.rnaSeqData
                                   FROM rnaSeqResult AS t1
                                   INNER JOIN rnaSeqLibrary AS t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId
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
        #data[0] = bgeeGeneId, data[1] = rnaSeqExperimentId, data[2] = rnaSeqLibraryId, data[3] = detectionFlag, data[4] = rnaSeqData
        $results{$data[0]}->{$data[1]}->{$data[2]}->{'call'}    = $data[3];
        $results{$data[0]}->{$data[1]}->{$data[2]}->{'quality'} = $data[4];
    }
    print "\t\tDone, ", scalar(keys %results), ' genes retrieved. ',
          "Generating expression summary...\n";

    # now iterating the genes to insert expression data
    # (one row for a gene-condition)
    for my $geneId ( keys %results ){
        # get the summary of expression and quality from the list of rnaSeqLibraryIds
        # where this gene is analyzed (for this condition)
        my ($reasonForExclusion, $summary) = Utils::summarizeExperimentCallAndQuality(\%{$results{$geneId}});

        my $expressionId = undef;

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

        # Now update the related rnaSeqResults
        if ( $debug ){
            my $printExprId = 'notRetrievedForLog';
            print "UPDATE rnaSeqResult AS t1 ",
                  "INNER JOIN rnaSeqLibrary AS t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId ",
                  "INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId ",
                  "SET expressionId=$printExprId, reasonForExclusion=$reasonForExclusion ",
                  "WHERE bgeeGeneId=$geneId and t3.exprMappedConditionId = $exprMappedConditionId ",
                  "and t1.reasonForExclusion != '".$Utils::EXCLUDED_FOR_PRE_FILTERED."'\n";
            if ( $reasonForExclusion eq $Utils::CALL_NOT_EXCLUDED ) {
                for my $expId ( keys %{ $summary } ) {
                    print "INSERT INTO rnaSeqExperimentExpression (expressionId, rnaSeqExperimentId,
                               presentHighRNASeqLibraryCount, presentLowRNASeqLibraryCount,
                               absentHighRNASeqLibraryCount, absentLowRNASeqLibraryCount,
                               rnaSeqExperimentCallDirection, rnaSeqExperimentCallQuality)
                           VALUES ($printExprId, $expId,
                               $summary->{$expId}->{'pstHighEvidenceCount'},
                               $summary->{$expId}->{'pstLowEvidenceCount'},
                               $summary->{$expId}->{'absHighEvidenceCount'},
                               $summary->{$expId}->{'absLowEvidenceCount'},
                               $summary->{$expId}->{'expCall'}, $summary->{$expId}->{'expCallQuality'})\n";
                }
            }
        } else {
            $updResult->execute($expressionId, $reasonForExclusion, $geneId, $exprMappedConditionId)
                or die $updResult->errstr;
            if ( $reasonForExclusion eq $Utils::CALL_NOT_EXCLUDED ){
                for my $expId ( keys %{ $summary } ){
                    $insExpSummary->execute($expressionId, $expId,
                        $summary->{$expId}->{'pstHighEvidenceCount'},
                        $summary->{$expId}->{'pstLowEvidenceCount'},
                        $summary->{$expId}->{'absHighEvidenceCount'},
                        $summary->{$expId}->{'absLowEvidenceCount'},
                        $summary->{$expId}->{'expCall'}, $summary->{$expId}->{'expCallQuality'})
                    or die $insExpSummary->errstr;
                }
            }
        }
    }
}

$bgee->disconnect;
print "Done\n"  if ( $debug );

exit 0;

