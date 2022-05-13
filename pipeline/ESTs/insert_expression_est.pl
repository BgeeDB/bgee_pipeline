#!/usr/bin/env perl

# last update Feb. 2017 Frederic Bastian: adapt to new call computations in Bgee 14

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

$| = 1; # stdout not out in memory buffer

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug)          = (0);
my %opts = ('bgee=s' => \$bgee_connector,     # Bgee connector string
            'debug'  => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=user=username__pass=mypass__host=127.0.0.1__port=3306__name=bgee_v15
\t-bgee      Bgee connector string
\t-debug     printing the update/insert SQL queries, not executing them
\n";
    exit 1;
}


# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);


##########################################
# GET CONDITIONS STUDIED WITH EST        #
##########################################
print "Retrieving conditions...\n";

my $queryConditions = $bgee->prepare('SELECT DISTINCT t2.exprMappedConditionId FROM estLibrary AS t1
                                      INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId');
$queryConditions->execute()  or die $queryConditions->errstr;
my @exprMappedConditions = ();
while ( my @data = $queryConditions->fetchrow_array ){
    push(@exprMappedConditions, $data[0]);
}

print 'Done, ', scalar(@exprMappedConditions), " conditions retrieved.\n";

##########################################
# PREPARE QUERIES                        #
##########################################

# insert/update expression
my $insUpExpr   = $bgee->prepare('INSERT INTO expression (bgeeGeneId, conditionId) VALUES (?, ?)
                                  ON DUPLICATE KEY UPDATE expressionId=LAST_INSERT_ID(expressionId)');


# Query to update expressedSequenceTag with the expressionId
my $updResult = $bgee->prepare('UPDATE expressedSequenceTag AS t1
                                INNER JOIN estLibrary AS t2 ON t1.estLibraryId = t2.estLibraryId
                                INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                SET expressionId = ?
                                WHERE bgeeGeneId = ? and t3.exprMappedConditionId = ?');

# Insertion into estLibraryExpression
my $insExpSummary = $bgee->prepare('INSERT INTO estLibraryExpression
                                    (expressionId, estLibraryId, estCount)
                                    VALUES (?, ?, ?)');

# Add p-value in expressedSequenceTag
my $updPVal = $bgee->prepare('UPDATE expressedSequenceTag SET pValue=? WHERE expressionId=? AND estLibraryId=?');

# query to get all EST results for a condition
my $queryResults = $bgee->prepare("SELECT t1.bgeeGeneId, t1.estLibraryId, count(distinct estId)
                                   FROM expressedSequenceTag AS t1
                                   INNER JOIN estLibrary AS t2 ON t1.estLibraryId = t2.estLibraryId
                                   INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                                   WHERE t3.exprMappedConditionId = ?
                                   GROUP BY t1.estLibraryId, t1.bgeeGeneId");

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
        #data[0] = bgeeGeneId, data[1] = estLibraryId, data[2] = EST count
        $results{$data[0]}->{$data[1]} = $data[2];
    }
    print "\t\tDone, ", scalar(keys %results), ' genes retrieved. ',
          "Generating expression summary...\n";

    # now iterating the genes to insert expression data
    # (one row for a gene-condition)
    for my $geneId ( keys %results ){
        # Because there is no EST "experiments" (EST libraries are considered independent),
        # and because EST data can only produced "present" expression calls,
        # we don't use the method Utils::summarizeExperimentCallAndQuality(\%{$results{$geneId}}

        my $expressionId = undef;

        # insert or update the expression table
        if ( $debug ){
            print "INSERT INTO expression (bgeeGeneId, conditionId) VALUES ($geneId, $exprMappedConditionId)...\n";
        } else {
            $expressionId = add_expression($geneId, $exprMappedConditionId);
        }

        # Now update the related expressedSequenceTag and estLibraryExpression tables
        if ( $debug ){
            my $printExprId = 'notRetrievedForLog';
            print "UPDATE expressedSequenceTag AS t1
                   INNER JOIN estLibrary AS t2 ON t1.estLibraryId = t2.estLibraryId
                   INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId
                   SET expressionId = $printExprId
                   WHERE bgeeGeneId = $geneId and t3.exprMappedConditionId = $exprMappedConditionId\n";
            for my $libId ( keys %{ $results{$geneId} } ){
                print "INSERT INTO estLibraryExpression
                           (expressionId, estLibraryId, estCount)
                       VALUES ($printExprId, $libId,
                       $results{$geneId}->{$libId})\n";

                my $pval = 2**(-($results{$geneId}->{$libId}+1));
                print "UPDATE expressedSequenceTag SET pValue=$pval WHERE expressionId=$printExprId AND estLibraryId=$libId\n";
            }
        } else {
            $updResult->execute($expressionId, $geneId, $exprMappedConditionId)  or die $updResult->errstr;
            for my $libId ( keys %{ $results{$geneId} } ){
                $insExpSummary->execute($expressionId, $libId, $results{$geneId}->{$libId})
                    or die $insExpSummary->errstr;

                #p = 2^-(Count+1)
                my $pval = 2**(-($results{$geneId}->{$libId}+1));
                $updPVal->execute($pval, $expressionId, $libId)  or die $updPVal->errstr;
            }
        }
    }
}

$bgee->disconnect;
print "Done\n"  if ( $debug );

exit 0;

