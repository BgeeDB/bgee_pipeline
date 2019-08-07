#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Frederic Bastian, Feb. 2017

# prefill the various condition and expression tables based on %UTILS::conditionCombinations
# DEPRECATED, NOT USED ANYMORE, CONDITION PARAMETER COMBINATIONS ARE DETERMINED IN THE JAVA PIPELINE (see InsertPropagatedCalls.java)

#############################################################

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
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


#############################################################
# PREFILL CONDITION TABLES                                  #
#############################################################
COMB: for my $combination ( keys %Utils::conditionCombinations ){
    # do nothing for the raw conditions and expression calls, already computed by pipeline
    if ( $combination eq $Utils::allFieldCombination ) {
        next COMB;
    }
    print "***** combination $combination *****\n";

    my $condTable = $Utils::conditionCombinations{$combination}->{'condTable'};
    my $condId    = $Utils::conditionCombinations{$combination}->{'condId'};
    my $exprTable = $Utils::conditionCombinations{$combination}->{'exprTable'};
    my @fields    = @{ $Utils::conditionCombinations{$combination}->{'condFields'} };

    ######################################
    # INSERT CONDITIONS                  #
    ######################################
    my $sql = 'INSERT INTO '.$condTable.' ('
        .Utils::get_fields_for_sql_select(\@fields)
        .') SELECT DISTINCT '
        .Utils::get_fields_for_sql_select(\@fields)
        .' FROM cond WHERE cond.conditionId = cond.exprMappedConditionId'; # we use only conditions used in expression table
    if ( $debug ) {
        print $sql."\n";
    } else {
        my $insCond   = $bgee->prepare($sql);
        $insCond->execute()  or die $insCond->errstr;
    }

    ######################################
    # INSERT CONDITIONS                  #
    ######################################
    $sql = 'INSERT INTO '.$exprTable.' (bgeeGeneId, '.$condId
        .') SELECT DISTINCT t1.bgeeGeneId, t3.'.$condId
        .' FROM expression AS t1 INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId'
        .' INNER JOIN '.$condTable.' AS t3 ON '.Utils::get_fields_for_sql_join(\@fields, 't2', 't3');

    if ( $debug ) {
        print $sql."\n";
    } else {
        my $insExpr   = $bgee->prepare($sql);
        $insExpr->execute()  or die $insExpr->errstr;
    }
}

$bgee->disconnect;
print "Done\n"  if ( $debug );

exit 0;
