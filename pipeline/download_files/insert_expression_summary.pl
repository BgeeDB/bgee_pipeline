#!/usr/bin/env perl

# Julien Wollbrett, April 2025

# This script is reponsible for inserting gene summary expression sentence in the database

#############################################################

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Data::Dumper;

$| = 1;

# Define arguments & their default value
my ($bgee_connector)                = ('');
my ($expression_summary_file)       = ('');
my ($debug)       = (0);

my %opts = ('bgee=s'                    => \$bgee_connector,     # Bgee connector string
            'expression_summary_file=s' => \$expression_summary_file,
            'debug'                     => \$debug 
           );

my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $expression_summary_file eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEE_CMD) -processed_dir=\$(PROCESSED_DIR) -calls_dir=\$(CALLS_DIR)
\t-bgee                         Bgee connector string
\t-expression_summary_file      The path to the file containing summary expression sentence for all genes in Bgee
\t-debug                        Debug mode, do not insert data in the database
\n";
    exit 1;
}


# read the expression summary file. It is a tsv file containing three columns without header.
# gene_id, species_id and expression_summary_sentence
# the expression summary sentence is a string containing the expression summary for the gene
# in the species.
my %expression_summary_hash;
open(my $expression_summary_fh, '<', $expression_summary_file) or die "Could not open file '$expression_summary_file' $!";
while (my $line = <$expression_summary_fh>) {
    chomp $line;
    my ($gene_id, $species_id, $expression_summary_sentence) = split(/\t/, $line);
    $expression_summary_hash{$gene_id}{$species_id} = $expression_summary_sentence;
}
close($expression_summary_fh);

# connect to the database
my $dbh = Utils::connect_bgee_db($bgee_connector);
# check if the connection is successful

# insert the expression summary sentence in the database
my $insert_expression_summary_query = "UPDATE gene as t1 SET t1.expressionSummary = ? WHERE t1.geneId = ? AND t1.speciesId = ?";
my $insert_stmt = $dbh->prepare($insert_expression_summary_query);

# loop through the expression summary hash and insert the values in the database
foreach my $gene_id (keys %expression_summary_hash) {
    foreach my $species_id (keys %{$expression_summary_hash{$gene_id}}) {
        if ($debug) {
            print "UPDATE gene as t1 SET t1.expressionSummary = $expression_summary_hash{$gene_id}{$species_id} WHERE t1.geneId = $gene_id AND t1.speciesId = $species_id\n";
        } else {
            $insert_stmt->execute($expression_summary_hash{$gene_id}{$species_id}, $gene_id, $species_id);
        }
    }
}
$insert_stmt->finish();
# close the database connection
$dbh->disconnect();
