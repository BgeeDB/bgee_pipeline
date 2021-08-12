#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

require("$FindBin::Bin/../rna_seq_utils.pl");
$| = 1;

# Julien Roux, created November 2016

# USAGE: perl insert_feature_length -length_info=<...> <OPTIONAL: -debug=1>
# This script inserts the transcript lengths into the "transcript" table
# as well as the mapping from genes to transcripts
# The file with tehse infos is generated from Kallisto output files by the
# export_feature_length.pl script.
#
# -debug=1: if provided, run in verbose mode (print the update/insert SQL queries, not executing them)
#############################################################


# Define arguments & their default value
my ($bgee_connector) = ('');
my ($length_info)    = ('');
my ($debug)          = (0);
my %opts = ('bgee=s'        => \$bgee_connector,   # Bgee connector string
            'length_info=s' => \$length_info,     # generated_files/RNA_Seq/rna_seq_length_info.txt file
            'debug'         => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $length_info eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD) -length_info=\$(RNASEQ_LENGTH_INFO_FILEPATH)
\t-bgee             Bgee connector string
\t-length_info      rna_seq_length_info.txt file
\t-debug            printing the update/insert SQL queries, not executing them
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);


print "Retrieving transcript mapping and length info...\n";
my %all_transcripts;
open(my $IN, '<', $length_info)  or die "Could not read file [$length_info]\n";
while ( defined (my $line = <$IN>) ){
    next  if ( $line =~ /^#/ );

    chomp $line;
    # file format: speciesId transcriptId, geneId, length
    my @tmp = map { bgeeTrim($_) } split(/\t/, $line);
    my $speciesId    = $tmp[0];
    my $transcriptId = $tmp[1];
    my $geneId       = $tmp[2];
    my $length       = $tmp[3];

    # record the length
    $all_transcripts{$speciesId}->{$transcriptId}->{'geneId'} = $geneId;
    $all_transcripts{$speciesId}->{$transcriptId}->{'length'} = $length;
}

print "retrieving mapping from bgeeGeneId to geneId...\n";
my %genes;
# Get hash of geneId to bgeeGeneId mapping per species
for my $speciesId ( keys %all_transcripts ){
    $genes{$speciesId} = Utils::query_bgeeGene($bgee, $speciesId);
}

print "Inserting data into the database...\n";
# query for feature length insertion
my $insLength = $bgee->prepare('INSERT INTO transcript (bgeeGeneId, transcriptId, transcriptLength) VALUES (?,?,?)');

for my $species ( keys %all_transcripts ){
    for my $transcript ( keys %{$all_transcripts{$species}} ){
        if ( $debug ){
            print 'INSERT INTO transcript (bgeeGeneId, transcriptId, transcriptLength) VALUES (', $genes{$species}->{ $all_transcripts{$species}->{$transcript}->{'geneId'}}, ',', $transcript, ',', $all_transcripts{$species}->{$transcript}->{'length'}, ")\n";
        }
        else {
            $insLength->execute($genes{$species}->{ $all_transcripts{$species}->{$transcript}->{'geneId'}}, $transcript, $all_transcripts{$species}->{$transcript}->{'length'})  or die $insLength->errstr;
        }
    }
}

$insLength->finish();
print "Done\n";
exit 0;

