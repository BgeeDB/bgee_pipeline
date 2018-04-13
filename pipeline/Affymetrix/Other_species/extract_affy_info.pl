#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use File::Slurp;


my $tsv_file = $ARGV[0]  or die "Cannot read tsv file.\n$0 tsv_file\n";

print "#Experiment\tTissue\tLife_stage\tSex\tTreatment\tGSM\tPlatform\n";
for my $line ( read_file("$tsv_file", chomp => 1) ){
    # N2MASample.csv: N2 means standard strain of C. elegans + MA means MicroArray
    next  if ( $line =~ /^Experiment\t/ || $line =~/^#/ ); # Header in N2MASample.csv
    # Experiment    Tissue    Life stage    Sex    Treatment    GEO Record    Platform

    my @field = split("\t", $line);
    next  if ( $field[6] !~ /^GPL200$/ );             # ONLY microarray data, from GPL200 Platform
    next  if ( $field[5] !~ /^GEO record: GSM\d+$/ ); # ONLY experiment with GEO id

    # Keep only single tissue and single stages
    if ( $field[1] =~ /[,\|]/ ){
        warn "\tSeveral tissues: [$field[1]]\n";
        next;
    }
    if ( $field[2] =~ /[,\|]/ ){
        warn "\tSeveral stages:  [$field[2]]\n";
        next;
    }

    # Replace organ N.A. by whole organism
    if ( $field[1] eq 'N.A.' || $field[1] eq 'whole animal' ){
        $field[1] = 'UBERON:0001062';
    }
    # Replace stage N.A. by whole development
    if ( $field[2] eq 'N.A.' ){
        $field[2] = 'UBERON:0000104';
    }

    $field[0] =~ s{:.+$}{};
    $field[1] =~ s{\([^\)]+?\)}{}g; # Remove text label, keep only ontology ids
    $field[2] =~ s{\([^\)]+?\)}{}g;
    $field[5] =~ s{^GEO record: }{};

    print join("\t", @field), "\n";
}


exit 0;

