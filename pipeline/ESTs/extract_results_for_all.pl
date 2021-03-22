#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use File::Slurp;
use Sort::Naturally;

if ( $#ARGV ne 0 ){
    print "\n\t$0 dmel_cdna.results > dmel_mapping\n\n"; exit 2;
}

# IN
my %species_genes;
for my $line ( read_file($ARGV[0], chomp => 1) ){
    #Fields: query id, subject id, % identity, alignment length, mismatches, gap opens,    q. start, q. end, s. start, s. end, evalue,  bit score
    if ( $line =~ m/^([A-Z][a-z][a-z]?\.\d+)\s+(FBgn\d+)\|.+\s+(\d+)\.\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+e?\-?.?\d+)\s+\d+/ ){
        my $unigene  = $1;
        my $ensembl  = $2;
        my $identity = $3;
        my $evalue   = $4;
        $species_genes{$unigene}->{$ensembl}++  if ( $identity > 90 );
    }
    elsif ( $line =~ m/^([A-Z][a-z][a-z]?\.\d+)\s+(ENSG\d+)\|.+\s+(\d+)\.\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+e?\-?.?\d+)\s+\d+/ ){
        my $unigene  = $1;
        my $ensembl  = $2;
        my $identity = $3;
        my $evalue   = $4;
        $species_genes{$unigene}->{$ensembl}++  if ( $identity > 90 );
    }
    elsif ( $line =~ m/^([A-Z][a-z][a-z]?\.\d+)\s+(ENS[A-Z][A-Z][A-Z]G\d+)\|.+\s+(\d+)\.\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+e?\-?.?\d+)\s+\d+/ ){
        my $unigene  = $1;
        my $ensembl  = $2;
        my $identity = $3;
        my $evalue   = $4;
        $species_genes{$unigene}->{$ensembl}++  if ( $identity > 90 );
    }
}


# OUT
for my $unigene ( nsort keys %species_genes ){
    my $count = keys %{$species_genes{$unigene}};
    print "$unigene\t", keys %{$species_genes{$unigene}}, "\n"  if ( $count == 1 );
}

exit 0;

