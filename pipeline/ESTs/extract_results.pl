#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use File::Slurp;


if ( $#ARGV ne 0 ){
    print "\n\t$0 dmel_cdna.results > dmel_mapping\n\n"; exit 2;
}

# IN
my %flybase_genes;
for my $line ( read_file($ARGV[0], chomp => 1) ){
    #Fields: query id, subject id, % identity, alignment length, mismatches, gap opens,    q. start, q. end, s. start, s. end, evalue,  bit score
    if ( $line =~ m/Dm\.(\d+)\s+(FBgn\d+)\|.+\s+(\d+)\.\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+e?\-?.?\d+)\s+\d+/ ){
        my $unigene  = $1;
        my $flybase  = $2;
        my $identity = $3;
        my $evalue   = $4;
        $flybase_genes{$unigene}->{$flybase}++  if ( $identity > 90 );
    }
}


# OUT
for my $unigene ( sort {$a <=> $b} keys %flybase_genes ){
    my $count = keys %{$flybase_genes{$unigene}};
    print "Dm.$unigene\t", keys %{$flybase_genes{$unigene}}, "\n"  if ( $count == 1 );
}

exit 0;

