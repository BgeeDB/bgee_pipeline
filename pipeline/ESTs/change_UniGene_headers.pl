#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use File::Slurp;

if ( $#ARGV ne 0 ){
    print "\n\tUSAGE: $0 Dm.seq.uniq\n\n"; exit 1;
}

# Read file in list context and directly chomp lines
for my $line ( read_file($ARGV[0], chomp => 1) ){
    if ( $line =~ m/\>gnl.+\/ug=([A-Z][a-z][a-z]?\..+)\s\// ){
        print "> $1\n";
    }
    else {
        print "$line\n";
    }
}

exit 0;

