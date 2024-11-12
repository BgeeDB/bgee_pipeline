#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use File::Slurp;

# Format example:
#ID   Y34F_DROME              Reviewed;        1820 AA.
#AC   Q9W5D0; A8JUT7; A8JUT8; D8FT41; M9PDF0; M9PG59; M9PGC2; M9PGC6; M9PGJ1;
#AC   M9PIN4; O77432; O77433; Q1RKT6;
#//
#ID   Y3500_AZOVD             Reviewed;          75 AA.
#AC   C1DQL3;
#//

my $dat = $ARGV[0]  or die "\n\t$0 <sprot/trembl-ac-id.dat>\n\n";

my $mapping;
my $id = '';
my @ac;
for my $line ( read_file("$dat", chomp => 1) ){
    if ( $line =~ /^ID   (\w+)/ ){ #NOTE assume 1 ID perl ID line
        $id = $1;
    }
    elsif ( $line =~ /^AC   (.+)$/ ){
        push @ac, split(/; */, $1);
    }
    elsif ( $line =~ /^\/\// ){
        if ( $id ne '' && scalar @ac > 0 ){
            map { print "$_\t$id\n"; } @ac;
        }
        $id = '';
        @ac = ();
    }
}

exit 0;

