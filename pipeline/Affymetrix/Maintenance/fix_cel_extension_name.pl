#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use File::Slurp;
use File::Basename;

my $infile   = $ARGV[0] || '';
my $cel_path = $ARGV[1] || '';
die "\n\tNo input file: $0 affymetrixChip /var/bgee/extra/pipeline/Affymetrix/cel_data\n\n"  if ( $infile eq '' || $cel_path eq '' );


CHIPID:
for my $line ( read_file("$infile", chomp => 1) ){
    my @tmp = split(/\t/, $line);

    # $tmp[0] == 'Chip ID'
    # $tmp[1] == 'Experiment ID'
    if ( $tmp[0] =~ /^(.*?).CEL$/i ){
        my $prefix = $1;
        my @chipId = grep { /\/$prefix\.CEL/i }
                     glob($cel_path.'/'.$tmp[1].'/'.$prefix.'.*');
        if ( scalar @chipId > 1 ){
            warn "\tSeveral chip ID with the same prefix name [$prefix]: ", join(';', @chipId), "\n";
        }
        elsif ( scalar @chipId == 0 ){
            warn "\tMissing chip ID with this prefix name: [$prefix] in [$tmp[1]]\n";
        }

        if ( !-e $chipId[0] || -z $chipId[0] ){
            warn "Missing or empty file: [$chipId[0]]\n";
        }
        $chipId[0] =~ s{\.gz$}{}i;
        print basename($chipId[0]), "\n";
    }
    else {
        # Do not check header & MAS5 files ( non-*.CEL files)
        print $tmp[0], "\n";
    }
}

exit 0;

