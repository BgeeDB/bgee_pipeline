#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Frederic Bastian, created July 2012
# USAGE: perl check_mas5_filtered.pl
# a script to check that MAS5 files were properly filtered
#############################################################

$| = 1;

use File::Spec;
use Getopt::Long;
require 'mas5_utils.pl';
require 'affy_utils.pl';

my $fileDir = 'processed_mas5';


# Define arguments & their default value
my ($affyChipFilesDir, $affymetrixChip) = ('', '');
my ($debug) = (0);
my %opts = ('affyChipFilesDir=s' => \$affyChipFilesDir,
            'affymetrixChip=s'   => \$affymetrixChip,
            'debug'              => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $affyChipFilesDir eq '' || $affymetrixChip eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -affyChipFilesDir=\$(AFFYDATAPATH)  -affymetrixChip=\$(AFFY_CHIP_FILEPATH)
\t-affyChipFilesDir     Main Affymetrix directory
\t-affymetrixChip       affymetrixChip.tsv annotation file
\n";
    exit 1;
}


print "Getting chips information...\n"  if ( $debug );
# Get chip annotations
my %affyChips = getAllChips($affymetrixChip);

print "Checking if all required information is present...\n"  if ( $debug );
# Now, check that we have all required info
for my $expId ( %affyChips ){
    for my $chipId ( keys %{$affyChips{$expId}} ){
        # if the chip is commented, skip
        next  if ( $affyChips{$expId}{$chipId}{'commented'} eq '1');
        # if it is not a MAS5 file
        next  if ( $affyChips{$expId}{$chipId}{'normalizationTypeId'} ne '1' );

        my ($volume, $directories, $file) = File::Spec->splitpath($affyChipFilesDir, 1);
        my @dirs = File::Spec->splitdir($directories);
        push(@dirs, $fileDir);
        push(@dirs, $expId);

        my $mas5File = File::Spec->catfile(@dirs, $chipId);
        if ( !is_a_valid_filtered_file($mas5File) ){
            warn "Warning, invalid filtered MAS5 file for chipId: $chipId - expId: $expId - path: $mas5File\n";
        }
    }
}

exit 0;

