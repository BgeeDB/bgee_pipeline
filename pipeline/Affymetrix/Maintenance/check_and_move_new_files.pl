#!/usr/bin/env perl

# Frederic Bastian, June 2012
# Script to check whether new experiments taken from annotators computers are already present
# USAGE: perl check_new_files.pl <Directory where new files are stored, e.g. extra/pipeline/curation/Affymetrix/new_files/> <Directory where already present files are stored, e.g. extra/pipeline/curation/Affymetrix/> <OPTIONAL: 1. If this argument is set, it will move to the directory where files are permanently stored, new_files not already present on the server>
#######################################################

# Perl core modules
use strict;
use warnings;
use diagnostics;

use File::Copy;

$| = 1;

my $newFilesDir  = $ARGV[0] || '';
my $oldFilesDir  = $ARGV[1] || '';
my $move         = $ARGV[2] || 0;

if ( !-e $newFilesDir || !-d $newFilesDir ){
    die "Could not find [$newFilesDir], or is not a directory\n";
}

# We do not want to move all old files to the new directory because of an error in providing arguments
# so we check that the new directory includes 'new_files/' in it
if ( $newFilesDir !~ /new_files/ ){
    die "Are you sure you provided the correct new files directory?\n";
}
if ( !-e $oldFilesDir || !-d $oldFilesDir ){
    die "Could not find [$oldFilesDir], or is not a directory\n";
}
if ( $newFilesDir eq $oldFilesDir ){
    die "Same directories provided...\n";
}


my $dirNewCelData = $newFilesDir.'cel_data/';
opendir(IMD, $dirNewCelData)  or die("Cannot open directory $dirNewCelData\n");
my @files = readdir(IMD);
closedir(IMD);
my @newCelDataExp;
for my $newFile ( @files ){
    next  if ( $newFile =~ /^\./ );

    if ( -d $dirNewCelData.$newFile ){
        push(@newCelDataExp, $newFile);
    }
}


my @celDataAlreadyPresent;
my @celDataNotPresent;
my $dirCelDataAlreadyPresent = $oldFilesDir.'cel_data/';
for my $newFile ( @newCelDataExp ){
    if ( -e $dirCelDataAlreadyPresent.$newFile && -d $dirCelDataAlreadyPresent.$newFile ){
        push(@celDataAlreadyPresent, $newFile);
    }
    else {
        push(@celDataNotPresent,     $newFile);
    }
}


my $dirNewProcessedData = $newFilesDir.'processed_mas5_original_files/';
opendir(IMD, $dirNewProcessedData)  or die("Cannot open directory $dirNewProcessedData\n");
@files = readdir(IMD);
closedir(IMD);
my @newProcessedDataExp;
for my $newFile ( @files ){
    next  if ( $newFile =~ /^\./ );

    if ( -d $dirNewProcessedData.$newFile ){
        push(@newProcessedDataExp, $newFile);
    }
}


my @processedDataAlreadyPresent;
my @processedDataNotPresent;
my $dirProcessedDataAlreadyPresent = $oldFilesDir.'processed_mas5_original_files/';
for my $newFile ( @newProcessedDataExp ){
    if ( -e $dirProcessedDataAlreadyPresent.$newFile && -d $dirProcessedDataAlreadyPresent.$newFile ){
        push(@processedDataAlreadyPresent, $newFile);
    }
    else {
        push(@processedDataNotPresent,     $newFile);
    }
}


my $dirNewProcessedFilteredData = $newFilesDir.'processed_mas5/';
opendir(IMD, $dirNewProcessedFilteredData)  or die("Cannot open directory $dirNewProcessedFilteredData\n");
@files = readdir(IMD);
closedir(IMD);
my @newProcessedFilteredDataExp;
for my $newFile ( @files ){
    next  if ( $newFile =~ /^\./ );

    if ( -d @newProcessedFilteredDataExp.$newFile ){
        push(@newProcessedFilteredDataExp, $newFile);
    }
}


my @processedFilteredDataAlreadyPresent;
my @processedFilteredDataNotPresent;
my $dirProcessedFilteredDataAlreadyPresent = $oldFilesDir.'processed_mas5/';
for my $newFile ( @newProcessedFilteredDataExp ){
    if ( -e $dirProcessedFilteredDataAlreadyPresent.$newFile && -d $dirProcessedFilteredDataAlreadyPresent.$newFile ){
        push(@processedFilteredDataAlreadyPresent, $newFile);
    }
    else {
        push(@processedFilteredDataNotPresent, $newFile);
    }
}



print "New CEL data already present (= problem):\n";
if ( scalar(@celDataAlreadyPresent) == 0 ){
    print "none\n";
}
for my $file ( @celDataAlreadyPresent ){
    print "$file\n";
}

print "\nNew processed original MAS5 data already present (= problem):\n";
if ( scalar(@processedDataAlreadyPresent) == 0 ){
    print "none\n";
}
for my $file ( @processedDataAlreadyPresent ){
    print "$file\n";
}

print "\nNew processed filtered MAS5 data already present (= problem):\n";
if ( scalar(@processedFilteredDataAlreadyPresent) == 0 ){
    print "none\n";
}
for my $file ( @processedFilteredDataAlreadyPresent ){
    print "$file\n";
}

print "\n\nNew CEL data not present (= OK): ", scalar(@celDataNotPresent), "\n";
if ( scalar(@celDataNotPresent) == 0 ){
    print "none\n";
}
for my $file ( @celDataNotPresent ){
    print "$file";
    if ( defined $move && $move eq "1" ){
        if ( move($dirNewCelData.$file, $dirCelDataAlreadyPresent.$file) == 1 ){
            print ' - move to the server successful';
        }
        else {
            print " - error while moving the file: $!";
        }
    }
    print "\n";
}

print "\nNew processed original MAS5 data not present (= OK): ", scalar(@processedDataNotPresent), "\n";
if ( scalar(@processedDataNotPresent) == 0 ){
    print "none\n";
}
for my $file ( @processedDataNotPresent ){
    print "$file";
    if ( defined $move && $move eq "1" ){
        if ( move($dirNewProcessedData.$file, $dirProcessedDataAlreadyPresent.$file) == 1 ){
            print ' - move to the server successful';
        }
        else {
            print " - error while moving the file: $!";
        }
    }
    print "\n";
}

print "\nNew processed filtered MAS5 data not present (= OK): ", scalar(@processedFilteredDataNotPresent)," \n";
if ( scalar(@processedFilteredDataNotPresent) == 0 ){
    print "none\n";
}
for my $file ( @processedFilteredDataNotPresent ){
    print "$file";
    if ( defined $move && $move eq "1" ){
        if ( move($dirNewProcessedFilteredData.$file, $dirProcessedFilteredDataAlreadyPresent.$file) == 1 ){
            print ' - move to the server successful';
        }
        else {
            print " - error while moving the file: $!";
        }
    }
    print "\n";
}

exit 0;

