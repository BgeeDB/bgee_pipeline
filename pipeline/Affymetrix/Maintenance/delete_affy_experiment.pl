#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;


# Frederic Bastian, May 2012
# Script to delete files related to an Affymetrix experiment
# USAGE: perl delete_affy_experiment.pl <Directory where Affymetrix files are stored (e.g. /var/bgee/extra/pipeline/Affymetrix/)> <ExpId to delete>
#######################################################

$| = 1;

my $affyDir        = $ARGV[0];
my $affyChipInfo   = $ARGV[1];
my $expIdToDelete  = $ARGV[2];

if ( !defined $affyDir || !defined $expIdToDelete ){
	die "Wrong parameters:\n\t$0 \$(AFFYDATAPATH) \$(AFFY_CHIPINFO_FILEPATH) <ExpId_to_delete>\n";
}

# bioconductor/out/: file names follow the pattern experimentId_chipTypeId.out
remove_files($affyDir.'/bioconductor/out/', $expIdToDelete.'_.+?\.out$');
# bioconductor/differential/: Names follow the pattern experimentId_chipTypeId.out
remove_files($affyDir.'/bioconductor/differential/', $expIdToDelete.'_.+?\.out$');
# bioconductor/targets/: Names of these files follow the pattern experimentId___chipTypeId.target
remove_files($affyDir.'/bioconductor/targets/', $expIdToDelete.'___.+?\.target$');

# cel_data/: CEL files, stored experiment by experiment (e.g., cel_data/expId1/, cel_data/expId2/)
removeDir($affyDir.'/cel_data/'.$expIdToDelete.'/');
# processed_differential/: results of diff analyses, stored experiment by experiment (e.g., processed_differential/expId1/, processed_differential/expId2/)
removeDir($affyDir.'/processed_differential/'.$expIdToDelete.'/');
# processed_schuster/: results of the Schuster analyses, stored experiment by experiment
removeDir($affyDir.'/processed_schuster/'.$expIdToDelete.'/');
# processed_mas5/: results of the MAS5 analyses, stored experiment by experiment
removeDir($affyDir.'/processed_mas5/'.$expIdToDelete.'/');
removeDir($affyDir.'/processed_mas5_original_files/'.$expIdToDelete.'/');


# Clean affymetrixChipInformation file
system("grep -v -P '\t$expIdToDelete\t' $affyChipInfo > $affyChipInfo.tmp; mv $affyChipInfo.tmp $affyChipInfo");

exit 0;


sub remove_files {
    my ($dir, $pattern) = @_;
    # dir where the files to delete are stored
    # pattern of the files names to delete

    if ( -e $dir && -d $dir ){
        opendir(my $IMD1, $dir)  or die("Cannot open directory [$dir]\n");
        my @files = readdir($IMD1);
        closedir($IMD1);
        for my $file ( @files ){
            next  if ( $file eq '.' || $file eq '..' );

            if ( $file =~ /$pattern/ ){
                unlink($dir.$file)  or die "Cannot remove file [$dir$file]: $!\n";
                print "File [$dir$file] removed\n";
            }
        }
    }
    else {
        warn "[$dir] is not a directory\n";
    }
    return;
}

sub removeDir {
    my ($dir) = @_;

    if ( -e $dir && -d $dir ){
        opendir(my $IMD2, $dir)  or die("Cannot open directory [$dir]\n");
        my @files = readdir($IMD2);
        closedir($IMD2);
        for my $storedFile ( @files ){
            next  if ( $storedFile eq '.' || $storedFile eq '..' );

            if ( -d $dir.$storedFile ){
                removeDir($dir.$storedFile.'/');
            }
            else {
                unlink($dir.$storedFile)  or die "Cannot remove file [$dir.$storedFile]: $!\n";
            }
        }
        rmdir($dir)  or die "Cannot remove dir $dir\n";
        print "Directory [$dir] removed\n";
    }
    else {
        warn "[$dir] is not a directory\n";
    }
}

