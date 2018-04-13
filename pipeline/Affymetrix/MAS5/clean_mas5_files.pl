#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;


# Frederic Bastian, created June 2012
#
# USAGE: perl clean_mas5_files.pl <path to mas5 original files> <path to store cleaned mas5 files>
#   path to mas5 original files: path where not-already cleaned mas5 files are stored, usually '../Affymetrix/processed_mas5_original_files/'
#   path to store cleaned mas5 files: path to the directory where to store cleaned mas5 files: usually '../Affymetrix/processed_mas5/'
# This script cleans all mas5 files present in <path to mas5 original files>, and put them in <path to store cleaned mas5 files>:
# sometimes the columns of the mas5 files are not in the proper order, or there are additional header lines, etc...
# All the files are converted to a proper format, where columns are always in the same order, and where there is only one header line
#
# USAGE: perl clean_mas5_files.pl <path to mas5 original files> <path to store cleaned mas5 files> <file path> <probeset_id column> <call column> <signal column>
# Clean only one file.
#   file path: path to an original mas5 file, RELATIVE TO <path to mas5 original files> (basically, 'expId/chipId')
#   probeset_id column: index of the column containing probeset IDs in the file.  Usually, it is 0.
#   call column: index of the column containing the mas5 calls in the file.       Usually, it is 1.
#   signal column: index of the column containing signal intensities in the file. Usually, it is 2.
#############################################################

$| = 1;

require 'mas5_utils.pl';
require 'bgee_utils.pl';
use File::Spec;
use File::Slurp;


my $processed_mas5_original_files_path = $ARGV[0] || '';
my $processed_mas5_path                = $ARGV[1] || '';
if ( !-e $processed_mas5_original_files_path || !-d $processed_mas5_original_files_path || !-e $processed_mas5_path || !-d $processed_mas5_path ){
    print "\n\tDirectories do not exist\n\t$0 \$(MAS5ORIPATH) \$(MAS5PATH)\n\n";
    exit 1;
}
if ( $processed_mas5_path =~ /original_files/ ){
    die "I don't think you want to do that.\n";
}


my $filePath           = $ARGV[2];
my $probesetIdColIndex = $ARGV[3];
my $callColIndex       = $ARGV[4];
my $signalColIndex     = $ARGV[5];

if ( defined $filePath && defined $probesetIdColIndex && defined $callColIndex && defined $signalColIndex ){
    if ( !File::Spec->file_name_is_absolute($filePath) ){
        my $absPath = File::Spec->rel2abs($filePath, $processed_mas5_original_files_path);
        # check for the existence of the file
        if ( !-e $absPath ){
            die "File [$filePath] doest not exist relatively to $processed_mas5_original_files_path\n";
        }

        # create a hash of the proper column indexes in the file
        my %column_indexes;
        $column_indexes{'probeset_id'} = $probesetIdColIndex;
        $column_indexes{'call'}        = $callColIndex;
        $column_indexes{'signal'}      = $signalColIndex;

        createCleanFile($absPath, \%column_indexes);
    }
    else {
        die "You must provide a path relative to $processed_mas5_original_files_path\n";
    }
}
else {
    # we regenerate all clean files, so we remove all the already generated files
    # removeSubDir($processed_mas5_path);
    # now we generate the clean files
    browseDirectory($processed_mas5_original_files_path);
}

exit 0;


sub createCleanFile {
    my @args = @_;
    my $file               = $args[0];
    my $column_indexes_ref = $args[1];

    # we need the absolute path to be able to write identical directories in processed_mas5/
    # it is safer if the script wasn't launched from the proper directory
    if ( !File::Spec->file_name_is_absolute($file) ){
        # should never happen, browseDirectory converts all path to absolute path
        print "Error, an absolute path must be provided\n";
        return;
    }

    # now, we get the relative directory from $processed_mas5_original_files_path,
    # in order to build identical directories in $processed_mas5_path
    my $fileRelPath = File::Spec->abs2rel($file, $processed_mas5_original_files_path);

    # now, get the directories
    my ($volume, $directories, $fileName) = File::Spec->splitpath($fileRelPath);
    my @dirs = File::Spec->splitdir($directories);
    my @newDirs;
    push(@newDirs, $processed_mas5_path);
    # recreate the same directories in $processed_mas5_path
    for my $dir ( @dirs ){
        next  if ( $dir eq '' );

        push(@newDirs, $dir);
        my $newDir = File::Spec->catdir(@newDirs);
        if ( !-e $newDir ){
            mkdir($newDir);
        }
    }

    # now, we create the clean file,
    # only if it does not already exist, and if the $column_indexes_ref was not provided
    # (in that case, it means it is a manual run of the script for one individual file)
    my $newFile = File::Spec->catfile(@newDirs, $fileName);
    if ( -e $newFile && !defined $column_indexes_ref ){
        return;
    }

    if ( !defined $column_indexes_ref ){
        my %column_indexes  = get_mas5_columns($file);
        $column_indexes_ref = \%column_indexes;
    }

    writeCleanFile($file, $newFile, $column_indexes_ref);
}

sub writeCleanFile {
    my @args = @_;
    my $originalFile       = $args[0];
    my $newFile            = $args[1];
    my $column_indexes_ref = $args[2];
    #print "$column_indexes_ref->{'call'}, $column_indexes_ref->{'signal'}, $column_indexes_ref->{'probeset_id'}\n";

    my %corresp_call = get_corresp_call();
    if ( $column_indexes_ref->{'call'} != -1 && $column_indexes_ref->{'signal'} != -1 && $column_indexes_ref->{'probeset_id'} != -1 ){
        open (my $OUT, '>', "$newFile")  or die "Cannot write [$newFile]\n";
        # a hash to store lines to be written, to sort them by probeset_id before writing them
        # {probesetId}{'call'}
        # {probesetId}{'signal'}
        my %storeOutputLines;


        my $lineCount = 0;
        my $firstLine = 1;
        for my $line ( read_file("$originalFile") ){
            $line =~ s/(\r\n|\r)$/\n/;
            chomp $line;
            my @tmp = split(/\t/, $line);

            if ( defined $tmp[$column_indexes_ref->{'call'}] ){
                $tmp[$column_indexes_ref->{'call'}]        = bgeeTrim($tmp[$column_indexes_ref->{'call'}]);
            }
            if ( defined $tmp[$column_indexes_ref->{'signal'}] ){
                $tmp[$column_indexes_ref->{'signal'}]      = bgeeTrim($tmp[$column_indexes_ref->{'signal'}]);
            }
            if ( defined $tmp[$column_indexes_ref->{'probeset_id'}] ){
                $tmp[$column_indexes_ref->{'probeset_id'}] = bgeeTrim($tmp[$column_indexes_ref->{'probeset_id'}]);
            }

            if ( $firstLine == 1 ){
                if ( $line =~ /^#/ ){
                    # it is an additional header line, no problem
                    # but it is still not the header, we do not set $firstLine to 0
                }
                else { # it is the header
                    $firstLine = 0;
                    # write the header
                    if ( defined $tmp[$column_indexes_ref->{'probeset_id'}] && defined $tmp[$column_indexes_ref->{'call'}] && defined $tmp[$column_indexes_ref->{'signal'}] ){
                        print {$OUT} $tmp[$column_indexes_ref->{'probeset_id'}]."\t".$tmp[$column_indexes_ref->{'call'}]."\t".$tmp[$column_indexes_ref->{'signal'}]."\n";
                    }
                    else {
                        print {$OUT} $line."\n";
                    }
                }
                next;
            }
            $lineCount++;

            if ( !is_valid_processed_mas5_line($line, $column_indexes_ref) ){
                print "Warning, unrecognized line, # $lineCount, in $originalFile: $line".
                    ' - probeset_id: ' . $tmp[$column_indexes_ref->{'probeset_id'}].
                    ' - call: '        . $tmp[$column_indexes_ref->{'call'}].
                    ' - signal: '      . $tmp[$column_indexes_ref->{'signal'}].
                    "\n";
                next;
            }
            $storeOutputLines{$tmp[$column_indexes_ref->{'probeset_id'}]}{'call'}   = $corresp_call{$tmp[$column_indexes_ref->{'call'}]};
            $storeOutputLines{$tmp[$column_indexes_ref->{'probeset_id'}]}{'signal'} = $tmp[$column_indexes_ref->{'signal'}];
        }

        for my $probesetId ( sort keys %storeOutputLines ){
            print {$OUT} $probesetId."\t".
                $storeOutputLines{$probesetId}{'call'}."\t".
                $storeOutputLines{$probesetId}{'signal'}."\n";
        }
        close $OUT;
    }
    else {
        print "Error, invalid mas5 file format: $originalFile\n";
    }
}

sub browseDirectory {
    my @args = @_;
    my $dir  = $args[0];

    if ( -e $dir && -d $dir ){
        opendir(my $IMD, $dir)  or die("Cannot open directory $dir\n");
        my @files = File::Spec->no_upwards(readdir($IMD));
        closedir($IMD);

        for my $storedFile ( @files ){
            next  if ( $storedFile =~ /\/\./ || $storedFile =~ /^\./ || $storedFile =~ /not_separated/ );

            my $absPath = $storedFile;
            if ( !File::Spec->file_name_is_absolute($storedFile) ){
                $absPath = File::Spec->rel2abs($storedFile, $dir) ;
            }
            if ( -d $absPath ){
                browseDirectory($absPath);
            }
            else {
                createCleanFile($absPath);
            }
        }
    }
    else {
        print "Warning: $dir is not a directory\n";
    }
}

sub removeSubDir {
    my @args  = @_;
    my $dir   = $args[0];
    my $rmDir = $args[1];

    if ( -e $dir && -d $dir ){
        opendir(my $IMD2, $dir)  or die("Cannot open directory $dir\n");
        my @files = File::Spec->no_upwards(readdir($IMD2));
        closedir($IMD2);

        for my $storedFile ( @files ){
            next  if ( $storedFile =~ /\/\./ || $storedFile =~ /^\./ );

            my $absPath = $storedFile;
            if ( !File::Spec->file_name_is_absolute($storedFile) ){
                $absPath = File::Spec->rel2abs($storedFile, $dir) ;
            }
            if ( -d $absPath ){
                removeSubDir($absPath, 1);
            }
            else {
                unlink($absPath)  or die "Cannot remove file $absPath: $!\n";
            }
        }

        if ( defined $rmDir && $rmDir == 1 ){
            rmdir($dir)  or die "Cannot remove dir $dir\n";
        }
    }
    else {
        warn "[$dir] is not a directory\n";
    }
}

