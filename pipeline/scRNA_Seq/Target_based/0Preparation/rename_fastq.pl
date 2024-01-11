#!/usr/bin/env perl

## This script allows to rename fastq files
## When bam files are downloaded from SRA we use the tool bamtofastq to generate fastq files.
## The name of those files does not follow our standard naming so we have to rename them
## More than just the name, it hapens that one file type (I1, R1, R2) exists for several lanes (L001, L002, ...) for the same run. We merge those lanes in one file.

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Copy qw(move);
use File::Basename;
use File::Slurp;

## Define arguments & their default value
my ($fastqPath, $runId, $renaming) = ('', '', '');
my %opts = ('fastqPath=s'           => \$fastqPath,
            'runId=s'               => \$runId,
            'renaming=s'            => \$renaming
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ($fastqPath eq '' || $runId eq '' || $renaming eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0 -fastqPath=... -runId=... >> $@.tmp 2> $@.warn
\t-fastqPath                directory where fastq files of the run have been downloaded.
\t-runId                    ID of the run
\t-renaming                 type of renaming to apply. For now only 2 possibilities. Either \"bamtofastq\" or \"fastqdump\"
\n";
    exit 1;
}

if ($renaming eq "bamtofastq") {
    renameFastqFilesFromBamToFastq("I1", $fastqPath, 0);
    renameFastqFilesFromBamToFastq("R1", $fastqPath, 1);
    renameFastqFilesFromBamToFastq("R2", $fastqPath, 1);
} elsif ($renaming eq "fastqdump") {
    ## SRA does not provide precise names. There are no information about I1, R1 or R2 in the file names.
    ## It is not easy to find back which file contains which info.
    ## It looks like, for 10x v2 and 10x v3:
    ## 3 files ==> _1.fastq corresponds to I1. _2.fastq corresponds to R1. _3.fastq corresponds to R2.
    ## 2 files ==> _1.fastq corresponds to R1. _2.fastq corresponds to R2.
    ## As we found counter examples to this asumption we decided to rename files based on length of sequences
    ## expected sequence length: I1 < R1 < R2. This asumption is always true for 10X.
    ## It is definitly not bullet proof but we did not find any other criteria.
    ##TODO think of a better solution to rename SRA fastq files
    opendir D, $fastqPath or die "Could not open dir: $!\n";
    my @filelist = grep(/.fastq$/i, readdir D);
    my %fileToSeqLength;
    ## we never expect to have 1 file or more than 3 files. Update implementation if it can happen
    if (scalar @filelist <= 1 || scalar @filelist > 3 ) {
        die "wrong number of files. Expected 2 or 3 but have ".scalar @filelist." for run $runId";
    }
    foreach my $file (@filelist) {
        if ($file =~ /.fastq$/i) {
            my $filePath = "$fastqPath/$file";
            $fileToSeqLength{$filePath} = length(`cat $filePath | tail -n +2 | head -1`);
            print "file ".$file." has a sequence length of $fileToSeqLength{$filePath} bp.\n";
        }
    }
    #sort name of the files per sequence length
    my @orderedFilelist = sort { $fileToSeqLength{$a} <=> $fileToSeqLength{$b} } keys(%fileToSeqLength);
    #rename the files depending on sequence length
    if (scalar @orderedFilelist == 2) {
        move $orderedFilelist[0], "$fastqPath/${runId}_R1.fastq";
        move $orderedFilelist[1], "$fastqPath/${runId}_R2.fastq";
    } elsif (scalar @orderedFilelist == 3) {
        move $orderedFilelist[0], "$fastqPath/${runId}_I1.fastq";
        move $orderedFilelist[1], "$fastqPath/${runId}_R1.fastq";
        move $orderedFilelist[2], "$fastqPath/${runId}_R2.fastq";
    }
} else {
    die "unknown renaming approach.";
}

#rename fastq files and combine them if several lanes exist (e.g L001, L002, ...)
sub renameFastqFilesFromBamToFastq {
    my ($fastqFileType, $fastqPath, $fileIsMandatory) = @_;
    opendir D, $fastqPath or die "Could not open dir: $!\n";
    my @filelist = grep(/_${fastqFileType}_|_$fastqFileType.f/i, readdir D);
    my $renamedFile = "$fastqPath/${runId}_$fastqFileType.fastq.gz";
    #job die if no R1 file
    if ($fileIsMandatory) {
        if (scalar @filelist == 0) {
            die "no $fastqFileType fastq file found for that run";
        }
    }
    # rename R1 file in only one file
    if (scalar @filelist == 1) { 
        if("$fastqPath/$filelist[0]" ne $renamedFile) {
            move $filelist[0], "$renamedFile";
        }
    # merge files into a new one called with the proper name if several lanes
    } elsif (scalar @filelist > 1) {
        my $mergedCommand = "cat ";
        my @removeCmds;
        foreach my $file (sort @filelist) {
            $mergedCommand .= "$fastqPath/$file ";
            push(@removeCmds, "rm $fastqPath/$file");
        }
        $mergedCommand .= "> $renamedFile";
        system($mergedCommand);
        foreach my $removeCmd (sort @removeCmds) {
            system($removeCmd);
        }
    }
}