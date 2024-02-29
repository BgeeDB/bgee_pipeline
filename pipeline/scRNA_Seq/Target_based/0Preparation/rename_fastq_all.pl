#!/usr/bin/env perl

# This script is not integrated in the pipeline is has been used to solve a renaming issue
# after downloading fastq files (see https://www.biostars.org/p/9588563/#9588575). This problem
# should not happen again so this script could potentially be deleted.
# Kept it for now...

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;

## Define arguments & their default value
my ($metadataFile, $excludedLibraries, $fastqDir, $whitelistFolder, $queue, $account) = ('', '', '', '', '', '');
my ($doNotDownload)          = (0);
my %opts = ('metadataFile=s'        => \$metadataFile,
            'excludedLibraries=s'   => \$excludedLibraries,
            'fastqDir=s'            => \$fastqDir,
            'whitelistFolder=s'     => \$whitelistFolder,
            'queue=s'               => \$queue,
            'account=s'             => \$account
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $fastqDir eq '' || $excludedLibraries eq '' ||
    $whitelistFolder eq '' || $queue eq '' || $account eq '') {
    print "\n\tInvalid or missing argument:
\t-metadataFile            path to the rna_seq_sample_info.txt file
\t-excludedLibraries       file containing the ID of all libraries not to download
\t-whitelistFolder         path to directory containing 10X whitelists
\t-fastqDir                parent directory containing all FASTQ files
\t-queue                   queue to use to run jobs on the cluster
\t-account                 account to use to run jobs on the cluster
\n";
    exit 1;
}

require("$FindBin::Bin/../../target_base_utils.pl");

# Info of processed libraries coming from the pipeline
my $isTargetBased = 1;
my %processedLibraries = get_processed_libraries_info($metadataFile, $isTargetBased);
system("module use /software/module/");
system("module load R/3.6.1");

foreach my $experimentId (keys %processedLibraries){
    foreach my $libraryId (keys %{$processedLibraries{$experimentId}}){
        my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
        my $libDirectory = "$fastqDir/$speciesId/$libraryId";
        foreach my $runId (keys %{$processedLibraries{$experimentId}{$libraryId}{'runIds'}}){
            my $fastqRunDir = "$libDirectory/$runId/FASTQ/";
            next if ! -e $fastqRunDir;
            print "Rscript rename_fastq.R fastqPath=${fastqRunDir} runId=$runId renaming=fastqdump whitelistFolder=$whitelistFolder\n";
        }
    }
}
