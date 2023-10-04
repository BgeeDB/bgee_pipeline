#!/usr/bin/env perl

## This script allows to download FASTQ files from SRA. It is highly inspired from the download script of the target based pipeline

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Try::Tiny;
use FindBin;
use Cwd;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;
use File::Path qw(make_path);
use File::Basename;

## Define arguments & their default value
my ($metadataFile, $parallelJobs, $downloadedLibraries, $outputDir, $queue, $account) = ('', '', '', '', '',  '');
my %opts = ('metadataFile=s'        => \$metadataFile,
            'parallelJobs=s'        => \$parallelJobs,
            'downloadedLibraries=s' => \$downloadedLibraries,
            'outputDir=s'           => \$outputDir,
            'queue=s'               => \$queue,
            'account=s'             => \$account
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $parallelJobs eq '' || $outputDir eq '' || $downloadedLibraries eq '' ||
    $queue eq '' || $account eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0 -metatadataFile=... -parallelJobs=50 -outputDir=...  >> $@.tmp 2> $@.warn
\t-metadataFile            file containing metadata necessary to download each run
\t-parallelJobs            maximum number of jobs to run in parallel
\t-downloadedLibraries     file containing the ID of all alreaydy downloaded libraries for single cell
\t-outputDir               directory where FASTQ files are downloaded/generated
\t-queue                   queue to use to run jobs on the cluster
\t-account                 account to use to run jobs on the cluster
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# Info of processed libraries coming from the pipeline
my $isTargetBased = 0;
my %processedLibraries = get_processed_libraries_info($metadataFile, $isTargetBased);

# Read already downloaded libraries
my %alreadyDownloaded = map { $_ => 1 } read_file("$downloadedLibraries", chomp=>1);

my %sbatchToRun = ();

## create directory necessary to store sbatch files
make_path("$outputDir/sbatch");
make_path("$outputDir/clusterOutput");

#store initial dir location to be able to move for symlink generation and then come back later
my $initialDir = getcwd;

my $jobPrefix = "download_";
my $jobs_created = 0;
my $experimentDirName = "EXPERIMENTS";
my $experimentOutputDir = "$outputDir/$experimentDirName";
## first create sbatch files and add them to an array of sbatch to run
foreach my $experimentId (keys %processedLibraries){
    foreach my $libraryId (keys %{$processedLibraries{$experimentId}}){
        next if ( exists $alreadyDownloaded{$libraryId} );
        foreach my $runId (keys %{$processedLibraries{$experimentId}{$libraryId}}){
            next if $runId eq "speciesId";
            $jobs_created++;
            #create sbatch file and
            my $source = $processedLibraries{$experimentId}{$libraryId}{$runId}{'downloadSource'};
            my $jobName = "$jobPrefix$runId";
            ## Use 4Gb of memory. Should maybe be increase depending on the run to download
            ## ask for 4 cpus as it is the number of threads used by the bamtofastq tool
            my $sbatchTemplate = Utils::sbatch_template($queue, $account, 1,
              4, "$outputDir/clusterOutput/$jobName.out", "$outputDir/clusterOutput/$jobName.err",
              $jobName);
            my $libDirectory = "$experimentOutputDir/$experimentId/$libraryId";
            my $fastq_ftp = $processedLibraries{$experimentId}{$libraryId}{$runId}{'fastqFTP'};
            my @fastq_files = split(";", $fastq_ftp);
            $sbatchTemplate .= "module load gcc/10.4.0;\nmodule load sratoolkit/3.0.0;\n\n";
            ## download fastq from SRA split into several FASTQ files if paired library
            if($processedLibraries{$experimentId}{$libraryId}{$runId}{'libraryType'} eq "single") {
                $sbatchTemplate .= "fastq-dump --outdir $libDirectory $runId &&\n";
            } else {
                $sbatchTemplate .= "fastq-dump --outdir $libDirectory --split-files $runId &&\n";
            }
            $sbatchTemplate .= "find $libDirectory \\( -name '*.fastq' \\) -exec gzip {} \\; &&\n";
            $sbatchTemplate .= "touch $libDirectory/$runId.done";
            ## create sbatch file and add its path to the hash of sbatch files
            my $sbatchFilePath = "$outputDir/sbatch/$jobName.sbatch";
            $sbatchToRun{$experimentId}{$libraryId}{$runId} = $sbatchFilePath;
            open(FH, '>', $sbatchFilePath) or die $!;
            print FH $sbatchTemplate;
            close(FH);
        }
    }
}
print "created $jobs_created sbatch files.\n";

my $numberJobRun = 0;
my $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
foreach my $experimentId (keys %sbatchToRun){
    foreach my $libraryId (keys %{$sbatchToRun{$experimentId}}){
        my $libDirectory = "$experimentOutputDir/$experimentId/$libraryId";
        foreach my $runId (keys %{$sbatchToRun{$experimentId}{$libraryId}}){
            next if (-f "$libDirectory/$runId.done");
            $numberJobRun++;
            $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
            while ($jobsRunning >= $parallelJobs) {
                sleep(15);
                $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
            }
            make_path("$libDirectory");
            chdir "$libDirectory";
            system("sbatch $sbatchToRun{$experimentId}{$libraryId}{$runId}");
        }
    }
}

while ($jobsRunning > 0) {
    sleep(15);
    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
}
print "all download finished properly. Run $numberJobRun jobs\n";
# wait 60 seconds to be sure .done files had enough time to be created
sleep(20);
print "now start to generate symlinks of libraries per species\n";
foreach my $experimentId (keys %sbatchToRun){
    foreach my $libraryId (keys %{$sbatchToRun{$experimentId}}){
        my $libDirectory = "$experimentOutputDir/$experimentId/$libraryId";
        my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
        my $done = 1;
        foreach my $runId (keys %{$sbatchToRun{$experimentId}{$libraryId}}){
            next if $runId eq "speciesId";
            if(! -e "$libDirectory/$runId.done") {
                $done = 0;
                print "error for run $runId file $libDirectory/$runId.done does not exist.\n";
            }
        }
        if ($done) {
            # if all the run of this library have been properly downloaded then
            # we generate a symlink
            print "done download of library $libraryId\n";
            my $speciesDir = "$outputDir/$speciesId";
            make_path("$speciesDir");
            chdir "$speciesDir";
            system("ln -s ../$experimentDirName/$experimentId/$libraryId $libraryId");
            # we also add this library to the list of already generated libraries
            $alreadyDownloaded{$libraryId} = 1;
        } else {
            warn "Did not properly download the library $libraryId";
            rmdir $libDirectory;
        }
    }
}

#now go back to original location to update file listing all downloaded libraries
chdir "$initialDir";
print "Finally update the file containing downloaded libraries";
open my $outFh, "> ", "$downloadedLibraries" or die "Cannot write: $! \n";
foreach my $libraryId (keys %alreadyDownloaded) {
    print $outFh "$libraryId\n";
}
close $outFh;