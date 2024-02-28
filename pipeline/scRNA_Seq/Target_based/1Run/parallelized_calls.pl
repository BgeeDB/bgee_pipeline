#!/usr/bin/env perl

## This script allows generate calls for all required libraries. It takes as input the metadata_info_10X.txt
## file and test for each library if calls have already been generated. If not, a job is run.
## Once all jobs finished, a tsv file and a plot summarizing present/absent calls are generated.

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Try::Tiny;
use FindBin;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;
use File::Path qw(make_path);
use File::Basename;
use Time::Piece;
use Data::Dumper;
## Define arguments & their default value
my ($metadataFile, $parallelJobs, $refIntergenicFolder, $cellTypeFolder, $outputDir,
    $queue, $account, $pathToCallsScript, $pathToSummaryScript, $pValueCutoff,
    $rLibs) = ('', '', '', '', '', '',  '', '', '' , '', '');
my %opts = ('metadataFile=s'                     => \$metadataFile,
            'parallelJobs=s'                     => \$parallelJobs,
            'refIntergenicFolder=s'              => \$refIntergenicFolder,
            'cellTypeFolder=s'                   => \$cellTypeFolder,
            'outputDir=s'                        => \$outputDir,
            'queue=s'                            => \$queue,
            'account=s'                          => \$account,
            'pathToCallsScript=s'                => \$pathToCallsScript,
            'pathToSummaryScript=s'              => \$pathToSummaryScript,
            'pValueCutoff=s'                     => \$pValueCutoff,
            'rLibs=s'                            => \$rLibs
           );

######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $parallelJobs eq '' || $refIntergenicFolder eq '' ||             
    $cellTypeFolder eq '' || $outputDir eq '' || $queue eq '' || $account eq '' || $pathToCallsScript eq '' ||
    $pathToSummaryScript eq '' || $pValueCutoff eq '' || $rLibs eq '') {
    print "\n\tInvalid or missing argument:
\t-metadataFile                 file metadata_info containing all run to process
\t-parallelJobs                  maximum number of jobs to run in parallel
\t-refIntergenicFolder           Path to the directory containing reference intergenic sequences
\t-cellTypeFolder                Path to the directory containing celltype expression quantification
\t-outputDir                     path where should be saved results of this script
\t-queue                         queue to use to run jobs on the cluster
\t-account                       account to use to run jobs on the cluster
\t-pathToCallsScript             path to the R script that will generate calls for each library
\t-pathToSummaryScript           path to the R script that will generate tsv and pdf summary
\t-rLibs                         path to the directory containing R packages
\t-pValueCutoff                  pValue cutoff used to consider a call present/absent
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# Info of processed libraries coming from the pipeline
my %processedLibraries = get_processed_libraries_info($metadataFile, 1);

print Dumper(\%opts);
my %sbatchToRun = ();

my $jobPrefix = "calls_";

## create directories
my $clusterOutput = "${outputDir}/clusterOutput/";
my $sbatchLocation = "${outputDir}/sbatch/";
if (! -d $clusterOutput) {
    make_path($clusterOutput) or die "Couldn't create $clusterOutput directory, $!";
}
if (! -d $sbatchLocation) {
    make_path($sbatchLocation) or die "Couldn't create $sbatchLocation directory, $!";
}
my $jobs_created = 0;
## first create sbatch files and add them to an array of sbatch to run
foreach my $experimentId (keys %processedLibraries){
    foreach my $libraryId (keys %{$processedLibraries{$experimentId}}){
        # do not consider libraries with calls already processed properly
	    next if -e "${outputDir}/${libraryId}/${libraryId}_stats.tsv";
        #create sbatch file and
        my $jobName = "${jobPrefix}${libraryId}";
        my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
        my $speciesName = $processedLibraries{$experimentId}{$libraryId}{'speciesName'};
	## use 5Gb of memory.
	    my $sbatchTemplate = Utils::sbatch_template($queue, $account, 1,
          2, "${clusterOutput}${jobName}.out", "${clusterOutput}/${jobName}.err",
          $jobName);

        #TODO: move modules management to a script attribute
        $sbatchTemplate .= "module use /software/module/\n";
	    $sbatchTemplate .= "module add R/3.6.1;\n";
        $sbatchTemplate .= "\nexport R_LIBS_USER=$rLibs\n\n";
	$speciesName =~ s/ /_/g;
        my $commandToRun = "Rscript $pathToCallsScript libraryId=\\\"$libraryId\\\"".
            " speciesId=\\\"$speciesId\\\" speciesName=\\\"$speciesName\\\"".
            " celltypeFolder=\\\"$cellTypeFolder\\\" refIntergenicFolder=\\\"$refIntergenicFolder\\\"".
            " pValueCutoff=\\\"$pValueCutoff\\\" callsOutputFolder=\\\"$outputDir\\\"";

	    $sbatchTemplate .= "$commandToRun\n";
        my $sbatchFilePath = "$sbatchLocation$jobName.sbatch";
        $sbatchToRun{$experimentId}{$libraryId} = $sbatchFilePath;
        open(FH, '>', $sbatchFilePath) or die $!;
        print FH $sbatchTemplate;
        close(FH);
        $jobs_created++;
    }
}

print "created $jobs_created sbatch files.\n";

my $numberJobRun = 0;
my $startTime = localtime->strftime('%Y-%m-%dT%H:%M:%S');
my $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
foreach my $experimentId (keys %sbatchToRun){
    foreach my $libraryId (keys %{$sbatchToRun{$experimentId}}){
	    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
	    while ($jobsRunning >= $parallelJobs) {
            sleep(15);
            $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
    	}
	    system("sbatch $sbatchToRun{$experimentId}{$libraryId} >/dev/null");
        $numberJobRun++;
    }
}

print "all jobs created properly. Run $numberJobRun jobs\n";

while ($jobsRunning > 0) {
    sleep(15);
    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
}

my %jobs_status = Utils::count_status_finished_jobs($jobPrefix, $startTime);
print $jobs_status{"completed"}." jobs completed, ".$jobs_status{"failed"}." jobs failed, ".$jobs_status{"out_of_memory"}." jobs failed with an out of memory issue and ".$jobs_status{"cancelled"}." jobs have been cancelled.\n";

