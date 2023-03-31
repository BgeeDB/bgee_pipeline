#!/usr/bin/env perl

## This script allows to run the knee plot filtering + celltype identification bus for all required libraries

##TODO: integrate this script in the pipeline.
## example of command :
##           perl 1Run/parallelized_qc.pl -metadataFile=/users/jwollbre/Documents/git/bgee_pipeline/generated_files/scRNA_Seq/Target_based/all_downloaded_metadata_info_10X.txt -parallelJobs=100 -fastqDir=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/FASTQ/target_base/ -gtfDir=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/rna_seq/GTF_15/ -scRNASeqInfoFile=/users/jwollbre/Documents/git/bgee_pipeline/generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt -kallistoResults=/scratch/beegfs/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/kallisto_bus_results_bgee_v15_0 -queue=normal -account=mrobinso_bgee_sensitive -pathToScript=$PWD/1Run/kallisto_bus_one_lib.R > kallisto.out 2>> kallisto.err

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

## Define arguments & their default value
my ($metadataFile, $parallelJobs, $gtfDir, $barcodeFolder, $kallistoResults, $outputDir,
    $queue, $account, $pathToScript) = ('', '', '', '',  '', '', '' , '', '');
my %opts = ('metadataFile=s'                     => \$metadataFile,
            'parallelJobs=s'                     => \$parallelJobs,
            'gtfDir=s'                           => \$gtfDir,
            'barcodeFolder=s'                    => \$barcodeFolder,
            'kallistoResults=s'                  => \$kallistoResults,
            'outputDir=s'                        => \$outputDir,
	        'queue=s'                            => \$queue,
	        'account=s'                          => \$account,
            'pathToScript=s'                     => \$pathToScript
           );

command_arg <- c("libraryId", "output")

######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $parallelJobs eq '' || $gtfDir eq '' || $barcodeFolder eq '' ||             
    $kallistoResults eq '' || $outputDir eq '' || $queue eq '' || $account eq '' || $pathToScript eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0 -metadataFile=... -parallelJobs=50 -outputDir=...  >> $@.tmp 2> $@.warn
\t-metadataFile                 file metadata_info containing all run to process
\t-parallelJobs                  maximum number of jobs to run in parallel
\t-gtfDir                        Folder where is placed the informative files as transcriptomes index + gtf_all
\t-barcodeFolder                 Path to the directory containing barcode to celltype annotation for all experiments
\t-kallistoResults               path where should be saved all kallisto bus results per library_ID
\t-outputDir                     path where should be saved results of this script
\t-queue                         queue to use to run jobs on the cluster
\t-account                       account to use to run jobs on the cluster
\t-pathToScript                  path to the R script that will run kallisto for each library
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../target_base_utils.pl");

# Info of processed libraries coming from the pipeline
my %processedLibraries = get_processed_libraries_info($metadataFile);

my %sbatchToRun = ();

## create directory necessary to store sbatch files
make_path("$outputDir/sbatch");
my $jobPrefix = "cellId_";
my $clusterOutput = "${$outputDir}/clusterOutput/";
my $sbatchLocation = "${$outputDir}/sbatch/";
if (! -d $clusterOutput) {
    mkdir($clusterOutput) or die "Couldn't create $clusterOutput directory, $!";
}
if (! -d $sbatchLocation) {
    mkdir($sbatchLocation) or die "Couldn't create $sbatchLocation directory, $!";
}
my $jobs_created = 0;
## first create sbatch files and add them to an array of sbatch to run
foreach my $experimentId (keys %processedLibraries){
    foreach my $libraryId (keys %{$processedLibraries{$experimentId}}){
        #create sbatch file and
	    my $jobName = "${jobPrefix}${libraryId}";
        my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
        ## use 50Gb of memory. Should maybe be increase depending on the run to download
	    my $sbatchTemplate = Utils::sbatch_template($queue, $account, 1,
          50, "${clusterOutput}${jobName}.out", "${clusterOutput}/${jobName}.err",
          $jobName);

        #TODO: move modules management to a script attribute
        $sbatchTemplate .= "module use /software/module/\n";
        $sbatchTemplate .= "export PATH=/software/bin:\$PATH\n";
	    $sbatchTemplate .= "module add R/3.6.1\n";
        $sbatchTemplate .= "module load UHTS/Analysis/kallisto/0.46.0\n";
        my $commandToRun = "R CMD BATCH --no-save --no-restore \'--args libraryId=\"$libraryId\"".
            " scRNASeq_Info=\"$scRNASeqInfoFile\" kallisto_bus_results=\"$kallistoResults\"".
            " folderSupport=\"$gtfDir\" infoFolder=\"$barcodeFolder\" output=\"$outputDir\"\'".
            " $pathToScript ${clusterOutput}/${jobName}.Rout";

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
my $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
foreach my $experimentId (keys %sbatchToRun){
    foreach my $libraryId (keys %{$sbatchToRun{$experimentId}}){
	    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
	    while ($jobsRunning >= $parallelJobs) {
             sleep(15);
             $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
    	}
	    system("sbatch $sbatchToRun{$experimentId}{$libraryId}");
        $numberJobRun++;
    }
}

print "all jobs created properly. Run $numberJobRun jobs\n";

while ($jobsRunning > 0) {
    sleep(15);
    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
}

print "all jobs finished\n";

