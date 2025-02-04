#!/usr/bin/env perl

## This script allows to run the knee plot filtering + celltype identification bus for all required libraries

##TODO: integrate this script in the pipeline.
## example of command :
##           perl 1Run/parallelized_qcCellIdentification.pl -metadataFile=/users/jwollbre/Documents/git/bgee_pipeline/generated_files/scRNA_Seq/Target_based/all_downloaded_metadata_info_10X.txt  -parallelJobs=100 -gtfDir=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/rna_seq/GTF_15/ -barcodeFolder=/users/jwollbre/Documents/git/bgee_pipeline/generated_files/scRNA_Seq/Target_based/cleaned_barcodes/ -kallistoResults=/scratch/beegfs/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/kallisto_bus_results_bgee_v15_0/ -outputDir=/scratch/beegfs/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/QC_CellType_identification/ -queue=normal -account=mrobinso_bgee_sensitive -pathToScript=$PWD/1Run/cellIdentification.R > cellIdentification.out 2>> cellIdentification.err

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

## Define arguments & their default value
my ($metadataFile, $parallelJobs, $gtfDir, $barcodeFolder, $kallistoResults, $outputDir,
    $queue, $account, $pathToScript, $rLibs) = ('', '', '', '', '',  '', '', '' , '', '');
my %opts = ('metadataFile=s'                     => \$metadataFile,
            'parallelJobs=s'                     => \$parallelJobs,
            'gtfDir=s'                           => \$gtfDir,
            'barcodeFolder=s'                    => \$barcodeFolder,
            'kallistoResults=s'                  => \$kallistoResults,
            'outputDir=s'                        => \$outputDir,
            'queue=s'                            => \$queue,
            'account=s'                          => \$account,
            'pathToScript=s'                     => \$pathToScript,
            'rLibs=s'                            => \$rLibs
           );

######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $parallelJobs eq '' || $gtfDir eq '' || $barcodeFolder eq '' ||             
    $kallistoResults eq '' || $outputDir eq '' || $queue eq '' || $account eq '' || $pathToScript eq '' ||
    $rLibs eq '') {
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
\t-rLibs                         path to the directory containing R packages
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# Info of processed libraries coming from the pipeline
my %processedLibraries = get_processed_libraries_info($metadataFile, 1);

my %sbatchToRun = ();

## create directory necessary to store sbatch files
make_path("$outputDir/sbatch");

## uncompress all gene2biotype files before running jobs. Redirect output to /dev/null as we don't want
## to have an error if no .gene2biotype.xz file exist (they are potentially already uncompressed)
system("unxz $gtfDir/*.gene2biotype.xz >/dev/null");
my $jobPrefix = "cellId_";
my $clusterOutput = "${outputDir}/clusterOutput/";
my $sbatchLocation = "${outputDir}/sbatch/";
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
        # do not consider libraries with celltype identification already processed properly
	next if -e "${outputDir}/${libraryId}/done";
        #create sbatch file and
        my $jobName = "${jobPrefix}${libraryId}";
        my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
        ## use 30Gb of memory. Should maybe be increase depending on the library to process
	    my $sbatchTemplate = Utils::sbatch_template($queue, $account, 1,
          30, "${clusterOutput}${jobName}.out", "${clusterOutput}/${jobName}.err",
          $jobName);

        #TODO: move modules management to a script attribute
        $sbatchTemplate .= "module use /software/module/\n";
	$sbatchTemplate .= "module add R/3.6.1;\n";
        $sbatchTemplate .= "\nexport R_LIBS_USER=$rLibs\n\n";
        my $commandToRun = "R CMD BATCH --no-save --no-restore \'--args libraryId=\"$libraryId\"".
            " scRNASeq_Info=\"$metadataFile\" kallisto_bus_results=\"$kallistoResults\"".
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

## compress back all .gene2biotype files
system("xz $gtfDir/*.gene2biotype >/dev/null");
