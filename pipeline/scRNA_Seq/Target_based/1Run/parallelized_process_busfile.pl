#!/usr/bin/env perl

## This script allows to run bustools bus for all required libraries

## example of command :
##           perl 1Run/parallelized_process_busFile.pl -metadataFile=/users/jwollbre/Documents/git/bgee_pipeline/generated_files/scRNA_Seq/Target_based/all_downloaded_metadata_info_10X.txt -parallelJobs=100 -fastqDir=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/FASTQ/target_base/ -gtfDir=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/rna_seq/GTF_15/ -scRNASeqInfoFile=/users/jwollbre/Documents/git/bgee_pipeline/generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt -whiteListPath=/users/jwollbre/Documents/git/bgee_pipeline/source_files/scRNA_Seq/Target_based/ -kallistoResults=/scratch/beegfs/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/kallisto_bus_results_bgee_v15_0 -queue=normal -account=mrobinso_bgee_sensitive -pathToScript=$PWD/1Run/process_busfile_one_lib.R > kbustools.out 2>> bustools.err

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
my ($metadataFile, $parallelJobs, $fastqDir, $gtfDir, $scRNASeqInfoFile, $whiteListPath,
    $kallistoResults, $queue, $account, $pathToScript) = ('', '', '', '',  '', '', '' , '', '', '');
my %opts = ('metadataFile=s'        => \$metadataFile,
            'parallelJobs=s'        => \$parallelJobs,
            'fastqDir=s'            => \$fastqDir,
            'gtfDir=s'              => \$gtfDir,
            'scRNASeqInfoFile=s'    => \$scRNASeqInfoFile,
            'whiteListPath=s'       => \$whiteListPath,
            'kallistoResults=s'     => \$kallistoResults,   
            'queue=s'               => \$queue,
            'account=s'             => \$account,
            'pathToScript=s'        => \$pathToScript
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $parallelJobs eq '' || $fastqDir eq '' || $gtfDir eq '' ||
    $scRNASeqInfoFile eq '' || $whiteListPath eq '' || $kallistoResults eq '' || $queue eq '' || $account eq '' || $pathToScript eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0 -metatadataFile=... -parallelJobs=50 -outputDir=...  >> $@.tmp 2> $@.warn
\t-metadataFile            file containing metadata necessary to download each run
\t-parallelJobs            maximum number of jobs to run in parallel
\t-fastqDir                directory where FASTQ files are downloaded/generated
\t-gtfDir                  Folder where is placed the informative files as transcriptomes index + gtf_all
\t-scRNASeqInfoFile        Path to the scRNA_Seq_info_TargetBased file
\t-whiteListPath           Path to the folder containing the barcode_whitelist files
\t-kallistoResults         path where should be saved all kallisto bus results per library_ID
\t-queue                   queue to use to run jobs on the cluster
\t-account                 account to use to run jobs on the cluster
\t-pathToScript            path to the R script that will run kallisto for each library
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# Info of processed libraries coming from the pipeline
my %processedLibraries = get_processed_libraries_info($metadataFile, 1);

my %sbatchToRun = ();

## create directory necessary to store sbatch files
make_path("$kallistoResults/sbatch_bustools");
my $jobPrefix = "bustools_";
my $clusterOutput = "${kallistoResults}/clusterOutput_bustools/";
my $sbatchLocation = "${kallistoResults}/sbatch_bustools/";
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
        next if -e "$kallistoResults/$libraryId/done_bustools";
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
        $sbatchTemplate .= "module add UHTS/Analysis/bustools/0.40.0;\n";
        my $commandToRun = "R CMD BATCH --no-save --no-restore \'--args libraryId=\"$libraryId\"".
            " whiteList_Path=\"$whiteListPath\" folder_gtf=\"$gtfDir\" scRNASeq_Info=\"$scRNASeqInfoFile\"".
            " kallisto_bus_results=\"$kallistoResults\"\' $pathToScript ${clusterOutput}/${jobName}.Rout";
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

# tx2gene files have to be used by the scripts. As some of those files could be compressed
# we uncompress all tx2gene files prior to running the jobs
system("unxz $gtfDir/*.tx2gene.xz");
my $startTime = localtime->strftime('%Y-%m-%dT%H:%M:%S');
my $numberJobRun = 0;
my $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
foreach my $experimentId (keys %sbatchToRun){
    foreach my $libraryId (keys %{$sbatchToRun{$experimentId}}){
        while ($jobsRunning >= $parallelJobs) {
             sleep(15);
             $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
        }
        system("sbatch $sbatchToRun{$experimentId}{$libraryId} >/dev/null");
        $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
	$numberJobRun++;
    }
}

print "all jobs created properly. Run $numberJobRun jobs\n";

while ($jobsRunning > 0) {
    sleep(15);
    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
}

# now that all jobs have been finished we compress back all tx2gene files
system("xz $gtfDir/*.tx2gene");

print "all jobs finished\n";

my %jobs_status = Utils::count_status_finished_jobs($jobPrefix, $startTime);
    print $jobs_status{"completed"}." jobs completed, ".$jobs_status{"failed"}." jobs failed, ".$jobs_status{"out_of_memory"}." jobs failed with an out of memory issue and ".$jobs_status{"cancelled"}." jobs have been cancelled.\n";

