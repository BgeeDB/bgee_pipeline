#!/usr/bin/env perl

## This script allows to download FASTQ files from SRA

##XXX What is the best approach to download FASTQ files from SRA (using sra-toolkit)?
##    - use fastq-dump or fasterq-dump to directly download FASTQ files
##    - use prefetch to download BAM files and then use bamtofastq from 10X Genomics to convert to FASTQ
##    - use a mix of these of these approach depending on metatadata retrieved from SRA
##         - For some libraries 2 FASTQ files are available but for some others (initially submitted on ENA),
##           only one FASTQ file is available. It is then necessary to use the BAM file to generate the FASTQs
##         - have to find a robust criteria to know which approach to use
##         TODO Find a robust criteria

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

## Define arguments & their default value
my ($metatadataFile, $parallelJobs, $outputDir) = ('', '',  '');
my %opts = ('metatadataFile=s'      => \$metatadataFile,
            'parallelJobs=s'        => \$parallelJobs,
            'outputDir=s'           => \$outputDir
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metatadataFile || $parallelJobs eq '' || $outputDir eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0 -metatadataFile=... -parallelJobs=50 -outputDir=...  >> $@.tmp 2> $@.warn
\t-metadataFile            file containing metadata necessary to download each run
\t-parallelJobs            maximum number of jobs to run in parallel
\t-outputDir               directory where FASTQ files are downloaded/generated
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../target_base_utils.pl");

$queue = ;
$account = ;

# Info of processed libraries coming from the pipeline
my %processedLibraries = get_processed_libraries_info($metadataFile);

my @sbatchToRun = ();

my $jobPrefix = "downloadSRA_";
## first create sbatch files and add them to an array of sbatch to run
foreach my $experimentId (keys %processedLibraries){
    foreach my $libraryId (keys %processedLibraries{$experimentId}){
        foreach my $runId (keys %processedLibraries{$experimentId}{$libraryId}){
            #create sbatch file and
            my $jobName = "$jobPrefix$runId";
            my $sbtachTemplate = sbatch_template($queue, $account, 1,
                1, "$outputDir/$jobName.out", "$outputDir/$jobName.err",
                $jobName);
            ## load SRA sra-toolkit
            $sbatchTemplate .= "module load gcc/10.4.0;\n
                                 module load sratoolkit/3.0.0;\n\n";
            my $submittedFtp = $processedLibraries{$experimentId}{$libraryId}{$runId}{$submittedFTP};
            ## download fastq from SRA
            if ($submittedFtp eq 'NA') {
                $sbatchTemplate .= "fastq-dump --outdir $outputDir --split-files $runId\n";
            }
            ## download BAM from SRA
            } else {
                $sbatchTemplate .= "prefetch --type bam --max-size 9999999999 -O $outputDir $runId\n";
            }
            ## create sbatch file and add its path to the hash of sbatch files
            my $sbatchFilePath = "$outputDir/sbatch/$jobName.sbatch"
            push (@sbatchToRun, $sbatchFilePath);
            open(FH, '>', $sbatchFile) or die $!;
            print FH $sbatchTemplate;
            close(FH);
            push (@sbatchToRun, $sbatchFile);
        }
    }
}

my jobsRunning = check_active_jobs_number($jobPrefix);
for (my $sbatchFile in @sbatchToRun) {
    jobsRunning = check_active_jobs_number($jobPrefix);
    while ($jobsRunning >= $parallelJobs) {
        sleep(15);
        jobsRunning = check_active_jobs_number($jobPrefix);
    }
    system("sbatch $sbatchFile");
}

while ($jobsRunning > 0) {
    sleep(15);
   jobsRunning = check_active_jobs_number($jobPrefix);
}

print "all download finished properly";



# attributes:
# - metatadata_file
# - number of parallel jobs
# - output dir where raw reads are downloaded

# steps
# - read metadata file
# - check if FASTQ or BAM files have to be downloaded
#    - run corresponding download
# - check number of job running
#    - if more than max of parallel download then wait
#    - wait until en of last job to end the script

