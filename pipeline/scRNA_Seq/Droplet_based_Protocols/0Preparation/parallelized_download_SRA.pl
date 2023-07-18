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
use Try::Tiny;
use FindBin;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;
use File::Path qw(make_path);
use File::Basename;

## Define arguments & their default value
my ($metadataFile, $parallelJobs, $outputDir, $bamtofastq, $queue, $account) = ('', '', '', '', '',  '');
my %opts = ('metadataFile=s'        => \$metadataFile,
            'parallelJobs=s'        => \$parallelJobs,
            'outputDir=s'           => \$outputDir,
        'bamtofastq=s'          => \$bamtofastq,
        'queue=s'               => \$queue,
        'account=s'             => \$account
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $parallelJobs eq '' || $outputDir eq '' || $bamtofastq eq '' ||
    $queue eq '' || $account eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0 -metatadataFile=... -parallelJobs=50 -outputDir=...  >> $@.tmp 2> $@.warn
\t-metadataFile            file containing metadata necessary to download each run
\t-parallelJobs            maximum number of jobs to run in parallel
\t-outputDir               directory where FASTQ files are downloaded/generated
\t-bamtofastq              directory where the bamtofastq tool from 10X is installed
\t-queue                   queue to use to run jobs on the cluster
\t-account                 account to use to run jobs on the cluster
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
make_path("$outputDir/clusterOutput");
my $jobPrefix = "download_";
my $jobs_created = 0;
## first create sbatch files and add them to an array of sbatch to run
foreach my $experimentId (keys %processedLibraries){
    foreach my $libraryId (keys %{$processedLibraries{$experimentId}}){
        foreach my $runId (keys %{$processedLibraries{$experimentId}{$libraryId}}){
            next if $runId eq "speciesId";
            $jobs_created++;
            #create sbatch file and
            my $source = $processedLibraries{$experimentId}{$libraryId}{$runId}{'downloadSource'};
            my $jobName = "$jobPrefix$runId";
            ## Use 4Gb of memory. Should maybe be increase depending on the run to download
            ## ask for 4 cpus as it is the number of threads used by the bamtofastq tool
            my $sbatchTemplate = Utils::sbatch_template($queue, $account, 4,
              4, "$outputDir/cluterOutput/$jobName.out", "$outputDir/clusterOutput/$jobName.err",
              $jobName);
            $sbatchTemplate .= "export bamtofastq=$bamtofastq\n";
            my $runDirectory = "$outputDir/$experimentId/$libraryId/$runId";
            my $submittedFtp = $processedLibraries{$experimentId}{$libraryId}{$runId}{'submittedFTP'};
            my @bamInfos = split(";", $submittedFtp);
            if ($source eq "SRA") {
                ## load SRA sra-toolkit
                $sbatchTemplate .= "module load gcc/10.4.0;\nmodule load sratoolkit/3.0.0;\n\n";
                ## download fastq from SRA
                if ($submittedFtp eq 'NA') {
                    $sbatchTemplate .= "fastq-dump --outdir $runDirectory/FASTQ --split-files $runId &&\n";
                   ## download BAM from SRA
                } else {
                    $sbatchTemplate .= "{\nprefetch --type bam --max-size 9999999999 -O $runDirectory $runId ||\n";
                    $sbatchTemplate .= "wget --no-verbose --directory-prefix=$runDirectory $bamInfos[0]\n} &&\n";
                    $sbatchTemplate .= "if [ -d $runDirectory/FASTQ ]; then rm -rf $runDirectory/FASTQ; fi &&\n";
                    $sbatchTemplate .= "$bamtofastq --reads-per-fastq=90000000000 --nthreads=4 $runDirectory/".basename($bamInfos[0])." $runDirectory/FASTQ &&\n";
                    $sbatchTemplate .= "find $runDirectory -name '*.fastq.gz' -o -name '*.fq.gz' -exec mv -t $runDirectory/FASTQ {} + &&\n";
                    $sbatchTemplate .= "rm $runDirectory/".basename($bamInfos[0])." &&\n";
                }
            } elsif ($source eq "EBI") {
                # if have to download fastq files
                if ($bamInfos[0] =~ m/\.fastq\.gz$/) {
                    foreach my $bamInfo (@bamInfos) {
                        #download
                        $sbatchTemplate .= "wget --no-verbose --directory-prefix=$runDirectory $bamInfo &&\n";
                        if ($bamInfo =~ m/_I1_/) {
                            $sbatchTemplate .= "mv $runDirectory/".basename($bamInfo)." $runDirectory/FASTQ/".$runId."_I1.fastq.gz &&\n";
                        } elsif ($bamInfo =~ m/_R1_/) {
                            $sbatchTemplate .= "mv $runDirectory/".basename($bamInfo)." $runDirectory/FASTQ/".$runId."_R1.fastq.gz &&\n";
                        } elsif ($bamInfo =~ m/_R2_/) {
                            $sbatchTemplate .= "mv $runDirectory/".basename($bamInfo)." $runDirectory/FASTQ/".$runId."_R2.fastq.gz &&\n";
                        } else {
                            warn "did not manage to rename the file $bamInfo";
                        }
                    }
                } elsif ($bamInfos[0] =~ m/\.bam$/) {
                    $sbatchTemplate .= "wget --no-verbose --directory-prefix=$runDirectory $bamInfos[0] &&\n";
                    $sbatchTemplate .= "if [ -d $runDirectory/FASTQ ]; then rm -rf $runDirectory/FASTQ; fi &&\n";
                    $sbatchTemplate .= "$bamtofastq --reads-per-fastq=90000000000 --nthreads=4 $runDirectory/".basename($bamInfos[0])." $runDirectory/FASTQ &&\n";
                    $sbatchTemplate .= "find $runDirectory -name '*.fastq.gz' -o -name '*.fq.gz' -exec mv -t $runDirectory/FASTQ {} + &&\n";
                    $sbatchTemplate .= "rm $runDirectory/".basename($bamInfos[0])." &&\n";
                } else {
                    warn "unrecognized file extension for $bamInfos[0]";
                }
            } elsif ($source eq "FCA") {
                $sbatchTemplate .= "wget --no-verbose --directory-prefix=$runDirectory $bamInfos[0] &&\n";
                $sbatchTemplate .= "mv $runDirectory/".basename($bamInfos[0])." $runDirectory/FASTQ/${runId}_R1.fastq.gz &&\n";
                $sbatchTemplate .= "wget --no-verbose --directory-prefix=$runDirectory $bamInfos[1] &&\n";
                $sbatchTemplate .= "mv $runDirectory/".basename($bamInfos[1])." $runDirectory/FASTQ/${runId}_R2.fastq.gz &&\n";
            }
            $sbatchTemplate .= "find $runDirectory \\( -name '*.fastq' \\) -exec gzip {} \\; &&\n";
            $sbatchTemplate .= "touch $runDirectory/done";
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
        foreach my $runId (keys %{$sbatchToRun{$experimentId}{$libraryId}}){
            my $runDirectory = "$outputDir/$experimentId/$libraryId/$runId";
            next if (-f "$runDirectory/done");
            $numberJobRun++;
            $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
            while ($jobsRunning >= $parallelJobs) {
                    sleep(15);
                    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
                }
            make_path("$runDirectory/FASTQ");
            chdir "$outputDir/$experimentId/$libraryId";
            system("sbatch $sbatchToRun{$experimentId}{$libraryId}{$runId}");
        }
    }
}

while ($jobsRunning > 0) {
    sleep(15);
    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
}

print "all download finished properly. Run $numberJobRun jobs\n";
##TODO add IDs of downloaded run/libraries to a file in order not to redownload them later