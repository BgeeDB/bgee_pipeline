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
use Cwd;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;
use File::Path qw(make_path remove_tree);
use File::Basename;

## Define arguments & their default value
my ($metadataFile, $parallelJobs, $outputDir, $bamtofastq, $queue, $account, $downloadedLibraries) = ('', '', '', '', '',  '', '');
my $doNotDownload = (0);
my %opts = ('metadataFile=s'        => \$metadataFile,
            'parallelJobs=s'        => \$parallelJobs,
            'downloadedLibraries=s' => \$downloadedLibraries,
            'outputDir=s'           => \$outputDir,
            'bamtofastq=s'          => \$bamtofastq,
            'queue=s'               => \$queue,
            'account=s'             => \$account,
            'doNotDownload'         => \$doNotDownload
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $parallelJobs eq '' || $outputDir eq '' || $bamtofastq eq '' ||
    $queue eq '' || $account eq '') {
    print "\n\tInvalid or missing argument:
\te.g. $0 -metatadataFile=... -parallelJobs=50 -outputDir=...  >> $@.tmp 2> $@.warn
\t-metadataFile            file containing metadata necessary to download each run
\t-parallelJobs            maximum number of jobs to run in parallel
\t-downloadedLibraries     file containing the ID of all alreaydy downloaded libraries for single cell
\t-outputDir               directory where FASTQ files are downloaded/generated
\t-bamtofastq              directory where the bamtofastq tool from 10X is installed
\t-queue                   queue to use to run jobs on the cluster
\t-account                 account to use to run jobs on the cluster
\t-doNotDownload           (optional) option used to check libraries that have to be downloaded without downloading them.It generates symlink for libraries
                           properly downloaded and add them to the file listing already downloaded libraries. This option is useful if the script has been killed
                           before generating symlink or updating the file listing download libraries. This can happen when the download is too long. The script can
                           then be killed by cluster admins.
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# Info of processed libraries coming from the pipeline
my $isTargetBased = 1;
my %processedLibraries = get_processed_libraries_info($metadataFile, $isTargetBased);

my $experimentDirName = "EXPERIMENTS";
my $experimentOutputDir = "$outputDir/$experimentDirName";

# Read already downloaded libraries
my %alreadyDownloaded = map { $_ => 1 } read_file("$downloadedLibraries", chomp=>1);

my %sbatchToRun = ();

## create directory
my $sbatchDir = "$outputDir/sbatch";
my $clusterOutputDir = "$outputDir/clusterOutput";
make_path("$sbatchDir");
make_path("$clusterOutputDir");

#store initial dir location to be able to move for symlink generation and then come back later
my $initialDir = getcwd;

my $jobPrefix = "download_";
my $jobs_created = 0;
## first create sbatch files and add them to an array of sbatch to run
foreach my $experimentId (keys %processedLibraries){
    foreach my $libraryId (keys %{$processedLibraries{$experimentId}}){
        my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
        my $libDirectory = "$outputDir/$speciesId/$libraryId";
        next if ( exists $alreadyDownloaded{$libraryId} );
        foreach my $runId (keys %{$processedLibraries{$experimentId}{$libraryId}}){
            next if $runId eq "speciesId";
            my $runDirectory = "$libDirectory/$runId";
            if(! $doNotDownload) {
                next if (-f "$runDirectory/done");
            }
            $jobs_created++;
            #create sbatch file and
            my $source = $processedLibraries{$experimentId}{$libraryId}{$runId}{'downloadSource'};
            my $jobName = "$jobPrefix$runId";
            ## Use 4Gb of memory. Should maybe be increase depending on the run to download
            ## ask for 4 cpus as it is the number of threads used by the bamtofastq tool
            my $sbatchTemplate = Utils::sbatch_template($queue, $account, 4,
              50, "$clusterOutputDir/$jobName.out", "$clusterOutputDir/$jobName.err",
              $jobName);
            $sbatchTemplate .= "export bamtofastq=$bamtofastq\n";
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
                        if ($bamInfo =~ m/_I1_/ || $bamInfo =~ m/_I1.f/) {
                            $sbatchTemplate .= "mv $runDirectory/".basename($bamInfo)." $runDirectory/FASTQ/".$runId."_I1.fastq.gz &&\n";
                        } elsif ($bamInfo =~ m/_R1_/ || $bamInfo =~ m/_R1.f/) {
                            $sbatchTemplate .= "mv $runDirectory/".basename($bamInfo)." $runDirectory/FASTQ/".$runId."_R1.fastq.gz &&\n";
                        } elsif ($bamInfo =~ m/_R2_/ || $bamInfo =~ m/_R2.f/) {
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
            my $sbatchFilePath = "$sbatchDir/$jobName.sbatch";
            $sbatchToRun{$experimentId}{$libraryId}{$runId} = $sbatchFilePath;
            $sbatchToRun{$experimentId}{$libraryId}{'speciesId'} = $speciesId;
            open(FH, '>', $sbatchFilePath) or die $!;
            print FH $sbatchTemplate;
            close(FH);
        }
    }
}

print "created $jobs_created sbatch files.\n";

my $numberJobRun = 0;

if(! $doNotDownload) {
    my $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
    foreach my $experimentId (keys %sbatchToRun){
        foreach my $libraryId (keys %{$sbatchToRun{$experimentId}}){
            my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
            my $libDirectory = "$outputDir/$speciesId/$libraryId";
            foreach my $runId (keys %{$sbatchToRun{$experimentId}{$libraryId}}){
                my $runDirectory = "$libDirectory/$runId";
                next if (-f "$runDirectory/done");
                next if $runId eq "speciesId";
                $numberJobRun++;
                $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
                while ($jobsRunning >= $parallelJobs) {
                        sleep(15);
                        $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
                    }
                make_path("$runDirectory/FASTQ");
                chdir "$libDirectory";
                system("sbatch $sbatchToRun{$experimentId}{$libraryId}{$runId}");
            }
        }
    }

    while ($jobsRunning > 0) {
        sleep(15);
        $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
    }
}

print "all download finished properly. Run $numberJobRun jobs\n";




# wait 20 seconds to be sure .done files had enough time to be created
sleep(20);
print "now start to generate symlinks of libraries per species\n";
foreach my $experimentId (keys %sbatchToRun){
    foreach my $libraryId (keys %{$sbatchToRun{$experimentId}}){
        my $speciesId = $sbatchToRun{$experimentId}{$libraryId}{'speciesId'};
        my $libDirectory = "$outputDir/$speciesId/$libraryId";
        my $done = 1;
        foreach my $runId (keys %{$sbatchToRun{$experimentId}{$libraryId}}){
            next if $runId eq "speciesId";
            if(! -e "$libDirectory/$runId/done") {
                $done = 0;
            }
        }
        if ($done) {
            # if all the run of this library have been properly downloaded then
            # we generate a symlink
            print "done download of library $libraryId\n";
            my $currentExpDir = "$experimentOutputDir/$experimentId";
            make_path("$currentExpDir");
            chdir "$currentExpDir";
            system("ln -s ../../$speciesId/$libraryId $libraryId");
            # we also add this library to the list of already generated libraries
            $alreadyDownloaded{$libraryId} = 1;
        } else {
            warn "Did not properly download the library $libraryId";
            remove_tree("$libDirectory");
        }
    }
}

#now go back to original location to update file listing all downloaded libraries
chdir "$initialDir";
print "Finally update the file containing downloaded libraries";
open my $outFh, "> ", "$downloadedLibraries" or die "Cannot write: $! \n";
foreach my $libraryId (sort keys %alreadyDownloaded) {
    print $outFh "$libraryId\n";
}
close $outFh;