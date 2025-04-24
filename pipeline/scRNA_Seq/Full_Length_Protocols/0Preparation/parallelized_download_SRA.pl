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
use File::Path qw(make_path remove_tree);
use File::Basename;

## Define arguments & their default value
my ($metadataFile, $parallelJobs, $downloadedLibraries, $outputDir, $queue, $account) = ('', '', '', '', '',  '');
my ($excludedLibraries, $loadApptainerModule, $containerCmd, $speciesIds, $doNotDownload) = ('', '', '', '', 0);
my %opts = ('metadataFile=s'        => \$metadataFile,
            'parallelJobs=s'        => \$parallelJobs,
            'downloadedLibraries=s' => \$downloadedLibraries,
            'excludedLibraries=s'   => \$excludedLibraries,
            'outputDir=s'           => \$outputDir,
            'loadApptainerModule=s' => \$loadApptainerModule,
            'containerCmd=s'        => \$containerCmd,
            'queue=s'               => \$queue,
            'account=s'             => \$account,
            'speciesIds=s'          => \$speciesIds,
            'doNotDownload'         => \$doNotDownload
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $parallelJobs eq '' || $outputDir eq '' || $downloadedLibraries eq '' ||
    $queue eq '' || $account eq '' || $loadApptainerModule eq '' || $containerCmd eq '' || $excludedLibraries eq '' ) {
    print "\n\tInvalid or missing argument:
\te.g. $0 -metatadataFile=... -parallelJobs=50 -outputDir=...  >> $@.tmp 2> $@.warn
\t-metadataFile            file containing metadata necessary to download each run
\t-parallelJobs            maximum number of jobs to run in parallel
\t-downloadedLibraries     file containing the ID of all alreaydy downloaded libraries for single cell
\t-excludedLibraries       file containing the ID of all libraries not to download
\t-outputDir               directory where FASTQ files are downloaded/generated
\t-loadApptainerModule     command loading all modules necessary to run apptainer
\t-containerCmd            command to load the Bgee container with apptainer and provide all necessary bindings
\t-queue                   queue to use to run jobs on the cluster
\t-account                 account to use to run jobs on the cluster
\t-speciesIds              (optional) list of species for which FASTQ files have to be downloaded. It has to be
                           formatted as speciesIds separated by a comma. If no speciesId is provided then all
                           libraries are downloaded.
\t-doNotDownload           (optional) option used to check libraries that have to be downloaded without downloading them.
                           It generates symlink for libraries properly downloaded and add them to the file listing already
                           downloaded libraries. This option is useful if the script has been killed before generating
                           symlink or updating the file listing download libraries. This can happen when the download is
                           too long. The script can then be killed by cluster admins.
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# Info of processed libraries coming from the pipeline
my $isTargetBased = 0;
my %processedLibraries = get_processed_libraries_info($metadataFile, $isTargetBased);

# retrieve speciesIds to download from script argument
my @speciesIdsToDownload = split(',', $speciesIds);

# Read already downloaded libraries
my %alreadyDownloaded = map { $_ => 1 } read_file("$downloadedLibraries", chomp=>1);

# Read excluded libraries
open(my $excluded, $excludedLibraries) || die "failed to read sample excluded file: $!";
my @excluded_libraries;
while (my $line = <$excluded>) {
    chomp $line;
     ## skip comment lines
    next  if ( ($line =~ m/^#/) or ($line =~ m/^\"#/) );
    my @line = split(/\t/, $line);
    if ($line[1] eq "TRUE") {
        push(@excluded_libraries, $line[0])
    }
}

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
        my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
        next if ( ! grep(/$speciesId/, @speciesIdsToDownload)  && $speciesIds ne '');
        next if ( exists $alreadyDownloaded{$libraryId} );
        next if ( grep( /^$libraryId$/, @excluded_libraries));
        my $libDirectory = "$outputDir/$speciesId/$libraryId";
        foreach my $runId (keys %{$processedLibraries{$experimentId}{$libraryId}{"runIds"}}){
            if ( $runId =~ /^[SEDC]RR\d+/ ){ #S: SRA/NCBI; E: EBI; D: DDBJ; C: GSA_China
                if(! $doNotDownload) {
                    next if (-f "$libDirectory/$runId.done");
                }
                $jobs_created++;
                #create sbatch file and
                my $jobName = "$jobPrefix$runId";
                ## Use 4Gb of memory. Should maybe be increase depending on the run to download
                ## ask for 4 cpus as it is the number of threads used by the bamtofastq tool
                my $sbatchTemplate = Utils::sbatch_template($queue, $account, 1,
                4, "$outputDir/clusterOutput/$jobName.out", "$outputDir/clusterOutput/$jobName.err",
                $jobName);
                $sbatchTemplate .= $loadApptainerModule =~ s/; /;\n/r;
                $sbatchTemplate .= "\n\n";
                my $libDirectory = "$experimentOutputDir/$experimentId/$libraryId";
                ## download fastq from SRA split into several FASTQ files if paired library
                my $prefix      = "$libDirectory/$runId";
                my $fastq_fastp = '';
                my $fastq_R     = '';
                if($processedLibraries{$experimentId}{$libraryId}{$runId}{'libraryType'} eq "single") {
                    $sbatchTemplate .= "$containerCmd fastq-dump --outdir $libDirectory --gzip $runId &&\n";
                    $fastq_fastp = "$prefix.fastq.gz";
                    $fastq_R     = $fastq_fastp;
                } elsif ($processedLibraries{$experimentId}{$libraryId}{$runId}{'libraryType'} eq "paired") {
                    $sbatchTemplate .= "$containerCmd fastq-dump --outdir $libDirectory --split-files --gzip $runId &&\n";
                    $fastq_fastp = "${prefix}_1.fastq.gz -I ${prefix}_2.fastq.gz";
                    $fastq_R     = "${prefix}_1.fastq.gz    ${prefix}_2.fastq.gz";
                } else {
                    warn "unrecognized type of library for $libraryId\n";
                }
                if ( !-e "$prefix.fastp.html.xz" || !-e "$prefix.fastp.json.xz" ){
                    $sbatchTemplate .= "$containerCmd fastp -i $fastq_fastp --json $prefix.fastp.json --html $prefix.fastp.html  > $prefix.fastp.log 2>$prefix.fastp.log &&\n";
                    $sbatchTemplate .= "xz -9 $prefix.fastp.html $prefix.fastp.json &&\n";
                }
                if ( !-e "$prefix.R.stat" ){
                    $sbatchTemplate .= "$containerCmd  /bin/echo \"#min\tmax\tmedian\tmean\" > $prefix.R.stat &&\n";
                    #NOTE for cases like SRX1372530 with paired-end files coming with a single-end file in the same run, use ${prefix}*.fastq.gz ???
                    $sbatchTemplate .= "zcat $fastq_R | sed -n '2~4p' | awk '{print length(\$0)}' | $containerCmd Rscript -e 'd<-scan(\"stdin\", quiet=TRUE);cat(min(d), max(d), median(d), mean(d), sep=\"\\t\");cat(\"\\n\")' >> $prefix.R.stat &&\n";
                }
                $sbatchTemplate .= "touch $libDirectory/$runId.done";

                ## create sbatch file and add its path to the hash of sbatch files
                my $sbatchFilePath = "$outputDir/sbatch/$jobName.sbatch";
                $sbatchToRun{$experimentId}{$libraryId}{$runId} = $sbatchFilePath;
                $sbatchToRun{$experimentId}{$libraryId}{'speciesId'} = $speciesId;
                open(FH, '>', $sbatchFilePath) or die $!;
                print FH $sbatchTemplate;
                close(FH);
            } else {
                print "runId $runId does not start with a valid prefix. It should start with [SEDC]RR.\n";
                system("rm -f $libDirectory/sra/$runId.sra*");
            }
        }
    }
}
print "created $jobs_created sbatch files.\n";

my $numberJobRun = 0;
if(! $doNotDownload) {
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
}
print "all download finished properly. Run $numberJobRun jobs\n";
# wait 20 seconds to be sure .done files had enough time to be created
sleep(20);
print "now start to generate symlinks of libraries per species\n";
foreach my $experimentId (keys %sbatchToRun){
    foreach my $libraryId (keys %{$sbatchToRun{$experimentId}}){
        my $speciesId = $sbatchToRun{$experimentId}{$libraryId}{'speciesId'};
        my $libDirectory = "$experimentOutputDir/$experimentId/$libraryId";
        my $done = 1;
        foreach my $runId (keys %{$sbatchToRun{$experimentId}{$libraryId}}){
            next if $runId eq "speciesId";
            if(! -e "$libDirectory/$runId.done") {
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
            remove_tree($libDirectory);
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