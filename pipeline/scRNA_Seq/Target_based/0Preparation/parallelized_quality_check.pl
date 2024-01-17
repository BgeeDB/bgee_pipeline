#!/usr/bin/env perl

## This script allows to download FASTQ files from SRA

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Try::Tiny;
use FindBin;
use Cwd;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;
use Time::Piece;
use File::Path qw(make_path remove_tree);
use File::Basename;
use Cpanel::JSON::XS;
use LWP::Simple;

## Define arguments & their default value
my ($metadataFile, $parallelJobs, $outputDir, $queue, $account) = ('', '', '', '', '');
my $qualityScriptPath = '';
my %opts = ('metadataFile=s'        => \$metadataFile,
            'parallelJobs=s'        => \$parallelJobs,
            'outputDir=s'           => \$outputDir,
            'queue=s'               => \$queue,
            'account=s'             => \$account,
	    'qualityScriptPath=s'   => \$qualityScriptPath
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $parallelJobs eq '' || $outputDir eq '' ||
    $queue eq '' || $account eq '' || $qualityScriptPath eq '') {
    print "\n\tInvalid or missing argument:
\t-metadataFile            file containing metadata necessary to download each run
\t-parallelJobs            maximum number of jobs to run in parallel
\t-outputDir               directory where FASTQ files are downloaded/generated
\t-queue                   queue to use to run jobs on the cluster
\t-account                 account to use to run jobs on the cluster
\t-qualityScriptPath       path to the script allowing to run fastp and to generate .R.stat file
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# Info of processed libraries coming from the pipeline
my $isTargetBased = 1;
my %processedLibraries = get_processed_libraries_info($metadataFile, $isTargetBased);

my %sbatchToRun = ();

## create directory
my $jobPrefix = "quality_";
my $sbatchDir = "$outputDir/${jobPrefix}sbatch";
my $clusterOutputDir = "$outputDir/${jobPrefix}clusterOutput";
make_path("$sbatchDir");
make_path("$clusterOutputDir");

my $jobs_created = 0;
## first create sbatch files and add them to an array of sbatch to run
foreach my $experimentId (keys %processedLibraries){
    foreach my $libraryId (keys %{$processedLibraries{$experimentId}}){
        my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
        my $libDirectory = "$outputDir/$speciesId/$libraryId";
        foreach my $runId (keys %{$processedLibraries{$experimentId}{$libraryId}}){
            next if $runId eq "speciesId";
            my $runDirectory = "$libDirectory/$runId";
            #create sbatch file and
            my $jobName = "$jobPrefix$runId";
            ## Use 30Gb of memory. Should maybe be increase depending on the run
            my $sbatchTemplate = Utils::sbatch_template($queue, $account, 1,
              30, "$clusterOutputDir/$jobName.out", "$clusterOutputDir/$jobName.err",
              $jobName);
	    $sbatchTemplate .= "module load gcc/10.4.0;\nmodule use /software/module/;\nmodule load R/3.6.1;\nmodule load fastp/0.23.2;\nexport PATH=/software/bin:\$PATH;\n";
	    $sbatchTemplate .= "perl $qualityScriptPath -run_path=$runDirectory -run_id=$runId\n";
            $jobs_created++;
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

# if jobs had to be run
if ($jobs_created > 0) {
    
    my $numberJobRun = 0;
    my $startTime = localtime->strftime('%Y-%m-%dT%H:%M:%S');

    my $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
    foreach my $experimentId (keys %sbatchToRun){
        foreach my $libraryId (keys %{$sbatchToRun{$experimentId}}){
            my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
            my $libDirectory = "$outputDir/$speciesId/$libraryId";
            foreach my $runId (keys %{$sbatchToRun{$experimentId}{$libraryId}}){
                next if $runId eq "speciesId";
		my $runDirectory = "$libDirectory/$runId";
                $numberJobRun++;
                while ($jobsRunning >= $parallelJobs) {
                    sleep(15);
                    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
                }
                system("sbatch $sbatchToRun{$experimentId}{$libraryId}{$runId}");
                $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
            }
        }
    }

    while ($jobsRunning > 0) {
        sleep(15);
        $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
    }
 
    print "all download jobs finished. Run $numberJobRun jobs.\n";

    my %jobs_status = Utils::count_status_finished_jobs($jobPrefix, $startTime);
    print $jobs_status{"completed"}." jobs completed, ".$jobs_status{"failed"}." jobs failed, ".$jobs_status{"out_of_memory"}." jobs failed with an out of memory issue and ".$jobs_status{"cancelled"}." jobs have been cancelled.\n";
}
