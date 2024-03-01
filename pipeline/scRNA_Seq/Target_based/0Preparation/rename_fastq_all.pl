#!/usr/bin/env perl

# This script is not integrated in the pipeline is has been used to solve a renaming issue
# after downloading fastq files (see https://www.biostars.org/p/9588563/#9588575). This problem
# should not happen again so this script could potentially be deleted.
# Kept it for now...

# example of command to run it :
# nohup perl 0Preparation/rename_fastq_all.pl -metadataFile=../../../generated_files/scRNA_Seq/Target_based/metadata_info_10X.txt -excludedLibraries=../../../generated_files/scRNA_Seq/Target_based/librariesExcluded.tsv -whitelistFolder=../../../source_files/scRNA_Seq/Target_based/ -fastqDir=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/FASTQ/scRNAseq/10X/ -queue=urblauna -account=mrobinso_bgee_sensitive &
# 

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;
use File::Path qw(make_path remove_tree);
use Time::Piece;

## Define arguments & their default value
my ($metadataFile, $excludedLibraries, $fastqDir, $whitelistFolder, $queue, $account) = ('', '', '', '', '', '');
my ($doNotDownload)          = (0);
my %opts = ('metadataFile=s'        => \$metadataFile,
            'excludedLibraries=s'   => \$excludedLibraries,
            'fastqDir=s'            => \$fastqDir,
            'whitelistFolder=s'     => \$whitelistFolder,
            'queue=s'               => \$queue,
            'account=s'             => \$account
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $fastqDir eq '' || $excludedLibraries eq '' ||
    $whitelistFolder eq '' || $queue eq '' || $account eq '') {
    print "\n\tInvalid or missing argument:
\t-metadataFile            path to the rna_seq_sample_info.txt file
\t-excludedLibraries       file containing the ID of all libraries not to download
\t-whitelistFolder         path to directory containing 10X whitelists
\t-fastqDir                parent directory containing all FASTQ files
\t-queue                   queue to use to run jobs on the cluster
\t-account                 account to use to run jobs on the cluster
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# Info of processed libraries coming from the pipeline
my $isTargetBased = 1;
my %processedLibraries = get_processed_libraries_info($metadataFile, $isTargetBased);

## create directory
my $jobPrefix = "rename_";
my $sbatchDir = "$fastqDir/${jobPrefix}sbatch";
my $clusterOutputDir = "$fastqDir/${jobPrefix}clusterOutput";
make_path("$sbatchDir");
make_path("$clusterOutputDir");

my $jobs_created = 0;
my %sbatchToRun = ();
foreach my $experimentId (keys %processedLibraries){
    foreach my $libraryId (keys %{$processedLibraries{$experimentId}}){
        my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
        my $libDirectory = "$fastqDir/$speciesId/$libraryId";
        foreach my $runId (keys %{$processedLibraries{$experimentId}{$libraryId}{'runIds'}}){
            my $fastqRunDir = "$libDirectory/$runId/FASTQ/";
            next if ! -e $fastqRunDir;
            my $jobName = "$jobPrefix$runId";
            my $sbatchTemplate = Utils::sbatch_template($queue, $account, 1,
              5, "$clusterOutputDir/$jobName.out", "$clusterOutputDir/$jobName.err", $jobName);
            $sbatchTemplate .= "module use /software/module/;\nmodule load R/3.6.1;\nexport PATH=/software/bin:\$PATH;\n\n";

            $sbatchTemplate .=  "Rscript /users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Target_based/0Preparation/rename_fastq.R fastqPath=\\\"${fastqRunDir}\\\" runId=\\\"$runId\\\" renaming=\\\"fastqdump\\\" whitelistFolder=\\\"$whitelistFolder\\\" || exit 1;";
            $jobs_created++;
            ## create sbatch file and add its path to the hash of sbatch files
            my $sbatchFilePath = "$sbatchDir/$jobName.sbatch";
            $sbatchToRun{$experimentId}{$libraryId}{"runIds"}{$runId} = $sbatchFilePath;
            $sbatchToRun{$experimentId}{$libraryId}{'speciesId'} = $speciesId;
            open(FH, '>', $sbatchFilePath) or die $!;
            print FH $sbatchTemplate;
            close(FH);
        }
    }
}

print "created $jobs_created sbatch files.\n";
my $parallelJobs = 100;
# if jobs had to be run
if ($jobs_created > 0) {
    
    my $numberJobRun = 0;
    my $startTime = localtime->strftime('%Y-%m-%dT%H:%M:%S');

    my $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
    foreach my $experimentId (keys %sbatchToRun){
        foreach my $libraryId (keys %{$sbatchToRun{$experimentId}}){
            my $speciesId = $processedLibraries{$experimentId}{$libraryId}{'speciesId'};
            my $libDirectory = "$fastqDir/$speciesId/$libraryId";
            foreach my $runId (keys %{$sbatchToRun{$experimentId}{$libraryId}{"runIds"}}){
                $numberJobRun++;
                while ($jobsRunning >= $parallelJobs) {
                    sleep(15);
                    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
                }
                system("sbatch ".$sbatchToRun{$experimentId}{$libraryId}{"runIds"}{$runId});
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
