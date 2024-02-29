#!/usr/bin/env perl

## This script allows to download FASTQ files from SRA. It also creates symlink in the EXPERIMENTS folder.

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Try::Tiny;
use FindBin;
use Cwd;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;
use File::Path qw(make_path remove_tree);
use File::Basename;

## Define arguments & their default value
my ($metadataFile, $parallelJobs, $excludedLibraries, $downloadedLibraries, $encryptFile, $outputDir, $queue, $account, $speciesIds) = ('', '', '', '', '', '',  '', '', '');
my ($doNotDownload)          = (0);
my %opts = ('metadataFile=s'        => \$metadataFile,
            'parallelJobs=s'        => \$parallelJobs,
            'downloadedLibraries=s' => \$downloadedLibraries,
            'excludedLibraries=s'   => \$excludedLibraries,
            'encryptFile=s'         => \$encryptFile,
            'outputDir=s'           => \$outputDir,
            'queue=s'               => \$queue,
            'account=s'             => \$account,
            'speciesIds=s'          => \$speciesIds,
            'doNotDownload'         => \$doNotDownload
           );


######################## Check arguments ########################
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$metadataFile || $parallelJobs eq '' || $outputDir eq '' || $downloadedLibraries eq '' ||
    $queue eq '' || $account eq '' || $encryptFile eq '') {
    print "\n\tInvalid or missing argument:
\t-metadataFile            path to the rna_seq_sample_info.txt file
\t-parallelJobs            maximum number of jobs to run in parallel
\t-downloadedLibraries     file containing the ID of all alreaydy downloaded libraries
\t-excludedLibraries       file containing the ID of all libraries not to download
\t-outputDir               directory where FASTQ files are downloaded/generated
\t-encryptFile             path to the encryption file for secure libraries (e.g GTEx)
\t-queue                   queue to use to run jobs on the cluster
\t-account                 account to use to run jobs on the cluster
\t-speciesIds              (optional) list of species for which FASTQ files have to be downloaded. It has to be formatted as speciesIds separated by a comma.
                           If no speciesId is provided then all libraries are downloaded
\t-doNotDownload           (optional) option used to check libraries that have to be downloaded without downloading them.It generates symlink for libraries
                           properly downloaded and add them to the file listing already downloaded libraries. This option is useful if the script has been killed
                           before generating symlink or updating the file listing download libraries. This can happen when the download is too long. The script can
                           then be killed by cluster admins.

\n";
    exit 1;
}

require("$FindBin::Bin/../rna_seq_utils.pl");

# Private experimentId to store encrypted
my @private_exp_id     = ('SRP012682'); # E.g. GTEx

my @speciesIdsToDownload = split(',', $speciesIds);
(my $speciesIdsDirSuffix = $speciesIds) =~ s/,/_/g;

print "speciesIds to filter on \"$speciesIds\"\n";

# Info of processed libraries coming from the pipeline
my $isTargetBased = 0;

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

my $sbatchDir = $speciesIds eq '' ? "$outputDir/sbatch_allSpecies" : "$outputDir/sbatch_$speciesIdsDirSuffix";
my $clusterOutputDir = $speciesIds eq '' ? "$outputDir/clusterOutput_allSpecies" : "$outputDir/clusterOutput_$speciesIdsDirSuffix";

## create directory necessary to store sbatch files
make_path($sbatchDir);
make_path($clusterOutputDir);

#store initial dir location to be able to move for symlink generation and then come back later
my $initialDir = getcwd;

my $jobPrefix = "download_";
my $jobs_created = 0;
my $experimentDirName = "EXPERIMENTS";
my $experimentOutputDir = "$outputDir/$experimentDirName";
make_path($experimentOutputDir);
## first create sbatch files and add them to an array of sbatch to run
open(my $ANNOTATION, '<', "$metadataFile")  or die "\n\tCannot read/open [$metadataFile]\n\n";
#libraryId   experimentId   speciesId   organism        genomeFilePath                        database   platform                       libraryType   libraryInfo   readLength   runIds
#SRX081869   GSE30352       9031        Gallus gallus   gallus_gallus/Gallus_gallus.Galgal4   Ensembl    Illumina Genome Analyzer IIx   SINGLE                      76           SRR306710
my $count = 0;
my @missing;
LIB:
while (<$ANNOTATION>){
    chomp $_;
    my @tmp = map { s{^"|"$}{}g; $_ } # remove quotes
              split(/\t/, $_);
    my $sra_list      = $tmp[10];
    my $libraryId     = $tmp[0];
    my $experimentId  = $tmp[1];
    my $speciesId     = $tmp[2];
    my $libraryType   = $tmp[7];
    my $libDirectory = "$outputDir/$speciesId/$libraryId";
    #if speciesIds to download have been provided then consider only libraries ccorresponding to those species
    next LIB if ( ! grep(/$speciesId/, @speciesIdsToDownload)  && $speciesIds ne '');
    # Header or commented library
    next LIB if ( $libraryId =~ /^#/ || $sra_list =~ /^#/ );
    #do not consider already downloaded libraries
    next LIB if ( exists $alreadyDownloaded{$libraryId} );
    #do not consider excluded libraries
    next LIB if ( grep( /^$libraryId$/, @excluded_libraries));
    SRA:
    for my $runId ( sort split(/,/, $sra_list) ){
        if ( $runId =~ /^[SEDC]RR\d+/ ){ #S: SRA/NCBI; E: EBI; D: DDBJ; C: GSA_China
            if(! $doNotDownload) {
                next if (-f "$libDirectory/$runId.done");
            }
            $jobs_created++;
            #create sbatch file and
            my $jobName = "$jobPrefix$runId";
            ## Use 30Gb of memory. Should maybe be increase depending on the run to download
            ## ask for 1 cpu
            my $sbatchTemplate = Utils::sbatch_template($queue, $account, 1,
              30, "$clusterOutputDir/$jobName.out", "$clusterOutputDir/$jobName.err",
              $jobName);
            $sbatchTemplate .= "module load gcc/10.4.0;\nmodule load sratoolkit/3.0.0;\nmodule load fastp/0.23.2\nmodule use /software/module;\nmodule load R/3.6.1;\n\n";
            ## download fastq from SRA split into several FASTQ files if paired library and manage input files of fastp and R.stat
            my $prefix      = "$libDirectory/$runId";
            my $fastq_fastp = '';
            my $fastq_R     = '';
            if($libraryType eq "SINGLE") {
                $sbatchTemplate .= "fastq-dump --outdir $libDirectory --gzip $runId &&\n";
                $fastq_fastp = "$prefix.fastq.gz";
                $fastq_R     = $fastq_fastp;
            } elsif ($libraryType eq "PAIRED"){
                $sbatchTemplate .= "fastq-dump --outdir $libDirectory --split-files --gzip $runId &&\n";
                $fastq_fastp = "${prefix}_1.fastq.gz -I ${prefix}_2.fastq.gz";
                $fastq_R     = "${prefix}_1.fastq.gz    ${prefix}_2.fastq.gz";
            } else {
                warn "unrecognized type of library $libraryType for library $libraryId\n";
            }
            # Run FastP (A quality control tool for high throughput sequence data) for ALL SRR (runs)
            # as well as basic read length statistics with R
            #NOTE Would be nice to have all basic stats from FastP (currently some are done in R)
            if ( !-e "$prefix.fastp.html.xz" || !-e "$prefix.fastp.json.xz" ){
                $sbatchTemplate .= "fastp -i $fastq_fastp --json $prefix.fastp.json --html $prefix.fastp.html  > $prefix.fastp.log 2>$prefix.fastp.log &&\n";
                $sbatchTemplate .= "xz -9 $prefix.fastp.html $prefix.fastp.json &&\n";
            }
            if ( !-e "$prefix.R.stat" ){
                $sbatchTemplate .= "/bin/echo \"#min\tmax\tmedian\tmean\" > $prefix.R.stat &&\n";
                #NOTE for cases like SRX1372530 with paired-end files coming with a single-end file in the same run, use ${prefix}*.fastq.gz ???
                $sbatchTemplate .= "zcat $fastq_R | sed -n '2~4p' | awk '{print length(\$0)}' | Rscript -e 'd<-scan(\"stdin\", quiet=TRUE);cat(min(d), max(d), median(d), mean(d), sep=\"\\t\");cat(\"\\n\")' >> $prefix.R.stat &&\n";
            }
            # If private (need encryption):
            if ( (scalar grep { /^$experimentId$/ } @private_exp_id) >= 1 ){
                for my $fastq ( glob("$libDirectory/*.gz") ){
                    #NOTE Replace -salt by -d for decrypting and gz.enc as input and gz as output
                    $sbatchTemplate .= "openssl enc -aes-128-cbc -salt -in $fastq -out $fastq.enc -pass file:$encryptFile  &&  rm -f $fastq &&\n";
                }
            }
            $sbatchTemplate .= "touch $libDirectory/$runId.done";

            ## create sbatch file and add its path to the hash of sbatch files
            my $sbatchFilePath = "$sbatchDir/$jobName.sbatch";
            $sbatchToRun{$experimentId}{$libraryId}{$runId} = $sbatchFilePath;
            $sbatchToRun{$experimentId}{$libraryId}{'speciesId'} = $speciesId;
            open(FH, '>', $sbatchFilePath) or die $!;
            print FH $sbatchTemplate;
            close(FH);
        } else {
            warn "\t[$runId] is not an SRA id\n";
            system("rm -f $libDirectory/sra/$runId.sra*");
        }
    }
}
print "created $jobs_created sbatch files.\n";

my $numberJobRun = 0;

if(! $doNotDownload) {
    my $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
    foreach my $experimentId (keys %sbatchToRun){
        foreach my $libraryId (keys %{$sbatchToRun{$experimentId}}){
            my $speciesId = $sbatchToRun{$experimentId}{$libraryId}{'speciesId'};
            my $libDirectory = "$outputDir/$speciesId/$libraryId";
            foreach my $runId (keys %{$sbatchToRun{$experimentId}{$libraryId}}){
                next if (-f "$libDirectory/$runId.done");
                next if $runId eq "speciesId";
                $numberJobRun++;
                $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
                while ($jobsRunning >= $parallelJobs) {
                    sleep(15);
                    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
                }
                make_path("$libDirectory");
                chdir "$libDirectory";
                system("sbatch $sbatchToRun{$experimentId}{$libraryId}{$runId}>/dev/null");
            }
        }
    }

    while ($jobsRunning > 0) {
        sleep(15);
        $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $jobPrefix);
    }
}
print "all download jobs finished. Run $numberJobRun jobs\n";
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
            remove_tree("$libDirectory");
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
