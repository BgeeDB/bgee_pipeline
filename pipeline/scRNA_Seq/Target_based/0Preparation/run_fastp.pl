#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Slurp;

# Julien Wollbrett, created March 2023
# Script used to run fastp and generate R.stat file for one run

#####################################################################

# Define arguments & their default value
my ($run_path, $run_id)  = ('', '');
my %opts = ('run_path=s'        => \$run_path,      # path to the directory containing run information
	    'run_id=s'		=> \$run_id         #ID of the run
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ($run_path eq '' || $run_id eq ''){
    print "\n\tInvalid or missing argument:
\te.g., $0  -run_path=\$(RUN_PATH) -run_id=\$(RUN_ID)> $@.tmp 2>warnings.$@
\t-run_path              path to the run directory \n
\t-run_id                ID of the run \n";
    exit 1;
}
my $fastq_path = "${run_path}/FASTQ/";
my $fastq_fastp = '';
my $fastq_R     = '';
if ( -e "${fastq_path}${run_id}_R1.fastq.gz" && -e "${fastq_path}${run_id}_R2.fastq.gz" ){
    $fastq_fastp = "${fastq_path}${run_id}_R1.fastq.gz -I ${fastq_path}${run_id}_R2.fastq.gz";
    $fastq_R     = "${fastq_path}${run_id}_R1.fastq.gz    ${fastq_path}${run_id}_R2.fastq.gz";
} else {
    die "FASTQ files ${fastq_path}${run_id}_R1.fastq.gz or ${fastq_path}${run_id}_R2.fastq.gz were not properly downloaded for run : ", $run_id, "\n";
}
# Run FastP (A quality control tool for high throughput sequence data) for ALL SRR (runs)
# as well as basic read length statistics with R
#NOTE Would be nice to have all basic stats from FastP (currently some are done in R)
if ( !-e "${run_path}/${run_id}.fastp.html.xz" || !-e "${run_path}/${run_id}.fastp.json.xz" ){
    system("fastp -i $fastq_fastp --json ${run_path}/${run_id}.fastp.json --html ${run_path}/${run_id}.fastp.html --thread 2  > ${run_path}/${run_id}.fastp.log 2>&1")==0
        or do { warn "\tfastp failed for [${run_path}${run_id}]\n" };
    system("xz -9 ${run_path}/${run_id}.fastp.html ${run_path}/${run_id}.fastp.json");
}
my $numberRstatLines = 0;
my $rStatFile = "${run_path}/${run_id}.R.stat";
if (-e "$rStatFile") {
    $numberRstatLines = `wc -l < $rStatFile`;
}
# test number of lines as an out of memory issue while generating stats will leave the run with a R.stat file containing one single line
if ($numberRstatLines != 2) {
    system("/bin/echo \"#min\tmax\tmedian\tmean\" > $rStatFile");
    system("zcat $fastq_R | sed -n '2~4p' | awk '{print length(\$0)}' | Rscript -e 'd<-scan(\"stdin\", quiet=TRUE);cat(min(d), max(d), median(d), mean(d), sep=\"\\t\");cat(\"\\n\")' >> $rStatFile");
}

