#!/usr/bin/env perl

## This script allows to limit the number of jobs to run on a slurm cluster
## based on a file listing all sbatch commands

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use File::Path qw(make_path);
use File::Basename;
use File::Slurp;

## Define arguments & their default value
my ($input_file, $parallel_jobs, $job_prefix, $account) = ('', '', '', '');
my %opts = ('input_file=s'        => \$input_file,
            'parallel_jobs=s'     => \$parallel_jobs,
            'job_prefix=s'        => \$job_prefix,
            'account=s'           => \$account,
           );

my $test_options = Getopt::Long::GetOptions(%opts);

print "$input_file";
my $numberJobRun = 0;
my $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $job_prefix);
for my $line ( read_file($input_file, chomp=>1) ){
    next  if ( $line =~ /^#!/);
    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $job_prefix);
	while ($jobsRunning >= $parallel_jobs) {
        sleep(15);
        $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $job_prefix);
    }
    system("$line");
    $numberJobRun++;
}

print "all jobs created properly. Run $numberJobRun jobs\n";

while ($jobsRunning > 0) {
    sleep(15);
    $jobsRunning = Utils::check_active_jobs_number_per_account_and_name($account, $job_prefix);
}

print "all jobs finished\n";

