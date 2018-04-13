#!/usr/bin/env perl

# Frederic Bastian, created December 2012
# USAGE: perl copy_files_to_ftp.pl
# copy affy data that have passed our quality controls and filtering for duplicates step, to the FTP server
# using ssh/scp
# this script assumes that the host where to put files is prd.vital-it.ch,
# that the username is bbgee, and that the files must be put in /db/ftp/pub/databases/Bgee/affymetrix_data/
#
#FIXME to redo completely: too many hard coded paths and mysql connection info
#FIXME FTP is now on ftpbgee and files are on bigbgee !!!
#
#############################################################

# Perl core modules
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../..";
use Utils;

$| = 1;

my $dir       = '/db/ftp/pub/databases/Bgee/affymetrix_data/';
my $loginHost = 'bbgee@prd.vital-it.ch';

# Don't check success of these commands, maybe the dir does not exist
print "Create new_cel_files and new_mas5_files directories on host...\n";
system('ssh', $loginHost, "rm -rf $dir/new_cel_files/");
system('ssh', $loginHost, "rm -rf $dir/new_mas5_files/");
# here we control the success
system('ssh', $loginHost, "mkdir $dir/new_cel_files/")==0   or do {die "error: $!\n"};
system('ssh', $loginHost, "mkdir $dir/new_mas5_files/")==0  or do {die "error: $!\n"};
print "Done\n";


print "Retrieving data from Bgee database...\n";
# Get the chips that passed our QC procedures
my $bgee_connector = 'user=root__pass=XXX__host=127.0.0.1__port=3306__name=bgee_v14';
my $bgee  = Utils::connect_bgee_db($bgee_connector);
my @chips = ();

my $ret = $bgee->prepare('SELECT microarrayExperimentId, affymetrixChipId, normalizationType, anatEntityId, stageId FROM affymetrixChip ORDER BY microarrayExperimentId, normalizationType');
$ret->execute();
while ( my @data = $ret->fetchrow_array ){
	my %chip = ();
	$chip{'expId'}         = $data[0];
	$chip{'chipId'}        = $data[1];
	$chip{'normalization'} = $data[2];
	$chip{'organId'}       = $data[3];
	$chip{'stageId'}       = $data[4];
    push(@chips, \%chip);
}
$ret->finish;
$bgee->disconnect;
print "Done\n";

my $ftpCelFilesDir    = 'new_cel_files/';
my $ftpMas5FilesDir   = 'new_mas5_files/';
my $localCelFilesDir  = '/var/bgee/extra/pipeline/Affymetrix/cel_data/';
my $localMas5FilesDir = '/var/bgee/extra/pipeline/Affymetrix/processed_mas5/';

my $previousExpId         = '';
my $previousNormalization = '';
print "Start uploading...\n";
# provide the annotations
open(my $OUT, '>', 'annotations.tsv');
print {$OUT} "Experiment ID\tChip ID\tfile type\tanatomical entity ID\tstage ID\n";
for my $chipRef ( @chips ){
	# try to create the directory for the experiment
	my $ftpDirToUse = '';
	my $filePath    = '';
	if ( $chipRef->{'normalization'} ne 'MAS5' ){
	    $ftpDirToUse = $dir.$ftpCelFilesDir.$chipRef->{'expId'}.'/';
	    $filePath    = $localCelFilesDir.$chipRef->{'expId'}.'/'.$chipRef->{'chipId'};

	    if ( -e $filePath.'.cel' ){
	    	$filePath = $filePath.'.cel';
	    }
        elsif ( -e $filePath.'.cel.gz' ){
	    	$filePath = $filePath.'.cel.gz';
	    }
        elsif ( -e $filePath.'.CEL' ){
	    	$filePath = $filePath.'.CEL';
	    }
        elsif ( -e $filePath.'.CEL.gz' ){
	    	$filePath = $filePath.'.CEL.gz';
	    }
        elsif ( -e $filePath.'.gz' ){
	    	$filePath = $filePath.'.gz';
	    }
	}
    else {
		$ftpDirToUse = $dir.$ftpMas5FilesDir.$chipRef->{'expId'}.'/';
	    $filePath    = $localMas5FilesDir.$chipRef->{'expId'}.'/'.$chipRef->{'chipId'};
	}

	if ( !-e $filePath ){
		print "Problem, could not find [$filePath]\n";
		next;
	}
	if ( $previousExpId ne $chipRef->{'expId'} || $previousNormalization ne $chipRef->{'normalization'} ){
		print "\t", $chipRef->{'expId'}, "\t", $chipRef->{'normalization'}, "\n";
        system('ssh', $loginHost, "mkdir $ftpDirToUse")==0  or do {die "error: $!\n"};
	}

	# copy the chip file
	system('scp', $filePath, $loginHost.':'.$ftpDirToUse.'.')==0  or do {die "error: $!\n"};
    print {$OUT} $chipRef->{'expId'}, "\t", $chipRef->{'chipId'}, "\t",
              ($chipRef->{'normalization'} eq 'gcRMA') ? 'CEL file' : 'MAS5 file', "\t",
              $chipRef->{'organId'}, "\t", $chipRef->{'stageId'}, "\n";


    $previousExpId         = $chipRef->{'expId'};
    $previousNormalization = $chipRef->{'normalization'};
}
close $OUT;

# transfering the annotation file
system('scp', 'annotations.tsv', $loginHost.':'.$dir.'.')==0  or do {die "error: $!\n"};
print "Finished\n";

exit;

