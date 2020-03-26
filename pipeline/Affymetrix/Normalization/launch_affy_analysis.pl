#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Roux, created 15/05/08
# USAGE: perl launch_affy_analysis.pl
#        perl launch_affy_analysis.pl <EXP_ID>

# Run this script on cluster:
# To do:
#
# rsync -Wav -essh ~/work/bgee/extra/pipeline/pipeline/*.pl jroux@Rserv.vital-it.ch:/scratch/temporary/jroux/pipeline/pipeline/
# rsync -Wav -essh ~/work/bgee/extra/pipeline/Affymetrix/cel_data/ jroux@Rserv.vital-it.ch:/scratch/temporary/jroux/pipeline/Affymetrix/cel_data/
# rsync -Wav -essh ~/work/bgee/extra/pipeline/Affymetrix/processed_schuster/ jroux@Rserv.vital-it.ch:/scratch/temporary/jroux/pipeline/Affymetrix/processed_schuster/
# rsync -Wav -essh ~/work/bgee/extra/pipeline/Affymetrix/bioconductor/ jroux@Rserv.vital-it.ch:/scratch/temporary/jroux/pipeline/Affymetrix/bioconductor/
# rsync -Wav -essh ~/work/bgee/extra/pipeline/Affymetrix/affymetrixChip jroux@Rserv.vital-it.ch:/scratch/temporary/jroux/pipeline/Affymetrix/
#
#
# ssh jroux@prd.vital-it.ch
# cd /scratch/temporary/jroux/pipeline/pipeline/
# perl launch_affy_analysis.pl exp_id
# perl launch_affy_analysis.pl (for all experiments not already normalized)
#
# rsync -Wav -essh /scratch/temporary/jroux/pipeline/Affymetrix/bioconductor/out/ admin@130.223.48.225:~/work/bgee/extra/pipeline/Affymetrix/bioconductor/out/
# rsync -Wav -essh /scratch/temporary/jroux/pipeline/Affymetrix/bioconductor/affinities/ admin@130.223.48.225:~/work/bgee/extra/pipeline/Affymetrix/bioconductor/affinities/
# rsync -Wav -essh /scratch/temporary/jroux/pipeline/Affymetrix/processed_schuster/ admin@130.223.48.225:~/work/bgee/extra/pipeline/Affymetrix/processed_schuster/
####################################

use Getopt::Long;
use FindBin;

require 'affy_utils.pl';
require 'bgee_utils.pl';

# Define arguments & their default value
my ($affymetrixChip) = ('');
my ($affyChipInformation, $chipTypeQual) = ('', '');
my ($cel_data, $bioconductorout, $bioconductoraffin, $processed_schuster) = ('', '', '', '');
my %opts = (
            'affyChipInformation=s'  => \$affyChipInformation,
            'chipTypeQual=s'         => \$chipTypeQual,
            'affymetrixChip=s'       => \$affymetrixChip,
            'cel_data=s'             => \$cel_data,
            'bioconductorout=s'      => \$bioconductorout,
            'bioconductoraffin=s'    => \$bioconductoraffin,
            'processed_schuster=s'   => \$processed_schuster,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $affymetrixChip eq '' || $affyChipInformation eq '' || $chipTypeQual eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -affyChipInformation=\$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=\$(AFFY_CHIPTYPEQUAL_FILEPATH) -affymetrixChip=\$(AFFY_CHIP_FILEPATH) -cel_data=\$(CELPATH) -processed_schuster=\$(SCHUSTERPATH) -bioconductorout=\$(BIOCONDUCTOROUT) -bioconductoraffin=\$(BIOCONDUCTORAFFIN)
\t-affyChipInformation        affymetrixChipInformation                     pipeline   file
\t-chipTypeQual               chipTypeCorrespondencesAndQualityThresholds   pipeline   file
\t-affymetrixChip             affymetrixChip                                annotation file
\t-cel_data                   cel_data                directory
\t-bioconductorout            bioconductor/out        directory
\t-bioconductoraffin          bioconductor/affinities directory
\t-processed_schuster         processed_schuster      directory
\n";
    exit 1;
}


$| = 1; # stdout not in memory buffer


################################################
# Launch R in batch for normalizing affy arrays
################################################

# Read affymetrixChip file
my %chip;
open(my $IN1, '<', "$affymetrixChip")  or die "Can't read file [$affymetrixChip]\n";
my $line = <$IN1>; #header
if ( defined $ARGV[0] ){
    while ( defined ($line = <$IN1>) ){
        chomp $line;
        # skip commentary line
        next  if (( $line =~ /^#/ ) or ($line =~ /^\"#/));

        my @tmp = map { bgeeTrim($_) }
                  split("\t", $line);

        # Add only the experiment in argument
        # Check that $tmp[0] is a cel file
        if ( $tmp[1] eq $ARGV[0] && $tmp[3] eq 2 ){
            # experiment -> array -> cel file
            $chip{$tmp[1]}->{$tmp[2]}->{$tmp[0]}++;
        }
    }
}
else {
    while ( defined ($line = <$IN1>) ){
        chomp $line;
        # skip commentary line
        next  if (( $line =~ /^#/ ) or ($line =~ /^\"#/));

        my @tmp = map { bgeeTrim($_) }
                  split("\t", $line);

        # experiment -> array -> cel file
        # Check that $tmp[0] is a cel file
        if ( $tmp[3] eq 2 ){
            $chip{$tmp[1]}->{$tmp[2]}->{$tmp[0]}++;
        }
    }
}
close $IN1;

# Get additional information (quality score, percent present, ...)
my %affyChipsInfo = getAllChipsInfo($affyChipInformation);
# Get the quality thresholds
my %chipTypeInfo  = getChipTypesInformation($chipTypeQual);

# #Printing
#for my $exp ( keys %chip ){
#   print "$exp\n";
#   for my $array ( keys %{$chip{$exp}} ){
#       print "\t$array\n";
#       for my $cel ( keys %{$chip{$exp}->{$array}} ){
#           print "\t\t$cel\n";
#       }
#       }
#}

my $path_cel       = $cel_data;
my $path_out       = $bioconductorout;
my $path_processed = $processed_schuster;

my %out_files;
# Check that the experiment/array is not already normalized
opendir(my $DIR1, $path_out)  or die "Can't opendir [$path_out]: $!\n";
while ( defined(my $file = readdir($DIR1)) ){
    $out_files{$file}++;
}
closedir($DIR1);

if ( (keys %chip) eq 0 ){
    print "\tProblem reading the file affymetrixChip: no .cel file found!\n"
}

for my $exp ( sort keys %chip ){
    for my $array ( keys %{$chip{$exp}} ){
        print "\n$exp\t$array...\n";
        my $file_name = $exp.'_'.$array.'.out';

        # If the exp_array was not normalized before
        if ( !exists $out_files{$file_name} ){
            # mkdir processed, if not exists
            mkdir $path_processed.$exp  if ( ! -e "$path_processed$exp" );

            # Give all cell files to R (instead of the directory)
            my $cel_files = undef;
            for my $cel ( keys %{$chip{$exp}->{$array}} ){
                # gunzip files if needed
                if ( -e "$path_cel$exp/$cel.gz" && !-e "$path_cel$exp/$cel" ){
                    if ( system("gunzip $path_cel$exp/$cel.gz") != 0 ){
                        print "\tError: could not gunzip $path_cel$exp/$cel.gz, file skipped\n";
                        next;
                    }
                }

                my $affyChipInfo = undef;
                if ( defined $affyChipsInfo{$exp}{$cel} ){
                    $affyChipInfo = \%{$affyChipsInfo{$exp}{$cel}};
                }
                if ( isChipIncompatibleOrLowQuality(2, $affyChipInfo, \%chipTypeInfo, $array) ){
                    print "\tIncompatible or low quality file skipped: expId: $exp - chipId $cel\n";
                    next;
                }
                if ( !defined $affyChipsInfo{$exp}{$cel} ){
                    print "\tFile skipped, information was not generated for expId: $exp - chipId $cel\n";
                    next;
                }
                if ( !defined $cel_files ){
                    $cel_files = 'c("';
                }
                $cel_files = $cel_files.$cel.'","';
            }

            if ( !defined $cel_files ){
                print "\tNo cel file selected, $exp - $array skipped\n";
                next;
            }
            $cel_files = substr($cel_files, 0, -2); # remove last 2 characters ,"
            $cel_files = $cel_files.')';

            # Launch R script
            ## R CMD BATCH --no-save --no-restore '--args celpath=$path_cel.$exp filenames=$cel_files array=$array processed=$path_processed.$exp' affy_analysis.R $exp_$array.out

            my @args = ("R CMD BATCH --no-save --no-restore '--args celpath=\"".$path_cel.$exp."\" filenames=".$cel_files." array=\"".$array."\" processed=\"".$path_processed.$exp."/\" affin=\"".$bioconductoraffin."\"' $FindBin::Bin/affy_analysis.R ".$path_out.$file_name);
            system(@args)==0  or do{warn "\tsystem @args failed: $? "; system("mv $path_out$file_name $path_out$file_name.PROB");};
        }
        else {
            print "\tAlready normalized!\n";
        }
    }
}
exit 0;

