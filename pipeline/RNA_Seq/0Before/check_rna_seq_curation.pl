#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Frederic Bastian, created November 2012
# Julien Roux, updated Mar 2016

# This script check various common mistakes in annotation
# (blanck spaces in IDs, wrong column order, missing files, etc)

#############################################################

use Getopt::Long;
use File::Spec;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

$| = 1;
require("$FindBin::Bin/../rna_seq_utils.pl");

# Define arguments & their default value
my ($bgee_connector)                = ('');
my ($RNAseqExperiment, $RNAseqLib)  = ('', '');
my ($allRes)                        = ('');
my ($debug)                         = (0);
my %opts = ('bgee=s'                => \$bgee_connector,     # Bgee connector string
            'RNAseqExperiment=s'    => \$RNAseqExperiment,
            'RNAseqLib=s'           => \$RNAseqLib,
            'allRes=s'              => \$allRes,
            'debug'                 => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $RNAseqExperiment eq '' || $RNAseqLib eq '' || $allRes eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD) -RNAseqExperiment=\$(RNASEQ_EXPERIMENT_FILEPATH) -RNAseqLib=\$(RNASEQ_LIB_FILEPATH) -allRes=\$(RNASEQALLRES) before|after
\t-bgee                 Bgee connector string
\t-RNAseqExperiment     RNAseq Experiment from annotation file
\t-RNAseqLib            RNAseq Libraries from annotation file
\t-allRes               gene_results directory path
\t-debug                more verbose output
\n";
    exit 1;
}
# before|after: check annotations before or after data generation;
# if after data generation, existence and coherence of data are tested.
# If undefined, "after" is assumed

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

my $afterDataGeneration = defined $ARGV[0] && $ARGV[0] eq 'before' ? 0
                        :                                            1;

#########################
# Check RNASeqExperiment
#########################
print "Checking RNASeqExperiment...\n"  if ( $debug );

my %tsv = %{ Utils::read_spreadsheet("$RNAseqExperiment", "\t", 'csv', '"', 1) };
my %experimentsCount;
# This check on experiments is partial, and is here only because it perfoms "historical" checks.
# For more complete checks, see getAllAnnotatedExperiments in rna_seq_utils.pl (below)
for my $line ( 0..$#{$tsv{'experimentId'}} ) {
    my $exp_id = $tsv{'experimentId'}[$line];
    if ( $exp_id !~ /\s/ ) {
        $experimentsCount{$exp_id}++;
        warn "Warning! [$exp_id] char length too long\n"  if ( length($exp_id) > 255 );;
    }
    else {
        warn "Problem! Space inserted for [$exp_id] Id\n";
    }
}
for my $exp ( keys %experimentsCount ) {
    die "Problem! [$exp] inserted many times\n"  if ( $experimentsCount{$exp} > 1 );
}

# More complete check for experiments in sub getAllAnnotatedExperiments.
# my %experiments = getAllAnnotatedExperiments('../../rna_seq/rnaSeqExperiment');
# RNAseqExperiment has already been parsed, so we can re-use it:
my %experiments = getAllAnnotatedExperiments2( \%tsv );
my %experimentSources = ();
for my $expId ( keys %experiments ) {
    # store the source to check that it exists in the database
    $experimentSources{$experiments{$expId}->{'source'}}++;
}

my $selSrc = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName = ?');
for my $source ( keys %experimentSources ) {
    $selSrc->execute($source)  or die $selSrc->errstr;
    unless ( my @data = $selSrc->fetchrow_array ) {
        warn "Problem! [$source] in RNASeqExperiment.tsv does not exist\n";
    }
}
$selSrc->finish;
print "Done\n"  if ( $debug );

#######################################
# retrieve all organs/stages from Bgee
#######################################
my %organs = %{ Utils::getBgeedbOrgans($dbh) };
# $organs{organId}->{'startStageId'} = startStageId
# $organs{organId}->{'endStageId'}   = endStageId

my %stages = %{ Utils::getBgeedbStages($dbh) };
# $stages{stageId}->{'leftBound'}  = leftBound
# $stages{stageId}->{'rightBound'} = rightBound
$dbh->disconnect;

######################
# Check RNASeqLibrary
######################
print "Checking RNASeqLibrary...\n"  if ( $debug );

my %tsv2 = %{ Utils::read_spreadsheet("$RNAseqLib", "\t", 'csv', '"', 1) };
# Store all library IDs
my %libraryIds;
# Store all experiment IDs
my %experimentIds;
my $previousExpId;
for my $line ( 0..$#{$tsv2{'libraryId'}} ) {
    # line commented, skipped
    next  if ( $tsv2{'libraryId'}[$line] =~ /^#/ );
    # Check wrong chars depending on the column
    if ( $tsv2{'libraryId'}[$line]    !~ /^[\w\-\.]+$/
      || $tsv2{'experimentId'}[$line] !~ /^[\w\-\.]+$/
      || $tsv2{'platform'}[$line]     !~ /^[ \w\-\.]+$/
      || $tsv2{'uberonId'}[$line]     !~ /^[\w]+:[\w]+$/
      || $tsv2{'stageId'}[$line]      !~ /^[\w]+:[\w]+$/ ) {
        warn "Problem! wrong character inserted in line $line \[$tsv2{'libraryId'}[$line]] [$tsv2{'experimentId'}[$line]] [$tsv2{'platform'}[$line]] [$tsv2{'uberonId'}[$line]] [$tsv2{'stageId'}[$line]]\n";
    }

    my $libraryId = $tsv2{'libraryId'}[$line];
    my $expId     = $tsv2{'experimentId'}[$line];
    my $organId   = $tsv2{'uberonId'}[$line];
    my $stageId   = $tsv2{'stageId'}[$line];
    warn "Warning, [$libraryId] used several times\n"    if ( exists $libraryIds{$libraryId} );
    warn "Warning, [$libraryId] char length too long\n"  if ( length($libraryId) > 255 );

    if ( !defined $previousExpId || $previousExpId ne $expId ) {
        if ( exists $experimentIds{$expId} ) {
            warn "Warning, not continuous experiment [$expId] (split in different part of the file)\n";
        }
        if ( !exists $experiments{$expId} ) {
            warn "Warning, [$expId] not existing in the experiment annotation file\n";
        }
    }
    $experimentIds{$expId}++;

    # Check existence of the directory for the processed results
    if ( $afterDataGeneration ) {
        my $libraryFile = $allRes.'/'.$libraryId.'/abundance_gene_level+new_tpm+new_fpkm+calls.tsv';
        warn "Warning, processed file for experiment [$expId], library [$libraryId] does not exist\n"  if ( !-e $libraryFile );
    }

    warn "Warning, [$organId] anatId  does not exist in the database\n"  if ( !exists $organs{$organId} );
    warn "Warning, [$stageId] stageId does not exist in the database\n"  if ( !exists $stages{$stageId} );

    # Check impossible annotation: organ not existing at the annotated stage
    if ( exists $organs{$organId} && exists $stages{$stageId} ) {
        my $startStageId = $organs{$organId}->{'startStageId'};
        my $endStageId   = $organs{$organId}->{'endStageId'};
        if ( $stages{$startStageId}->{'leftBound'} <= $stages{$stageId}->{'rightBound'} && $stages{$endStageId}->{'rightBound'} >= $stages{$stageId}->{'leftBound'} ) {
            #OK
        }
        else {
            warn "Warning, impossible annotation, [$organId] not existing at [$stageId]\n";
        }
    }

    $previousExpId = $expId;
    $libraryIds{$libraryId}++;
}
print "Done.\n"  if ( $debug );

print "Checking experiments with no libraries annotated...\n"  if ( $debug );
my %annotations = getAllRnaSeqAnnotations2( \%tsv2 );
for my $expId ( keys %experiments ) {
    if ($experiments{$expId}->{'commented'} eq 1){
        warn "Warning, experiment [$expId] with no annotated libraries\n"  if ( !exists $annotations{$expId} );
    }
}
print "Done\n"  if ( $debug );

if ( !$afterDataGeneration ) {
    print "All done.\n"  if ( $debug );
    exit 0;
}

################################################
# Check existing files with missing annotations
################################################
print "Checking existing data file with missing annotations...\n"  if ( $debug );

# Read all processed data folder and check if these data are annotated (they don't necessarily have to be)
opendir(my $IMD, $allRes)  or die("Cannot open directory [$allRes]\n");
for my $storedFile ( File::Spec->no_upwards(readdir($IMD)) ) {
    # next if hidden file/dir
    next  if ( $storedFile =~ /^\./ );

    if ( -d $allRes.$storedFile ) {
        # Check that the library was present in the annotation file
        my $libraryId = $storedFile;
        if ( !exists $annotations{$storedFile} ) {
        warn "Warning: library [$storedFile] present on the server in [$allRes], not in the rnaSeqLibrary annotation file\n";
        }
    }
}
closedir $IMD;
print "All done\n"  if ( $debug );

exit 0;

