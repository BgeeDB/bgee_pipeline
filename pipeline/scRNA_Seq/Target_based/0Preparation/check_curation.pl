#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
# This script check various common mistakes in annotation and absence of annotated ontology terms in the database

#############################################################

use Getopt::Long;
use File::Spec;
use FindBin;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;

$| = 1;
require("$FindBin::Bin/../../rna_seq_utils.pl");

# Define arguments & their default value
my ($bgee_connector)                = ('');
my ($experiments, $libraries, $barcodeDir)  = ('', '', '');
my ($extraMapping)                  = ('');
my %opts = ('bgee=s'                => \$bgee_connector,     # Bgee connector string
            'experiments=s'         => \$experiments,
            'libraries=s'           => \$libraries,
            'barcodeDir=s'          => \$barcodeDir,
            'extraMapping=s'        => \$extraMapping,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $experiments eq '' || $libraries eq '' || $barcodeDir eq '' || $extraMapping eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD) -experiments=\$(RNASEQ_EXPERIMENT_FILEPATH) -libraries=\$(RNASEQ_LIB_FILEPATH) -extraMapping=\$(EXTRAMAPPING_FILE) before|after
\t-bgee                 Bgee connector string
\t-experiments          Droplet based experiments from annotation file
\t-libraries            Droplet based libraries from annotation file
\t-barcodeDir           Path to the directory containing the mapping barcode to celltype
\t-extraMapping         file containing remapping for ontology terms not yet inserted in the database
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

#########################
# Check experiments
#########################
print "Checking experiments...\n";

my %tsv = %{ Utils::read_spreadsheet("$experiments", "\t", 'csv', '', 1) };
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

# experiments have already been parsed, so we can re-use it:
my %experiments = getAllAnnotatedExperiments2( \%tsv );
my %experimentSources = ();
for my $expId ( keys %experiments ) {
    # store the source to check that it exists in the database
    $experimentSources{$experiments{$expId}->{'source'}}++;
}

my $selSrc = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName = ? AND category = "Single-cell RNA-Seq data source"');
for my $source ( keys %experimentSources ) {
    $selSrc->execute($source)  or die $selSrc->errstr;
    unless ( my @data = $selSrc->fetchrow_array ) {
        warn "Problem! [$source] in experiment file does not exist for single cell\n";
    }
}
$selSrc->finish;
print "Done checking experiments\n";

#######################################
# retrieve all organs/stages from Bgee
#######################################
my %organs = %{ Utils::getBgeedbOrgans($dbh) };
my %stages = %{ Utils::getBgeedbStages($dbh) };
$dbh->disconnect;

# Parse extra mapping info for currently too up-to-date annotations
##UnmappedId    UnmappedName    UberonID    UberonName    Comment
my %extra = map  { my @tmp = split(/\t/, $_, -1); if ( $tmp[2] ne '' && $tmp[0] ne '' ){ $tmp[0] => $tmp[2] } else { 'nonono' => 'nonono' } }
            grep { !/^#/ }
            read_file("$extraMapping", chomp => 1);


######################
# Check libraries
######################
print "Checking libraries annotation...\n";

my %tsv2 = %{ Utils::read_spreadsheet("$libraries", "\t", 'csv', '', 1) };
# Store all library IDs
my %libraryIds;
# Store all experiment IDs
my %experimentIds;
my $previousExpId;
my %missingOrganIdFromLib;
my %missingStageIdFromLib;
my %missingCellTypeIdFromLib;
for my $line ( 0..$#{$tsv2{'libraryId'}} ) {
    # line commented, skipped
    next  if ( $tsv2{'libraryId'}[$line] =~ /^#/ );
    # Check wrong chars depending on the column
    if ( $tsv2{'libraryId'}[$line]    !~ /^[\w\-\.]+$/
      || $tsv2{'experimentId'}[$line] !~ /^[\w\-\.]+$/
      || $tsv2{'platform'}[$line]     !~ /^[ \w\-\.]+$/
      || $tsv2{'anatId'}[$line]     !~ /^[\w]+:[\w]+$/
      || $tsv2{'cellTypeId'}[$line]     !~ /^[\w]+:[\w]+$/
      || $tsv2{'stageId'}[$line]      !~ /^[\w]+:[\w]+$/ ) {
        warn "Problem! wrong character inserted in line $line \[$tsv2{'libraryId'}[$line]] [$tsv2{'experimentId'}[$line]] [$tsv2{'platform'}[$line]] [$tsv2{'anatId'}[$line]] [$tsv2{'cellTypeId'}[$line]] [$tsv2{'stageId'}[$line]]\n";
    }

    my $libraryId    = $tsv2{'libraryId'}[$line];
    my $expId        = $tsv2{'experimentId'}[$line];
    my $organId      = $tsv2{'anatId'}[$line];
    my $cellTypeId   = $tsv2{'cellTypeId'}[$line];
    my $stageId      = $tsv2{'stageId'}[$line];

    warn "Warning, [$libraryId] used several times\n"    if ( exists $libraryIds{$expId}{$libraryId} );
    warn "Warning, [$libraryId] char length too long\n"  if ( length($libraryId) > 255 );

    if ( !defined $previousExpId || $previousExpId ne $expId ) {
        if ( exists $experimentIds{$expId} ) {
            warn "Warning, not continuous experiment [$expId] (split in different part of the file)\n";
        }
        if ( !exists $experiments{$expId} ) {
            warn "Warning, [$expId] not existing in the experiment annotation file\n";
        }
    }
    $experimentIds{$expId} = 1;
    $missingOrganIdFromLib{$organId} = 1 if (!exists $organs{$organId} && ! exists $extra{$organId});
    # Do not control celltype yet as the celltype at library level is not inserted in the database
    #$missingCellTypeIdFromLib{$cellTypeId} = 1 if (!exists $organs{$cellTypeId} && ! exists $extra{$cellTypeId});
    $missingStageIdFromLib{$stageId} = 1 if (!exists $stages{$stageId} && ! exists $extra{$stageId});

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
    $libraryIds{$expId}{$libraryId} = 1;
}
for my $organId (sort keys %missingOrganIdFromLib) {
    warn "Warning, [$organId] organId does not exist in the database or in the extra mapping file\n";
}
## commented for now as we do not use the celltype annotation coming from the scRNASeqLibrary_merged.tsv file
#for my $cellTypeId (sort keys %missingCellTypeIdFromLib) {
#    warn "Warning, [$cellTypeId] cellTypeId does not exist in the database or in the extra mapping file\n";
#}
for my $stageId (sort keys %missingStageIdFromLib) {
    warn "Warning, [$stageId] stageId does not exist in the database or in the extra mapping file\n";
}
print "Done checking libraries\n";

print "Checking celltypes in barcode annotation files\n";

for my $expId (keys %libraryIds) {
    my %missingCellTypesFromBarcode = ();
    my $experimentBarcodeFile = "$barcodeDir/scRNASeq_barcode_$expId.tsv";
    print "check barcode file $experimentBarcodeFile for experiment $expId\n";
    my %librariesOfExperiment = map { $_ => 1 } keys %{$libraryIds{$expId}};
    my %barcodeTsv = %{ Utils::read_spreadsheet("$experimentBarcodeFile", "\t", 'csv', '', 1) };
    for my $line ( 0..$#{$barcodeTsv{'library'}} ) {
        my $libraryId = $barcodeTsv{'library'}[$line];
        next if (! exists $librariesOfExperiment{$libraryId});
        my $cellTypeId = $barcodeTsv{'cellTypeId'}[$line];
        $missingCellTypesFromBarcode{$cellTypeId} = 1 if (!exists $organs{$cellTypeId} && ! exists $extra{$cellTypeId});
    }
    for my $cellTypeId (sort keys %missingCellTypesFromBarcode) {
        warn "Warning, [$cellTypeId] cellTypeId  used in barcode annotation for experiment $expId does not exist in the database or in the extra mapping file\n";
    }
}
print "Done checking celltypes in barcode annotation files\n";
print "Checking experiments with no libraries annotated...\n";
my %annotations = getAllRnaSeqAnnotations2( \%tsv2 );
for my $expId ( keys %experiments ) {
    if ($experiments{$expId}->{'commented'} eq 1){
        warn "Warning, experiment [$expId] with no annotated libraries\n"  if ( !exists $annotations{$expId} );
    }
}
print "Done checking experiments with no libraries annotated.\n";

exit 0;

