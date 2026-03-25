#!/usr/bin/env perl

# Frederic Bastian, May 2015
# Julien Wollbrett, April 2024

# This script is reponsible for generating the download files
# containing our processed expression data (not the calls,
# but the CPM/TPM, read counts, ranks, log signal intensities, etc).
# See opt debug message below for information about parameters.
#
# Since Bgee 15.1 all RNASeq datatypes (bulk, full length and droplet based)
# are stored in the same tables in the database. The code has been refactored
# to be able to use the same function for each of these datatypes

#############################################################

use strict;
use warnings;
use diagnostics;

use Archive::Tar;
use Getopt::Long;
use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Copy qw(move);
use File::Basename qw(basename);
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Parallel::ForkManager;

$| = 1;

# Define arguments & their default value
my ($bgee_connector)                = ('');
my ($speciesArg)                    = ('');
my ($rnaSeqDir)                     = ('');
my ($flScRnaSeqDir)                 = ('');
my ($dbScRnaSeqDir)                 = ('');
my ($bgeeVersion)                   = ('');
# by default generate 2 files in parallel
my ($parallelJobs)                  = (1);
my ($debug)                         = (0);
my %opts = ('bgee=s'                => \$bgee_connector,     # Bgee connector string
            'speciesArg=s'         => \$speciesArg,
            'rnaSeqDir=s'          => \$rnaSeqDir,
            'flScRnaSeqDir=s'      => \$flScRnaSeqDir,
            'dbScRnaSeqDir=s'      => \$dbScRnaSeqDir,
            'bgeeVersion=s'        => \$bgeeVersion,
            'parallelJobs=i'       => \$parallelJobs,
            'debug'                => \$debug
           );

# Check arguments
my $emptyArg = '-';
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $speciesArg eq '' ||
    $rnaSeqDir eq '' || $flScRnaSeqDir eq '' || $dbScRnaSeqDir eq '' || $bgeeVersion eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEE_CMD) -speciesArg=\$(SPECIES_ARG) -rnaSeqDir=\$(RNA_SEQ_DIR) -flScRnaSeqDir=\$(SC_RNA_SEQ_FL_DIR)
\t-bgee                 Bgee    connector string
\t-speciesArg           A comma-separated list of species IDs to generate files for, or '-' to generate files for all species
\t-rnaSeqDir            The path to generate RNA-Seq related files, or '-' for not generating data for RNA-Seq
\t flScRnaSeqDir        The path to generate full length single cell RNA-Seq related files, or '-' for not generating data for full length single cell RNA-Seq
\t dbScRnaSeqDir        The path to generate droplet based single cell RNA-Seq related files, or '-' for not generating data for droplet based single cell RNA-Seq
\t-bgeeVersion          The Bgee release for which the files are being generated, to generate FTP links to correct version, e.g. 'bgee_v15'.
\t-parallelJobs         The number of species for which files are generated in parallel
\t-debug                more verbose output
\n";
    exit 1;
}

#############################################################
# Retrieve information needed for all species / all data types


# -----------------------
# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);
# species list
my @speciesList = ();
if ($speciesArg ne $emptyArg) {
    @speciesList = split(',', $speciesArg);
}

# -----------------------
# retrieve species information (notably to generate file names)
# $species{speciesId}{'name'} = genus . species
my %species;

my $speciesSql = 'SELECT speciesId, genus, species FROM species';
if (@speciesList) {
    $speciesSql .= ' WHERE speciesId IN (';
    for my $i (0 .. $#speciesList) {
        if ($i > 0) {
            $speciesSql .= ', ';
        }
        $speciesSql .= $speciesList[$i];
    }
    $speciesSql .= ')';
}
my $speciesStmt = $dbh->prepare($speciesSql);
$speciesStmt->execute()  or die $speciesStmt->errstr;
while ( my @data = $speciesStmt->fetchrow_array ){
    $species{$data[0]}{'name'} = $data[1].' '.$data[2];
}


# -----------------------
# Retrieve information about data sources
# $dataSources{dataSourceId}{'name'} = data source name
# $dataSources{dataSourceId}{'experimentBaseUrl'} = experiment base URL (replace '[experiment_id]' with exp ID)
# $dataSources{dataSourceId}{'evidenceBaseUrl'} = evidence base URL (replace '[evidence_id]' with evidence ID, such as chip ID or library ID)
my %dataSources;
my $expIdTag      = quotemeta '[experiment_id]';
my $evidenceIdTag = quotemeta '[evidence_id]';

my $sourceStmt = $dbh->prepare('SELECT dataSourceId, dataSourceName, experimentUrl, evidenceUrl FROM dataSource');
$sourceStmt->execute()  or die $sourceStmt->errstr;
while ( my @data = $sourceStmt->fetchrow_array ){
    $dataSources{$data[0]}{'name'}              = $data[1];
    $dataSources{$data[0]}{'experimentBaseUrl'} = $data[2];
    $dataSources{$data[0]}{'evidenceBaseUrl'}   = $data[3];
}

#disconnect the connection to avoid issues with the ForkManager
$dbh->disconnect;

#############################################################
# Now, for each species, generate the files

# Link to FTP storing our files
my $ftpFilePath = 'https://www.bgee.org/ftp/';

my $pm = new Parallel::ForkManager($parallelJobs);
for my $speId ( keys %species ){

    # Forks and returns the pid for the child
    my $pid = $pm->start and next;

    # -----------------------
    # bulk RNA-Seq data
    if ($rnaSeqDir ne $emptyArg) {
        print "Generating Bulk RNA-Seq files for $speId\n";
        my $isSingleCell = 0;
        my $isSampleMultiplexing = 0;
        my $absDir = $rnaSeqDir;
        if ( !File::Spec->file_name_is_absolute($rnaSeqDir) ){
            $absDir = File::Spec->rel2abs($rnaSeqDir) ;
        }
        generateRnaSeqFiles($speId, $species{$speId}{'name'}, \%dataSources, $absDir, $isSingleCell, $isSampleMultiplexing,
            $bgee_connector);
        print "Done generating Bulk RNA-Seq files for $speId...\n";
    }
    # -----------------------
    # Full Length single cell RNA-Seq data
    if ($flScRnaSeqDir ne $emptyArg) {
        print "Generating Full Length single cell RNA-Seq files for $speId...\n";
        my $isSingleCell = 1;
        my $isSampleMultiplexing = 0;
        my $absDir = $flScRnaSeqDir;
        if ( !File::Spec->file_name_is_absolute($flScRnaSeqDir) ){
            $absDir = File::Spec->rel2abs($flScRnaSeqDir) ;
        }
        generateRnaSeqFiles($speId, $species{$speId}{'name'}, \%dataSources, $absDir, $isSingleCell, $isSampleMultiplexing,
            $bgee_connector);
        print "Done generating Full Length single cell RNA-Seq files for $speId...\n";
    }
    # -----------------------
    # Droplet based single cell RNA-Seq data
    if ($dbScRnaSeqDir ne $emptyArg) {
        print "Generating Droplet-based single-cell RNA-Seq files for $speId...\n";
        my $isSingleCell = 1;
        my $isSampleMultiplexing = 1;
        my $absDir = $dbScRnaSeqDir;
        if ( !File::Spec->file_name_is_absolute($dbScRnaSeqDir) ){
            $absDir = File::Spec->rel2abs($dbScRnaSeqDir) ;
        }
        generateRnaSeqFiles($speId, $species{$speId}{'name'}, \%dataSources, $absDir, $isSingleCell, $isSampleMultiplexing,
            $bgee_connector);
        print "Done generating Droplet-based single-cell RNA-Seq files for $speId...\n";
    }
    print "Done generating files for species $speId.\n";

    $pm->finish;

}

$pm->wait_all_children;

exit 0;

#############################################################
# SUBS

sub generateRnaSeqFiles {
    my @args = @_;
    my $speciesId             = $args[0];
    my $speciesName           = $args[1];
    my $dataSourcesRef        = $args[2];
    my $filesDir              = $args[3];
    my $isSingleCell          = $args[4];
    my $isSampleMultiplexing  = $args[5];
    my $bgee_connector        = $args[6];


    my $bgee_thread = Utils::connect_bgee_db($bgee_connector);

    my $speciesNameForFile = $speciesName;
    $speciesNameForFile =~ s/ /_/g;
    my $tmpDirName = $speciesNameForFile.'_tmp';
    my $speciesDirName = $speciesNameForFile;
    my $expFileName = '';
    my $libFileName = '';
    my $fileNamePattern = '';
    my $expLibTarballName = '';
    my $exprFilePath = $bgeeVersion.'/download/processed_expr_values/';
    if (!$isSingleCell) {
        $expFileName = $speciesNameForFile.'_RNA-Seq_experiments.tsv';
        $libFileName = $speciesNameForFile.'_RNA-Seq_libraries.tsv';
        $fileNamePattern = '_RNA-Seq_read_counts_TPM';
        $expLibTarballName = $speciesNameForFile.'_RNA-Seq_experiments_libraries.tar.gz';
        $exprFilePath .= 'rna_seq/'.$speciesNameForFile.'/';
    } elsif ($isSingleCell) {
        if ($isSampleMultiplexing) {
            $expFileName = $speciesNameForFile.'_Droplet-Based_SC_RNA-Seq_experiments.tsv';
            $libFileName = $speciesNameForFile.'_Droplet-Based_SC_RNA-Seq_libraries.tsv';
            $fileNamePattern = '_Droplet-Based_SC_RNA-Seq_read_counts_CPM';
            $expLibTarballName = $speciesNameForFile.'_Droplet-Based_SC_RNA-Seq_experiments_libraries.tar.gz';
            $exprFilePath .= 'droplet_based/'.$speciesNameForFile.'/';
        } elsif (!$isSampleMultiplexing) {
            $expFileName = $speciesNameForFile.'_Full-Length_SC_RNA-Seq_experiments.tsv';
            $libFileName = $speciesNameForFile.'_Full-Length_SC_RNA-Seq_libraries.tsv';
            $fileNamePattern = '_Full-Length_SC_RNA-Seq_read_counts_TPM';
            $expLibTarballName = $speciesNameForFile.'_Full-Length_SC_RNA-Seq_experiments_libraries.tar.gz';
            $exprFilePath .= 'full_length/'.$speciesNameForFile.'/';

        }
    }
    if ($expFileName eq '' || $libFileName eq '') {
        die "unrecognized RNASeq datatype";
    }
    # recreate temp directory for experiment and library information
    remove_tree(File::Spec->catdir( $filesDir, $tmpDirName));
    make_path(File::Spec->catdir( $filesDir, $tmpDirName));
    
    # -------------------------------------------------
    # First, we retrieve information about experiments
    # $exp{'expId'}                = exp ID
    # $exp{'name'}                 = exp name
    # $exp{'desc'}                 = exp description
    # $exp{'sourceId'}             = data source ID
    # $exp{'libCount'}             = library count
    # $exp{'conditionCount'}       = informative condition count (expression-mapped condition IDs)
    # $exp{'anatEntityStageCount'} = count of couples organ/stage
    # $exp{'anatEntityCount'}      = anatomical entity count
    # $exp{'stageCount'}           = stage count
    # sex types that are not 'not annotated' and 'NA'
    # $exp{'sexCount'}             = informative sex information count
    # strain count that are not 'not annotated', 'NA', 'confidential_restricted_data'
    # $exp{'strainCount'}          = informative strain count
    my @experiments = ();

    my $sqlExpPart =
    'SELECT t1.rnaSeqExperimentId, t1.rnaSeqExperimentName, t1.rnaSeqExperimentDescription, '
              .'t1.dataSourceId, COUNT(DISTINCT t2.rnaSeqLibraryId) AS libCount, '
              # to count the number of conditions, we need to used the expression-mapped condition IDs, 
              # which merge non-informative annotations such as 'not annotated' and 'NA'
              .'COUNT(DISTINCT t4.exprMappedConditionId) AS conditionCount, '
              .'COUNT(DISTINCT t4.anatEntityId, stageId) AS anatEntityStageCount, '
              .'COUNT(DISTINCT t4.anatEntityId) AS anatEntityCount, '
              # distinct cellTypeIds are also exported to the bulk RNASeq files even if the
              # number should always be 1 (root of the ontology)
              .'COUNT(DISTINCT t4.cellTypeId) AS cellTypeCount, '
              .'COUNT(DISTINCT t4.stageId) AS stageCount, '
              .'COUNT(DISTINCT t4.sex)  AS sexCount, '
              .'COUNT(DISTINCT t4.strain) AS strainCount '
              .'FROM rnaSeqExperiment AS t1 '
              .'INNER JOIN rnaSeqLibrary AS t2 ON t1.rnaSeqExperimentId = t2.rnaSeqExperimentId '
              .'INNER JOIN rnaSeqLibraryAnnotatedSample AS t3 ON t2.rnaSeqLibraryId = t3.rnaSeqLibraryId '
              .'INNER JOIN cond AS t4 ON t3.conditionId = t4.conditionId '
              .'WHERE speciesId = ? '
              .'AND t2.rnaSeqTechnologyIsSingleCell = ? '
              .'AND t2.sampleMultiplexing = ? '
              .'GROUP BY t1.rnaSeqExperimentId '
              .'ORDER BY libCount DESC, conditionCount DESC, anatEntityStageCount DESC, '
              .'anatEntityCount DESC, cellTypeCount DESC, stageCount DESC, sexCount DESC, '
              .'strainCount DESC, t1.rnaSeqExperimentId';

    my $stmt = $bgee_thread->prepare($sqlExpPart);
    $stmt->execute($speciesId, $isSingleCell, $isSampleMultiplexing)  or die $stmt->errstr;

    while ( my @data = $stmt->fetchrow_array ){
        my %exp;
        $exp{'expId'}                = $data[0];
        $exp{'name'}                 = $data[1];
        $exp{'desc'}                 = $data[2];
        $exp{'sourceId'}             = $data[3];
        $exp{'libCount'}             = $data[4];
        $exp{'conditionCount'}       = $data[5];
        $exp{'anatEntityStageCount'} = $data[6];
        $exp{'anatEntityCount'}      = $data[7];
        $exp{'cellTypeCount'}        = $data[8];
        $exp{'stageCount'}           = $data[9];
        $exp{'sexCount'}             = $data[10];
        $exp{'strainCount'}          = $data[11];

        push @experiments, \%exp;
    }
    
    # recreate species directory for experiment and library information
    if( @experiments > 0){
    	remove_tree(File::Spec->catdir( $filesDir, $speciesDirName));
    	make_path(File::Spec->catdir( $filesDir, $speciesDirName));
    }

    # Print experiment information into file
    my $expFile = File::Spec->catfile($filesDir, $tmpDirName, $expFileName);
    open(my $fh, '>', $expFile) or die "Could not open file '$expFile' $!";
    print $fh "Experiment ID\tExperiment name\t"
              ."Library count\tCondition count\tOrgan-stage count\t"
              ."Organ count\tCelltype Count\tStage count\tSex count\tStrain count\t"
              ."Data source\tData source URL\tBgee normalized data URL\t"
              ."Experiment description\n";
    for my $exp ( @experiments ){
        print $fh $exp->{'expId'}."\t";
        
        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        my $name = $exp->{'name'};
        $name =~ s/"/'/g;
        print $fh '"'.$name.'"'."\t";
        
        print $fh $exp->{'libCount'}."\t"
                  .$exp->{'conditionCount'}."\t"
                  .$exp->{'anatEntityStageCount'}."\t"
                  .$exp->{'anatEntityCount'}."\t"
                  .$exp->{'cellTypeCount'}."\t"
                  .$exp->{'stageCount'}."\t"
                  .$exp->{'sexCount'}."\t"
                  .$exp->{'strainCount'}."\t";
        my $sourceUrl  = 'NA';
        my $sourceName = '';
        if (defined $dataSourcesRef->{$exp->{'sourceId'}}) {
            if ($dataSourcesRef->{$exp->{'sourceId'}}->{'experimentBaseUrl'}) {
                $sourceUrl = $dataSourcesRef->{$exp->{'sourceId'}}->{'experimentBaseUrl'};
                $sourceUrl =~ s/$expIdTag/$exp->{'expId'}/g;
            }
            $sourceName = $dataSourcesRef->{$exp->{'sourceId'}}->{'name'};
        }
        print $fh $sourceName."\t".$sourceUrl."\t";

        my $resultsFileName = "${speciesNameForFile}${fileNamePattern}_".$exp->{'expId'}.'.tsv.gz';
        $resultsFileName =~ s/ /_/g;
        print $fh $ftpFilePath.$exprFilePath.$resultsFileName."\t";

        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        my $desc = $exp->{'desc'};
        $desc =~ s/"/'/g;
        print $fh '"'.$desc.'"';

        print $fh "\n";

    }
    close $fh;

    # -------------------------------------------------
    # Now, information about libraries.
    # TODO: add field allMappedReadsCount when database is updated
    # $lib{'expId'}                            = experiment ID
    # $lib{'libId'}                            = library ID
    # $lib{'platformId'}                       = RNA-Seq  platform ID
    # $lib{'tmmFactor'}                        = TMM normalization factor
    # $lib{'tpmThreshold'}                     = Threshold of TPM value to consider a gene as expressed
    # $lib{'fpkmThreshold'}                    = Threshold of FPKM value to consider a gene as expressed
    # $lib{'allGenesPercentPresent'}           = percentage of genes called present
    # $lib{'proteinCodingGenesPercentPresent'} = percentage of protein coding genes called present
    # $lib{'intergenicRegionsPercentPresent'}  = percentage of intergenic regions called present
    # $lib{'allReadsCount'}                    = total number of reads in library, including those not mapped.
    # $lib{'mappedReadsCount'}                 = number of reads mapped by pseudo alignement.
    # $lib{'minReadLength'}                    = min. length of the reads
    # $lib{'maxReadLength'}                    = max. length of the reads
    # $lib{'libraryType'}                      = enum('single','paired')
    # $lib{'libraryOrientation'}               = enum('forward','reverse','unstranded')
    # $lib{'anatEntityId'}                     = anat. entity ID annotated for this library
    # $lib{'anatEntityName'}                   = anat. entity name annotated for this library
    # $lib{'stageId'}                          = stage ID annotated for this library
    # $lib{'stageName'}                        = stage name annotated for this library
    # $lib{'sex'}                              = annotated sex info (not mapped for expression table)
    # $lib{'strain'}                           = annotated strain info (not mapped for expression table)
    # $lib{'exprMappedAnatEntityId'}           = anat. entity ID remapped for expression table for this library
    # $lib{'exprMappedAnatEntityName'}         = anat. entity name remapped for expression table for this library
    # $lib{'exprMappedStageId'}                = stage ID remapped for expression table for this library
    # $lib{'exprMappedStageName'}              = stage name remapped for expression table for this library
    # $lib{'exprMappedSex'}                    = sex info remapped for expression table
    # $lib{'exprMappedStrain'}                 = strain info remapped for expression table
    # $lib{'rnaSeqLibraryAnnotatedSampleDistinctRankCount'}         = count of distinct ranks in the library
    # $lib{'maxRank'}                          = maximum rank in the corresponding global condition
    # $lib{'sourceId'}                         = data source ID
    # $lib{'runIds'}                           = IDs of runs used, separated by '|'
    my @libs = ();

    my $sql = 'SELECT t1.rnaSeqExperimentId, t1.rnaSeqLibraryId, t3.anatEntityId, t4.anatEntityName, ';
    if ($isSingleCell) {
        $sql .= 
               't3.cellTypeId, t5.anatEntityName, t2.cellTypeAuthorAnnotation, '
              .'t20.cellTypeId AS exprMappedCellTypeId, t40.anatEntityName AS exprMappedCellTypeName, ';
    }
    $sql .=    't3.stageId, t6.stageName, t3.sex, '
              .'t3.strain, t2.anatEntityAuthorAnnotation, t2.stageAuthorAnnotation, '
              .'t1.rnaSeqSequencerName, t1.strandSelection, t1.cellCompartment, '
              .'t1.sequencedTranscriptPart, t1.rnaSeqPopulationCaptureId, t1.genotype, t1.allReadsCount, '
              .'t1.mappedReadsCount, t1.minReadLength, t1.maxReadLength, t1.libraryType, t2.abundanceThreshold, '
              .'t2.abundanceUnit, t2.tmmFactor, t2.allGenesPercentPresent, t2.proteinCodingGenesPercentPresent, '
              .'t2.intergenicRegionsPercentPresent, '
              .'t20.anatEntityId AS exprMappedAnatEntityId, t30.anatEntityName AS exprMappedAnatEntityName, '
              .'t20.stageId AS exprMappedStageId, t50.stageName AS exprMappedStageName, '
              .'t20.sex AS exprMappedSex, t20.strain AS exprMappedStrain, '
              .'t2.rnaSeqLibraryAnnotatedSampleDistinctRankCount, '
              # TODO to change if we ever use globalMaxRank instead of maxRank?
              # But then we would have ranks not only for conditions with data,
              # so I guess it would not be present in this file. To rethink in this case.
              .'t2.rnaSeqLibraryAnnotatedSampleMaxRank, '
              .'t7.dataSourceId, '
              .'GROUP_CONCAT(DISTINCT t8.rnaSeqRunId ORDER BY t8.rnaSeqRunId SEPARATOR "|") AS runIds, '
              # this column corresponds to an internal ID. It is used here to be able to group by annotated sample and
              # then to retrieve proper information for libraries containing several annotated samples (e.g droplet based libraries). It is not exported in the file, but it is needed for the query to work properly.
              .'t2.rnaSeqLibraryAnnotatedSampleId '
              .'FROM rnaSeqLibrary AS t1 '
              .'INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId '
              .'INNER JOIN cond AS t3 ON t3.conditionId = t2.conditionId '
              .'INNER JOIN anatEntity AS t4 ON t4.anatEntityId = t3.anatEntityId '
              .'INNER JOIN cond AS t20 ON t3.exprMappedConditionId = t20.conditionId ';
    if ($isSingleCell) {
        $sql .= 
               'INNER JOIN anatEntity AS t5 ON t3.cellTypeId = t5.anatEntityId '
              .'INNER JOIN anatEntity AS t40 ON t20.cellTypeId = t40.anatEntityId ';
    }
    $sql .=    'INNER JOIN stage AS t6 ON t3.stageId = t6.stageId '
              .'INNER JOIN anatEntity AS t30 ON t20.anatEntityId = t30.anatEntityId '
              .'INNER JOIN stage AS t50 ON t20.stageId = t50.stageId '
              .'INNER JOIN ('.$sqlExpPart.') AS t7 ON t1.rnaSeqExperimentId =   t7.rnaSeqExperimentId '
              .'LEFT OUTER JOIN rnaSeqRun AS t8 ON t1.rnaSeqLibraryId = t8.rnaSeqLibraryId '
              .'WHERE t3.speciesId = ? '
              .'AND t1.rnaSeqTechnologyIsSingleCell = ? '
              .'AND t1.sampleMultiplexing = ? '
              .'GROUP BY t2.rnaSeqLibraryAnnotatedSampleId '
              .'ORDER BY libCount DESC, conditionCount DESC, anatEntityStageCount DESC, '
              .'anatEntityCount DESC, cellTypeCount DESC, '
              .'stageCount DESC, sexCount DESC, strainCount DESC, '
              .'t1.rnaSeqExperimentId, t3.anatEntityId, t3.stageId, t3.cellTypeId, '
              .'t3.sex, t3.strain, t1.rnaSeqLibraryId';
    print "query : $sql\n";
    $stmt = $bgee_thread->prepare($sql);
    $stmt->execute($speciesId, $isSingleCell, $isSampleMultiplexing, $speciesId, $isSingleCell, $isSampleMultiplexing)  or die $stmt->errstr;

    while ( my @data = $stmt->fetchrow_array ){
        my %lib;
        my $i = 0;

        $lib{'expId'}                            = $data[$i++];
        $lib{'libId'}                            = $data[$i++];
        $lib{'anatEntityId'}                     = $data[$i++];
        $lib{'anatEntityName'}                   = $data[$i++];
        if ($isSingleCell) {
            $lib{'cellTypeId'}                   = $data[$i++];
            $lib{'cellTypeName'}                 = $data[$i++];
            $lib{'cellTypeAuthorAnnotation'}     = $data[$i++];
            $lib{'exprMappedCellTypeId'}         = $data[$i++];
            $lib{'exprMappedCellTypeName'}       = $data[$i++];
        }

        $lib{'stageId'}                          = $data[$i++];
        $lib{'stageName'}                        = $data[$i++];
        $lib{'sex'}                              = $data[$i++];
        $lib{'strain'}                           = $data[$i++];
        $lib{'anatEntityAuthorAnnotation'}       = $data[$i++];
        $lib{'stageAuthorAnnotation'}            = $data[$i++];
        $lib{'sequencerName'}                    = $data[$i++];
        $lib{'strandSelection'}                  = $data[$i++];
        $lib{'cellCompartment'}                  = $data[$i++];
        $lib{'sequencedTranscriptPart'}          = $data[$i++];
        $lib{'rnaSeqPopulationCapture'}          = $data[$i++];
        $lib{'genotype'}                         = $data[$i++];
        $lib{'allReadsCount'}                    = $data[$i++];
        $lib{'mappedReadsCount'}                 = $data[$i++];
        $lib{'minReadLength'}                    = $data[$i++];
        $lib{'maxReadLength'}                    = $data[$i++];
        $lib{'libraryType'}                      = $data[$i++];
        $lib{'abundanceThreshold'}               = $data[$i++];
        $lib{'abundanceUnit'}                    = $data[$i++];
        $lib{'tmmFactor'}                        = $data[$i++];
        $lib{'allGenesPercentPresent'}           = $data[$i++];
        $lib{'proteinCodingGenesPercentPresent'} = $data[$i++];
        $lib{'intergenicRegionsPercentPresent'}  = $data[$i++];
        $lib{'exprMappedAnatEntityId'}           = $data[$i++];
        $lib{'exprMappedAnatEntityName'}         = $data[$i++];
        $lib{'exprMappedStageId'}                = $data[$i++];
        $lib{'exprMappedStageName'}              = $data[$i++];
        $lib{'exprMappedSex'}                    = $data[$i++];
        $lib{'exprMappedStrain'}                 = $data[$i++];
        $lib{'rnaSeqLibraryAnnotatedSampleDistinctRankCount'}  = $data[$i++];
        $lib{'rnaSeqLibraryAnnotatedSampleMaxRank'}            = $data[$i++];
        $lib{'sourceId'}                         = $data[$i++];
        $lib{'runIds'}                           = $data[$i++];
        push @libs, \%lib;
    }
    # Print library information into file
    my $libFile = File::Spec->catfile($filesDir, $tmpDirName, $libFileName);
    open($fh, '>', $libFile) or die "Could not open file '$libFile' $!";
    print $fh  "Experiment ID\tLibrary ID\tAnatomical entity ID\tAnatomical entity name\t"
              ."Anatomical entity author annotation\t";
    if ($isSingleCell) {
        print $fh "Celltype ID\tCelltype name\tCelltype author annotation\t";
    }
    print $fh  "Stage ID\tStage name\tStage author annotation\tSex\tStrain\t"
              ."Expression mapped anatomical entity ID\tExpression mapped anatomical entity name\t";
    if ($isSingleCell) {
        print $fh "Expression mapped celltype ID\tExpression mapped celltype name\t";
    }
    print $fh  "Expression mapped stage ID\tExpression mapped stage name\t"
              ."Expression mapped sex\tExpression mapped strain\t"
              ."Platform ID\tProtocol\tLibrary type\tLibrary orientation\t"
              ."TMM normalization factor\tExpression threshold\tExpression unit\t"
              ."Cell compartment\tSequenced transcript part\tGenotype\t"
              ."Read count\tMapped read count\t"
              ."Min. read length\tMax. read length\tAll genes percent present\t"
              ."Protein coding genes percent present\tIntergenic regions percent present\t"
              ."Distinct rank count\tMax rank in the expression mapped condition\tRun IDs\t"
              ."Data source\tData source URL\tBgee normalized data URL\tRaw file URL\n";

    for my $lib ( @libs ){
        print $fh $lib->{'expId'}."\t"
                  .$lib->{'libId'}."\t";

        print $fh $lib->{'anatEntityId'}."\t";
        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        my $toPrint = $lib->{'anatEntityName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        $toPrint = $lib->{'anatEntityAuthorAnnotation'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        if ($isSingleCell) {
            print $fh $lib->{'cellTypeId'}."\t";
            $toPrint = $lib->{'cellTypeName'};
            $toPrint =~ s/"/'/g;
            print $fh '"'.$toPrint.'"'."\t";
            $toPrint = '';
            if (defined $lib->{'cellTypeAuthorAnnotation'}) {
                $toPrint = $lib->{'cellTypeAuthorAnnotation'};
                $toPrint =~ s/"/'/g;
            }
            print $fh '"'.$toPrint.'"'."\t";

        }

        print $fh $lib->{'stageId'}."\t";
        $toPrint = $lib->{'stageName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";
        $toPrint = $lib->{'stageAuthorAnnotation'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";
        
        print $fh $lib->{'sex'}."\t";
        
        $toPrint = $lib->{'strain'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        print $fh $lib->{'exprMappedAnatEntityId'}."\t";
        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        $toPrint = $lib->{'exprMappedAnatEntityName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";
        if ($isSingleCell) {
            print $fh $lib->{'exprMappedCellTypeId'}."\t";
            # we replace double quotes with simple quotes, and we surround with double quotes
            # the values to escape potential special characters
            $toPrint = $lib->{'exprMappedCellTypeName'};
            $toPrint =~ s/"/'/g;
            print $fh '"'.$toPrint.'"'."\t";
        }

        print $fh $lib->{'exprMappedStageId'}."\t";
        $toPrint = $lib->{'exprMappedStageName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        print $fh $lib->{'exprMappedSex'}."\t";

        $toPrint = $lib->{'exprMappedStrain'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        my $genotype = '';
        if(defined $lib->{'genotype'}) {
            $genotype = $lib->{'genotype'};
        }

        print $fh $lib->{'sequencerName'}."\t".$lib->{'rnaSeqPopulationCapture'}."\t".$lib->{'libraryType'}."\t".$lib->{'strandSelection'}."\t"
            .$lib->{'tmmFactor'}."\t".$lib->{'abundanceThreshold'}."\t".$lib->{'abundanceUnit'}."\t"
            .$lib->{'cellCompartment'}."\t".$lib->{'sequencedTranscriptPart'}."\t".$genotype."\t"
            .$lib->{'allReadsCount'}."\t".$lib->{'mappedReadsCount'}."\t";

        print $fh $lib->{'minReadLength'}."\t".$lib->{'maxReadLength'}."\t"
            .$lib->{'allGenesPercentPresent'}."\t".$lib->{'proteinCodingGenesPercentPresent'}."\t"
            .$lib->{'intergenicRegionsPercentPresent'}."\t";

        my $maxRank = "NA";
        if(defined $lib->{'rnaSeqLibraryAnnotatedSampleMaxRank'}) {
            $maxRank = $lib->{'rnaSeqLibraryAnnotatedSampleMaxRank'};
        }
        my $distinctRank = "NA";
        if(defined $lib->{'rnaSeqLibraryAnnotatedSampleDistinctRankCount'}) {
            $distinctRank = $lib->{'rnaSeqLibraryAnnotatedSampleDistinctRankCount'};
        }
        print $fh $distinctRank."\t".$maxRank."\t";

        if (defined $lib->{'runIds'} && $lib->{'runIds'}) {
            print $fh $lib->{'runIds'};
        } else {
            print $fh 'NA';
        }
        print $fh "\t";

        my $sourceUrl  = 'NA';
        my $sourceName = '';
        if (defined $dataSourcesRef->{$lib->{'sourceId'}}) {
            if ($dataSourcesRef->{$lib->{'sourceId'}}->{'evidenceBaseUrl'}) {
                $sourceUrl = $dataSourcesRef->{$lib->{'sourceId'}}->{'evidenceBaseUrl'};
                $sourceUrl =~ s/$expIdTag/$lib->{'expId'}/g;
                $sourceUrl =~ s/$evidenceIdTag/$lib->{'libId'}/g;
            }
            $sourceName = $dataSourcesRef->{$lib->{'sourceId'}}->{'name'};
        }
        print $fh $sourceName."\t".$sourceUrl."\t";

        my $resultsFileName = "${speciesNameForFile}${fileNamePattern}_".$lib->{'expId'}.'.tsv.gz';
        $resultsFileName =~ s/ /_/g;
        print $fh $ftpFilePath.$exprFilePath.$resultsFileName."\t";

        print $fh 'https://trace.ncbi.nlm.nih.gov/Traces/study/?acc='.$lib->{'libId'}."\n";
    }
    close $fh;

    # -------------------------------------------------
    # Now, we generate the RNA-Seq results files,
    # one file, per experiment.
    # But with GTEx, we can't load all rnaSeqResults at once in memory, 
    # so we query one library at a time.
    # XXX: Do we normalize libraries of different types independently?
    # First, we get the list of exp (and not the library types for now)
    # $fileParams{expId} = 1;
    my %fileParams;
    for my $lib ( @libs ){
    	if (!exists $fileParams{$lib->{'expId'}}) {
    		$fileParams{$lib->{'expId'}} = ();
    	}
    	push @{$fileParams{$lib->{'expId'}}}, $lib->{'libId'};
    }

    # Now, generate the files. We'll store the path to all gz files generated,
    # to make one giant tar.gz at the end.
    my @tarFileNames = ();
    for my $expId ( keys %fileParams ){
        # we will store all names of files generated for an experiment, to pack all of them together
        my @resultsFileNames = ();
    
        # UPDATES Bgee 15.2
        # file header to update
        # - read count -> can be "read count" or "UMI count"
        # - TPM -> can be "TPM" or "CPM"
        # - FPKM -> removed
        # - cellTypeId -> new column
        # - cellTypeName -> new column

        # Because we reopen the connection in the expId loop we cannot prepare the query outside the loop.
        # XXX: left outer join to expression to retrieve the global call quality?
        $sql =           'SELECT t4.rnaSeqExperimentId, t3.rnaSeqLibraryId, t4.libraryType, t2.geneId, '
                        .'t5.anatEntityId, t6.anatEntityName, ';
        if ($isSingleCell) {
            $sql .=      't5.cellTypeId, t7.anatEntityName, t3.cellTypeAuthorAnnotation, ';
        }
        $sql .=          't5.stageId, t8.stageName, t5.sex, t5.strain, ';
        if ($isSingleCell && $isSampleMultiplexing) {
            $sql .=      't1.UMIsCount, ';
        } else {
            $sql .=      't1.ReadsCount, ';
        }
        $sql .=          't1.abundance, t1.rawRank, t1.pValue, t4.usedInPropagatedCalls '
                        .'FROM rnaSeqLibraryAnnotatedSampleGeneResult AS t1 '
                        .'INNER JOIN gene AS t2 ON t2.bgeeGeneId = t1.bgeeGeneId '
                        .'INNER JOIN rnaSeqLibraryAnnotatedSample AS t3 ON t3.rnaSeqLibraryAnnotatedSampleId = t1.rnaSeqLibraryAnnotatedSampleId '
                        .'INNER JOIN rnaSeqLibrary AS t4 ON t4.rnaSeqLibraryId = t3.rnaSeqLibraryId '
                        .'INNER JOIN cond AS t5 ON t5.conditionId = t3.conditionId '
                        .'INNER JOIN anatEntity AS t6 ON t6.anatEntityId = t5.anatEntityId ';
        if ($isSingleCell) {
            $sql .=      'INNER JOIN anatEntity AS t7 ON t7.anatEntityId = t5.cellTypeId ';
        }
        $sql .=          'INNER JOIN stage AS t8 ON t8.stageId = t5.stageId '
                        # removed that useless join as part of Bgee 15.2
                        #.'LEFT OUTER JOIN expression AS t7 ON t1.expressionId = t7.expressionId '
                        .'WHERE t4.rnaSeqLibraryId = ? '
                        .'AND t1.expressionId IS NOT NULL '
                        .'AND t5.speciesId = ? '
                        .'AND t4.rnaSeqTechnologyIsSingleCell = ? '
                        .'AND t4.sampleMultiplexing = ? '
                        .'ORDER BY t2.geneId';
        $stmt = $bgee_thread->prepare($sql);
                
        # Note: actually, it is not enough to simply sort the libraries based on their ID, 
        # we need to apply the same sorting as when generating the library info file. 
        # A simple solution is to redo a SQL query with the same ORDER BY clause. 
        my $selectOrderedLibraryIds =   'SELECT distinct t1.rnaSeqLibraryId FROM rnaSeqLibrary AS t1 '
                                       .'INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId '
                                       .'INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId '
                                       .'WHERE t1.rnaSeqExperimentId = ? '
                                       .'AND t3.speciesId = ? '
                                       .'ORDER BY t3.anatEntityId, ';
        if ($isSingleCell) {
            $selectOrderedLibraryIds .= 't3.cellTypeId, ';
        }
        $selectOrderedLibraryIds .=  't3.stageId, t3.sex, t3.strain, t1.rnaSeqLibraryId';
        
        my $getExpLibs = $bgee_thread->prepare($selectOrderedLibraryIds);
        my $resultsFileName = "${speciesNameForFile}${fileNamePattern}_".$expId.'.tsv';
        $resultsFileName =~ s/ /_/g;
        push @resultsFileNames, $resultsFileName;
        my $resultsFile = File::Spec->catfile($filesDir, $tmpDirName,
            $resultsFileName);
        open(my $fh, '>', $resultsFile) or die "Could not open file '$resultsFile' $!";
        print $fh "Experiment ID\tLibrary ID\tLibrary type\tGene ID\t"
                  ."Anatomical entity ID\tAnatomical entity name\t";
        if ($isSingleCell) {
            print $fh "Celltype ID\tCelltype name\tCelltype author annotation\t";
        }
        print $fh "Stage ID\tStage name\tSex\tStrain\t";
        if ($isSampleMultiplexing) {
            print $fh "UMI count\tCPM\t";
        } else {
            print $fh "Read count\tTPM\t";
        }
        print $fh "Rank\tDetection flag\tpValue\tState in Bgee\n";
        
        $getExpLibs->execute($expId, $speciesId) or die $getExpLibs->errstr;

        while ( my @libs = $getExpLibs->fetchrow_array ){
            my $libId = $libs[0];
            # Retrieve data from database.
            $stmt->execute($libId, $speciesId, $isSingleCell, $isSampleMultiplexing) or die $stmt->errstr;
            while ( my @data = $stmt->fetchrow_array) {
                my $i = 0;

                # we write the data directly to not store them in memory
                print $fh $data[$i++]."\t".$data[$i++]."\t".$data[$i++]."\t".$data[$i++]."\t"
                    .$data[$i++]."\t";
                # we replace double quotes with simple quotes, and we surround with double quotes
                # the values to escape potential special characters
                my $toPrint = $data[$i++];
                $toPrint =~ s/"/'/g;
                print $fh '"'.$toPrint.'"'."\t";
                #print celltypeId, celltypeName and cellTypeAuthorAnnotation
                if ($isSingleCell) {
                    print $fh $data[$i++]."\t";
                    $toPrint = $data[$i++];
                    $toPrint =~ s/"/'/g;
                    print $fh '"'.$toPrint.'"'."\t";
                    # print celltypeAuthorAnnotation that can be null
                    $toPrint = $data[$i++];
                    if (defined $toPrint) {
                        $toPrint =~ s/"/'/g;
                        print $fh '"'.$toPrint.'"'."\t";
                    } else {
                        print $fh '""'."\t";
                    }
                    

                }
                print $fh $data[$i++]."\t";
                $toPrint = $data[$i++];
                $toPrint =~ s/"/'/g;
                print $fh '"'.$toPrint.'"'."\t";

                print $fh $data[$i++]."\t";
                # print strain
                $toPrint = $data[$i++];
                $toPrint =~ s/"/'/g;
                print $fh '"'.$toPrint.'"'."\t";

                # Read count and abundance (TPM or CPM), rank, detection

                print $fh $data[$i++]."\t".$data[$i++]."\t";
                # rank, detection flag and pValue
                my $rank = "NA";
                if (defined $data[$i]) {
                    $rank = $data[$i];
                }
                $i++;
                print $fh $rank."\t";
                my $pValue = "NA";
                if(defined $data[$i]) {
                    $pValue = $data[$i];
                }
                $i++;
                if (defined $pValue && $pValue <= 0.05) {
                    print $fh "present\t";
                } else {
                    print $fh "absent\t";
                }
                print $fh $pValue."\t";
                if ($data[$i] == 1) {
                	print $fh 'Part of a call';
                } elsif ($data[$i] == 0) {
                    print $fh 'Not part of a call';
                } else {
                    die "usedInPropagatedCalls should either be equal to 1 (used to generate propagated calls)".
                        " or 0 (not used to generate propagated calls) but was ".$data[$i];
                }
                $i++;
                print $fh "\n";
            }
        }
        close $fh;
        
        # We close the connection before compressing the files, it can take quite some time
        $bgee_thread->disconnect;
        # we compress all files for an experiment together, and delete the uncompressed files
        # to save disk space
        my $fileName = "${speciesNameForFile}${fileNamePattern}_${expId}.tsv.gz";
        $fileName =~ s/ /_/g;
        unlink File::Spec->catfile($filesDir, $speciesDirName, $fileName);
        # use shell command to compress files because perl compression modules load data in memory
        # to compress. This is not possible because GTeX experiment is too big.
        # https://metacpan.org/pod/distribution/Archive-Tar/lib/Archive/Tar.pm#FAQ
        # https://www.perlmonks.org/?node_id=201013
        system("gzip -cvf ".File::Spec->catfile($filesDir, $tmpDirName, $resultsFileName)
            ." > ".File::Spec->catfile($filesDir, $speciesDirName, $fileName));
        #it is now safe to remove uncompressed file to save disk space
        unlink File::Spec->catfile($filesDir, $tmpDirName, $resultsFileName); 
        # Store the name of this experiment gz file, to make one giant tar.gz of all experiments
        # at the end.
        push @tarFileNames, $fileName;

        # reopen the connection for further use 
        $bgee_thread = Utils::connect_bgee_db($bgee_connector);
    }


    # -------------------------------------------------
    #everything went fine, we move and tar.gz the tmp files
    
    # We close the connection before compressing the files, it can take quite some time
    $bgee_thread->disconnect;

	if (scalar(keys %fileParams) > 0){
    	unlink File::Spec->catfile($filesDir, $speciesDirName, $expLibTarballName);
	    my $tar = Archive::Tar->new();
	    $tar->add_files(File::Spec->catfile($filesDir, $tmpDirName, $expFileName),
	    	File::Spec->catfile($filesDir, $tmpDirName, $libFileName));
	    for my $file_objs ($tar->get_files){
            $tar->rename( $file_objs, basename($file_objs->full_path));
	    }
	    $tar->write(
	        File::Spec->catfile($filesDir, $speciesDirName, $expLibTarballName), 
		        COMPRESS_GZIP);

 	   	# # Compress each file with RNA-Seq results independently, and add them to a global tar.gz file
 	   	# unlink File::Spec->catfile($filesDir, $speciesDirName, "${speciesNameForFile}${fileNamePattern}.tar.gz");
   		# $tar = Archive::Tar->new();
    
    	# my @tarFilePaths = ();
	    # for my $tarFileName ( @tarFileNames ) {
	    #    	push @tarFilePaths, File::Spec->catfile($filesDir, $speciesDirName, $tarFileName);
	    # }
	    # $tar->add_files(@tarFilePaths);
	    # for my $file_objs ($tar->get_files){
        #     $tar->rename( $file_objs, basename($file_objs->full_path));
	    # }
	    # $tar->write(
	    #     File::Spec->catfile($filesDir, $speciesDirName, "${speciesNameForFile}${fileNamePattern}.tar.gz"),
	    #     COMPRESS_GZIP);
	}

    remove_tree(File::Spec->catdir( $filesDir, $tmpDirName));
}