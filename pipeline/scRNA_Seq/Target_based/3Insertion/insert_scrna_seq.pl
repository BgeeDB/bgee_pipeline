#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use File::Slurp;
use lib "$FindBin::Bin/../../../"; # Get lib path for Utils.pm
use Utils;
use Parallel::ForkManager;

$| = 1; # no buffering of output

# Julien Wollbrett, created December 2022
# This script insert target base experiments, libraries, annotatedSamples, individualSamples, 
# and individualSamples results
#  * the individualSamples result are at the cell (=barcode) level. They correspond to the raw count as they
#    are provided by kallisto bus.
#  * the annotatedSamples result are at the celltype level. They correspond to each annotation of a library.
#    To obtain this result we sum up all rawCounts of cell corresponding to the same celltype. It is this
#    information that is used to generate the propagated calls.

# Julien Wollbrett, updated March 2023
# The script now also insert annotatedSample info coming from the Bgee pipeline (e.g proportion coding present,
# cpm threshold, ...), and annotatedSampleGeneResult which Correspond to Bgee calls.

#####################################################################

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($targetBaseLibrary, $kallistoResults, $callsResults, $sexInfo)  = ('', '', '', '');
my ($singleCellExperiment, $bgeeLibraryInfo) = ('', '');

my ($pipelineCallsSummary, $pipelineReportFile, $filteredBarcodeDir) = ('', '', '');
my $extraMapping = '';
my $numberCore = 1;
#my ($library_stats, $report_info = ('', '');
my ($debug) = (0);
my %opts = ('bgee=s'                 => \$bgee_connector,       # Bgee connector string
            'targetBaseLibrary=s'    => \$targetBaseLibrary,    # target base RNAseq library annotations
            'singleCellExperiment=s' => \$singleCellExperiment, # single cell RNASeq experiment annotations
            'bgeeLibraryInfo=s'      => \$bgeeLibraryInfo,      # metadata_info_10X.txt file
            'extraMapping=s'         => \$extraMapping,         # file containing remapping of ontology terms to Uberon when initial term is not in the database yet
            'filteredBarcodeDir=s'   => \$filteredBarcodeDir,   # path to the directory containing all filtered barcode to celltype files
            'pipelineCallsSummary=s' => \$pipelineCallsSummary, # path to the file containing a summary of processing calls info at library/celltype level (e.g percentage protein coding present, ...)
            'pipelineReportFile=s'   => \$pipelineReportFile,   # path to the file containing summary of kallisto (reads, reads mapped, ..,) and stats about read length
            'kallistoResults=s'      => \$kallistoResults,      # path to dir containing kallisto/bustools results for all libraries
            'callsResults=s'         => \$callsResults,         # path to dir containing calls results for all libraries
            'sexInfo=s'              => \$sexInfo,              # generated_files/uberon/uberon_sex_info.tsv
            'numberCore=s'           => \$numberCore,           # number of cores corresponding to number of threads used to insert data in the database
            'debug'                  => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $targetBaseLibrary eq '' || $singleCellExperiment eq '' ||
    $bgeeLibraryInfo eq '' || $filteredBarcodeDir eq '' || $pipelineCallsSummary eq '' || $pipelineReportFile eq '' || $kallistoResults eq '' ||
    $sexInfo eq '' || $numberCore eq '' || $extraMapping eq '' || $callsResults eq ''){
    print "\n\tInvalid or missing argument:
\te.g., $0  -bgee=\$(BGEECMD) -targetBaseLibrary=RNASeqLibrary_full.tsv -singleCellExperiment=RNASeqExperiment_full.tsv -bgeeLibraryInfo=metadata_info_10X.txt -sexInfo=\$(UBERON_SEX_INFO_FILE_PATH) > $@.tmp 2>warnings.$@
\t-bgee                    Bgee connector string
\t-targetBaseLibrary       targetBaseLibrary annotation file
\t-singleCellExperiment    singleCellExperiment file
\t-bgeeLibraryInfo         metadata_info_10X.txt containing libraries to process
\t-extraMapping            file containing remapping for ontology terms not yet inserted in the database
\t-filteredBarcodeDir      path to the directory containing all filtered barcode to celltype files
\t-pipelineCallsSummary    path to the file containing a summary of processing calls info at library/celltype level
\t-pipelineReportFile      path to the file containing summary of kallisto (reads, reads mapped, ..,) and stats about read length
\t-kallistoResults         path to dir containing kallisto/bustools results for all libraries
\t-callsResults            path to dir containing calls results for all libraries
\t-sexInfo                 file containing sex-related info about anatomical terms
\t-numberCore              number of threads used to insert data in the database.
\t-debug                   (optional) insertions are not made, just printed
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# initialize variables

my $barcodeToCelltypeFilePattern = "${filteredBarcodeDir}/scRNASeq_barcode_EXP_ID.tsv";
my $ABSENT_CALL_NOT_RELIABLE = "absent call not reliable";
my $NOT_EXCLUDED = "not excluded";

########################################################################
############################ Queries ###################################
########################################################################

# SELECT QUERIES
my $select_datasource = "SELECT dataSourceName, dataSourceId FROM dataSource ".
                        "WHERE category =\'Single-cell RNA-Seq data source\'";

my $select_populationCapture = "SELECT rnaSeqPopulationCaptureId FROM rnaSeqPopulationCapture";

my $select_biotype = "SELECT geneBioTypeId, geneBioTypeName FROM geneBioType";

# retrieve experiment already inserted to not reinsert already existing ones (e.g in case
# of experiment containing a mix of full length and target based RNA-Seq)
my $select_experiments = 'select rnaSeqExperimentId from rnaSeqExperiment';

# INSERT QUERIES
my $insert_experiment = 'INSERT INTO rnaSeqExperiment (rnaSeqExperimentId,'.
                        'rnaSeqExperimentName, rnaSeqExperimentDescription, dataSourceId)'.
                        ' VALUES (?, ?, ?, ?)';

my $insert_libraries =  'INSERT INTO rnaSeqLibrary (rnaSeqLibraryId, rnaSeqExperimentId,'.
                        'rnaSeqSequencerName, rnaSeqTechnologyName, rnaSeqTechnologyIsSIngleCell,'.
                        'sampleMultiplexing, libraryMultiplexing, strandSelection,'.
                        'cellCompartment, sequencedTranscriptPart, fragmentation,'.
                        'rnaSeqPopulationCaptureId, genotype, allReadsCount, mappedReadsCount,'.
                        ' minReadLength, maxReadLength, libraryType)'.
                        ' VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)';

my $insert_annotatedSamples =   'INSERT INTO rnaSeqLibraryAnnotatedSample (rnaSeqLibraryId, conditionId,'.
                                'cellTypeAuthorAnnotation, anatEntityAuthorAnnotation, stageAuthorAnnotation,'.
                                'abundanceUnit, meanAbundanceReferenceIntergenicDistribution,'.
                                'sdAbundanceReferenceIntergenicDistribution, abundanceThreshold,'.
                                'allGenesPercentPresent, proteinCodingGenesPercentPresent,'.
                                'intergenicRegionsPercentPresent, pValueThreshold, allUMIsCount, mappedUMIsCount,'.
                                'multipleLibraryIndividualSample, barcode, time, timeUnit, freeTextAnnotation)'.
                                ' VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)';

my $update_sumUMIs_annotatedSamples =  'UPDATE rnaSeqLibraryAnnotatedSample set mappedUMIsCount = ? where '.
                                        'rnaSeqLibraryAnnotatedSampleId = ?';

my $select_annotatedSampleId =  'SELECT rnaSeqLibraryAnnotatedSampleId FROM '.
                                'rnaSeqLibraryAnnotatedSample WHERE rnaSeqLibraryId = ? AND '.
                                'conditionId = ? AND cellTypeAuthorAnnotation = ?';

my $insert_individualSamples =  'INSERT INTO rnaSeqLibraryIndividualSample (rnaSeqLibraryAnnotatedSampleId,'.
                                'barcode, sampleName) VALUES (?, ?, ?)';

my $update_sumUMIs_individualSamples =  'UPDATE rnaSeqLibraryIndividualSample set mappedUMIsCount = ? where '.
                                        'rnaSeqLibraryIndividualSampleId = ?';

my $select_individualSampleId = 'SELECT rnaSeqLibraryIndividualSampleId FROM '.
                                'rnaSeqLibraryIndividualSample WHERE rnaSeqLibraryAnnotatedSampleId = ? AND '.
                                'barcode = ? and sampleName = ?';

my $insert_run = 'INSERT INTO rnaSeqRun (rnaSeqRunId, rnaSeqLibraryId) VALUES (?, ?)';

my $insert_annotatedSampleGeneResult =  'INSERT INTO rnaSeqLibraryAnnotatedSampleGeneResult ('.
                                        'rnaSeqLibraryAnnotatedSampleId, bgeeGeneId, abundanceUnit, abundance,'.
                                        'readsCount, UMIsCount, zScore, pValue,'.
                                        'reasonForExclusion) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)';

my $insert_individualSampleGeneResult = 'INSERT INTO rnaSeqLibraryIndividualSampleGeneResult ('.
                                        'rnaSeqLibraryIndividualSampleId, bgeeGeneId, abundanceUnit, abundance,'.
                                        'readsCount, UMIsCount, rnaSeqData, reasonForExclusion)'.
                                        'VALUES (?, ?, ?, ?, ?, ?, ?, ?)';

########################################################################
########################### Main Script ################################
########################################################################

## Initialize Bgee db connections
# connection used to insert data in rnaSeqExperiment, rnaSeqLibrary, rnaSeqAnnotatedSample
# and rnaSeqIndividualSample tables. Keep AutoCommit as each insert has to be done directly
# to retrieve autoincrement IDs of conditionId, rnaSeqLibraryAnnotatedSampleId and
# rnaSeqLibraryIndividualSampleId
my $bgee_metadata = Utils::connect_bgee_db($bgee_connector);
$bgee_metadata->{AutoCommit} = 1;
# AutoInactiveDestroy = 1 avoid the connection to be automatically destroyed at the end of a connection
# in a ForkManager. Then have to manually close the connection at the end of the script.
$bgee_metadata->{AutoInactiveDestroy} = 1;
# automatically reconnect when connection is lost
$bgee_metadata->{mysql_auto_reconnect} = 1;

## Load 
# Library info from manual curation used to launch the pipeline
my %libraries         = getTargetBaseCuratedLibrariesAnnotation($targetBaseLibrary);
print "\t", scalar keys %libraries, " experiments with libraries mapped.\n";
# Info of processed libraries coming from the pipeline
my %processedLibraries = get_processed_libraries_info($bgeeLibraryInfo, 1);
# Experiment annotation coming from flat files
my @experimentType = ('3\'end', 'Full-length and 3\'end', '5\'end');
my %experiments       = getSingleCellExperiments($singleCellExperiment,
    @experimentType);
# Stats of pipeline processing for each library/celltype (libraryAnnotatedSample level)
my %callsPipelineSummary = getCallsSummaryAtLibraryAnnotatedLevel($pipelineCallsSummary);
my %pipelineReport = getAllRnaSeqReportInfo($pipelineReportFile);

# sex-related information needed for sub 'insert_get_condition'
my $anatSexInfo    = Utils::get_anat_sex_info($sexInfo);
my $speciesSexInfo = Utils::get_species_sex_info($bgee_metadata);

# Parse extra mapping info for currently too up-to-date annotations
##UnmappedId    UnmappedName    UberonID    UberonName    Comment
my %extra = map  { my @tmp = split(/\t/, $_, -1); if ( $tmp[2] ne '' && $tmp[0] ne '' ){ $tmp[0] => $tmp[2] } else { 'nonono' => 'nonono' } }
            grep { !/^#/ }
            read_file("$extraMapping", chomp => 1);

################
# DATA SOURCES #
################
my %bgeeDataSources = ();
my $selSrc = $bgee_metadata->prepare($select_datasource);
$selSrc->execute()  or die $selSrc->errstr;
while ( my @data = $selSrc->fetchrow_array ){
    $bgeeDataSources{$data[0]} = $data[1];
}
$selSrc->finish;

######################
# POPULATION CAPTURE #
######################
my @populationCapture = ();
my $selPopCapture = $bgee_metadata->prepare($select_populationCapture);
$selPopCapture->execute()  or die $selSrc->errstr;
while ( my @data = $selPopCapture->fetchrow_array ){
    push(@populationCapture, $data[0])
}
$selPopCapture->finish;

############
# BIOTYPES #
############

my %biotypeNameToBiotypeId = ();
my $selBiotypes = $bgee_metadata->prepare($select_biotype);
$selBiotypes->execute()  or die $selBiotypes->errstr;
while ( my @data = $selBiotypes->fetchrow_array ){
    $biotypeNameToBiotypeId{$data[1]} = $data[0];
}
$selBiotypes->finish;

## retrieve already inserted experiment Ids to not try to reinsert them
my $selExp = $bgee_metadata->prepare($select_experiments);
$selExp->execute() or die $selExp->errstr;
my @insertedExpIds = ();
while ( my @data = $selExp->fetchrow_array ){
    push(@insertedExpIds, $data[0]);
}


################################################
# GET GENE INTERNAL IDS                        #
# GET ORGAN, STAGE, AND CONDITIONS INFORMATION #
################################################

my %genes;
# go over all libraries to check all species with data
for my $expId ( keys %processedLibraries ){
    foreach my $libraryId ( keys %{$processedLibraries{$expId}} ){
        $genes{$libraries{$expId}->{$libraryId}->{'speciesId'}} = ();
    }
}
# Get hash of geneId to bgeeGeneId mapping per species
for my $speciesId ( keys %genes ){
    $genes{$speciesId} = Utils::query_bgeeGene($bgee_metadata, $speciesId);
}

# Get already known conditions
my $conditions = Utils::query_conditions($bgee_metadata);

# Get simpler (upper level) stage equivalences
my $stage_equivalences = Utils::get_stage_equivalences($bgee_metadata);


####################
# START INSERTION  #
####################

my $inserted = 0;

my $insExp = $bgee_metadata->prepare($insert_experiment);
my $insLib = $bgee_metadata->prepare($insert_libraries);


for my $expId ( sort keys %processedLibraries ){
    if (grep( /^$expId$/, @insertedExpIds)) {
        print "\t$expId already inserted. Will anyway insert not yet inserted libraries from this experiment\n";
    } else {
        print "\tinsert $expId\n";
        if ( $debug ){
            binmode(STDOUT, ':utf8');

            print 'INSERT INTO rnaSeqExperiment: ',
                $expId, ' - ', $experiments{$expId}->{'name'}, ' - ',
                $experiments{$expId}->{'description'}, ' - ',
                $bgeeDataSources{$experiments{$expId}->{'source'}}, "\n";
            if(!exists($bgeeDataSources{$experiments{$expId}->{'source'}})) {
                print "databasource ".$experiments{$expId}->{'source'}." used for annotation of experiment ".
                $expId." does not exist in the database";
            }
        }
        else {
            $insExp->execute($expId, $experiments{$expId}->{'name'},
                $experiments{$expId}->{'description'},
                $bgeeDataSources{$experiments{$expId}->{'source'}}) or die $insExp->errstr;
        }
    }
    # Load barcode to cell type mapping for this experiment
    my $barcodeToCellTypeFile = $barcodeToCelltypeFilePattern =~ s/EXP_ID/$expId/r;
    #my %barcodesTsv = %{ Utils::read_spreadsheet("$barcodeToCellTypeFile", "\t", 'csv', '"', 1) };
    my %barcodesToCellType = getBarcodeToCellType($barcodeToCellTypeFile);
    print "Start inserting libraries for $expId...\n";

    ## parse libraries a first time in order to be able to insert libraries, annotated samples, individual samples and conditions
    ## it is done to avoid key collision if creating them in parallel
    LIBRARY:
    for my $libraryId ( sort keys %{$processedLibraries{$expId}} ){

        if (!exists $pipelineReport{$libraryId}) {
            warn "pipeline report file not generated for library $libraryId. The library is not inserted.";
            next;
        }
        if (!exists $callsPipelineSummary{$libraryId}) {
            warn "calls summary file not generated for library $libraryId. The library is not inserted.";
            next;
        }
        # prepare queries
        my $insAnnotatedSample = $bgee_metadata->prepare($insert_annotatedSamples);
        my $selectAnnotatedSampleId = $bgee_metadata->prepare($select_annotatedSampleId);

        #No need to check already inserted libraries here as the pipeline removed them from the list of libraries to process
        print "\tInsert $libraryId and conditions from $expId\n";
        # For now all target base are polyA. Check that population capture polyA is already present in the database
        if ( !grep(/^polyA$/, @populationCapture ) ) {
            die "polyA population capture not found in the database";
        }
        # insert libraries
        # retrieve one random runId of this library in order to access to libraryType info
        my $runId = ();
        for my $key (sort keys %{ $processedLibraries{$expId}->{$libraryId}->{'runIds'} }) {
            if ($key ne "speciesId") {
                $runId = $key;
                last;
            }
        }
        my $libraryType = $processedLibraries{$expId}->{$libraryId}->{'runIds'}->{$runId}->{'libraryType'};
        if ( $debug ){
            print 'INSERT INTO rnaSeqLibrary: ', $libraryId,                              ' - ',
                  $expId, ' - ', $libraries{$expId}->{$libraryId}->{'platform'},          ' - ',
                  #technologyName
                  $libraries{$expId}->{$libraryId}->{'protocol'},                         ' - ',
                  # isSingleCell (always true for target base RNASeq)
                  '1',                                                                    ' - ',
                  # sampleMultiplexing (always true for target base)
                  '1',                                                                    ' - ',
                  # libraryMultiplexing (for now always false for target base)
                  '0',                                                                    ' - ',
                  # strandSelection. No information in annotation file.
                  # Could maybe be detected from the platform or the technology
                  'NA',                                                                   ' - ',
                  # cellCompartment. No information in the annotation file.
                  # Could maybe be detected from the technology name
                  $libraries{$expId}->{$libraryId}->{'cellCompartment'},                  ' - ',
                  $libraries{$expId}->{$libraryId}->{'sequencedTranscriptPart'},          ' - ',
                  # fragmentation.
                  # TODO. Left empty for now but should be updated to take read length info
                  '0',                                                                    ' - ',
                  # rnaSeqPopulationCaptureId. For now all target base are polyA
                  'polyA',                                                                ' - ',
                  $libraries{$expId}->{$libraryId}->{'genotype'},                         ' - ',
                  $pipelineReport{$libraryId}->{'allReadsCount'},                         ' - ',
                  $pipelineReport{$libraryId}->{'mappedReadsCount'},                      ' - ',
                  $pipelineReport{$libraryId}->{'minReadLength'},                         ' - ',
                  $pipelineReport{$libraryId}->{'maxReadLength'},                         ' - ',
                  $libraryType,"\n";
        } else {
            $insLib->execute($libraryId,
                            $expId,
                            $libraries{$expId}->{$libraryId}->{'platform'},
                            $libraries{$expId}->{$libraryId}->{'protocol'},
                            1,
                            1,
                            0,
                            'NA',
                            $libraries{$expId}->{$libraryId}->{'cellCompartment'},
                            $libraries{$expId}->{$libraryId}->{'sequencedTranscriptPart'},
                            0,
                            'polyA',
                            $libraries{$expId}->{$libraryId}->{'genotype'},
                            $pipelineReport{$libraryId}->{'allReadsCount'},
                            $pipelineReport{$libraryId}->{'mappedReadsCount'},
                            $pipelineReport{$libraryId}->{'minReadLength'},
                            $pipelineReport{$libraryId}->{'maxReadLength'},
                            $libraryType
                        )  or die $insLib->errstr;
        }

        # Now start to insert conditions
        # first remap ontology terms if required
        my $anatEntityId = $libraries{$expId}->{$libraryId}->{'anatEntityId'};
        if (exists $extra{$anatEntityId}) {

            $anatEntityId = $extra{$anatEntityId};
            $libraries{$expId}->{$libraryId}->{'anatEntityId'} = $anatEntityId;
        }
        my $stageId = $libraries{$expId}->{$libraryId}->{'stageId'};
        if (exists $extra{$stageId}) {
            $stageId = $extra{$stageId};
            $libraries{$expId}->{$libraryId}->{'stageId'} = $stageId;
        }
    
        my %clusterToAnnotatedSampleId;
        for my $clusterId (sort keys %{$barcodesToCellType{$libraryId}{'clusters'}} ) {
            my $cellTypeId = $barcodesToCellType{$libraryId}{'clusters'}{$clusterId}{'cellTypeId'};
            my $authorCellTypeAnnotation = $barcodesToCellType{$libraryId}{'clusters'}{$clusterId}{'authorCellTypeAnnotation'};

            # Get conditionId/exprMappedConditionId for this library
            # Updates also the hash of existing conditions
            my $condKeyMap;
            if ($debug) {
                print 'If condition does not already exist run insert into cond',
                $anatEntityId,                                         ' - ',
                $stageId,                                              ' - ',
                $cellTypeId,                                           ' - ',
                $libraries{$expId}->{$libraryId}->{'sex'},             ' - ',
                $libraries{$expId}->{$libraryId}->{'strain'},          ' - ',
                $libraries{$expId}->{$libraryId}->{'speciesId'}, "\n";
            } else {
                ($condKeyMap, $conditions) = Utils::insert_get_condition($bgee_metadata,
                                 $conditions,
                                 $stage_equivalences,
                                 $anatEntityId,
                                 $stageId,
                                 $libraries{$expId}->{$libraryId}->{'speciesId'},
                                 $libraries{$expId}->{$libraryId}->{'sex'},
                                 $libraries{$expId}->{$libraryId}->{'strain'},
                                 $anatSexInfo, $speciesSexInfo,
                                 $libraryId, '',
                                 $cellTypeId
                                );
    
            }

            # now insert annotatedSamples
            my $annotatedSampleId = ();

            # generate path to file containing calls for this library/celusterId
            my $pathToCallFile = "${callsResults}/${libraryId}/Calls_cellPop_".
                "${libraryId}_${clusterId}_genic.tsv";
            if (! -e $pathToCallFile) {
                # only insert annotated sample Id and condition Id
                $annotatedSampleId = insert_get_annotated_sample($libraryId, $condKeyMap->{'conditionId'},
                    '', '', '', 'cpm', 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, '' , undef, '', '', $insAnnotatedSample,
                    $selectAnnotatedSampleId, $debug);
            } else {
                # first insert annotated sample with metrics from our pipeline
                $annotatedSampleId = insert_get_annotated_sample($libraryId, $condKeyMap->{'conditionId'},
                    $authorCellTypeAnnotation, $libraries{$expId}->{$libraryId}->{'authorAnatEntityAnnotation'},
                    $libraries{$expId}->{$libraryId}->{'authorStageAnnotation'}, 'cpm',
                    $callsPipelineSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'meanRefIntergenic'},
                    $callsPipelineSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'sdRefIntergenic'},
                    $callsPipelineSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'abundanceThreshold'},
                    $callsPipelineSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'allGenesPercentPresent'},
                    $callsPipelineSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'proteinCodingGenesPercentPresent'},
                    $callsPipelineSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'intergenicRegionsPercentPresent'},
                    $callsPipelineSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'pValueThreshold'},
                    # allUMIsCount -> 0, mappedUMIsCount -> 0 as they will be updated once all gene result have been processed
                    0, 0,
                    ## multipleLibraryIndividualSample -> 1, barcode -> '', time -> undef, timeUnit -> '', freeTextAnnotation -> '' 
                    1, '', undef, '', '', $insAnnotatedSample, $selectAnnotatedSampleId, $debug);

                # then load file containing all calls for this library/celltypeId
                my %callsOneAnnotatedSample = getCallsInfoPerLibrary($pathToCallFile);
                
                # used to sum all UMI and then update the table annotatedSample to populate the mappedUMIsCount column
                my $mappedUMIsAnnotatedSample = 0;
            }

            # map the celltypeId to the corresponding annotated sample 
            $clusterToAnnotatedSampleId{$authorCellTypeAnnotation}{$cellTypeId} = $annotatedSampleId;
        }

        #prepare query to update sumUMIs per rnaSeqLibraryAnnotatedSample
        $selectAnnotatedSampleId->finish;
        $insAnnotatedSample->finish;

        # finished insertion on anotated sample. Now start to insert individual samples
        my $insIndividualSample = $bgee_metadata->prepare($insert_individualSamples);
        my $selectIndividualSampleId = $bgee_metadata->prepare($select_individualSampleId);
        ## Now start to insert individual samples if needed
        # for now we only insert barcodes mapped to a cell type. It could be possible
        # to insert other cell types in a different table
        # (e.g rnaSeqLibraryIndevidualSampleNotAnnotated(rnaSeqLibraryId, barcode, ...))
        my %barcodeToIndividualSampleId;
        my @barcodesArray;
        for my $barcode (sort keys %{$barcodesToCellType{$libraryId}{'barcodes'}}) {
            my $cellTypeId = $barcodesToCellType{$libraryId}{'barcodes'}{$barcode}{'cellTypeId'};
            my $authorCellTypeAnnotation = $barcodesToCellType{$libraryId}{'barcodes'}{$barcode}{'authorCellTypeAnnotation'};
            my $annotatedSampleId = $clusterToAnnotatedSampleId{$authorCellTypeAnnotation}{$cellTypeId};
            my $individualSampleId = insert_get_individual_sample($insIndividualSample, $selectIndividualSampleId,
                $annotatedSampleId, $barcode, "", $debug);
            $barcodeToIndividualSampleId{$barcode} = $individualSampleId;
            push (@barcodesArray, $barcode);
        }

        ## close prepared statements.
        $selectIndividualSampleId->finish;
        $insIndividualSample->finish;

    }

    # parse libraries a second time in parallel to fasten insertion
    my $pm = new Parallel::ForkManager($numberCore);
    for my $libraryId ( sort keys %{$processedLibraries{$expId}} ){
        #start parallelization
        my $pid = $pm->start and next;
        if (!exists $pipelineReport{$libraryId} || !exists $callsPipelineSummary{$libraryId}) {
            $pm->finish;
            next;
        }
        #thread specific connection to the DB with autocommit
        my $bgee_thread = Utils::connect_bgee_db($bgee_connector);
        $bgee_thread->{AutoCommit} = 1;
        #TODO remove these 2 connection and manually commit after each insertion
        # bgee connection with autocommit disabled to be able to manually commit.
        my $bgee_data = Utils::connect_bgee_db($bgee_connector);
        $bgee_data->{AutoCommit} = 0;

        # prepare queries
        my $selectAnnotatedSampleId = $bgee_thread->prepare($select_annotatedSampleId);
        my $insAnnotatedSampleGeneResult = $bgee_data->prepare($insert_annotatedSampleGeneResult);


        # read count sparse matrix for all barcodes and genes of the library. It
        # corresponds to raw data per cell coming from kallisto/bustools. There was
        # no postprocessing filtering based on barcodes or celltype
        my %cpmMatrix =read_sparse_matrix("$kallistoResults/$libraryId/gene_counts", "gene");
        my %countMatrix = read_sparse_matrix("$kallistoResults/$libraryId/cpm_counts", "cpm_counts");

        # Now start to insert annotated samples
        my %clusterToAnnotatedSampleId;
        # init counter to update sumUMIs per annotated sample
        my %updateSumUMIs;

        for my $clusterId (sort keys %{$barcodesToCellType{$libraryId}{'clusters'}} ) {
            my $cellTypeId = $barcodesToCellType{$libraryId}{'clusters'}{$clusterId}{'cellTypeId'};
            my $authorCellTypeAnnotation = $barcodesToCellType{$libraryId}{'clusters'}{$clusterId}{'authorCellTypeAnnotation'};
            # condition were all inserted previously. Now we only retrieve $condKeyMap
            my $condKeyMap = ();
            if(! $debug) {
                ($condKeyMap, $conditions) = Utils::insert_get_condition($bgee_metadata,
                                 $conditions,
                                 $stage_equivalences,
                                 $libraries{$expId}->{$libraryId}->{'anatEntityId'},
                                 $libraries{$expId}->{$libraryId}->{'stageId'},
                                 $libraries{$expId}->{$libraryId}->{'speciesId'},
                                 $libraries{$expId}->{$libraryId}->{'sex'},
                                 $libraries{$expId}->{$libraryId}->{'strain'},
                                 $anatSexInfo, $speciesSexInfo,
                                 $libraryId, '',
                                 $cellTypeId
                                );
    
            }
            my $conditionId = $condKeyMap->{'conditionId'};
            # boolean used to define if calls can be inserted for this library/celltype
            my $insertCalls = 1;


            my $annotatedSampleId = ();

            # generate path to file containing calls for this library/celusterId
            my $pathToCallFile = "${callsResults}/${libraryId}/Calls_cellPop_".
                "${libraryId}_${clusterId}_genic.tsv";
            if (! -e $pathToCallFile) {
                warn "calls file does not exist for library $libraryId and celltype $cellTypeId. ",
                "Condition and annotated sample have been inserted in order to be able to save info at ",
                "cell level but no calls will be generated.";
                $insertCalls = 0;
            } elsif ($callsPipelineSummary{$libraryId}{$cellTypeId}{$authorCellTypeAnnotation}{'abundanceThreshold'} eq "NA") {
                warn "Abundance threshold not available for library $libraryId, cellType $cellTypeId and cell-type author annotation",
                " $authorCellTypeAnnotation. Condition and annotated sample have been inserted in order to be able to save info at ",
                "cell level but no calls will be generated.";
                $insertCalls = 0;
            }
            if ($debug) {
                print "SELECT rnaSeqLibraryAnnotatedSampleId FROM ".
                      "rnaSeqLibraryAnnotatedSample WHERE rnaSeqLibraryId = $libraryId AND ".
                      "conditionId = conditionId AND cellTypeAuthorAnnotation = $authorCellTypeAnnotation";
                $annotatedSampleId = "annotatedSampleId";
            } else {
                $selectAnnotatedSampleId->execute($libraryId, $conditionId, $authorCellTypeAnnotation) or die $selectAnnotatedSampleId->errstr;
                my @insertedExpIds = ();
                my @data = $selectAnnotatedSampleId->fetchrow_array;
                if (scalar @data == 0) {
                    die "Did not find any annotated sample for library $libraryId, conditionId $conditionId and author annotation $authorCellTypeAnnotation";
                } elsif (scalar @data > 1) {
                    die "Found more than one annotated sample for library $libraryId, conditionId $conditionId and author annotation $authorCellTypeAnnotation";
                }
                $annotatedSampleId = $data[0];
            }
            if ($insertCalls) {
                # load file containing all calls for this library/celltypeId
                my %callsOneAnnotatedSample = getCallsInfoPerLibrary($pathToCallFile);
                
                # sum all UMI and then update the table annotatedSample to populate the mappedUMIsCount column
                my $mappedUMIsAnnotatedSample = 0;

                #Start transaction with TRANSACTION ISOLATION LEVEL READ UNCOMMITTED
                for my $geneId (sort keys %callsOneAnnotatedSample) {
                    $mappedUMIsAnnotatedSample += $callsOneAnnotatedSample{$geneId}->{'sumUMI'};
                    my $bgeeGeneId = $genes{$libraries{$expId}->{$libraryId}->{'speciesId'}}{$geneId};
                    my $reasonForExclusion = $NOT_EXCLUDED;
                    if ($callsOneAnnotatedSample{$geneId}{'pValue'} > 0.05) {
                        $reasonForExclusion = $ABSENT_CALL_NOT_RELIABLE;
                    }
                    if ($debug) {
                        my $zScoreNullable = $callsOneAnnotatedSample{$geneId}{'zScore'};
                        if (!defined $zScoreNullable) {
                            $zScoreNullable = 'undef'
                        }
                        print 'INSERT INTO rnaSeqLibraryAnnotatedSampleGeneResult: ',
                        'annotatedSampleId', ' - ', $bgeeGeneId, ' - ', "cpm", ' - ',
                        ## the 0 is for the readsCount. The only info we have for target base
                        ## is the UMICounts so weleave readCounts at 0
                        $callsOneAnnotatedSample{$geneId}{'cpm'}, ' - ', 0, ' - ',
                        $callsOneAnnotatedSample{$geneId}{'sumUMI'}, ' - ',
                        $zScoreNullable, ' - ',
                        $callsOneAnnotatedSample{$geneId}{'pValue'}, ' - ',
                        $reasonForExclusion, "\n";
                    } else {
                        $insAnnotatedSampleGeneResult->execute($annotatedSampleId, $bgeeGeneId,
                        ## the 0 is for the readsCount. The only info we have for target base
                        ## is the UMICounts so weleave readCounts at 0
                        'cpm', $callsOneAnnotatedSample{$geneId}{'cpm'}, 0,
                        $callsOneAnnotatedSample{$geneId}{'sumUMI'},
                        $callsOneAnnotatedSample{$geneId}{'zScore'},
                        $callsOneAnnotatedSample{$geneId}{'pValue'},
                        $reasonForExclusion)
                            or die $insAnnotatedSampleGeneResult->errstr;
                    }
                }
                #now update annotated sample to insert sumUMIs
                $updateSumUMIs{$annotatedSampleId} = $mappedUMIsAnnotatedSample;
                #commit after each barcode
                $bgee_data->commit;
            }

            # map the celltypeId to the corresponding annotated sample 
            $clusterToAnnotatedSampleId{$authorCellTypeAnnotation}{$cellTypeId} = $annotatedSampleId;
        }

        #prepare query to update sumUMIs per rnaSeqLibraryAnnotatedSample
        my $updateUMIAnnotatedSample = $bgee_thread->prepare($update_sumUMIs_annotatedSamples);

        for my $annotatedSampleIdToUpdate (sort keys %updateSumUMIs) {
            if ($debug) {
                print "UDPATE rnaSeqLibraryAnnotatedSample SET mappedUMIsCount = $updateSumUMIs{$annotatedSampleIdToUpdate}".
                " WHERE rnaSeqLibraryAnnotatedSampleId = 'annotatedSampleId';\n";
            } else {
                $updateUMIAnnotatedSample->execute($updateSumUMIs{$annotatedSampleIdToUpdate}, $annotatedSampleIdToUpdate);
            }            
        }
        $selectAnnotatedSampleId->finish;
        $updateUMIAnnotatedSample->finish;
        $insAnnotatedSampleGeneResult->finish;

        # finished insertion on anotated sample. Now start to insert individual samples
        my $selectIndividualSampleId = $bgee_thread->prepare($select_individualSampleId);
        my $insIndividualSampleGeneResult = $bgee_data->prepare($insert_individualSampleGeneResult);
        ## Now start to insert individual samples if needed
        # for now we only insert barcodes mapped to a cell type. It could be possible
        # to insert other cell types in a different table
        # (e.g rnaSeqLibraryIndevidualSampleNotAnnotated(rnaSeqLibraryId, barcode, ...))
        my %barcodeToIndividualSampleId;
        my @barcodesArray;
        for my $barcode (sort keys %{$barcodesToCellType{$libraryId}{'barcodes'}}) {
            my $cellTypeId = $barcodesToCellType{$libraryId}{'barcodes'}{$barcode}{'cellTypeId'};
            my $authorCellTypeAnnotation = $barcodesToCellType{$libraryId}{'barcodes'}{$barcode}{'authorCellTypeAnnotation'};
            my $annotatedSampleId = $clusterToAnnotatedSampleId{$authorCellTypeAnnotation}{$cellTypeId};
            my $individualSampleName = '';
            if ($debug) {
                print "SELECT rnaSeqLibraryIndividualSampleId FROM ".
                      "rnaSeqLibraryIndividualSample WHERE rnaSeqLibraryAnnotatedSampleId = annotatedSampleId AND ".
                      "barcode = $barcode and sampleName = $individualSampleName";
            }
            $selectIndividualSampleId->execute($annotatedSampleId, $barcode, '') or die $selectIndividualSampleId->errstr;
            my @data = $selectIndividualSampleId->fetchrow_array;
            if (scalar @data == 0) {
                die "Did not find any individual sample for annotated sample $annotatedSampleId, barcode $barcode and name $individualSampleName";
            } elsif (scalar @data > 1) {
                die "Found more than one individual sample for annotated sample $annotatedSampleId, barcode $barcode and name $individualSampleName";
            }
            my $individualSampleId = $data[0];
            $barcodeToIndividualSampleId{$barcode} = $individualSampleId;
            push (@barcodesArray, $barcode);
        }

        # then start to insert abundance per cell. parallelized per celltype
        my %mappedUMIsIndividualSamples;
        foreach my $barcode (@barcodesArray) {

            # For each barcode start a transaction with TRANSACTION ISOLATION LEVEL READ UNCOMMITTED
            my $individualSampleId = $barcodeToIndividualSampleId{$barcode};

            my %subsetCpmMatrix = %cpmMatrix{$barcode};
            my %subsetCountMatrix = %countMatrix{$barcode};

            for my $geneId ( keys %{$subsetCountMatrix{$barcode}} ) {
                # check that the gene is present in the database. It is both a
                # security check and a way to remove intergenic regions
                #TODO: using a for loop on %subsetCpmMatrix would avoid parsing intergenic. To test in a future release
                next if (! exists $genes{$libraries{$expId}->{$libraryId}->{'speciesId'}}{$geneId} ||
                    ! exists $subsetCountMatrix{$barcode}{$geneId});
                my $bgeeGeneId = $genes{$libraries{$expId}->{$libraryId}->{'speciesId'}}{$geneId};
                if (! exists $subsetCpmMatrix{$barcode}{$geneId}) {
                    warn "Warning, gene $geneId has count for barcode $barcode but",
                        " no abundance was generated";
                    next;
                }
                if (exists $mappedUMIsIndividualSamples{$individualSampleId}) {
                    $mappedUMIsIndividualSamples{$individualSampleId} += $subsetCountMatrix{$barcode}{$geneId};
                } else {
                    $mappedUMIsIndividualSamples{$individualSampleId} = $subsetCountMatrix{$barcode}{$geneId};
                }
                if ($debug) {
                    print 'INSERT INTO rnaSeqLibraryIndividualSampleGeneResult: ',
                    $individualSampleId, ' - ', $bgeeGeneId, ' - ', "cpm", ' - ',
                    $subsetCpmMatrix{$barcode}{$geneId}, ' - ', 0, ' - ',
                    $subsetCountMatrix{$barcode}{$geneId}, ' - ',
                    "high quality", ' - ', 'not excluded', "\n";
                } else {
                    $insIndividualSampleGeneResult->execute($individualSampleId, $bgeeGeneId,
                        'cpm', $subsetCpmMatrix{$barcode}{$geneId}, 0,
                        $subsetCountMatrix{$barcode}{$geneId}, 'high quality', 'not excluded')
                            or die $insIndividualSampleGeneResult->errstr;
                }
            }
            #commit after each barcode
            $bgee_data->commit;
        }

        #now update the table rnaSeqIndividualSample to add mappedUMIsCount
        my $updateUMIIndividualSample = $bgee_data->prepare($update_sumUMIs_individualSamples);
        foreach my $individualSampleId (sort keys %mappedUMIsIndividualSamples) {
            $updateUMIIndividualSample->execute($individualSampleId, $mappedUMIsIndividualSamples{$individualSampleId});
        }
        #commit only once update of all indifivudalsamples
        $bgee_data->commit;
        ## close prepared statements.
        $selectIndividualSampleId->finish;
        $updateUMIIndividualSample->finish;
        # $insAnnotatedSample->finish;
        # $selectAnnotatedSampleId->finish;
        $insIndividualSampleGeneResult->finish;
        $bgee_data->disconnect;
        $bgee_thread->disconnect;
        $pm->finish;




    }
    $pm->wait_all_children
}

$bgee_metadata->disconnect;
exit 0;

