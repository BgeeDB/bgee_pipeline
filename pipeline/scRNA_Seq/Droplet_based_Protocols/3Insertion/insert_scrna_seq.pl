#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use File::Slurp;
use lib "$FindBin::Bin/../../../"; # Get lib path for Utils.pm
##use Utils;
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
my ($singleCellExperiment, $bgeeLibraryInfo, $sourceDir) = ('', '', '');
my ($pipelineCallsSummary, $pipelineReportFile) = ('', '');
my $numberCore = 1;
#my ($library_stats, $report_info = ('', '');
my ($debug)                      = (0);
my %opts = ('bgee=s'                 => \$bgee_connector,       # Bgee connector string
            'targetBaseLibrary=s'    => \$targetBaseLibrary,    # target base RNAseq library annotations
            'singleCellExperiment=s' => \$singleCellExperiment, # single cell RNASeq experiment annotations
            'bgeeLibraryInfo=s'      => \$bgeeLibraryInfo,      # metadata_info_10X.txt file
            'pipelineCallsSummary=s' => \$pipelineCallsSummary, # path to the file containing a summary of processing calls info at library/celltype level (e.g percentage protein coding present, ...)
            'pipelineReportFile=s'   => \$pipelineReportFile    # path to the file containing summary of kallisto (reads, reads mapped, ..,) and stats about read length
            'kallistoResults=s'      => \$kallistoResults,      # path to dir containing kallisto/bustools results for all libraries
            'callsResults=s'         => \$callsResults,         # path to dir containing calls results for all libraries
            'sourceDir=s'            => \$sourceDir,            # path to the directory containing source files of the target base pipeline
            'sexInfo=s'              => \$sexInfo,              # generated_files/uberon/uberon_sex_info.tsv
            'numberCore=s'           => \$numberCore,           # number of cores corresponding to number of threads used to insert data in the database
            'debug'                  => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $targetBaseLibrary eq '' || $singleCellExperiment eq '' ||
    $bgeeLibraryInfo eq '' || $pipelineCallsSummary eq '' || $pipelineReportFile eq '' || $kallistoResults eq '' || $sourceDir eq '' ||
    $sexInfo eq '' || $numberCore eq '' || $callsResults eq ''){
    print "\n\tInvalid or missing argument:
\te.g., $0  -bgee=\$(BGEECMD) -targetBaseLibrary=RNASeqLibrary_full.tsv -singleCellExperiment=RNASeqExperiment_full.tsv -bgeeLibraryInfo=metadata_info_10X.txt -sexInfo=\$(UBERON_SEX_INFO_FILE_PATH) > $@.tmp 2>warnings.$@
\t-bgee                    Bgee connector string
\t-targetBaseLibrary       targetBaseLibrary annotation file
\t-singleCellExperiment    singleCellExperiment file
\t-bgeeLibraryInfo         metadata_info_10X.txt file
\t-pipelineCallsSummary    path to the file containing a summary of processing calls info at library/celltype level
\t-pipelineReportFile      path to the file containing summary of kallisto (reads, reads mapped, ..,) and stats about read length
\t-kallistoResults         path to dir containing kallisto/bustools results for all libraries
\t-callsResults            path to dir containing calls results for all libraries
\t-sourceDir               path to the directory containing sources files of the target bas pipeline
\t-sexInfo                 file containing sex-related info about anatomical terms
\t-numberCore              number of threads used to insert data in the database.
\t-debug                   (optional) insertions are not made, just printed
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../target_base_utils.pl");

# initialize variables

my $barcodeToCelltypeFilePattern = "$sourceDir/scRNASeq_barcode_EXP_ID.tsv";

########################################################################
############################ Queries ###################################
########################################################################

# SELECT QUERIES
my $select_datasource = "SELECT dataSourceName, dataSourceId FROM dataSource ".
                        "WHERE category =\'Single-cell RNA-Seq data source\'";

my $select_populationCapture = "SELECT rnaSeqPopulationCaptureId FROM rnaSeqPopulationCapture";

my $select_biotype = "SELECT geneBioTypeId, geneBioTypeName FROM geneBioType";

# INSERT QUERIES
my $insert_experiment = 'INSERT INTO rnaSeqExperimentDev (rnaSeqExperimentId,'.
                        'rnaSeqExperimentName, rnaSeqExperimentDescription, dataSourceId)'.
                        ' VALUES (?, ?, ?, ?)';

my $insert_libraries =  'INSERT INTO rnaSeqLibraryDev (rnaSeqLibraryId, rnaSeqExperimentId,'.
                        'rnaSeqSequencerName, rnaSeqTechnologyName, rnaSeqTechnologyIsSIngleCell,'.
                        'sampleMultiplexing, libraryMultiplexing, strandSelection,'.
                        'cellCompartment, sequencedTranscriptPart, fragmentation,'.
                        'rnaSeqPopulationCaptureId, genotype, allReadsCount, mappedReadsCount,'.
                        ' minReadLength, maxReadLength, libraryType)'.
                        ' VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)';

my $insert_annotatedSamples =   'INSERT INTO rnaSeqLibraryAnnotatedSampleId (rnaSeqLibraryId,'.
                                'conditionId, abundanceUnit,'.
                                'meanAbundanceReferenceIntergenicDistribution,'.
                                'sdAbundanceReferenceIntergenicDistribution, tmmFactor,'.
                                'abundanceThreshold, allGenesPercentPresent,'.
                                'proteinCodingGenesPercentPresent, intergenicRegionsPercentPresent,'.
                                'pValueThreshold, mappedUMIsCount, multipleLibraryIndividualSample)'.
                                ' VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)';

my $select_annotatedSampleId =  'SELECT rnaSeqLibraryAnnotatedSampleId FROM '.
                                'rnaSeqLibraryAnnotatedSampleDev WHERE conditionId = ? AND '.
                                'rnaSeqLibraryId = ?';

my $insert_individualSamples =  'INSERT INTO rnaSeqLibraryIndividualSampleDev (rnaSeqLibraryAnnotatedSampleId,'.
                                'barcode, sampleName) VALUES (?, ?, ?)';

my $select_individualSampleId = 'SELECT rnaSeqLibraryIndividualSampleId FROM '.
                                'rnaSeqLibraryIndividualSampleDev WHERE rnaSeqLibraryAnnotatedSampleId = ? AND '.
                                'barcode = ? and sampleName = ?';

my $insert_run = 'INSERT INTO rnaSeqRun (rnaSeqRunId, rnaSeqLibraryId) VALUES (?, ?)';

my $insert_annotatedSampleGeneResult =  'INSERT INTO rnaSeqLibraryAnnotatedSampleGeneResultDev ('.
                                        'rnaSeqLibraryAnnotatedSampleId, bgeeGeneId, abundanceUnit, abundance,'.
                                        'readsCount, UMIsCount, zScore, pValue, rnaSeqData,'.
                                        'reasonForExclusion) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)';

my $insert_individualSampleGeneResult = 'INSERT INTO rnaSeqLibraryIndividualSampleGeneResultDev ('.
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

## Load 
# Library info from manual curation used to launch the pipeline
my %libraries         = getTargetBaseCuratedLibrariesAnnotation($targetBaseLibrary);
print "\t", scalar keys %libraries, " experiments with libraries mapped.\n";
# Info of processed libraries coming from the pipeline
my %processedLibraries = get_processed_libraries_info($bgeeLibraryInfo);
# Experiment annotation coming from flat files
my @experimentType = ('3\'end', 'Full-length and 3\'end');
my %experiments       = getSingleCellExperiments($singleCellExperiment,
    @experimentType);
# Stats of pipeline processing for each library/celltype (libraryAnnotatedSample level)
my %callsPipelineSummary = getCallsSummaryAtLibraryAnnotatedLevel($pipelineCallsSummary);
my %pipelineReport = getAllRnaSeqReportInfo($pipelineReportFile);

# sex-related information needed for sub 'insert_get_condition'
my $anatSexInfo    = Utils::get_anat_sex_info($sexInfo);
my $speciesSexInfo = Utils::get_species_sex_info($bgee_metadata);

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

######################
# INSERT EXPERIMENTS #
######################
print "Inserting experiments...\n";
my $insExp = $bgee_metadata->prepare($insert_experiment);
for my $expId ( sort keys %processedLibraries ){
    print "\t$expId\n";
    if ( $debug ){
        binmode(STDOUT, ':utf8');
        print 'INSERT INTO rnaSeqExperiment: ',
            $expId, ' - ', $experiments{$expId}->{'name'}, ' - ',
            $experiments{$expId}->{'description'}, ' - ',
            $bgeeDataSources{$experiments{$expId}->{'source'}}, "\n";
    }
    else {
        $insExp->execute($expId, $experiments{$expId}->{'name'},
            $experiments{$expId}->{'description'},
            $bgeeDataSources{$experiments{$expId}->{'source'}}) or die $insExp->errstr;
    }
}
$insExp->finish();
print "Done inserting experiments\n\n";

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


###################################################
# INSERT LIBRARIES, CONDITION, ANNOTATED SAMPLES  #
###################################################
print "Inserting libraries condition and annotated samples...\n";
# prepare queries
my $insLib = $bgee_metadata->prepare($insert_libraries);
my $insAnnotatedSample = $bgee_metadata->prepare($insert_annotatedSamples);
my $selectAnnotatedSampleId = $bgee_metadata->prepare($select_annotatedSampleId);
my $insIndividualSample = $bgee_metadata->prepare($insert_individualSamples);
my $selectIndividualSampleId = $bgee_metadata->prepare($select_individualSampleId);
my $inserted = 0;

for my $expId ( sort keys %processedLibraries ){

    # Load barcode to cell type mapping
    my $barcodeToCellTypeFile = $barcodeToCelltypeFilePattern =~ s/EXP_ID/$expId/r;
    #my %barcodesTsv = %{ Utils::read_spreadsheet("$barcodeToCellTypeFile", "\t", 'csv', '"', 1) };
    my %barcodesToCellType = getBarcodeToCellType($barcodeToCellTypeFile);

    print "Start inserting libraries...\n";
    LIBRARY:
    for my $libraryId ( sort keys %{$processedLibraries{$expId}} ){
        # read count sparse matrix for all barcodes and genes of the library. It
        # corresponds to raw data per cell coming from kallisto/bustools. There was
        # no postprocessing filtering based on barcodes or celltype
        my %sparseMatrixCount = read_sparse_matrix("$kallistoResults/$libraryId/gene_counts", "gene");
        my %sparseMatrixCpm = read_sparse_matrix("$kallistoResults/$libraryId/cpm_counts", "cpm_counts");
        my %callsPerLibrary = getCallsInfoPerLibrary("$callsResults/$libraryId/Calls_$libraryId.tsv");

        print "\tInsert $libraryId from $expId\n";
        # For now all target base are polyA. Check that population capture polyA is already present in the database
        if ( !grep(/^polyA$/, @populationCapture ) ) {
            die "polyA population capture not found in the database";
        }

        # insert libraries
        if ( $debug ){
            print 'INSERT INTO rnaSeqLibrary: ', $libraryId,                                     ' - ',
                  $expId, ' - ', $libraries{$expId}->{$libraryId}->{'platform'},                 ' - ',
                  #technologyName
                  $libraries{$expId}->{$libraryId}->{'protocol'},                                ' - ',
                  # isSingleCell (always true for target base RNASeq)
                  '1',                                                                           ' - ',
                  # sampleMultiplexing (always true for target base)
                  '1',                                                                           ' - ',
                  # libraryMultiplexing (for now always false for target base)
                  '0',                                                                           ' - ',
                  # strandSelection. No information in annotation file.
                  # Could maybe be detected from the platform or the technology
                  'NA',                                                                          ' - ',
                  # cellCompartment. No information in the annotation file.
                  # Could maybe be detected from the technology name
                  $libraries{$expId}->{$libraryId}->{'cellCompartment'},                         ' - ',
                  $libraries{$expId}->{$libraryId}->{'sequencedTranscriptPart'},                 ' - ',
                  # fragmentation.
                  # TODO. Left empty for now but should be updated to take read length info
                  '0',                                                                           ' - ',
                  # rnaSeqPopulationCaptureId. For now all target base are polyA
                  'polyA',                                                                       ' - ',
                  $libraries{$expId}->{$libraryId}->{'genotype'},                                ' - ',
                  $pipelineReport{$libraryId}->{'allReadsCount'},
                  $pipelineReport{$libraryId}->{'mappedReadsCount'},
                  $pipelineReport{$libraryId}->{'minReadLength'},
                  $pipelineReport{$libraryId}->{'maxReadLength'},
                  $processedLibraries{$expId}->{$libraryId}->{'libraryType'},"\n";
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
                            $processedLibraries{$expId}->{$libraryId}->{'libraryType'}
                        )  or die $insLib->errstr;
        }

        # Now start to insert annotated samples
        # This script does not manage insertion of calls per cell type but only
        # counts at cell (= barcode) level. In the database schema, counts at cell level
        # ( table rnaSeqLibraryIndividualSample) are linked to a celltype (table rnaSeqLibraryAnnotatedSample).
        # In order to be able to insert the counts at cell level we first insert all
        # metadata at cell type level (table rnaSeqLibraryAnnotatedSample). However only mandatory
        # information are inserted (libraryId, rnaSeqLibraryAnnotatedSAmpleId, conditionId)
        # all information at celltype level that results from running our pipeline will
        # be inserted in an other script.
        my %celltypeToAnnotatedSampleId;
        for my $cellTypeId (sort keys %{$barcodesToCellType{$libraryId}{'cellTypes'}} ) {
            print "\t\tInsert celltype $cellTypeId for library $libraryId\n";
            # Get conditionId/exprMappedConditionId for this library
            # Updates also the hash of existing conditions
            my $condKeyMap;
            if ($debug) {
                print 'If condition does not already exist run insert into cond',
                $libraries{$expId}->{$libraryId}->{'anatEntityId'},    ' - ',
                $libraries{$expId}->{$libraryId}->{'stageId'},         ' - ',
                $cellTypeId,                                           ' - ',
                $libraries{$expId}->{$libraryId}->{'sex'},             ' - ',
                $libraries{$expId}->{$libraryId}->{'strain'},          ' - ',
                $libraries{$expId}->{$libraryId}->{'speciesId'}, "\n";
            } else {
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
            # We consider the fine-grained (low-level) conditionId for insertion of annotated sample: $condKeyMap->{'conditionId'}
            my $annotatedSampleId = insert_get_annotated_sample($insAnnotatedSample,
                $selectAnnotatedSampleId, 
                $callsPipelineSummary{libraryId}{cellTypeId}->{'abundanceThreshold'},
                $callsPipelineSummary{libraryId}{cellTypeId}->{'allGenesPercentPresent'},
                $callsPipelineSummary{libraryId}{cellTypeId}->{'proteinCodingGenesPercentPresent'},
                $callsPipelineSummary{libraryId}{cellTypeId}->{'intergenicRegionsPercentPresent'},
                $callsPipelineSummary{libraryId}{cellTypeId}->{'pValueThreshold'},
                $callsPipelineSummary{libraryId}{cellTypeId}->{'meanRefIntergenic'},
                $callsPipelineSummary{libraryId}{cellTypeId}->{'sdRefIntergenic'},
                $callsPipelineSummary{libraryId}{cellTypeId}->{'mappedUMIs'},
                ## 1 here means true for isSingleCell
                1, $condKeyMap->{'conditionId'}, $libraryId, $debug);


            # could create a hash celltype -> annotatedSampleId
            $celltypeToAnnotatedSampleId{$cellTypeId} = $annotatedSampleId;

            # Then insert rnaSeqAnnotatedSampleGeneResult
            # It is parallelized. Each thread will insert all calls for one annotatedSample
            my $pm = new Parallel::ForkManager($numberCore);
            for my $cellType (sort keys %callsPerLibrary) {
                my $pid = $pm->start and next;
                my $bgee_data = Utils::connect_bgee_db($bgee_connector);
                # disable autocommit for $bgee_data. Allows to manually commit after each library.
                $bgee_data->{AutoCommit} = 0;
                my $insAnnotatedSampleGeneResult = $bgee_data->prepare($insert_annotatedSampleGeneResult);
                for my $geneId (sort keys %{$callsPerLibrary{$cellType}} ) {
                    my $bgeeGeneId = $genes{$libraries{$expId}->{$libraryId}->{'speciesId'}}{$geneId};
                    if ($debug) {
                        print 'INSERT INTO rnaSeqLibraryAnnotatedSampleGeneResult: ',
                        $annotatedSampleId, ' - ', $bgeeGeneId, ' - ', "cpm", ' - ',
                        $callsPerLibrary{$cellType}{$geneId}{'cpm'}, ' - ', 0, ' - ',
                        $callsPerLibrary{$cellType}{$geneId}{'sumUMI'}, ' - ',
                        $callsPerLibrary{$cellType}{$geneId}{'zScore'}, ' - ',
                        $callsPerLibrary{$cellType}{$geneId}{'pValue'}, ' - ',
                        "high quality", ' - ', 'not excluded', "\n";
                    } else {
                        $insAnnotatedSampleGeneResult->execute($annotatedSampleId, $bgeeGeneId,
                        'cpm', $callsPerLibrary{$cellType}{$geneId}{'cpm'}, 0,
                        $callsPerLibrary{$cellType}{$geneId}{'sumUMI'},
                        $callsPerLibrary{$cellType}{$geneId}{'zScore'},
                        $callsPerLibrary{$cellType}{$geneId}{'pValue'},
                        'high quality', 'not excluded')
                            or die $insAnnotatedSampleGeneResult->errstr;
                    }
                }
                #commit after each barcode
                $bgee_data->commit;
                $insAnnotatedSampleGeneResult->finish;
                $bgee_data->disconnect;
                $pm->finish;
            }
            $pm->wait_all_children
        }

        ## Now start to insert individual samples
        # for now we only insert barcodes mapped to a cell type. It could be possible
        # to insert other cell types in a different table
        # (e.g rnaSeqLibraryIndevidualSampleNotAnnotated(rnaSeqLibraryId, barcode, ...))
        my %barcodeToIndividualSampleId;
        my @barcodesArray;
        for my $barcode (sort keys %{$barcodesToCellType{$libraryId}{'barcodes'}}) {
            my $annotatedSampleId = $celltypeToAnnotatedSampleId{$barcodesToCellType{$libraryId}{'barcodes'}{$barcode}{'cellTypeId'}};
            my $individualSampleId = insert_get_individual_sample($insIndividualSample, $selectIndividualSampleId,
                $annotatedSampleId, $barcode, "", $debug);
            $barcodeToIndividualSampleId{$barcode} = $individualSampleId;
            push (@barcodesArray, $barcode);
        }

        # now start to insert abundance per cell. parallelized per celltype
        my $pm = new Parallel::ForkManager($numberCore);
        foreach (@barcodesArray) {
            my $pid = $pm->start and next;
            my $bgee_data = Utils::connect_bgee_db($bgee_connector);
            # disable autocommit for $bgee_data. Allows to manually commit after each library.
            $bgee_data->{AutoCommit} = 0;
            my $insIndividualSampleGeneResult = $bgee_data->prepare($insert_individualSampleGeneResult);
            
            my $barcode = $_;
            my $individualSampleId = $barcodeToIndividualSampleId{$barcode};

            for my $geneId (sort keys %{$sparseMatrixCount{$barcode}} ) {
                # check that the gene is present in the database. It is both a
                # security check and a way to remove intergenic regions
                next if (! exists $genes{$libraries{$expId}->{$libraryId}->{'speciesId'}}{$geneId} ||
                    ! exists $sparseMatrixCount{$barcode}{$geneId});
                my $bgeeGeneId = $genes{$libraries{$expId}->{$libraryId}->{'speciesId'}}{$geneId};
                if (! exists $sparseMatrixCpm{$barcode}{$geneId}) {
                    warn "Warning, gene $geneId has count for barcode $barcode but",
                        " no abundance was generated";
                    next;
                }
                if ($debug) {
                    print 'INSERT INTO rnaSeqLibraryIndividualSampleGeneResult: ',
                    $individualSampleId, ' - ', $bgeeGeneId, ' - ', "cpm", ' - ',
                    $sparseMatrixCpm{$barcode}{$geneId}, ' - ', 0, ' - ',
                    $sparseMatrixCount{$barcode}{$geneId}, ' - ',
                    "high quality", ' - ', 'not excluded', "\n";
                } else {
                    $insIndividualSampleGeneResult->execute($individualSampleId, $bgeeGeneId,
                        'cpm', $sparseMatrixCpm{$barcode}{$geneId}, 0,
                        $sparseMatrixCount{$barcode}{$geneId}, 'high quality', 'not excluded')
                            or die $insIndividualSampleGeneResult->errstr;
                }
            }
            #commit after each barcode
            $bgee_data->commit;
            $insIndividualSampleGeneResult->finish;
            $bgee_data->disconnect;
            $pm->finish;
        }
        $pm->wait_all_children
    }
}

## close prepared statements.
$insLib->finish;
$insAnnotatedSample->finish;
$selectAnnotatedSampleId->finish;
$insIndividualSample->finish;
$selectIndividualSampleId->finish;
$bgee_metadata->disconnect;
exit 0;

