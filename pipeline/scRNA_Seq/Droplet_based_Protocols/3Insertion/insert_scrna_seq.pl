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
use Parallel::Loops;

$| = 1; # no buffering of output

# Julien Wollbrett, created December 2022
# This file insert target base experiments, libraries, annotatedSamples, individualSamples, 
# and individualSamples results
#  * the individualSamples result are at the cell (=barcode) level. They correspond to the raw count as they
#    are provided by kallisto bus.
#  * the annotatedSamples result are at the celltype level. They correspond to each annotation of a library.
#    To obtain this result we sum up all rawCounts of cell corresponding to the same celltype. It is this
#    information that is used to generate the propagated calls.
#XXX the script should probably first insert cond then annoated samples, then individual samples for all target base data.
#    A 2nd iteration of this script will then insert all the annotated/individual results. The advantage of this approach
#    is that it is possible to avoid conflicts of autoincrement primary keys (conditionId, annotatedSampleId,
#    individualSampleId) without potential colision by running the script linearily. Then, it could be possible to insert
#    each experiment in parallel.

#####################################################################

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($targetBaseLibrary, $allResults, $sexInfo)  = ('', '', '');
my ($singleCellExperiment, $libraryInfo, $sourceDir) = ('', '', '');
#my ($library_stats, $report_info = ('', '');
my ($debug)                      = (0);
my %opts = ('bgee=s'                => \$bgee_connector,       # Bgee connector string
            'targetBaseLibrary=s'   => \$targetBaseLibrary,    # target base RNAseq library annotations
            'singleCellExperiment=s'=> \$singleCellExperiment, # single cell RNASeq experiment annotations
            'libraryInfo=s'         => \$libraryInfo,          # metadata_info_info_10X.txt file
            'allResults=s'          => \$allResults,           # path to dir containing results for all libraries
            'sourceDir=s'           => \$sourceDir,            # path to the directory containing source files of the target base pipeline
            'sexInfo=s'             => \$sexInfo,              # generated_files/uberon/uberon_sex_info.tsv
            'debug'                 => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $targetBaseLibrary eq '' || $singleCellExperiment eq '' ||
    $libraryInfo eq '' || $allResults eq '' || $sourceDir eq '' ||$sexInfo eq ''){
    print "\n\tInvalid or missing argument:
\te.g., $0  -bgee=\$(BGEECMD) -targetBaseLibrary=RNASeqLibrary_full.tsv -singleCellExperiment=RNASeqExperiment_full.tsv -libraryInfo=metadata_info_10X.txt -sexInfo=\$(UBERON_SEX_INFO_FILE_PATH) > $@.tmp 2>warnings.$@
\t-bgee                    Bgee connector string
\t-targetBaseLibrary       targetBaseLibrary annotation file
\t-singleCellExperiment    singleCellExperiment file
\t-bgeelibraryInfo             metadata_info_10X.txt file
\t-allResults              allResults directory
\t-sourceDir               path to the directory containing sources files of the target bas pipeline
\t-sexInfo                 file containing sex-related info about anatomical terms
\t-debug                   insertions are not made, just printed
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../target_base_utils.pl");

# initialize variables

my $barcodeToCelltypeFilePattern = "$sourceDir/scRNASeq_barcode_EXP_ID.tsv";

#TODO create a script argument for this variable
my $maxProcs = 15;
my $pl = Parallel::Loops->new($maxProcs);

####################### FUNCTIONS ###############################
# for now all the functions are in this script but some of them could have to move
# to a generic scrna_seq_utils.pl script if more technologies are added
#################################################################


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
                        'rnaSeqPopulationCaptureId, libraryType)'.
                        ' VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)';

# Do not anymore insert information related to calls creation in this script (done in an other script)
# my $insert_annotatedSamples =   'INSERT INTO rnaSeqLibraryAnnotatedSampleId (rnaSeqLibraryAnnotatedSampleId,'.
#                               'rnaSeqLibraryId, conditionId, abundanceUnit,'.
#                                'meanAbundanceReferenceIntergenicDistribution,'.
#                                'sdAbundanceReferenceIntergenicDistribution, tmmFactor, abundanceThreshold,'.
#                                'allGenesPercentPresent, proteinCodingGenesPercentPresent,'.
#                                'intergenicRegionsPercentPresent, pValueThreshold, allReadsCount, allUMIsCount,'.
#                                'mappedReadsCount, mappedUMIsCount, minReadLength, maxReadLength,'.
#                                'libraryMaxRank, libraryDistinctRankCount, multipleLibraryIndividualSample,'.
#                                'barcode, genotype) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,'.
#                                '?, ?, ?, ?, ?, ?)';

my $insert_annotatedSamples =   'INSERT INTO rnaSeqLibraryAnnotatedSampleDev (rnaSeqLibraryId,'.
                                'conditionId) VALUES (?, ?)';

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
                                        'readsCount, UMIsCount, zScore, pValue, detectionFlag, rnaSeqData,'.
                                        'reasonForExclusion)'.
                                        'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)';

my $insert_individualSampleGeneResult =  'INSERT INTO rnaSeqLibraryIndividualSampleGeneResultDev ('.
                                        'rnaSeqLibraryIndividualSampleId, bgeeGeneId, abundanceUnit, abundance,'.
                                        'readsCount, UMIsCount, rnaSeqData, reasonForExclusion)'.
                                        'VALUES (?, ?, ?, ?, ?, ?, ?, ?)';
########################################################################
########################## Main Script #################################
########################################################################

## Initialize Bgee db connections
# connection used to insert data in rnaSeqExperiment, rnaSeqLibrary, rnaSeqAnnotatedSample
# and rnaSeqIndividualSample tables. Keep AutoCommit as each insert has to be done directly
# to retrieve autoincrement IDs of conditionId, rnaSeqLibraryAnnotatedSampleId and
# rnaSeqLibraryIndividualSampleId
my $bgee_metadata = Utils::connect_bgee_db($bgee_connector);
$bgee_metadata->{AutoCommit} = 1;

## Load 
# Library info from manual curation used to launch the pipeline
my %libraries         = getTargetBaseCuratedLibrariesAnnotation($targetBaseLibrary);
print "\t", scalar keys %libraries, " experiments with libraries mapped.\n";
# Info of processed libraries coming from the pipeline
my %processedLibraries = get_processed_libraries_info($libraryInfo);
# Experiment annotation coming from flat files
my @experimentType = ('3\'end', 'Full-length and 3\'end');
my %experiments       = getSingleCellExperiments($singleCellExperiment,
    @experimentType);

# Load sex-related information needed for sub 'insert_get_condition'
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
for my $expId ( sort keys %experiments ){
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
        my %sparseMatrixCount = read_sparse_matrix("$allResults/$libraryId/gene_counts", "gene");
        my %sparseMatrixCpm = read_sparse_matrix("$allResults/$libraryId/cpm_counts", "cpm_counts");
        #TODO: in R remove intergenic regions from count matrix and then calculate cpm using
        #NormalizeCount function from Seurat R package
        #my %sparseMatrixCpm   = read_sparse_matrix("$allResults/$libraryId/gene_cpm", "gene");

        print "\tInsert $libraryId from $expId\n";
        # For now all target base are polyA. Check that population capture polyA is already present in the database
        if ( !grep(/^polyA$/, @populationCapture ) ) {
            die "polyA population capture not found in the database";
        }

        # insert libraries
        if ( $debug ){
            print 'INSERT INTO rnaSeqLibrary: ', $libraryId,                                       ' - ',
                  $expId, ' - ', $libraries{$expId}->{$libraryId}->{'platform'},                 ' - ',
                  #technologyName
                  $libraries{$expId}->{$libraryId}->{'protocol'},                                ' - ',
                  # isSingleCell (always true for target base RNASeq)
                  '1',                                                                             ' - ',
                  # sampleMultiplexing (always true for target base)
                  '1',                                                                             ' - ',
                  # libraryMultiplexing (for now always false for target base)
                  '0',                                                                             ' - ',
                  # strandSelection. No information in annotation file.
                  # Could maybe be detected from the platform or the technology
                  'NA',                                                                            ' - ',
                  # cellCompartment. No information in the annotation file.
                  # Could maybe be detected from the technology name
                  $libraries{$expId}->{$libraryId}->{'cellCompartment'},                          ' - ',
                  $libraries{$expId}->{$libraryId}->{'sequencedTranscriptPart'},                  ' - ',
                  # fragmentation.
                  # TODO. Left empty for now but should be updated to take read length info
                  '0',                                                                            ' - ',
                  # rnaSeqPopulationCaptureId. For now all target base are polyA
                  'polyA',                                                                        ' - ',
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
                $selectAnnotatedSampleId,$condKeyMap->{'conditionId'}, $libraryId, $debug);
            # could create a hash celltype -> annotatedSampleId
            $celltypeToAnnotatedSampleId{$cellTypeId} = $annotatedSampleId;
        }

        ## Now start to insert individual samples
        # for now we only insert barcodes mapped to a cell type. It coult be possible
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

        $pl->foreach( \@barcodesArray, sub {
            my $bgee_data = Utils::connect_bgee_db($bgee_connector);
            # disable autocommit for $bgee_data. Allows to manually commit after each library.
            $bgee_data->{AutoCommit} = 0;
            my $insIndividualSampleGeneResult = $bgee_data->prepare($insert_individualSampleGeneResult);
            
            my $barcode = $_;
            my $individualSampleId = $barcodeToIndividualSampleId{$barcode};

            # before insertion in rnaSeqLibraryIndividualSampleGeneResult we need to
            # process bustools output that contain geneId, barcode and counts to:
            # - remove intergenic regions
            # - calculate cpm for each gene (should probably be done before in R using sparse matrix)
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
            #            insert_individual_sample_gene_result($insIndividualSampleGeneResult, \%sparseMatrixCount,
#                \%sparseMatrixCpm, \%{$genes{$libraries{$expId}->{$libraryId}->{'speciesId'}}},
#                $individualSampleId, $barcode, $debug);
        });
        

    }
}

## close prepared statements.
$insLib->finish;
$insAnnotatedSample->finish;
$selectAnnotatedSampleId->finish;
$insIndividualSample->finish;
$selectIndividualSampleId->finish;

exit 0;

