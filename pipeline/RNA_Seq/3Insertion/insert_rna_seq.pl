#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use File::Slurp;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;
$| = 1; # no buffering of output

# Frederic Bastian, created November 2012
# Julien Roux, updated Oct 2016
#####################################################################


my $abundance_file = 'gene_level_abundance+calls.tsv';

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($extraMapping)   = ('');
my ($rnaSeqLibrary, $all_results, $sex_info)  = ('', '', '');
my ($rnaSeqExperiment, $library_info, $excluded_libraries, $excluded_biotypes, $library_stats, $report_info) = ('', '', '', '', '', '');
my ($debug)                      = (0);
my ($Aport, $Sport)              = (0, 0);
my %opts = ('bgee=s'                => \$bgee_connector,     # Bgee connector string
            'rnaSeqLibrary=s'       => \$rnaSeqLibrary,      # RNAseqLibrary
            'rnaSeqExperiment=s'    => \$rnaSeqExperiment,   # RNAseqExperiment
            'library_info=s'        => \$library_info,       # rna_seq_sample_info.txt file
            'excluded_libraries=s'  => \$excluded_libraries, # rna_seq_sample_excluded.txt file
            'excluded_biotypes=s'   => \$excluded_biotypes,  # biotypes_excluded_absent_calls.tsv file
            'library_stats=s'       => \$library_stats,      # presence_absence_all_samples.txt
            'report_info=s'         => \$report_info,        # reports_info_all_samples.txt
            'all_results=s'         => \$all_results,        # /var/bgee/extra/pipeline/rna_seq/all_results_bgee_v15/
            'sex_info=s'            => \$sex_info,           # generated_files/uberon/uberon_sex_info.tsv
            'extraMapping=s'        => \$extraMapping,       # Extra mapping for too up-to-date ontology terms
            'debug'                 => \$debug,
            'Aport=i'               => \$Aport,              # ID MAPPING anatomy port socket
            'Sport=i'               => \$Sport,              # ID MAPPING stage   port socket
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $rnaSeqLibrary eq '' || $rnaSeqExperiment eq '' || $library_info eq ''  || $excluded_libraries eq '' || $excluded_biotypes eq '' || $library_stats eq '' || $report_info eq '' || $all_results eq '' || $sex_info eq '' || $Aport == 0 || $Sport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g., $0  -bgee=\$(BGEECMD) -rnaSeqLibrary=RNASeqLibrary_full.tsv -rnaSeqExperiment=RNASeqExperiment_full.tsv -library_info=\$(RNASEQ_SAMPINFO_FILEPATH) -excluded_libraries=\$(RNASEQ_SAMPEXCLUDED_FILEPATH) -excluded_biotypes=\$(RNASEQ_BIOTYPEEXCLUDED_FILEPATH) -library_stats=\$(RNASEQSAMPSTATS) -report_info=\$(RNASEQREPORTINFO) -all_results=\$(RNASEQALLRES) -sex_info=\$(UBERON_SEX_INFO_FILE_PATH) -extraMapping=\$(EXTRAMAPPING_FILEPATH) -Aport=\$(IDMAPPINGPORT) -Sport=\$(STGMAPPINGPORT)    > $@.tmp 2>warnings.$@
\t-bgee                Bgee connector string
\t-rnaSeqLibrary       RNAseqLibrary annotation file
\t-rnaSeqExperiment    RNAseqExperiment file
\t-library_info        rna_seq_sample_info.txt file
\t-excluded_libraries  rna_seq_sample_excluded.txt file
\t-excluded_biotypes   file containing the mapping between protocol and biotypes not used to generate absent calls
\t-library_stats       presence_absence__all_samples.txt
\t-report_info         reports_info_all_samples.txt
\t-all_results         all_results directory
\t-sex_info            file containing sex-related info about anatomical terms
\t-extraMapping        Extra mapping file
\t-debug               insertions are not made, just printed
\t-Aport               ID MAPPING anatomy port socket
\t-Sport               ID MAPPING stage   port socket
\n";
    exit 1;
}

require("$FindBin::Bin/../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

# Library info used to launch the pipeline
my %libraries         = getAllRnaSeqLibrariesInfo($library_info);
print "\t", scalar keys %libraries, " experiments with libraries mapped.\n";

# Excluded libraries (after mapping step)
my %excludedLibraries = getExcludedLibraries($excluded_libraries);
print "\t", scalar keys %excludedLibraries, " libraries excluded.\n";

my $count_libs = 0;
my %all_species; # record all species
for my $expId ( sort keys %libraries ){
    for my $libraryId ( sort keys %{$libraries{$expId}} ){
        if ( exists($excludedLibraries{$libraryId}) ){
            delete $libraries{$expId}->{$libraryId};
        } else {
            $all_species{$libraries{$expId}->{$libraryId}->{'speciesId'}}++;
            $count_libs++;
            unless ( -s "$all_results/$libraryId/$abundance_file" ){
                die "Missing or empty processed data file $all_results/$libraryId/$abundance_file for library $libraryId! Please check that the transfer from cluster was successful. Otherwise this library should maybe be added to the file of excluded libraries?\n";
            }
        }
    }
}
print "\t", $count_libs, " libraries mapped and to be inserted.\n";

# Library info generated by the pipeline
my %librariesStats    = getAllRnaSeqLibrariesStats($library_stats);
# Library info generated by the pipeline
my %reportInfo        = getAllRnaSeqReportInfo($report_info);

# Library annotation coming from flat files
my %tsv = %{ Utils::read_spreadsheet("$rnaSeqExperiment", "\t", 'csv', '"', 1) };
my %experiments       = getAllAnnotatedExperiments2( \%tsv );

# Using getAllRnaSeqAnnotations2 doesn't really work because header and commented libraries start with "#"
# %tsv = %{ Utils::read_spreadsheet("$rnaSeqLibrary", "\t", 'csv', '"', 1) };
# my %annotations     = getAllRnaSeqAnnotations2( \%tsv );
my %annotations       = getAllRnaSeqAnnotations($rnaSeqLibrary);
my $commented = 0;
my $species_not_included = 0;
$count_libs = 0;
for my $expId ( sort keys %annotations ){
    for my $libraryId ( sort keys %{$annotations{$expId}} ){
        $count_libs++;
        if ( $annotations{$expId}->{$libraryId}->{'commented'} ){
            $commented++;
        } elsif ( !exists($all_species{ $annotations{$expId}->{$libraryId}->{'speciesId'} }) ){
            $species_not_included++;
        } elsif ( (!exists $libraries{$expId}->{$libraryId}) and (!exists $excludedLibraries{$libraryId}) ){
            print "\t", $libraryId, " library annotated but not mapped. This may be for different reasons: is it in the RNASeqLibrary_worm_exclusion.tsv file? Is it a CAGE/RACE/weird experiment? Is it in SRA format (SRX/ERX)?\n";
        }
    }
}
print "\t", scalar keys %annotations, " experiments annotated.\n";
print "\t", $count_libs, " libraries annotated.\n";
print "\t", $commented, " libraries commented.\n";
print "\t", $species_not_included, " libraries from species not included into Bgee.\n";

# Load sex-related information needed for sub 'insert_get_condition'
my $anatSexInfo    = Utils::get_anat_sex_info($sex_info);
my $speciesSexInfo = Utils::get_species_sex_info($bgee);

################
# DATA SOURCES #
################
my %bgeeDataSources = ();
my $selSrc = $bgee->prepare("SELECT dataSourceName, dataSourceId FROM dataSource WHERE category =\'RNA-Seq data source\'");
$selSrc->execute()  or die $selSrc->errstr;
while ( my @data = $selSrc->fetchrow_array ){
    $bgeeDataSources{$data[0]} = $data[1];
}
$selSrc->finish;

############
# BIOTYPES #
############

my %biotypeNameToBiotypId = ();
my $selBiotypes = $bgee->prepare("SELECT geneBioTypeId, geneBioTypeName FROM geneBioType");
$selBiotypes->execute()  or die $selBiotypes->errstr;
while ( my @data = $selBiotypes->fetchrow_array ){
    $biotypeNameToBiotypId{$data[1]} = $data[0];
}
$selBiotypes->finish;

######################
# INSERT EXPERIMENTS #
######################
#first need to retrieve already inserted experimentIds
my $selExp = $bgee->prepare("select rnaSeqExperimentId from rnaSeqExperiment");
$selExp->execute()  or die $selExp->errstr;
my %insertedExp = ();
while ( my @data = $selExp->fetchrow_array ){
    $insertedExp{$data[0]} = 1;
}
print "Inserting experiments...\n";
my $insExp = $bgee->prepare('INSERT INTO rnaSeqExperiment (rnaSeqExperimentId, rnaSeqExperimentName, rnaSeqExperimentDescription, dataSourceId) VALUES (?, ?, ?, ?)');
for my $expId ( sort keys %experiments ){
    #do not insert the experiment if it is already in the database
    next if exists $insertedExp{$expId};
    print "\t$expId\n";
    if ( $debug ){
        binmode(STDOUT, ':utf8');
        print 'INSERT INTO rnaSeqExperiment: ',
            $expId, ' - ', $experiments{$expId}->{'name'}, ' - ',
            $experiments{$expId}->{'description'}, ' - ',
            $bgeeDataSources{$experiments{$expId}->{'source'}}, "\n";
    }
    else {
        $insExp->execute($expId, $experiments{$expId}->{'name'}, $experiments{$expId}->{'description'}, $bgeeDataSources{$experiments{$expId}->{'source'}})  or die $insExp->errstr;
    }
}
$insExp->finish();
print "Done\n\n";

##########################################
# INSERT POPULATION CAPTURE AND MAPPING  #
# POPULATION CAPTURE TO BIOTYPE NOT USED #
#        TO GENERATE ABSENT CALLS        #
##########################################

#first retrieve already inserted population captured
my $selPopCapture = $bgee->prepare('SELECT rnaSeqPopulationCaptureId  from rnaSeqPopulationCapture');
$selPopCapture->execute()  or die $selPopCapture->errstr;
my %insertedPopCaptured = ();
while ( my @data = $selPopCapture->fetchrow_array ){
    $insertedPopCaptured{$data[0]} = 1;
}

#retrieve population capture from population capture file
my %popCaptureToBiotypes = retrieveProtocolsToBiotypeExcludeAbsentCalls($excluded_biotypes);

# insert the new population captured
my $insPopCapture = $bgee->prepare('INSERT INTO rnaSeqPopulationCapture (rnaSeqPopulationCaptureId) VALUES (?)');
my $insPopCaptureToBiotype = $bgee->prepare(
    'INSERT INTO rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls (rnaSeqPopulationCaptureId, geneBioTypeId) VALUES (?, ?)');
for my $populationCaptureId ( keys %popCaptureToBiotypes ){
    # skip insertion for already inserted population capture
    next if exists $insertedPopCaptured{$populationCaptureId};
    # insert population capture
    if($debug) {
        print 'INSERT INTO rnaSeqPopulationCapture: ', $populationCaptureId, "\n";  
    } else {
        $insPopCapture->execute($populationCaptureId);
    }
    # insert pop. capture to excluded absent calls
    for my $biotypeName (@{$popCaptureToBiotypes{$populationCaptureId}}) {
        my $biotypeId = $biotypeNameToBiotypId{$biotypeName};
        if($debug) {
            print 'INSERT INTO rnaSeqProtocolToBiotypeExcludedAbsentCalls: ', $populationCaptureId,  ' - ', 
                                                                              $biotypeId,   "\n";
        } else {
            $insPopCaptureToBiotype->execute($populationCaptureId, $biotypeId) or die $insPopCaptureToBiotype->errstr;
        }
    }
    # add the population capture to the list of inserted ones
    $insertedPopCaptured{$populationCaptureId} = 1;
}
$selPopCapture->finish();
$insPopCapture->finish();
$insPopCaptureToBiotype->finish();

################################################
# GET GENE INTERNAL IDS                        #
# GET ORGAN, STAGE, AND CONDITIONS INFORMATION #
################################################

my %genes;
# go over all libraries to check all species with data
for my $expId ( keys %libraries ){
    foreach my $libraryId ( keys %{$libraries{$expId}} ){
        $genes{$libraries{$expId}->{$libraryId}->{'speciesId'}} = ();
    }
}
# Get hash of geneId to bgeeGeneId mapping per species
for my $speciesId ( keys %genes ){
    $genes{$speciesId} = Utils::query_bgeeGene($bgee, $speciesId);
}

# Parse extra mapping info for currently too up-to-date annotations
##UnmappedId    UnmappedName    UberonID    UberonName    Comment
my %extra = map  { my @tmp = split(/\t/, $_, -1); if ( $tmp[2] ne '' && $tmp[0] ne '' ){ $tmp[0] => $tmp[2] } else { 'nonono' => 'nonono' } }
            grep { !/^#/ }
            read_file("$extraMapping", chomp => 1);

# Get used stages & anatEntityId from annotation sheet
%tsv = %{ Utils::read_spreadsheet("$rnaSeqLibrary", "\t", 'csv', '"', 1) };
my @Stg  = @{ $tsv{'stageId'} };
my @Anat = @{ $tsv{'anatId'} };

# Fix mapping with extra mapping file
@Stg  = map { $extra{$_} || $_ } @Stg;
@Anat = map { $extra{$_} || $_ } @Anat;

my $doneAnat = Utils::get_anatomy_mapping(\@Anat, $Aport, 0);
my $doneStg  = Utils::get_anatomy_mapping(\@Stg,  $Sport, 0);

# Get already known conditions
my $conditions = Utils::query_conditions($bgee);

# Get simpler (upper level) stage equivalences
my $stage_equivalences = Utils::get_stage_equivalences($bgee);


################################
# INSERT LIBRARIES AND RESULTS #
################################
print "Inserting libraries and all results...\n";
# query for samples insertion
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
                                'multipleLibraryIndividualSample, barcode, time, timeUnit, physiologicalStatus)'.
                                ' VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)';

my $select_annotatedSampleId =  'SELECT rnaSeqLibraryAnnotatedSampleId FROM '.
                                'rnaSeqLibraryAnnotatedSample WHERE rnaSeqLibraryId = ? AND '.
                                'conditionId = ?';

my $insert_annotatedSampleGeneResult =  'INSERT INTO rnaSeqLibraryAnnotatedSampleGeneResult ('.
                                        'rnaSeqLibraryAnnotatedSampleId, bgeeGeneId, abundanceUnit, abundance,'.
                                        'readsCount, UMIsCount, zScore, pValue,'.
                                        'reasonForExclusion) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)';

# Excluded libraries
my $selectDiscardedLib = $bgee->prepare('SELECT rnaSeqLibraryId FROM rnaSeqLibraryDiscarded');
$selectDiscardedLib->execute()  or die $selectDiscardedLib->errstr;
my %insertedDiscardedLib = ();
while ( my @data = $selectDiscardedLib->fetchrow_array ){
    $insertedDiscardedLib{$data[0]} = 1;
}
my $insExcludedLib = $bgee->prepare('INSERT INTO rnaSeqLibraryDiscarded (rnaSeqLibraryId, rnaSeqLibraryDiscardReason) VALUES (?, ?)');
for my $libraryId ( sort keys %excludedLibraries ){
    next if exists $insertedDiscardedLib{$libraryId};
    if ( $debug ){
        print 'INSERT INTO rnaSeqLibraryDiscarded: ', $libraryId, "\n";
    }
    else {
        $insExcludedLib->execute($libraryId, $excludedLibraries{$libraryId})  or die $insExcludedLib->errstr;
    }
}

# query for runs insertion
my $insRun = $bgee->prepare('INSERT INTO rnaSeqRun (rnaSeqRunId, rnaSeqLibraryId) VALUES (?, ?)');

my $inserted = 0;

# used to commit after each library when condition and libraries were not inserted
print "disable autocommit. Manually commit for each library\n";
$bgee->{AutoCommit} = 0;

# prepare queries
my $insLib = $bgee->prepare($insert_libraries);
my $insAnnotatedSample = $bgee->prepare($insert_annotatedSamples);
my $selectAnnotatedSampleId = $bgee->prepare($select_annotatedSampleId);
my $insAnnotatedSampleGeneResult = $bgee->prepare($insert_annotatedSampleGeneResult);

for my $expId ( sort keys %libraries ){
    LIBRARY:
    for my $libraryId ( sort keys %{$libraries{$expId}} ){
        if ( !exists $reportInfo{$libraryId} ){
            warn "Report file does not contain this library [$libraryId]\n";
            next LIBRARY;
        }
        print "\t$expId $libraryId\n";

        # Remap to extra mapping if any
        $annotations{$expId}->{$libraryId}->{'anatId'} = $extra{ $annotations{$expId}->{$libraryId}->{'anatId'} } || $annotations{$expId}->{$libraryId}->{'anatId'};
        $annotations{$expId}->{$libraryId}->{'stageId'}  = $extra{ $annotations{$expId}->{$libraryId}->{'stageId'} }  || $annotations{$expId}->{$libraryId}->{'stageId'};

        if ( !exists $doneAnat->{$annotations{$expId}->{$libraryId}->{'anatId'}} || $doneAnat->{$annotations{$expId}->{$libraryId}->{'anatId'}} eq '' ){
            warn "[$annotations{$expId}->{$libraryId}->{'anatId'}] unmapped organ id for [$libraryId]\n";
            next LIBRARY;
        }
        # if ( !exists $doneStg->{$annotations{$expId}->{$libraryId}->{'stageId'}}   || $doneStg->{$annotations{$expId}->{$libraryId}->{'stageId'}}   eq '' ){
        #     warn "[$annotations{$expId}->{$libraryId}->{'stageId'}] unmapped stage id for [$libraryId]\n";
        #     next LIBRARY;
        # }

        # Get conditionId/exprMappedConditionId for this library
        # Updates also the hash of existing conditions
        my $condKeyMap;
        ($condKeyMap, $conditions) = Utils::insert_get_condition($bgee,
                                                                 $conditions,
                                                                 $stage_equivalences,
                                                                 $doneAnat->{$annotations{$expId}->{$libraryId}->{'anatId'}},
                                                                 $annotations{$expId}->{$libraryId}->{'stageId'},
                                                                 $annotations{$expId}->{$libraryId}->{'speciesId'},
                                                                 $annotations{$expId}->{$libraryId}->{'sex'},
                                                                 $annotations{$expId}->{$libraryId}->{'strain'},
                                                                 $anatSexInfo, $speciesSexInfo,
                                                                 $libraryId, '',
                                                                );
        # We consider the fine-grained (low-level) conditionId for insertion: $condKeyMap->{'conditionId'}

        # insert library
        if ( $debug ){
            print 'INSERT INTO rnaSeqLibrary: ', $libraryId,                              ' - ',
                  $expId, ' - ', $libraries{$expId}->{$libraryId}->{'platform'},          ' - ',
                  #technologyName always Illumina for now
                  #TODO remove that hardcoded insertion when refactoring insertion scripts
                  #Actually why this is information is required????
                  "Illumina",                         ' - ',
                  # isSingleCell (always true for target base RNASeq)
                  '0',                                                                    ' - ',
                  # sampleMultiplexing (always true for target base)
                  '0',                                                                    ' - ',
                  # libraryMultiplexing (for now always false for target base)
                  '0',                                                                    ' - ',
                  # strandSelection. No information in annotation file.
                  # Could maybe be detected from the platform or the technology
                  'NA',                                                                   ' - ',
                  # cellCompartment. No information in the annotation file
                  'cell',                  ' - ',
                  # sequencedTranscriptPart is always full length for bulk RNA-Seq
                  'full length',          ' - ',
                  # fragmentation.
                  $libraries{$expId}->{$libraryId}->{'readLength'},                       ' - ',
                  # rnaSeqPopulationCaptureId. For now all target base are polyA
                  $annotations{$expId}->{$libraryId}->{'populationCapture'},                                                                ' - ',
                  $annotations{$expId}->{$libraryId}->{'genotype'},                       ' - ',
                  $reportInfo{$libraryId}->{'allReadsCount'},                         ' - ',
                  $reportInfo{$libraryId}->{'mappedReadsCount'},                      ' - ',
                  $reportInfo{$libraryId}->{'minReadLength'},                         ' - ',
                  $reportInfo{$libraryId}->{'maxReadLength'},                         ' - ',
                  lc $libraries{$expId}->{$libraryId}->{'libraryType'},"\n";

        } else {
            $insLib->execute($libraryId,
                            $expId,
                            $libraries{$expId}->{$libraryId}->{'platform'},
                            #TODO remove that hardcoded insertion when refactoring insertion scripts
                            #Actually why this information is required????
                            'Illumina',
                            0,
                            0,
                            0,
                            # strandSelection. No information in annotation file.
                            # Could maybe be detected from the platform or the technology
                            'NA',
                            # cellCompartment. No information in the annotation file
                            'cell',
                            # sequencedTranscriptPart is always full length for bulk RNA-Seq
                            'full length',
                            $libraries{$expId}->{$libraryId}->{'readLength'},
                            $annotations{$expId}->{$libraryId}->{'populationCapture'},
                            $annotations{$expId}->{$libraryId}->{'genotype'},
                            $reportInfo{$libraryId}->{'allReadsCount'},
                            $reportInfo{$libraryId}->{'mappedReadsCount'},
                            $reportInfo{$libraryId}->{'minReadLength'},
                            $reportInfo{$libraryId}->{'maxReadLength'},
                            lc $libraries{$expId}->{$libraryId}->{'libraryType'}
                        )  or die $insLib->errstr;
        }

        # Now insert rnaSeqLibraryAnnotatedSample
        my $annotatedSampleId = insert_get_annotated_sample($libraryId, $condKeyMap->{'conditionId'},
                    #no celltype author annotation for bulk RNA-Seq
                    '',
                    $annotations{$expId}->{$libraryId}->{'freeTextOrgan'},
                    $annotations{$expId}->{$libraryId}->{'freeTextStage'},
                    #for now abundace unit is always tpm for bulk.
                    'tpm', 
                    $librariesStats{$libraryId}->{'meanIntergenic'},
                    $librariesStats{$libraryId}->{'sdIntergenic'},
                    $librariesStats{$libraryId}->{'cutoffTPM'},
                    $librariesStats{$libraryId}->{'allGenesPercentPresent'},
                    $librariesStats{$libraryId}->{'proteinCodingPercentPresent'},
                    $librariesStats{$libraryId}->{'intergenicRegionsPercentPresent'},
                    $librariesStats{$libraryId}->{'pValueThreshold'},
                    # no UMIs info for bulk. So allUMIsCount and mappedUMIsCount are 0
                    0, 0,
                    # multipleLibraryIndividualSample and barcode are false and empty for bulk
                    0, '',
                    # As of Bgee 15.2 time and timeUnit were not present in the annotation
                    # files. Those fields have been added for the SalmoBase project.
                    #TODO: implement the insertion of those metadata once they are available.
                    undef , undef , $annotations{$expId}->{$libraryId}->{'physiologicalStatus'},
                    $insAnnotatedSample, $selectAnnotatedSampleId, $debug);

        # insert runs
        for my $runId ( keys %{$libraries{$expId}->{$libraryId}->{'runIds'}} ){
            if ( $debug ){
                print 'INSERT INTO rnaSeqRun: ', $runId, ' - ', $libraryId, "\n";
            }
            else {
                $insRun->execute($runId, $libraryId)  or die $insRun->errstr;
            }
        }

        # insert genes results
        my %genesResults = getGenesResults("$all_results/$libraryId/$abundance_file");
        for my $geneId ( keys %genesResults ){
            $inserted++;
            # Note: pre-filtering exclusion is now managed in the script insert_rna_seq_expression.pl,
            # it used to be managed here.
            my $exclusion = $Utils::CALL_NOT_EXCLUDED;

            if ( $debug ){
                print 'INSERT INTO rnaSeqLibraryAnnotatedSampleGeneResult: ', $annotatedSampleId,   ' - ',
                      $genes{ $libraries{$expId}->{$libraryId}->{'speciesId'}}->{ $geneId }, ' - ',
                      'tpm',                                      ' - ',
                      $genesResults{$geneId}->{'TPM'},            ' - ',
                      $genesResults{$geneId}->{'estimatedCount'}, ' - ',
                      # no UMIs for bulk RNA-Seq
                      0,                                          ' - ',
                      $genesResults{$geneId}->{'zscore'},         ' - ',
                      $genesResults{$geneId}->{'pValue'},         ' - ',
                      $exclusion, "\n";
            }
            else {
                # pvalue and zscore can be null (if no read mapped). In this case BgeeCall retrieve "NA".
                # DBI use undef value to insert null in the database. That's why we modify "NA" to undef for zScore.
                # For the pvalue we decided replace NA with 1 in order to use this value as a datapoint to generate propagated calls 
                if ($genesResults{$geneId}->{'pValue'} eq "NA") {
                  $genesResults{$geneId}->{'pValue'} = 1;
                }
                if ($genesResults{$geneId}->{'zscore'} eq "NA") {
                  $genesResults{$geneId}->{'zscore'} = undef;
                }
                $insAnnotatedSampleGeneResult->execute($annotatedSampleId,
                                    # geneId is an ensembl ID, we need to get the bgeeGeneId
                                    $genes{ $libraries{$expId}->{$libraryId}->{'speciesId'}}->{ $geneId },
                                    'tpm',
                                    $genesResults{$geneId}->{'TPM'},
                                    $genesResults{$geneId}->{'estimatedCount'},
                                    # no UMIs for bulk RNA-Seq
                                    0,
                                    $genesResults{$geneId}->{'pValue'},
                                    $genesResults{$geneId}->{'pValue'},
                                    $exclusion,
                                   )  or die $insAnnotatedSampleGeneResult->errstr;
            }

        }
        # used to commit after each library when condition and libraries were not inserted
        $bgee->commit;
    }
}
# used to commit after each library when condition and libraries were not inserted
print "reactivate autocommit\n";
$bgee->{AutoCommit} = 1;

$insLib->finish();
$insRun->finish();
$insAnnotatedSample->finish();
$selectAnnotatedSampleId->finish();
$insAnnotatedSampleGeneResult->finish();

print "Done. You should have $inserted rows in the rnaSeqResult table.\nExiting\n";

exit 0;

