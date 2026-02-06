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
use Data::Dumper;

$| = 1; # no buffering of output

# TODO: 90% of this script is similar to insert_rna_seq.pl for bulk RNA-Seq.
#       We should refactor the code to avoid redundancy.
#       Differences: - populationCatpure not inserted
#                    - libraryType always 'full-length'
#                    - do not insert runIds


#####################################################################


my $abundance_file = 'gene_level_abundance+calls.tsv';

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($extraMapping)   = ('');
my ($rnaSeqLibrary, $all_results, $sex_info)  = ('', '', '');
my ($rnaSeqExperiment, $library_info, $excluded_libraries, $library_stats, $report_info) = ('', '', '', '', '');
my ($debug)                      = (0);
my ($Aport, $Sport)              = (0, 0);
my $threads = 2; # default number of parallel threads
my %opts = ('bgee=s'                => \$bgee_connector,     # Bgee connector string
            'rnaSeqLibrary=s'       => \$rnaSeqLibrary,      # scRNAseqLibrary_merged.tsv
            'rnaSeqExperiment=s'    => \$rnaSeqExperiment,   # scRNAseqExperiment.tsv
            'library_info=s'        => \$library_info,       # scrna_seq_sample_info.txt file
            'excluded_libraries=s'  => \$excluded_libraries, # scrna_seq_sample_excluded.txt file
            'library_stats=s'       => \$library_stats,      # presence_absence_all_samples.txt
            'report_info=s'         => \$report_info,        # reports_info_all_samples.txt
            'all_results=s'         => \$all_results,        # /var/bgee/extra/pipeline/rna_seq/all_results_bgee_v15/
            'sex_info=s'            => \$sex_info,           # generated_files/uberon/uberon_sex_info.tsv
            'extraMapping=s'        => \$extraMapping,       # Extra mapping for too up-to-date ontology terms
            'debug'                 => \$debug,
            'Aport=i'               => \$Aport,              # ID MAPPING anatomy port socket
            'Sport=i'               => \$Sport,              # ID MAPPING stage   port socket
            'threads=i'             => \$threads,            # number of parallel threads
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $rnaSeqLibrary eq '' || $rnaSeqExperiment eq '' || $library_info eq ''  || $excluded_libraries eq '' || $library_stats eq '' || $report_info eq '' || $all_results eq '' || $sex_info eq '' || $Aport == 0 || $Sport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g., $0  -bgee=\$(BGEECMD) -rnaSeqLibrary=RNASeqLibrary_full.tsv -rnaSeqExperiment=RNASeqExperiment_full.tsv -library_info=\$(RNASEQ_SAMPINFO_FILEPATH) -excluded_libraries=\$(RNASEQ_SAMPEXCLUDED_FILEPATH) -library_stats=\$(RNASEQSAMPSTATS) -report_info=\$(RNASEQREPORTINFO) -all_results=\$(RNASEQALLRES) -sex_info=\$(UBERON_SEX_INFO_FILE_PATH) -extraMapping=\$(EXTRAMAPPING_FILEPATH) -Aport=\$(IDMAPPINGPORT) -Sport=\$(STGMAPPINGPORT)    > $@.tmp 2>warnings.$@
\t-bgee                Bgee connector string
\t-rnaSeqLibrary       RNAseqLibrary annotation file
\t-rnaSeqExperiment    RNAseqExperiment file
\t-library_info        rna_seq_sample_info.txt file
\t-excluded_libraries  rna_seq_sample_excluded.txt file
\t-library_stats       presence_absence__all_samples.txt
\t-report_info         reports_info_all_samples.txt
\t-all_results         all_results directory
\t-sex_info            file containing sex-related info about anatomical terms
\t-extraMapping        Extra mapping file
\t-debug               insertions are not made, just printed
\t-Aport               ID MAPPING anatomy port socket
\t-Sport               ID MAPPING stage   port socket
\t-threads             number of parallel threads (default: 2)
\n";
    exit 1;
}

require("$FindBin::Bin/../../rna_seq_utils.pl");
require("$FindBin::Bin/../../target_base_utils.pl");

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

# Library info used to launch the pipeline
my %libraries = getAllFullLengthScRnaSeqLibrariesInfo($library_info);
print "\t", scalar keys %libraries, " experiments with libraries mapped.\n";

# # Excluded libraries (after mapping step)
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
my @experimentType = ('Full-length and 3\'end', 'Full-length', 'Full-length, 3\'end');
my %experiments       = getSingleCellExperiments($rnaSeqExperiment,
    @experimentType);

my %annotations       = loadFullLengthLibrariesAnnotation($rnaSeqLibrary);

# Load sex-related information needed for sub 'insert_get_condition'
my $anatSexInfo    = Utils::get_anat_sex_info($sex_info);
my $speciesSexInfo = Utils::get_species_sex_info($bgee);

################
# DATA SOURCES #
################
my %bgeeDataSources = ();
my $selSrc = $bgee->prepare("SELECT dataSourceName, dataSourceId FROM dataSource WHERE category =\'Single-cell RNA-Seq data source\'");
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

# Get used stages & anatEntityId from libraries already loaded
my (@Stg, @Anat);
for my $expId ( keys %libraries ){
    for my $libraryId ( keys %{$libraries{$expId}} ){
        push @Stg, $libraries{$expId}->{$libraryId}->{'stageId'} if $libraries{$expId}->{$libraryId}->{'stageId'};
        push @Anat, $libraries{$expId}->{$libraryId}->{'anatId'} if $libraries{$expId}->{$libraryId}->{'anatId'};
    }
}

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
                        'rnaSeqSequencerName, rnaSeqTechnologyName, rnaSeqTechnologyIsSingleCell,'.
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
                                'conditionId = ? AND cellTypeAuthorAnnotation = ?';

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

# used to commit after each library when condition and libraries were not inserted
print "disable autocommit. Manually commit for each library\n";
$bgee->{AutoCommit} = 0;

# prepare queries
my $insLib = $bgee->prepare($insert_libraries);
my $insAnnotatedSample = $bgee->prepare($insert_annotatedSamples);
my $selectAnnotatedSampleId = $bgee->prepare($select_annotatedSampleId);

# Store annotated sample IDs for later parallel insertion
my %annotatedSampleIds;

for my $expId ( sort keys %libraries ){
    LIBRARY:
    for my $libraryId ( sort keys %{$libraries{$expId}} ){
        if ( !exists $reportInfo{$libraryId} ){
            warn "Report file does not contain this library [$libraryId]\n";
            next LIBRARY;
        }
        print "\t$expId $libraryId\n";

        # Remap to extra mapping if any
        $libraries{$expId}->{$libraryId}->{'anatId'} = $extra{ $libraries{$expId}->{$libraryId}->{'anatId'} } || $libraries{$expId}->{$libraryId}->{'anatId'};
        $libraries{$expId}->{$libraryId}->{'stageId'}  = $extra{ $libraries{$expId}->{$libraryId}->{'stageId'} }  || $libraries{$expId}->{$libraryId}->{'stageId'};
        $libraries{$expId}->{$libraryId}->{'cellTypeId'} = $extra{ $libraries{$expId}->{$libraryId}->{'cellTypeId'} } || $libraries{$expId}->{$libraryId}->{'cellTypeId'};

        if ( !exists $doneAnat->{$libraries{$expId}->{$libraryId}->{'anatId'}} || $doneAnat->{$libraries{$expId}->{$libraryId}->{'anatId'}} eq '' ){
            warn "[$libraries{$expId}->{$libraryId}->{'anatId'}] unmapped organ id for [$libraryId]\n";
            next LIBRARY;
        }
        if ( !exists $doneStg->{$libraries{$expId}->{$libraryId}->{'stageId'}}   || $doneStg->{$libraries{$expId}->{$libraryId}->{'stageId'}}   eq '' ){
             warn "[$libraries{$expId}->{$libraryId}->{'stageId'}] unmapped stage id for [$libraryId]\n";
             next LIBRARY;
        }

        # Get conditionId/exprMappedConditionId for this library
        # Updates also the hash of existing conditions
        my $condKeyMap;
        ($condKeyMap, $conditions) = Utils::insert_get_condition($bgee,
                                                                 $conditions,
                                                                 $stage_equivalences,
                                                                 $doneAnat->{$libraries{$expId}->{$libraryId}->{'anatId'}},
                                                                 $libraries{$expId}->{$libraryId}->{'stageId'},
                                                                 $libraries{$expId}->{$libraryId}->{'speciesId'},
                                                                 $libraries{$expId}->{$libraryId}->{'sex'},
                                                                 $libraries{$expId}->{$libraryId}->{'strain'},
                                                                 $anatSexInfo, $speciesSexInfo,
                                                                 $libraryId, '',
                                                                 $libraries{$expId}->{$libraryId}->{'cellTypeId'}
                                                                );
        # We consider the fine-grained (low-level) conditionId for insertion: $condKeyMap->{'conditionId'}

        # insert library
        if ( $debug ){
            print 'INSERT INTO rnaSeqLibrary: ', $libraryId,                              ' - ',
                  $expId, ' - ', $libraries{$expId}->{$libraryId}->{'platform'},          ' - ',
                  $libraries{$expId}->{$libraryId}->{'protocol'},                         ' - ',
                  # isSingleCell (always true for full length single-cell RNASeq)
                  '1',                                                                    ' - ',
                  # sampleMultiplexing (always false for full length single-cell RNASeq)
                  '0',                                                                    ' - ',
                  # libraryMultiplexing (for now always false for full length single-cell RNASeq)
                  '0',                                                                    ' - ',
                  # strandSelection. No information in annotation file.
                  # Could maybe be detected from the platform or the technology
                  'NA',                                                                   ' - ',
                  # cellCompartment
                  $annotations{$expId}->{$libraryId}->{'cellCompartment'},                  ' - ',
                  # sequencedTranscriptPart is always full length for full-length scRNA-Seq
                  'full length',          ' - ',
                  # fragmentation.
                  $libraries{$expId}->{$libraryId}->{'readLength'},                       ' - ',
                  # rnaSeqPopulationCaptureId. For now all full-length scRNA-Seq annotated are polyA
                  'polyA',                                                                ' - ',
                  $libraries{$expId}->{$libraryId}->{'genotype'},                       ' - ',
                  $reportInfo{$libraryId}->{'allReadsCount'},                         ' - ',
                  $reportInfo{$libraryId}->{'mappedReadsCount'},                      ' - ',
                  $reportInfo{$libraryId}->{'minReadLength'},                         ' - ',
                  $reportInfo{$libraryId}->{'maxReadLength'},                         ' - ',
                  lc $libraries{$expId}->{$libraryId}->{'libraryType'},"\n";

        } else {
            $insLib->execute($libraryId,
                            $expId,
                            $libraries{$expId}->{$libraryId}->{'platform'},
                            $libraries{$expId}->{$libraryId}->{'protocol'},
                            1,
                            0,
                            0,
                            # strandSelection. No information in annotation file.
                            # Could maybe be detected from the platform or the technology
                            'NA',
                            # cellCompartment.
                            $annotations{$expId}->{$libraryId}->{'cellCompartment'},
                            # sequencedTranscriptPart is always full length for full length scRNA-Seq
                            'full length',
                            $libraries{$expId}->{$libraryId}->{'readLength'},
                            'polyA',
                            $libraries{$expId}->{$libraryId}->{'genotype'},
                            $reportInfo{$libraryId}->{'allReadsCount'},
                            $reportInfo{$libraryId}->{'mappedReadsCount'},
                            $reportInfo{$libraryId}->{'minReadLength'},
                            $reportInfo{$libraryId}->{'maxReadLength'},
                            lc $libraries{$expId}->{$libraryId}->{'libraryType'}
                        )  or die $insLib->errstr;
        }

        # Now insert rnaSeqLibraryAnnotatedSample
        my $annotatedSampleId = insert_get_annotated_sample($libraryId, $condKeyMap->{'conditionId'},
                    $annotations{$expId}->{$libraryId}->{'authorCellTypeAnnotation'},
                    $annotations{$expId}->{$libraryId}->{'authorAnatEntityAnnotation'},
                    $annotations{$expId}->{$libraryId}->{'authorStageAnnotation'},
                    #for now abundace unit is always tpm for full length scRNA-Seq.
                    'tpm', 
                    $librariesStats{$libraryId}->{'meanIntergenic'},
                    $librariesStats{$libraryId}->{'sdIntergenic'},
                    $librariesStats{$libraryId}->{'cutoffTPM'},
                    $librariesStats{$libraryId}->{'allGenesPercentPresent'},
                    $librariesStats{$libraryId}->{'proteinCodingPercentPresent'},
                    $librariesStats{$libraryId}->{'intergenicRegionsPercentPresent'},
                    $librariesStats{$libraryId}->{'pValueThreshold'},
                    # no UMIs info for full length scRNA-Seq. So allUMIsCount and mappedUMIsCount are 0
                    0, 0,
                    # multipleLibraryIndividualSample and barcode are false and empty for full length scRNA-Seq
                    0, '',
                    # As of Bgee 16.0, time, timeUnit and physiologicalStatus are not available for full length scRNA-Seq samples.
                    #TODO: implement the insertion of those metadata once they are available.
                    undef , undef , undef,
                    $insAnnotatedSample, $selectAnnotatedSampleId, $debug);

        # Store for parallel processing
        $annotatedSampleIds{$libraryId} = {
            annotatedSampleId => $annotatedSampleId,
            speciesId => $libraries{$expId}->{$libraryId}->{'speciesId'}
        };

        $bgee->commit;

    }
}
$insLib->finish();
$insAnnotatedSample->finish();
$selectAnnotatedSampleId->finish();

print "Reactivate autocommit\n";
$bgee->{AutoCommit} = 1;
$bgee->disconnect();

print "Done inserting libraries and samples.\n\n";

#######################################
# INSERT GENE RESULTS IN PARALLEL     #
#######################################
print "Inserting gene results in parallel...\n";
my $insertedGeneResults = 0;

# Determine number of processes
my $pm = Parallel::ForkManager->new($threads);

# Callback to retrieve inserted gene results count from each thread
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
    if (defined $data_ref) {
        $insertedGeneResults += $data_ref->{inserted};
    }
});

for my $libraryId ( sort keys %annotatedSampleIds ){
    $pm->start and next;

    # create new DB connection for each thread
    my $child_bgee = Utils::connect_bgee_db($bgee_connector);
    $child_bgee->{AutoCommit} = 0;

    my $child_insAnnotatedSampleGeneResult = $child_bgee->prepare($insert_annotatedSampleGeneResult);

    my $annotatedSampleId = $annotatedSampleIds{$libraryId}->{'annotatedSampleId'};
    my $speciesId = $annotatedSampleIds{$libraryId}->{'speciesId'};

    print "\tProcessing gene results for library: $libraryId\n";
    
    # insert genes results
    my %genesResults = getGenesResults("$all_results/$libraryId/$abundance_file");
    my $inserted = 0;
    for my $geneId ( keys %genesResults ){
        # Check if gene exists in database
        if ( !exists $genes{$speciesId}->{$geneId} ){
            die "Gene [$geneId] not found in Bgee database for species [$speciesId] in library [$libraryId].\n";
        }

        # Note: pre-filtering exclusion is now managed in the script insert_rna_seq_expression.pl,
        # it used to be managed here.
        my $exclusion = $Utils::CALL_NOT_EXCLUDED;

        if ( $debug ){
            print 'INSERT INTO rnaSeqLibraryAnnotatedSampleGeneResult: ', $annotatedSampleId,   ' - ',
                  $genes{ $speciesId}->{ $geneId }, ' - ',
                  'tpm',                                      ' - ',
                  $genesResults{$geneId}->{'TPM'},            ' - ',
                  $genesResults{$geneId}->{'estimatedCount'}, ' - ',
                  # no UMIs for full length scRNA-Seq
                  0,                                          ' - ',
                  $genesResults{$geneId}->{'zscore'},         ' - ',
                  $genesResults{$geneId}->{'pValue'},         ' - ',
                  $exclusion, "\n";
        }
        else {
            # pvalue and zscore can be null (if no read mapped). In this case BgeeCall retrieve "NA".
            # DBI use undef value to insert null in the database. That's why we modify "NA" to undef for zScore.
            # For the pvalue we decided replace NA with 1 in order to use this value as a datapoint to generate propagated calls 
            my $pValue = $genesResults{$geneId}->{'pValue'} eq "NA" ? 1 : $genesResults{$geneId}->{'pValue'};
            my $zscore = $genesResults{$geneId}->{'zscore'} eq "NA" ? undef : $genesResults{$geneId}->{'zscore'};
            $child_insAnnotatedSampleGeneResult->execute($annotatedSampleId,
                                # geneId is an ensembl ID, we need to get the bgeeGeneId
                                $genes{ $speciesId}->{ $geneId },
                                'tpm',
                                $genesResults{$geneId}->{'TPM'},
                                $genesResults{$geneId}->{'estimatedCount'},
                                # no UMIs for full length scRNA-Seq
                                0,
                                $zscore,
                                $pValue,
                                $exclusion,
                            )  or die $child_insAnnotatedSampleGeneResult->errstr;
        }
        $inserted++;

    }
    # used to commit after each library when condition and libraries were not inserted
    $child_bgee->commit;
    $child_insAnnotatedSampleGeneResult->finish();
    $child_bgee->disconnect();

    $pm->finish(0, { inserted => $inserted });

}

$pm->wait_all_children;

print "Done. $insertedGeneResults rows have been inserted in the rnaSeqLibraryAnnotatedSampleGeneResult table.\n";

exit 0;

