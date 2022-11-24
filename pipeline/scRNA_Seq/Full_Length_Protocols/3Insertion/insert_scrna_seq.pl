#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use File::Slurp;
use lib "$FindBin::Bin/../../.."; # Get lib path for Utils.pm
use Utils;
$| = 1; # no buffering of output

# Julien Wollbrett, Feb 2021
##
## This script is highly inspired from the insertion one of bulk RNASeq created by Frederic Bastian.
##
## Some parts are identical with bulk RNASeq as both pipelines use BgeeCall to generate calls
## NOTE : Depending on how close this script will be with the bulk RNASeq one it will require a refactoring once bgee 15.0 is ready
## differences with bulk RNASeq are listed below in order to more easily refactor:
## - no sample excluded file
## - One more annotation columns for condition (Cell ID/Cell name) 
## - sample info file not the same (no libraryInfo, genomePath, database)
## - for the moment does not use a generated RNASeq experiment but directly the file manually created by annotators (does not contain info about worm experiments). This file has more columns than the one generated during bulk RNASeq pipeline
## - 


#####################################################################


my $abundance_file = 'gene_level_abundance+calls.tsv';

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($extraMapping)   = ('');
my ($scRnaSeqLibrary, $all_results, $sex_info)  = ('', '', '');
my ($scRnaSeqExperiment, $library_info, $excluded_libraries, $library_stats, $report_info) = ('', '', '', '', '');
my ($debug)                      = (0);
my ($Aport, $Sport)              = (0, 0);
my %opts = ('bgee=s'                => \$bgee_connector,           # Bgee connector string
            'scRnaSeqExperiment=s'  => \$scRnaSeqExperiment,       # scRNAseqExperiment
            'library_info=s'        => \$library_info,             # NEW_scRNASeq_sample_info.txt file
            'excluded_libraries=s'  => \$excluded_libraries,       # Discard_scRNASeq_sample_info.txt file
            'library_stats=s'       => \$library_stats,            # presence_absence_all_samples.txt
            'report_info=s'         => \$report_info,              # reports_info_all_samples.txt
            'all_results=s'         => \$all_results,              # path to Bgeecall calls directory
            'sex_info=s'            => \$sex_info,                 # generated_files/uberon/uberon_sex_info.tsv
            'extraMapping=s'        => \$extraMapping,             # Extra mapping for too up-to-date ontology terms
            'debug'                 => \$debug,
            'Aport=i'               => \$Aport,                    # ID MAPPING anatomy port socket
            'Sport=i'               => \$Sport,                    # ID MAPPING stage   port socket
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $scRnaSeqExperiment eq '' || $library_info eq ''  || $excluded_libraries eq '' || $library_stats eq '' || $report_info eq '' || $all_results eq '' || $sex_info eq '' || $extraMapping eq '' || $Aport == 0 || $Sport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g., $0  -bgee=\$(BGEECMD) -scRnaSeqExperiment=scRNAseqExperiment.tsv -library_info=\$(SC_RNASEQ_SAMPINFO_PASS_FILEPATH) -excluded_libraries=\$(SC_RNASEQ_SAMPINFO_NOT_PASS_FILEPATH) -library_stats=\$(SC_RNASEQ_SAMP_STATS_FL) -report_info=\$(SC_RNASEQ_REPORT_INFO_FL) -all_results=\$(SC_RNASEQ_ALL_RES_FL) -sex_info=\$(UBERON_SEX_INFO_FILE_PATH) -extraMapping=\$(EXTRAMAPPING_FILEPATH) -Aport=\$(IDMAPPINGPORT) -Sport=\$(STGMAPPINGPORT)    > $@.tmp 2>warnings.$@
\t-bgee                Bgee connector string
\t-scRnaSeqExperiment  single cell RNAseq experiment file
\t-library_info        NEW_scRNASeq_sample_info.txt file
\t-excluded_libraries  Discard_scRNASeq_sample_info.txt file
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

require("$FindBin::Bin/../../rna_seq_utils.pl");

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

# Library info generated after QC. Contains only library that passed the QCs
my %libraries         = getAllFullLengthScRnaSeqLibrariesInfo($library_info);
print "\t", scalar keys %libraries, " experiments with libraries mapped.\n";
my %excludedLibraries = getAllFullLengthScRnaSeqLibrariesInfo($excluded_libraries);
print "\t", scalar keys %excludedLibraries, " experiments with libraries discarded after quality control.\n";

my $count_libs = 0;
my %all_species; # record all species
for my $expId ( sort keys %libraries ){
    for my $libraryId ( sort keys %{$libraries{$expId}} ){
        $all_species{$libraries{$expId}->{$libraryId}->{'speciesId'}}++;
        $count_libs++;
        unless ( -s "$all_results/$libraryId/$abundance_file" ){
            die "Missing or empty processed data file for library $libraryId! Please check that the transfer from cluster was successful. Otherwise this library should maybe be added to the file of excluded libraries?\n";
        }
    }
}
print "\t", $count_libs, " libraries mapped and to be inserted.\n";
print "\t", scalar keys %libraries, " experiment mapped and to be inserted.\n";

# Library info generated by the pipeline. Same function than for bulk RNASeq
my %librariesStats    = getAllRnaSeqLibrariesStats($library_stats);
# Library info generated by the pipeline. Same function than for bulk RNASeq
my %reportInfo        = getAllRnaSeqReportInfo($report_info);

# Experiment annotation coming from flat files
my %tsv = %{ Utils::read_spreadsheet("$scRnaSeqExperiment", "\t", 'csv', '"', 1) };
# use same function than bulk RNASeq.
# XXX: Not sure this function works if Experiments are commented (i.e line starts with #)
my %experiments       = getAllAnnotatedExperiments2( \%tsv );

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

######################
# INSERT EXPERIMENTS #
######################
print "Inserting experiments...\n";
my $insExp = $bgee->prepare('INSERT INTO scRnaSeqFullLengthExperiment (scRnaSeqFullLengthExperimentId, scRnaSeqFullLengthExperimentName, scRnaSeqFullLengthExperimentDescription, dataSourceId) VALUES (?, ?, ?, ?)');
# experiment annotation file contains a mix of full length and target base single cell rnaseq. 
# the list of experiments to insert is then defined from experiments of libraries to insert or libraries to discard.
my %experimentIds = ();
for my $expId ( sort keys %libraries ){
  $experimentIds{$expId} = 1;
}
for my $expId ( sort keys %excludedLibraries ){
  $experimentIds{$expId} = 1;
}

for my $expId ( sort keys %experimentIds ){
    print "\t$expId\n";
    # sanity check that the experiment is present in experiment annotations
    if(!exists $experiments{$expId}->{'name'}) {
      die "libraries of experiment [$expId] has to be inserted but the experiment is not described in the annotation file";
    }
    if ( $debug ){
        binmode(STDOUT, ':utf8');
        print 'INSERT INTO scRnaSeqFullLengthExperiment: ',
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


######################
# INSERT PLATFORMS   #
######################
print "Inserting platforms...\n";
# Retrieve all platforms used for the analyzed libraries
my %platforms = ();
for my $expId ( keys %libraries ){
    foreach my $libraryId ( keys %{$libraries{$expId}} ){
        if ( $libraries{$expId}->{$libraryId}->{'platform'} ne '' ){
            $platforms{$libraries{$expId}->{$libraryId}->{'platform'}}++;
        }
    }
}

my $insPlatform = $bgee->prepare('INSERT INTO scRnaSeqFullLengthPlatform (scRnaSeqFullLengthPlatformId, scRnaSeqFullLengthPlatformDescription) VALUES (?, ?)');
# now insert the platform(s)
for my $platformId ( sort keys %platforms ){
    print "\t$platformId\n";
    if ( $debug ){
        print 'INSERT INTO scRnaSeqFullLengthPlatform: ', $platformId, ' - ', "\n";
    }
    else {
        $insPlatform->execute($platformId, '')  or die $insPlatform->errstr;
    }
}
$insPlatform->finish();
print "Done\n\n";

######################
# INSERT ĜENOTYPES   #
######################
print "Inserting genotypes...\n";
#Retrieve already inserted genotypes (useful for incremental updates)
my %insertedGenotypes = ();
my $retrieveGenotype = $bgee->prepare('SELECT genotypeId, genotypeName from genotype');
$retrieveGenotype->execute()  or die $retrieveGenotype->errstr;
while ( my @data = $retrieveGenotype->fetchrow_array ){
    $insertedGenotypes{$data[1]} = $data[0];
}
# Retrieve all genotypes used for the analyzed libraries and not already inserted
my %genotypes = ();
for my $expId ( keys %libraries ){
    foreach my $libraryId ( keys %{$libraries{$expId}} ){
        if ( $libraries{$expId}->{$libraryId}->{'genotype'} ne '' &&
            !exists $insertedGenotypes{$libraries{$expId}->{$libraryId}->{'genotype'}}){
            $genotypes{$libraries{$expId}->{$libraryId}->{'genotype'}}++;
        }
    }
}
my $insGenotype = $bgee->prepare('INSERT INTO genotype (genotypeName) VALUES (?)');
# now insert the genotype(s)
for my $genotypeName ( sort keys %genotypes ){
    if ( $debug ){
        print 'INSERT INTO genotype: ', $genotypeName, "\n";
    }
    else {
        $insGenotype->execute($genotypeName) or die $insGenotype->errstr;
    }
}
$insGenotype->finish();
#retrieve again genotypes present in the database now that all of them are inserted
$retrieveGenotype->execute()  or die $retrieveGenotype->errstr;
while ( my @data = $retrieveGenotype->fetchrow_array ){
    $insertedGenotypes{$data[1]} = $data[0];
}
$retrieveGenotype->finish();
print "Done\n\n";


################################################
# GET GENE INTERNAL IDS                        #
# GET ORGAN, STAGE, AND CONDITIONS INFORMATION #
################################################

my %genes;
my @Stg;
my @Anat;

# go over all libraries to check all species with data and retrieve all
# used stages and anat. entities
for my $expId ( keys %libraries ){
    foreach my $libraryId ( keys %{$libraries{$expId}} ){
        $genes{$libraries{$expId}->{$libraryId}->{'speciesId'}} = ();
        if(defined($libraries{$expId}->{$libraryId}->{'stageId'})) {
            push(@Stg, $libraries{$expId}->{$libraryId}->{'stageId'});
        }
        if(defined($libraries{$expId}->{$libraryId}->{'uberonId'})) {
            push(@Anat, $libraries{$expId}->{$libraryId}->{'uberonId'});
        }
        if(defined($libraries{$expId}->{$libraryId}->{'cellTypeId'})) {
            push(@Anat, $libraries{$expId}->{$libraryId}->{'cellTypeId'});
        }
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
my $insLib = $bgee->prepare('INSERT INTO scRnaSeqFullLengthLibrary (scRnaSeqFullLengthLibraryId, scRnaSeqFullLengthExperimentId,
                             scRnaSeqFullLengthPlatformId, genotypeId, conditionId, tpmThreshold,
                             allGenesPercentPresent, proteinCodingGenesPercentPresent,
                             intergenicRegionsPercentPresent, meanTpmReferenceIntergenicDistribution, 
                             sdTpmReferenceIntergenicDistribution, pValueThreshold, allReadsCount, 
                             mappedReadsCount, minReadLength, maxReadLength, libraryType, libraryOrientation)
                             VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)');


# Excluded libraries

my $discarded_reason = "discarded after quality control";
my $insExcludedLib = $bgee->prepare('INSERT INTO scRnaSeqFullLengthLibraryDiscarded (scRnaSeqFullLengthLibraryId, scRnaSeqFullLengthLibraryDiscardReason) VALUES (?, ?)');
for my $exp ( sort keys %excludedLibraries ){
    for my $libraryId ( sort keys %{$excludedLibraries{$exp}}){
        if ( $debug ){
            print 'INSERT INTO scRnaSeqFullLengthLibraryDiscarded: ', $libraryId, ' - ', "$discarded_reason\n";
        }
        else {
            $insExcludedLib->execute($libraryId, $discarded_reason)  or die $insExcludedLib->errstr;
        }
    }
}

# query for genes results insertion
my $insResult = $bgee->prepare('INSERT INTO scRnaSeqFullLengthResult (scRnaSeqFullLengthLibraryId, bgeeGeneId, fpkm, tpm,
                                readsCount, pValue, zScore, detectionFlag, reasonForExclusion)
                                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)');

my $inserted = 0;

# used to commit after each library when condition and libraries were not inserted
#print "disable autocommit. Manually commit for each library\n";
#$bgee->{AutoCommit} = 0;

for my $expId ( sort keys %libraries ){
    LIBRARY:
    for my $libraryId ( sort keys %{$libraries{$expId}} ){
        if ( !exists $reportInfo{$libraryId} ){
            warn "Report file does not contain this library [$libraryId]\n";
            next LIBRARY;
        }
        print "\t$expId $libraryId\n";

        # Remap to extra mapping if any
        $libraries{$expId}->{$libraryId}->{'uberonId'}   = $extra{ $libraries{$expId}->{$libraryId}->{'uberonId'} }   || $libraries{$expId}->{$libraryId}->{'uberonId'};
        $libraries{$expId}->{$libraryId}->{'cellTypeId'} = $extra{ $libraries{$expId}->{$libraryId}->{'cellTypeId'} } || $libraries{$expId}->{$libraryId}->{'cellTypeId'};
        $libraries{$expId}->{$libraryId}->{'stageId'}    = $extra{ $libraries{$expId}->{$libraryId}->{'stageId'} }    || $libraries{$expId}->{$libraryId}->{'stageId'};

        if ( !exists $doneAnat->{$libraries{$expId}->{$libraryId}->{'uberonId'}} || $doneAnat->{$libraries{$expId}->{$libraryId}->{'uberonId'}} eq '' ){
            warn "[$libraries{$expId}->{$libraryId}->{'uberonId'}] unmapped organ id for [$libraryId]\n";
            next LIBRARY;
        }
        if ( !exists $doneAnat->{$libraries{$expId}->{$libraryId}->{'cellTypeId'}} || $doneAnat->{$libraries{$expId}->{$libraryId}->{'cellTypeId'}} eq '' ){
            warn "[$libraries{$expId}->{$libraryId}->{'cellTypeId'}] unmapped cell type id for [$libraryId]\n";
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
                                                                 $doneAnat->{$libraries{$expId}->{$libraryId}->{'uberonId'}},
                                                                 $doneStg->{$libraries{$expId}->{$libraryId}->{'stageId'}},
                                                                 $libraries{$expId}->{$libraryId}->{'speciesId'},
                                                                 $libraries{$expId}->{$libraryId}->{'sex'},
                                                                 $libraries{$expId}->{$libraryId}->{'strain'},
                                                                 $anatSexInfo, $speciesSexInfo,
                                                                 $libraryId, '',
                                                                 $doneAnat->{$libraries{$expId}->{$libraryId}->{'cellTypeId'}},
                                                                );
        # We consider the fine-grained (low-level) conditionId for insertion: $condKeyMap->{'conditionId'}

        # insert sample
        if ( $debug ){
            print 'INSERT INTO scRnaSeqFullLengthLibrary: ', $libraryId,                        ' - ',
                  $expId, ' - ', $libraries{$expId}->{$libraryId}->{'platform'},    ' - ',
                  $insertedGenotypes{$libraries{$expId}->{$libraryId}->{'genotype'}}, ' - ',
                  $condKeyMap->{'conditionId'},                                     ' - ',
                  $librariesStats{$libraryId}->{'cutoffTPM'},                       ' - ',
                  $librariesStats{$libraryId}->{'allGenesPercentPresent'},          ' - ',
                  $librariesStats{$libraryId}->{'proteinCodingPercentPresent'},     ' - ',
                  $librariesStats{$libraryId}->{'intergenicRegionsPercentPresent'}, ' - ',
                  $librariesStats{$libraryId}->{'meanIntergenic'},                  ' - ',
                  $librariesStats{$libraryId}->{'sdIntergenic'},                    ' - ',
                  $librariesStats{$libraryId}->{'pValueThreshold'},                 ' - ',
                  $reportInfo{$libraryId}->{'allReadsCount'},                       ' - ',
                  $reportInfo{$libraryId}->{'mappedReadsCount'},                    ' - ',
                  $reportInfo{$libraryId}->{'minReadLength'},                       ' - ',
                  $reportInfo{$libraryId}->{'maxReadLength'},                       ' - ',
                  $libraries{$expId}->{$libraryId}->{'libraryType'},                ' - ',
                  "NA\n";
        }
        else {
            $insLib->execute($libraryId,
                             $expId,
                             $libraries{$expId}->{$libraryId}->{'platform'},
                             $insertedGenotypes{$libraries{$expId}->{$libraryId}->{'genotype'}},
                             $condKeyMap->{'conditionId'},
                             $librariesStats{$libraryId}->{'cutoffTPM'},
                             $librariesStats{$libraryId}->{'allGenesPercentPresent'},
                             $librariesStats{$libraryId}->{'proteinCodingPercentPresent'},
                             $librariesStats{$libraryId}->{'intergenicRegionsPercentPresent'},
                             $librariesStats{$libraryId}->{'meanIntergenic'},
                             $librariesStats{$libraryId}->{'sdIntergenic'},
                             $librariesStats{$libraryId}->{'pValueThreshold'},
                             $reportInfo{$libraryId}->{'allReadsCount'},
                             $reportInfo{$libraryId}->{'mappedReadsCount'},
                             $reportInfo{$libraryId}->{'minReadLength'},
                             $reportInfo{$libraryId}->{'maxReadLength'},
                             $libraries{$expId}->{$libraryId}->{'libraryType'},
                            'NA',
                           )  or die $insLib->errstr;
        }

        # TODO: implement insertion of runs?
        

        # insert genes results
        my %genesResults = getGenesResults("$all_results/$libraryId/$abundance_file");
        for my $geneId ( keys %genesResults ){

            my $exclusion = $Utils::CALL_NOT_EXCLUDED;
            # for full-length single cell RNA-Seq absent calls are inserted but excluded
            if($genesResults{$geneId}->{'expressionCall'} eq "absent") {
                $exclusion = $Utils::EXCLUDED_FOR_ABSENT_CALLS;
            }
            
            $inserted++;
            
            if ( $debug ){
                print 'INSERT INTO scRnaSeqFullLengthResult: ', $libraryId,   ' - ', $genes{ $libraries{$expId}->{$libraryId}->{'speciesId'}}->{ $geneId }, ' - ',
                      $genesResults{$geneId}->{'FPKM'},           ' - ',
                      $genesResults{$geneId}->{'TPM'},            ' - ',
                      $genesResults{$geneId}->{'estimatedCount'}, ' - ',
                      $genesResults{$geneId}->{'pValue'},         ' - ',
                      $genesResults{$geneId}->{'zscore'},         ' - ',
                      $genesResults{$geneId}->{'expressionCall'}, ' - ',
                      $exclusion, "\n";
            }
            else {

                # pvalue and zscore can be null (if no read mapped). In this case BgeeCall retrieve "NA".
                # DBI use undef value to insert null in the database. That's why we modify "NA" to undef.
                if ($genesResults{$geneId}->{'pValue'} eq "NA") {
                  $genesResults{$geneId}->{'pValue'} = undef;
                }
                if ($genesResults{$geneId}->{'zscore'} eq "NA") {
                  $genesResults{$geneId}->{'zscore'} = undef;
                }
                
                $insResult->execute($libraryId,
                                    # geneId is an ensembl ID, we need to get the bgeeGeneId
                                    $genes{ $libraries{$expId}->{$libraryId}->{'speciesId'}}->{ $geneId },
                                    $genesResults{$geneId}->{'FPKM'},
                                    $genesResults{$geneId}->{'TPM'},
                                    $genesResults{$geneId}->{'estimatedCount'},
                                    $genesResults{$geneId}->{'pValue'},
                                    $genesResults{$geneId}->{'zscore'},
                                    $genesResults{$geneId}->{'expressionCall'},
                                    $exclusion,
                                   )  or die $insResult->errstr;
            }

        }
        #used to commit after each library when condition and libraries were not inserted
        #$bgee->commit;
    }
}
#used to commit after each library when condition and libraries were not inserted
#print "reactivate autocommit\n";
#$bgee->{AutoCommit} = 1;

$insLib->finish();
$insResult->finish();
print "Done. You should have $inserted rows in the scRnaSeqFullLengthResult table.\nExiting\n";

exit 0;

