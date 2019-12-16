#!/usr/bin/env perl
# Julien Roux, created 24/07/09
# modified 15/08/16, JR: revamped whole script for bgee v14 with condition Ids
# last modified Dec 2016, Frederic Bastian: split script between insert_affy.pl and insert_affy_expression.pl.
# and adapt to new global vars in Utils.pm

# Perl core modules
use strict;
use warnings;
use diagnostics;
use File::Slurp;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

$| = 1; # no buffering of output
require 'mas5_utils.pl';
require 'affy_utils.pl';
require 'bgee_utils.pl';


# Define arguments & their default value
my ($bgee_connector, $ensembl_connector)                = ('', '');
my ($chipTypeQual, $affymetrixChipInformation)          = ('', '');
my ($chipType, $microarrayExperiment, $affymetrixChip)  = ('', '', '');
my ($annotations, $processed_mas5, $processed_schuster) = ('', '', '');
my ($expType)                                           = ('');
my ($sex_info)                                          = ('');
my ($Aport, $Sport)                                     = (0, 0);
my ($debug, $test, $resume)                             = (0, 0, 0);
my %opts = ('bgee=s'                    => \$bgee_connector,     # Bgee connector string
            'ensembl=s'                 => \$ensembl_connector,  # Ensembl connector string
            'chipType=s'                => \$chipType,
            'chipTypeQual=s'            => \$chipTypeQual,
            'microarrayExperiment=s'    => \$microarrayExperiment,
            'affymetrixChipInfo=s'      => \$affymetrixChipInformation,
            'affymetrixChip=s'          => \$affymetrixChip,
            'annotations=s'             => \$annotations,
            'processed_mas5=s'          => \$processed_mas5,
            'processed_schuster=s'      => \$processed_schuster,
            'sex_info=s'                => \$sex_info,           # generated_files/uberon/uberon_sex_info.tsv
            'exp=s'                     => \$expType,
            'Aport=i'                   => \$Aport,              # ID MAPPING anatomy socket port
            'Sport=i'                   => \$Sport,              # ID MAPPING stage   socket port
            'debug|v'                   => \$debug,              # more verbose output
            'test'                      => \$test,               # test: no insertion of chips and experiments, but script stops at line 344 because no chip is present in the database
            'resume'                    => \$resume,             # option allowing to resume the script if it was interrupted. Does not insert again chips and experiments, as well as probsets of chips that were already inserted. Be careful when using this option, you need to know what you are doing (see comments in code below).
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $ensembl_connector eq '' || $sex_info eq '' || $Aport == 0 || $Sport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD) -ensembl=\$(ENSCMD)  -chipType=\$(AFFY_CHIPTYPE_FILEPATH) -chipTypeQual=\$(AFFY_CHIPTYPEQUAL_FILEPATH) -microarrayExperiment=\$(MICROARRAY_EXPERIMENT_FILEPATH) -affymetrixChipInfo=\$(AFFY_CHIPINFO_FILEPATH) -affymetrixChip=\$(AFFY_CHIP_FILEPATH) -annotations=\$(ANNOTATIONPATH) -processed_mas5=\$(MAS5PATH) -processed_schuster=\$(SCHUSTERPATH) -sex_info=\$(UBERON_SEX_INFO_FILE_PATH) -Aport=\$(IDMAPPINGPORT) -Sport=\$(STGMAPPINGPORT  -exp=<expr/no_expr/both>
\t-bgee                 Bgee    connector string
\t-ensembl              Ensembl connector string
\t-chipTypeQual         chipTypeCorrespondencesAndQualityThresholds   pipeline   file
\t-affymetrixChipInfo   affymetrixChipInformation                     pipeline   file
\t-chipType             chipType                                      annotation file
\t-microarrayExperiment microarrayExperiment                          annotation file
\t-affymetrixChip       affymetrixChip                                annotation file
\t-annotations          Affy annotation path
\t-processed_mas5       processed_mas5     dir
\t-processed_schuster   processed_schuster dir
\t-sex_info             file containing sex-related info about anatomical terms
\t-exp                  analyse expression for expr|no_expr|both      (default: both)
\t-Aport                ID MAPPING anatomy socket port
\t-Sport                ID MAPPING stage   socket port
\t-debug                more verbose output
\t-test                 test mode, without any insertions/updates in the database
\n";
    exit 1;
}
die "Invalid exp option: [$expType]\n"  if ( $expType ne 'both' && $expType ne 'expr' && $expType ne 'no_expr' );

# Ensembl connection via Ensembl API/Registry
my $reg = Utils::connect_ensembl_registry($ensembl_connector, 0);
# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

# Load sex-related information needed for sub 'insert_get_condition'
my $anatSexInfo    = Utils::get_anat_sex_info($sex_info);
my $speciesSexInfo = Utils::get_species_sex_info($dbh);

##############
# dataSource #
##############
# Check AE not already in dataSource
my %bgeeDataSources = ();
my $selSource = $dbh->prepare('SELECT dataSourceName, dataSourceId FROM dataSource WHERE category =\'Affymetrix data source\';');
$selSource->execute()  or die $selSource->errstr;
while ( my @data = $selSource->fetchrow_array ){
    $bgeeDataSources{$data[0]} = $data[1];
}
$selSource->finish;

######################################
# Read and insert chipType new lines #
######################################
open(my $IN3, '<', "$chipType")  or die "Can't read file [$chipType]\n";
my %chipTypes;
my $line = <$IN3>;   #header
while ( defined ($line = <$IN3>) ){
    chomp $line;
    next if (( $line =~ /^#(.+)/ ) or ($line =~ /^\"#(.+)/ )); # commented chips
    my @tmp = map { bgeeTrim($_) }
              split(/\t/, $line);
    # chipTypeId  chipTypeName  speciesId  Ensembl_Xref_table  Comments
    $chipTypes{$tmp[0]}->{'name'} = $tmp[1];
    $chipTypes{$tmp[0]}->{'speciesId'}  = $tmp[2];
    $chipTypes{$tmp[0]}->{'xref'} = $tmp[3];
}
close $IN3;

my %chipTypeInfo = getChipTypesInformation($chipTypeQual);

my $selChipType = $dbh->prepare('SELECT * FROM chipType WHERE chipTypeId = ?');
my $insChipType = $dbh->prepare('INSERT INTO chipType (chipTypeId, chipTypeName, cdfName, isCompatible, qualityScoreThreshold, percentPresentThreshold) VALUES (?, ?, ?, ?, ?, ?)');

for my $chip ( keys %chipTypes ){
    $selChipType->execute($chip)  or die $selChipType->errstr;
    my @tmp = $selChipType->fetchrow_array;

    # if chip type not inserted before
    if ( $#tmp == -1 ){
        print "Inserting chipType [$chipTypes{$chip}->{'name'}]...\n" if ( $debug );
        if ( !defined $chipTypeInfo{$chip} ){
            #die "Error, no info in chipTypeCorrespondencesAndQualityThresholds for chip $chip\n";
            print "\t[$chip] skipped, no info on it in chipTypeCorrespondencesAndQualityThresholds, should never be used.\n" if ( $debug );
            # if no info on this chip, should never be used...
            next;
        }
        my $status = 1;
        if ( $chipTypeInfo{$chip}->{'status'} ne 'compatible' ){
            $status = 0;
        }
        my $qualityScore = 0;
        if ( $chipTypeInfo{$chip}->{'qualityScore'} ne '' ){
            $qualityScore = $chipTypeInfo{$chip}->{'qualityScore'} + 0;
        }
        my $percentPresent = 0;
        if ( $chipTypeInfo{$chip}->{'percentPresent'} ne '' ){
            $percentPresent = $chipTypeInfo{$chip}->{'percentPresent'} + 0;
        }
        if ( !$test and !$resume){
            $insChipType->execute($chip, $chipTypes{$chip}->{'name'}, $chipTypeInfo{$chip}->{'correspondence'}, $status, $qualityScore, $percentPresent)
                or warn $insChipType->errstr, " for [$chip] [$chipTypes{$chip}->{'name'}]\n";
        }
    }
}
$selChipType->finish;
$insChipType->finish;

###################################
# Read microarrayExperiment infos #
###################################
my %exp_to_insert;
my %dataSources;
my %tsv = %{ Utils::read_spreadsheet("$microarrayExperiment", "\t", 'csv', '"', 1) };
# For all lines in experimentId column
for my $line ( 0..$#{$tsv{'experimentId'}} ){
    $exp_to_insert{ $tsv{'experimentId'}[$line] }->{'name'}   = $tsv{'experimentName'}[$line];
    $exp_to_insert{ $tsv{'experimentId'}[$line] }->{'desc'}   = $tsv{'experimentDescription'}[$line];
    $exp_to_insert{ $tsv{'experimentId'}[$line] }->{'source'} = $tsv{'experimentSource'}[$line];
    $dataSources{ $tsv{'experimentSource'}[$line] }           = ();
}

# Check that all data sources present in the file were present in the database
for my $key ( keys %dataSources ){
    if ( !exists $bgeeDataSources{$key} || !defined $bgeeDataSources{$key} ){
        die "Data source [$key] not found in the Bgee database\n";
    }
}

#######################################################
# Get gene, organ, stages and conditions information #
#######################################################
# Get stages & anatEntityId from annotation sheets
my %tsv2 = %{ Utils::read_spreadsheet("$affymetrixChip", "\t", 'csv', '"', 1) };
# #chipId	experimentId	chipTypeId	normalizationTypeId	detectionTypeId	anatEntityId	organName	uberonId	uberonName	stageId	stageName	infoOrgan	infoStage	sampleTitle	sampleSource	sampleDescription	sampleCharacteristics	organAnnotationStatus	organBiologicalStatus	stageAnnotationStatus	stageBiologicalStatus	sex	strain	speciesId	comment	annotatorId	lastModificationDate	replicate	infoReplicate
my @Stg  = @{ $tsv2{'stageId'} };
my @Anat = @{ $tsv2{'uberonId'} };
my $doneAnat = Utils::get_anatomy_mapping(\@Anat, $Aport, 0);
my $doneStg  = Utils::get_anatomy_mapping(\@Stg,  $Sport, 0);

# Get already known conditions
my $conditions = Utils::query_conditions($dbh);

# Get simpler (upper level) stage equivalences
my $stage_equivalences = Utils::get_stage_equivalences($dbh);


# Get a hash of mapping from Ensembl gene ids to bgeeGeneId for all species
# First get an array of all (distinct) species with affymetrix data
my @all_species = do { my %seen; grep { !$seen{$_}++ } @{ $tsv2{'speciesId'} } };

my %genes;
# Get hash of geneId to bgeeGeneId mapping per species
foreach my $speciesId (@all_species) {
  $genes{$speciesId} = Utils::query_bgeeGene($dbh, $speciesId);
}

############################################################
# Insert in microarrayExperiment and affymetrixChip tables #
############################################################
# get additional information (quality score, percent present, ...) for all chips
my %affyChipsInfo = getAllChipsInfo($affymetrixChipInformation);

my $insExp      = $dbh->prepare('INSERT INTO microarrayExperiment (microarrayExperimentId, microarrayExperimentName, microarrayExperimentDescription, dataSourceId) VALUES (?, ?, ?, ?)');
my $insAffyChip = $dbh->prepare('INSERT INTO affymetrixChip (affymetrixChipId, microarrayExperimentId, chipTypeId, normalizationType, detectionType, conditionId, qualityScore, percentPresent, scanDate)
                                 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)');

for my $experiment ( sort keys %exp_to_insert ){
    print "Inserting microarrayExperiment [$experiment]...\n" if ( $debug );
    if ( !$test and !$resume){
        $insExp->execute($experiment, $exp_to_insert{$experiment}->{'name'}, $exp_to_insert{$experiment}->{'desc'}, $bgeeDataSources{$exp_to_insert{$experiment}->{'source'}})  or warn $insExp->errstr;
    }

    my %affymetrixChip;
    # go through chips of this experiment one by one
    for my $line ( 0..$#{$tsv2{'chipId'}} ){

        if ( $tsv2{'experimentId'}[$line] eq $experiment ){
            # skip commented chips
            next  if (( $tsv2{'chipId'}[$line] =~ /^#/ ) or ($tsv2{'chipId'}[$line] =~ /^\"#/));

            # filenames have spaces substituted by '_' (see separate_affy_processed_mas5.pl)
            $tsv2{'chipId'}[$line] =~ s{ }{_}g;

            # additional information (quality score, percent present, ...) for this chip
            my $affyChipInfo = undef;
            if ( defined $affyChipsInfo{ $tsv2{'experimentId'}[$line] }->{ $tsv2{'chipId'}[$line] } ){
                $affyChipInfo = \%{$affyChipsInfo{ $tsv2{'experimentId'}[$line] }->{ $tsv2{'chipId'}[$line] }};
            }
            # can this chip be inserted?
            if ( isChipIncompatibleOrLowQuality($tsv2{'normalizationTypeId'}[$line], $affyChipInfo, \%chipTypeInfo, $tsv2{'chipTypeId'}[$line]) ){
                warn "\tIncompatible or low quality file skipped: expId: [$tsv2{'experimentId'}[$line]] - chipId [$tsv2{'chipId'}[$line]]\n";
                next;
            }
            if ( !defined $affyChipsInfo{ $tsv2{'experimentId'}[$line] }->{ $tsv2{'chipId'}[$line] } ){
                warn "\tWarning, file skipped, information was not generated for expId: [$tsv2{'experimentId'}[$line]] - chipId [$tsv2{'chipId'}[$line]]\n";
                next;
            }

            # Are annotated organ stages present in the database?
            if ( !exists $doneAnat->{ $tsv2{'uberonId'}[$line] } || $doneAnat->{ $tsv2{'uberonId'}[$line] } eq '' ){
                warn "[$tsv2{'uberonId'}[$line]] unmapped organ id\n";
                next;
            }
            if ( !exists $doneStg->{ $tsv2{'stageId'}[$line] }  || $doneStg->{ $tsv2{'stageId'}[$line] }  eq '' ){
                warn "[$tsv2{'stageId'}[$line]] unmapped stage id\n";
                next;
            }

            # chip Id used as hash key
            my $affy_id;
            if ( $tsv2{'chipId'}[$line] =~ /(.+)\.cel$/i ){ # case insensitive
              $affy_id = $1;
            } else {
              $affy_id = $tsv2{'chipId'}[$line];
            }

            $affymetrixChip{ $affy_id }->{'chipTypeId'}     = $tsv2{'chipTypeId'}[$line];
            $affymetrixChip{ $affy_id }->{'filename'}       = $tsv2{'chipId'}[$line];
            $affymetrixChip{ $affy_id }->{'norm'}           = $tsv2{'normalizationTypeId'}[$line] == 2 ? 'gcRMA'
                                                                           : $tsv2{'normalizationTypeId'}[$line] == 1 ? 'MAS5'
                                                                           : $tsv2{'normalizationTypeId'}[$line] == 3 ? 'RMA'
                                                                           : undef;
            $affymetrixChip{ $affy_id }->{'detect'}         = $tsv2{'detectionTypeId'}[$line] == 2 ? 'Schuster'
                                                                           : $tsv2{'detectionTypeId'}[$line] == 1 ? 'MAS5'
                                                                           : undef;
            # Normalize the encoding of sex
            my $sex = $tsv2{'sex'}[$line] eq 'F'     ? $Utils::FEMALE_SEX
                    : $tsv2{'sex'}[$line] eq 'M'     ? $Utils::MALE_SEX
                    : $tsv2{'sex'}[$line] eq 'H'     ? $Utils::HERMAPHRODITE_SEX
                    : $tsv2{'sex'}[$line] eq 'U'     ? $Utils::NOT_ANNOTATED_SEX
                    : $tsv2{'sex'}[$line] eq 'mixed' ? $Utils::MIXED_SEXES
                    : $tsv2{'sex'}[$line] eq 'NA'    ? $Utils::NA_SEX
                    : $Utils::NA_SEX ;

            # Insert a new conditionId if not already existing
            # Update the hash of existing conditions
            # Get conditionId/exprMappedConditionId for this chip
            my $condKeyMap;
            ($condKeyMap, $conditions) = Utils::insert_get_condition($dbh,
                                                                     $conditions,
                                                                     $stage_equivalences,
                                                                     $doneAnat->{ $tsv2{'uberonId'}[$line] },
                                                                     $doneStg->{  $tsv2{'stageId'}[$line] },
                                                                     $tsv2{'speciesId'}[$line],
                                                                     $sex,
                                                                     $tsv2{'strain'}[$line],
                                                                     $anatSexInfo, $speciesSexInfo,
                                                                     $experiment, '');

            # record the conditions in hash
            $affymetrixChip{ $affy_id }->{'conditionId'}           = $condKeyMap->{'conditionId'};

            # add last info to the hash
            if((defined $affyChipsInfo{$experiment}->{ $tsv2{'chipId'}[$line] }->{'qualityScore'}) and ($affyChipsInfo{$experiment}->{ $tsv2{'chipId'}[$line] }->{'qualityScore'} ne '')){
              $affymetrixChip{ $affy_id }->{'qualityScore'} = $affyChipsInfo{$experiment}->{ $tsv2{'chipId'}[$line] }->{'qualityScore'};
            } else {
              $affymetrixChip{ $affy_id }->{'qualityScore'} = 0;
            }
            $affymetrixChip{ $affy_id }->{'percentPresent'} = $affyChipsInfo{$experiment}->{ $tsv2{'chipId'}[$line] }->{'percentPresent'};
            if (defined $affyChipsInfo{$experiment}->{ $tsv2{'chipId'}[$line] }->{'scanDate'}){
              $affymetrixChip{ $affy_id }->{'scanDate'}     = $affyChipsInfo{$experiment}->{ $tsv2{'chipId'}[$line] }->{'scanDate'};
            } else {
              $affymetrixChip{ $affy_id }->{'scanDate'}     = ''
            }
        }
    }

    ## All chips to be inserted from this experiment
    for my $affy_id ( sort keys %affymetrixChip ){
        if ( !$test and !$resume){
            # bgeeAffymetrixChipId will be auto incrementd
            $insAffyChip->execute($affy_id, $experiment,
                                  $affymetrixChip{$affy_id}->{'chipTypeId'},
                                  $affymetrixChip{$affy_id}->{'norm'},
                                  $affymetrixChip{$affy_id}->{'detect'},
                                  $affymetrixChip{$affy_id}->{'conditionId'}, # we insert the low-level condition Id here
                                  $affymetrixChip{$affy_id}->{'qualityScore'},
                                  $affymetrixChip{$affy_id}->{'percentPresent'},
                                  $affymetrixChip{$affy_id}->{'scanDate'})
                or warn $insAffyChip->errstr;
        }
    }
}
$insExp->finish;
$insAffyChip->finish;

##################################
# Query the table affymetrixChip #
##################################
# store all chips in a hash:
# bgeeAffymetrixChipId -> experimentId/chipId/chipTypeId/norm
my %all_chips;
# store all chip types in a hash where keys are chip type IDs.
my %all_chipTypes;
# all experiments
my $selAffyChip = $dbh->prepare('SELECT bgeeAffymetrixChipId, microarrayExperimentId, affymetrixChipId, chipTypeId, normalizationType FROM affymetrixChip');
$selAffyChip->execute()  or die $selAffyChip->errstr;
die "Nothing in affymetrixChip table\n"  if ( !$test && $selAffyChip->rows() == 0 );

## We query conditionId, but for insertion in expression table, we nee dthe high-level conditions (exprMappedConditionId). We use a correspondance hash to map conditionId to exprMappedConditionId:
my $condition_equivalences = Utils::get_condition_equivalences($dbh);

while ( my @data = $selAffyChip->fetchrow_array ){
    $all_chips{ $data[0] }->{'experimentId'} = $data[1];
    $all_chips{ $data[0] }->{'chipId'}       = $data[2];
    $all_chips{ $data[0] }->{'chipTypeId'}   = $data[3];
    $all_chips{ $data[0] }->{'norm'}         = $data[4];
    $all_chipTypes{$data[3]}++;
}
$selAffyChip->finish;

## Bgee connection because not needed for a while
$dbh->disconnect;

########################
# Retrieve the mapping #
########################
my %annotation;

for my $chip ( sort keys %all_chipTypes ){
    print "Retrieving mapping of [$chip]...\n" if ( $debug );
    # No mapping
    if ( !exists $chipTypes{$chip}->{'xref'} || $chipTypes{$chip}->{'xref'} =~ /^\.out$/ ){
        warn "No mapping for [$chip] [$chipTypes{$chip}->{'xref'}]\;";
        next;
    }

    my %annot_pbsets;
    # if the annotation is in a local file
    if ( $chipTypes{$chip}->{'xref'} =~ /^.+\.out$/ ){
        # parse the annotation file
        my $problems = 0;
        for my $line ( read_file("$annotations/".$chipTypes{$chip}->{'xref'}, chomp => 1) ){
            my @tmp = split(/\t/, $line);
            # check that the gene is in Bgee
            if ( exists $genes{$chipTypes{$chip}->{'speciesId'}}->{$tmp[1]} ){
                # probeset -> bgeeGeneId
                $annot_pbsets{$tmp[0]}->{$genes{$chipTypes{$chip}->{'speciesId'}}->{$tmp[1]}}++;
            }
            else {
                $problems++;
            }
        }
        print "\tWarning! The mapping in the .out file is outdated ($problems probesets are no longer mapped) [$chip] [$chipTypes{$chip}->{'xref'}]\n"  if ( $problems > 0 );
    }
    # if the annotation is in an Ensembl table
    else {
        # Get adaptors for this organism (Array & Transcript)
        my $array_adaptor = $reg->get_adaptor($chipTypes{$chip}->{'speciesId'}, 'funcgen', 'Array');
        my $tx_adaptor    = $reg->get_adaptor($chipTypes{$chip}->{'speciesId'}, 'core',    'Transcript');

        # Get array object according to its name in chipType sheet (Xref field)
        my $array = $array_adaptor->fetch_by_name_vendor( $chipTypes{$chip}->{'xref'} );
        for my $pset ( @{ $array->get_all_ProbeSets() } ){
            my %distinct;
            #NOTE try with fetch_all_by_linked_transcript_Gene() when it will be repaired in the API
            for my $trans ( @{ $pset->get_all_Transcript_DBEntries } ){
                my $transcript = $tx_adaptor->fetch_by_stable_id( $trans->primary_id );
                if (defined $transcript){
                  $distinct{ $transcript->get_Gene->stable_id } = $pset->name;
                }
            }
            # probeset -> bgeeGeneId
            map { $annot_pbsets{ $distinct{$_} }->{ $genes{$chipTypes{$chip}->{'speciesId'}}->{$_} }++ } keys %distinct;
        }
        print "\tWarning! There seems to be a bad mapping for this chip [$chip] [$chipTypes{$chip}->{'xref'}]: ". keys(%annot_pbsets) ." probesets\n"  if ( keys(%annot_pbsets) < 100 );
    }

    # keep only probesets mapping to a unique gene id
    for my $pbset ( keys %annot_pbsets ){
        if ( !(scalar(keys %{$annot_pbsets{$pbset}}) eq 1) ){
            delete $annot_pbsets{$pbset};
        }
        else {
            for my $gene ( keys %{$annot_pbsets{$pbset}} ){
                $annot_pbsets{$pbset} = $gene;
                last;
            }
        }
    }
    # Keep the mapping to use it later
    $annotation{$chip} = \%annot_pbsets;
}
print "\n\n" if ( $debug );


###########################
# restart Bgee connection #
###########################
$dbh = Utils::connect_bgee_db($bgee_connector);

##################
# -resume option #
##################
# If this option is enabled, it means that some probesets have already been inserted
# and the script was stopped before completion. We want to start from where it has stopped,
# so we query the table affymetrixProbeset to get all chips with probesets already inserted.
# This assumes that each chip was inserted in full.

my %insertedChips;
my $selInsertedChip = $dbh->prepare('SELECT DISTINCT bgeeAffymetrixChipId FROM affymetrixProbeset;');
$selInsertedChip->execute()  or die $selInsertedChip->errstr;
while ( my @data = $selInsertedChip->fetchrow_array ){
  $insertedChips{$data[0]} = ();
}
$selInsertedChip->finish;

###############################
# Insert probesets            #
###############################
print("Inserting probesets...\n") if ( $debug );

# affymetrixProbeset
my $insAffyPset = $dbh->prepare('INSERT INTO affymetrixProbeset
                                (affymetrixProbesetId, bgeeAffymetrixChipId, bgeeGeneId,
                                normalizedSignalIntensity, detectionFlag, affymetrixData, reasonForExclusion)
                                VALUES (?, ?, ?, ?, ?, ?, ?)');

CHIP: for my $bgeeAffymetrixChipId (sort { $a <=> $b } keys %all_chips) {
  print("\tChip $bgeeAffymetrixChipId\n") if ( $debug );
  # see comments above for the use of this option
  next CHIP if ( $resume and exists($insertedChips{$bgeeAffymetrixChipId}) );

  my $chipTypeId = $all_chips{$bgeeAffymetrixChipId}->{'chipTypeId'};
  # return a hash probesetId -> call/signal/quality
  my $chip_pbsets = extractProbesetsFromFile($all_chips{$bgeeAffymetrixChipId}->{'experimentId'},
                                             $all_chips{$bgeeAffymetrixChipId}->{'chipId'},
                                             $all_chips{$bgeeAffymetrixChipId}->{'norm'},
                                             $processed_mas5, $processed_schuster, 0);

  PROBESET: for my $pbsetId ( keys %{$chip_pbsets} ) {
    if ( !exists $annotation{$chipTypeId}->{$pbsetId} || !defined $annotation{$chipTypeId}->{$pbsetId} ) {
      next PROBESET;
    }
    # Note: pre-filtering exclusion is now managed in the script insert_affy_expression.pl,
    # it used to be managed here.
    my $reasonForExclusion = $Utils::CALL_NOT_EXCLUDED;
    if ( !$test ) {
      $insAffyPset->execute($pbsetId, $bgeeAffymetrixChipId,
                            $annotation{$chipTypeId}->{$pbsetId},  #bgeeGeneId
                            $chip_pbsets->{$pbsetId}->{'signal'},
                            $chip_pbsets->{$pbsetId}->{'call'},
                            $chip_pbsets->{$pbsetId}->{'quality'},
                            $reasonForExclusion)
      or warn $insAffyPset->errstr, " for [$bgeeAffymetrixChipId][$pbsetId]\n";
    }
  }
}
$insAffyPset->finish;

$dbh->disconnect;
exit 0;

# Remark: we don't check that all processed files have the same length
# for some reasons, this is not always the case (bad submission in AE)

# TODO: add real debug mode where queries are printed, not executed
