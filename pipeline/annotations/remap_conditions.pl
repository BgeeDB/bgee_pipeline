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

#####################################################################
# THIS SCRIPT ASSUMES THAT:
# * A TABLE 'remapStage' HAS BEEN CREATED AND POPULATED,
# * AND A TABLE 'remapCond' EXISTS
# (see beginning of 'pipeline/post_processing/remap_conditions.sql').
# WE NEED THIS SCRIPT AND NOT A SIMPLE SQL FILE BECAUSE WE NEED TO USE
# THE FUNCTION insert_get_condition TO INSERT/REMAP CONDITIONS
#
# Frederic Bastian, created August 2019
#####################################################################

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($extraMapping)   = ('');
my ($rnaSeqLibrary, $all_results, $sex_info)  = ('', '', '');
my ($rnaSeqExperiment, $library_info, $excluded_libraries, $library_stats, $report_info) = ('', '', '', '', '');
my ($debug)                      = (0);
my ($Aport, $Sport)              = (0, 0);
my %opts = ('bgee=s'                => \$bgee_connector,     # Bgee connector string
            'sex_info=s'            => \$sex_info,           # generated_files/uberon/uberon_sex_info.tsv
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $sex_info eq '' ){
    print "\n\tInvalid or missing argument:
\te.g., $0  -bgee=\$(BGEECMD) -sex_info=\$(UBERON_SEX_INFO_FILE_PATH) > $@.tmp 2>warnings.$@
\t-bgee                Bgee connector string
\t-sex_info            file containing sex-related info about anatomical terms
\n";
    exit 1;
}


# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

##############################
# CONDITIONS TO UPDATE/REMAP #
##############################
# incorrect condition ID -> existing matching condition ID (0) if no matching exiting conditions
my %condToRemap = ();
my $query = $bgee->prepare("SELECT DISTINCT t1.conditionId AS incorrectConditionId, t3.conditionId AS remappedConditionId,
t1.anatEntityId, t1.stageId, t1.speciesId, t1.sex, t1.sexInferred, t1.strain, t2.remappedStageId
FROM cond AS t1
INNER JOIN remapStage AS t2 ON t1.stageId = t2.incorrectStageId AND t1.speciesId = t2.speciesId
LEFT OUTER JOIN cond AS t3 ON t1.speciesId = t3.speciesId AND t1.anatEntityId = t3.anatEntityId
    AND t2.remappedStageId = t3.stageId AND t1.sex = t3.sex AND t1.sexInferred = t3.sexInferred
    AND t1.strain = t3.strain;");
$query->execute()  or die $query->errstr;
while ( my @data = $query->fetchrow_array ){
	my $condId = $data[0];
    $condToRemap{$condId}->{'remappedConditionId'} = $data[1];
    $condToRemap{$condId}->{'anatEntityId'}        = $data[2];
    $condToRemap{$condId}->{'stageId'}             = $data[3];
    $condToRemap{$condId}->{'speciesId'}           = $data[4];
    $condToRemap{$condId}->{'sex'}                 = $data[5];
    $condToRemap{$condId}->{'sexInferred'}         = $data[6];
    $condToRemap{$condId}->{'strain'}              = $data[7];
    $condToRemap{$condId}->{'remappedStageId'}     = $data[8];
    # If the sex was inferred, we'll let the function infer it again,
    # so that the correct value of sexInferred will be inserted
    if ($condToRemap{$condId}->{'sexInferred'}) {
    	$condToRemap{$condId}->{'sex'} = $Utils::NOT_ANNOTATED_SEX;
    }
}
$query->finish;


# Get already known conditions
my $conditions = Utils::query_conditions($bgee);
# Load sex-related information needed for sub 'insert_get_condition'
my $anatSexInfo    = Utils::get_anat_sex_info($sex_info);
my $speciesSexInfo = Utils::get_species_sex_info($bgee);
# Get simpler (upper level) stage equivalences
my $stage_equivalences = Utils::get_stage_equivalences($bgee);

my $insertMapping = $bgee->prepare("INSERT INTO remapCond (incorrectConditionId, remappedConditionId) VALUES (?, ?)");
for my $condId (keys %condToRemap) {
	my $mappedCondId = $condToRemap{$condId}->{'remappedConditionId'};
	# we skip conditions for which there is a simple remapping to an existing condition to do
	next if $mappedCondId;
	
    # Get conditionId/exprMappedConditionId for this library
    # Updates also the hash of existing conditions
    my $condKeyMap;
    ($condKeyMap, $conditions) = Utils::insert_get_condition($bgee,
                                                             $conditions,
                                                             $stage_equivalences,
                                                             $condToRemap{$condId}->{'anatEntityId'},
                                                             $condToRemap{$condId}->{'remappedStageId'},
                                                             $condToRemap{$condId}->{'speciesId'},
                                                             $condToRemap{$condId}->{'sex'},
                                                             $condToRemap{$condId}->{'strain'},
                                                             $anatSexInfo, $speciesSexInfo,
                                                             '', '',
                                                            );
    $condToRemap{$condId}->{'remappedConditionId'} = $condKeyMap->{'conditionId'};
    $insertMapping->execute($condId, $condToRemap{$condId}->{'remappedConditionId'})  or die $insertMapping->errstr;
}

exit 0;