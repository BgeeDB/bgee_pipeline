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
# This script takes notably as entry a TSV file used for remapping,
# with the following columns:
# Condition ID to remap - New anat. entity ID - New dev. stage ID - New sex - New strain.
# The species ID is assumed to stay the same and will not be changed.
# This script will not update directly the conditions, rather, it will populate a table
# in the database, named 'remapCond', storing the condition ID to remap,
# and the ID of the condition it is remapped to. It is then manually that you should
# update/delete the conditions. This allows more flexibility and controls.
#
# Frederic Bastian, created August 2019
#####################################################################

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($remapping_file, $sex_info)   = ('', '');
my ($debug)                      = (0);
my %opts = ('bgee=s'                => \$bgee_connector,     # Bgee connector string
            'remapping_file=s'      => \$remapping_file,     # source_files/annotations/condition_remapping.tsv
            'sex_info=s'            => \$sex_info,           # generated_files/uberon/uberon_sex_info.tsv
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $sex_info eq ''  || $remapping_file eq '' ){
    print "\n\tInvalid or missing argument:
\te.g., $0  -bgee=\$(BGEECMD) -sex_info=\$(UBERON_SEX_INFO_FILE_PATH) > $@.tmp 2>warnings.$@
\t-bgee                Bgee connector string
\t-remapping_file      file containing the conditions to be remapped
\t-sex_info            file containing sex-related info about anatomical terms
\n";
    exit 1;
}


# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

##############################
# CONDITIONS TO UPDATE/REMAP #
##############################

# Get already known conditions
my $conditions = Utils::query_conditions($bgee);
# Load sex-related information needed for sub 'insert_get_condition'
my $anatSexInfo    = Utils::get_anat_sex_info($sex_info);
my $speciesSexInfo = Utils::get_species_sex_info($bgee);
# Get simpler (upper level) stage equivalences
my $stage_equivalences = Utils::get_stage_equivalences($bgee);

my $insertMapping = $bgee->prepare("INSERT INTO remapCond (incorrectConditionId, remappedConditionId) VALUES (?, ?)");
# TODO use column names rather than column order
REMAP: for my $line ( grep { !/^#/ && !/^Condition ID/} read_file("$remapping_file", chomp => 1) ){
    #Condition ID to remap  New anat. entity ID New dev. stage ID   New sex New strain
    #40558  CL:0000094  UBERON:0000066  male    NA
    my @tmp = split(/\t/, $line);

    my $condIdToRemap = $tmp[0];
    my $speciesId = 0;
    my $sexInferred = 0;
    # Get the species ID and sex inferrence of the existing condition to remap
    COND: for my $condHashKey (keys %{ $conditions }) {
        if ($condIdToRemap == $conditions{$condHashKey}->{'conditionId'}) {
            $speciesId   = $conditions{$condHashKey}->{'speciesId'};
            $sexInferred = $conditions{$condHashKey}->{'sexInference'};
            last COND;
        }
    }
    if (!$speciesId) {
        print "Condition with ID $condIdToRemap could not be found and is skipped\n";
        next REMAP;
    }
    # If the sex was inferred, we'll let the function infer it again,
    # so that the correct value of sexInferred will be inserted
    my $sex = $tmp[3];
    if ($sexInferred) {
        $sex = $Utils::NOT_ANNOTATED_SEX;
    }


    # Try to get or insert the new mapped condition
    my $condKeyMap;
    ($condKeyMap, $conditions) = Utils::insert_get_condition($bgee,
                                                             $conditions,
                                                             $stage_equivalences,
                                                             $tmp[1],
                                                             $tmp[2],
                                                             $speciesId,
                                                             $sex,
                                                             $tmp[4],
                                                             $anatSexInfo, $speciesSexInfo,
                                                             '', '',
                                                            );
    # If the sex was originally inferred, check that it is still the case
    if ($sexInferred && !$condKeyMap->{'sexInference'}) {
        die "Incorrect new sex inference for new cond ID ".$condKeyMap->{'conditionId'}.
            " - old condition ID: ".$condIdToRemap;
    }
    my $newCondId = $condKeyMap->{'conditionId'};

    $insertMapping->execute($condIdToRemap, $newCondId)  or die $insertMapping->errstr;
    print "Old condition ID $condIdToRemap to remap to condition ID $newCondId\n";
}

exit 0;