#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Roux, created 29/08/08
# USAGE: perl check_affy.pl before/after normalization
#
# TO DO: before normalization, check that the expression is
# annotated on organs that exist at the given stage
#############################################################

$| = 1; # stdout not put in memory buffer

use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

require 'affy_utils.pl';
require 'bgee_utils.pl';


# Define arguments & their default value
my ($bgee_connector) = ('');
my ($normalizationType, $detectionType)  = ('', '');
my ($chipType, $microarrayExperiment, $affymetrixChip) = ('', '', '');
my ($cel_data, $processed_mas5, $processed_schuster)   = ('', '', '');
my ($affyChipInformation, $chipTypeQual) = ('', '');
my %opts = ('bgee=s'                 => \$bgee_connector,     # Bgee connector string
            'normalizationType=s'    => \$normalizationType,
            'detectionType=s'        => \$detectionType,
            'chipType=s'             => \$chipType,
            'microarrayExperiment=s' => \$microarrayExperiment,
            'cel_data=s'             => \$cel_data,
            'processed_mas5=s'       => \$processed_mas5,
            'affyChipInformation=s'  => \$affyChipInformation,
            'chipTypeQual=s'         => \$chipTypeQual,
            'affymetrixChip=s'       => \$affymetrixChip,
            'processed_schuster=s'   => \$processed_schuster,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\S(BGEECMD) -normalizationType=\$(AFFY_NORMTYPE_FILEPATH) -detectionType=\$(AFFY_DETCTYPE_FILEPATH) -chipType=\$(AFFY_CHIPTYPE_FILEPATH) -microarrayExperiment=\$(MICROARRAY_EXPERIMENT_FILEPATH) -cel_data=\$(CELPATH) -processed_mas5=\$(MAS5PATH) -affyChipInformation=\$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=\$(AFFY_CHIPTYPEQUAL_FILEPATH) -affymetrixChip=\$(AFFY_CHIP_FILEPATH) -processed_schuster=\$(SCHUSTERPATH)   --  before/after
\t-bgee                       Bgee connector string
\t-normalizationType          normalizationType                             pipeline   file
\t-detectionType              detectionType                                 pipeline   file
\t-affyChipInformation        affymetrixChipInformation                     pipeline   file
\t-chipTypeQual               chipTypeCorrespondencesAndQualityThresholds   pipeline   file
\t-chipType                   chipType                                      annotation file
\t-microarrayExperiment       microarrayExperiment                          annotation file
\t-affymetrixChip             affymetrixChip                                annotation file
\t-cel_data                   cel_data           directory
\t-processed_mas5             processed_mas5     directory
\t-processed_schuster         processed_schuster directory
\n";
    exit 1;
}

if ( $ARGV[0] !~ /^(before|after)$/ ){
    print "\n\tInvalid or missing argument:
        \te.g. $0  -bgee=\$(BGEECMD) -normalizationType=\$(AFFY_NORMTYPE_FILEPATH) -detectionType=\$(AFFY_DETCTYPE_FILEPATH) -chipType=\$(AFFY_CHIPTYPE_FILEPATH) -microarrayExperiment=\$(MICROARRAY_EXPERIMENT_FILEPATH) -cel_data=\$(CELPATH) -processed_mas5=\$(MAS5PATH) -affyChipInformation=\$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=\$(AFFY_CHIPTYPEQUAL_FILEPATH) -affymetrixChip=\$(AFFY_CHIP_FILEPATH) -processed_schuster=\$(SCHUSTERPATH)   --  before/after\n";
    exit 2;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


##############################
# Retrieving species in BGee
##############################
my $selSpecies = $dbh->prepare('SELECT speciesId FROM species');
$selSpecies->execute()  or die $selSpecies->errstr;
my %all_species;
while ( my @data = $selSpecies->fetchrow_array ){
    $all_species{$data[0]}++;
}
$selSpecies->finish;


##########################
# Read normalizationType
##########################
open(my $IN0, '<', "$normalizationType")  or die "Can't read file [$normalizationType]\n";
my %norms;
my $line = <$IN0>; #header
while ( defined ($line = <$IN0>) ){
    chomp $line;
    my @tmp = map { bgeeTrim($_) }
              split(/\t/, $line);
    $norms{$tmp[0]} = $tmp[1];
}
close $IN0;


######################
# Read detectionType
######################
open(my $IN1, '<', "$detectionType")  or die "Can't read file [$detectionType]\n";
my %detects;
$line = <$IN1>; #header
while ( defined ($line = <$IN1>) ){
    chomp $line;
    my @tmp = map { bgeeTrim($_) }
              split(/\t/, $line);
    $detects{$tmp[0]} = $tmp[1];
}
close $IN1;


##################
# Check chipType
##################
open(my $IN2, '<', "$chipType")  or die "Can't read file [$chipType]\n";
my %chip_id;
$line = <$IN2>; #header
while ( defined ($line = <$IN2>) ){
    chomp $line;
    my @tmp = map { bgeeTrim($_) }
              split(/\t/, $line);

    if ( $tmp[0] =~ /\s/ ){
        die "Problem! Space inserted in chipTypeId for [$tmp[0]]\n";
    }
    else {
        $chip_id{$tmp[0]}++;
    }

    if ( $tmp[2] =~ /\s/ && exists $all_species{$tmp[2]} ){
        die "Problem! (Space or not in Bgee?) Species Id for [$tmp[0]]\n";
    }

    if ( $tmp[3] =~ /\s/ ){
        die "Problem! Space inserted in xref for [$tmp[0]]\n";
    }
}
close $IN2;

for my $chip ( keys %chip_id ){
    die "Problem! Two chipTypeIds are identical [$chip_id{$chip} | $chip]\n" if ( $chip_id{$chip} > 1 );
}


##############################
# Check microarrayExperiment
##############################
my %experiments;
my %tsv = %{ Utils::read_spreadsheet("$microarrayExperiment", "\t", 'csv', '"', 1) };
map { warn "Problem! Space inserted for $_ Id\n"  if ( $_ =~ /\s/ );
      warn "Missing experimentId\n"              if ( $_ eq '' );
      $experiments{$_}++;
    } @{ $tsv{'experimentId'} };

#NOTE Seb: replace this block with the read_spreadsheet function & the new microarrayExperiment sheet, in column
#open(my $IN3, '<', "$microarrayExperiment")  or die "Can't read file [$microarrayExperiment]\n";
#while ( defined (my $line = <$IN3>) ){
#    chomp $line;
#
#    if ( $line =~ /^EXPERIMENT_ID\t(.+)/ ){
#        my $exp_id = bgeeTrim($1);
#        if ( $exp_id !~ /\s/ ){
#            $experiments{$exp_id}++
#        }
#        else {
#            print "Problem! Space inserted for $exp_id Id\n";
#        }
#    }
#    if ( $line =~ /^EXPERIMENT_DESCRIPTION\:/ ){
#        print "Problem! remove semi columns after EXPERIMENT_DESCRIPTION [$line]\n";
#    }
#
#    if ( $line =~ /^$/ ){
#        die "Problem! Blank line in microarrayExperiment (instead of \"//\"?)\n";
#    }
#}
#close $IN3;

map { warn "Missing EXPERIMENT_NAME\n"         if ($_ eq ''); } @{ $tsv{'EXPERIMENT_NAME'} };
map { warn "Missing EXPERIMENT_DESCRIPTION\n"  if ($_ eq ''); } @{ $tsv{'EXPERIMENT_DESCRIPTION'} };
map { warn "Missing EXPERIMENT_STATUS\n"       if ($_ eq ''); } @{ $tsv{'EXPERIMENT_STATUS'} };

my %experimentSources;
map { warn "Missing EXPERIMENT_SOURCE\n"  if ( $_ eq '' );
      $experimentSources{$_}++;
    } @{ $tsv{'EXPERIMENT_SOURCE'} };

#NOTE Seb: replace this block with the read_spreadsheet function & the new microarrayExperiment sheet, in column
#open(my $IN4, '<', "$microarrayExperiment")  or die "Can't read file (2nd times) [$microarrayExperiment]\n";
#my $lineCount       = 0;
#my $totalLineCount  = 0;
#my $error           = '';
#my %experimentSources = ();
#while ( defined ($line = <$IN4>) ){
#    chomp $line;
#    my @tmp = map { bgeeTrim($_) }
#              split(/\t/, $line);
#    # First line describes an experiment,
#    # it can be either the experiment ID, or an URL linking to the experiment
#    if ( $lineCount == 0 ){
#        $error = 'missing EXPERIMENT_ID'           if ( $tmp[0] ne 'EXPERIMENT_ID' );
#    }
#    elsif ( $lineCount == 1 ){
#        $error = 'missing EXPERIMENT_NAME'         if ( $tmp[0] ne 'EXPERIMENT_NAME' );
#    }
#    elsif ( $lineCount == 2 ){
#        $error = 'missing EXPERIMENT_DESCRIPTION'  if ( $tmp[0] ne 'EXPERIMENT_DESCRIPTION' );
#    }
#    elsif ( $lineCount == 3 ){
#        if ( $tmp[0] eq 'EXPERIMENT_SOURCE' ){
#            # everything's OK
#            # store the soure to check that it exists in the database
#            $experimentSources{$tmp[1]}++;
#        }
#        else {
#            $error = 'missing EXPERIMENT_SOURCE';
#        }
#    }
#    elsif ( $lineCount == 4 ){
#        $error = 'missing EXPERIMENT_STATUS'       if ( $tmp[0] ne 'EXPERIMENT_STATUS' );
#    }
#    elsif ( $lineCount == 5 ){
#        if ( $tmp[0] eq 'COMMENT' ){
#            # optionnal comment line, everything's OK
#            # $line-- to match "//" during next iteration
#            $line--;
#        }
#        elsif ( $line =~ /^\/\/\s*$/ ){
#            # everything's OK
#        }
#        else {
#            $error = 'missing COMMENT or end of experiment //';
#        }
#    }
#    $totalLineCount++;
#    if ( $error ne '' ){
#        print "Badly formatted microarrayExperiment file at line $totalLineCount: $error\n";
#        $error = '';
#    }
#    if ( $line =~ /^\/\/\s*$/ ){
#        $lineCount = -1;
#        $error     = '';
#    }
#
#    $lineCount++;
#}
#close $IN4;


my $selDataSrc = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName = ?');
for my $source ( keys %experimentSources ){
    $selDataSrc->execute($source)  or die $selDataSrc->errstr;
    if ( my @data = $selDataSrc->fetchrow_array ){
    }
    else {
        die "Problem! [$source] in microarrayExperiment does not exist\n";
    }
}
$selDataSrc->finish;

for my $exp ( keys %experiments ){
    die "Problem! [$exp] inserted many times\n"  if ( $experiments{$exp} > 1 );
}


# Checking experiment directories not present in the file
opendir(my $IMD0, "$cel_data")  or die("Cannot open directory [$cel_data]\n");
for my $storedFile ( readdir($IMD0) ){
    next  if ( $storedFile =~ /^\./ );

    if ( -d $cel_data.$storedFile ){
        # Check that the experiment was present in the annotation file
        warn sprintf("Warning: experiment %-17s present on the server in %s, not in the microarrayExperiment annotation file\n", "[$storedFile]", $cel_data)
            if ( !exists $experiments{$storedFile} );
    }
}
closedir($IMD0);

opendir(my $IMD1, "$processed_mas5")  or die("Cannot open directory [$processed_mas5]\n");
for my $storedFile ( readdir($IMD1) ){
    next if ( $storedFile =~ /^\./ );

    if ( -d $processed_mas5.$storedFile ){
        # Check that the experiment was present in the annotation file
        warn sprintf("Warning: experiment %-17s present on the server in %s, not in the microarrayExperiment annotation file\n", "[$storedFile]", $processed_mas5)
            if ( !exists $experiments{$storedFile} );
    }
}
closedir($IMD1);



########################################
# retrieve all organs/stages from Bgee
########################################

# $organs{anatEntityId}{'startStageId'} = startStageId
# $organs{anatEntityId}{'endStageId'}   = endStageId
my %organs;
# $stages{stageId}{'leftBound'}  = leftBound
# $stages{stageId}{'rightBound'} = rightBound
my %stages;

my $selOrgan = $dbh->prepare('SELECT anatEntityId, startStageId, endStageId FROM anatEntity');
$selOrgan->execute()  or die $selOrgan->errstr;
while ( my @data = $selOrgan->fetchrow_array ){
    $organs{$data[0]}{'startStageId'} = $data[1];
    $organs{$data[0]}{'endStageId'}   = $data[2];
}
$selOrgan->finish;

my $selStage = $dbh->prepare('SELECT stageId, stageLeftBound, stageRightBound FROM stage');
$selStage->execute()  or die $selStage->errstr;
while ( my @data = $selStage->fetchrow_array ){
  $stages{$data[0]}{'leftBound'}  = $data[1];
  $stages{$data[0]}{'rightBound'} = $data[2];
}
$selStage->finish;


########################
# Check affymetrixChip
########################
my %all_filenames;

# get chip info
my %affyChipsInfo = getAllChipsInfo($affyChipInformation);
# get the incompatible chip type
my %chipTypeInfo  = getChipTypesInformation($chipTypeQual);

open(my $IN5, '<', "$affymetrixChip")  or die "Can't read file [$affymetrixChip]\n";
my %affy_filename;
my %affy_chip;
my %affy_organ;
my %affy_stage;
my $affy_experiment;
my $start = 0;
# Store experiments including some chips incompatible or commented
my %commentedOrIncompatibleExperiments = ();
# Store experiments including some chips compatible and uncommented.
# Used with the previous var, it will allow to identify experiments that include
# ONLY chips incompatible or commented.
my %nonCommentedAndCompatibleExperiments = ();
$line = <$IN5>;    #header
while ( defined ($line = <$IN5>) ){
    chomp $line;
    # ChipID  ExperimentID  chipTypeId  normalizationTypeId  detectionTypeId  organId  organName  UberonId  UberonName  stageId  stageName  infoOrgan  infoStage  sampleTitle  sampleSource  SampleDescription  SampleCharacteristics  organAnnotationStatus  organBiologicalStatus  stageAnnotationStatus  stageBiologicalStatus  sex  strain  comment  annotatorId  lastModificationDate
    # Check if it is a low quality or incompatible chip
    my @tmp = map { bgeeTrim($_) }
              split(/\t/, $line);

    # Skip commentary line
    # but store them to remove experiments with only commented or incompatible chips
    if (( $line =~ /^#/ ) or ( $line =~ /^\"#/ )){
        $commentedOrIncompatibleExperiments{$tmp[1]}++;
        next;
    }
    if ( !defined $affyChipsInfo{$tmp[1]}{$tmp[0]} || isChipIncompatibleOrLowQuality($tmp[3], $affyChipsInfo{$tmp[1]}{$tmp[0]}, \%chipTypeInfo, $tmp[2]) ){
        $commentedOrIncompatibleExperiments{$tmp[1]}++;
        next;
    }
    # At this point, it is a chip compatible and non commented, store the experiment
    $nonCommentedAndCompatibleExperiments{$tmp[1]}++;

    my $countColumn = 0;
    for my $element ( @tmp ){
        # After the 7th column, it's commentary that can include spaces
        #   columns 5 and 6 are organ and stage ID with : in it
        if ( ($countColumn <= 4 && $element !~ /^[\w\-\.]+$/) || (($countColumn == 7 || $countColumn == 9) && $element !~ /^[\w]+:[\w]+$/) ){
            print "Problem! wrong char inserted in line: $line ==>$element<== column: $countColumn\n";
        }
        $countColumn++;
    }

    my $anatEntityId = $tmp[7];
    my $stageId      = $tmp[9];
    # Check impossible annotation: organ not existing at the annotated stage
    if ( exists $organs{$anatEntityId} && exists $stages{$stageId} ){
        my $startStageId = $organs{$anatEntityId}{'startStageId'};
        my $endStageId   = $organs{$anatEntityId}{'endStageId'};
        if ( $stages{$startStageId}{'leftBound'} <= $stages{$stageId}{'rightBound'} && $stages{$endStageId}{'rightBound'} >= $stages{$stageId}{'leftBound'} ){
            #OK
        }
        else {
            print "Warning, impossible annotation, $anatEntityId not existing at $stageId\n";
        }
    }

    if ( $start ne 0 && $tmp[1] ne $affy_experiment ){
        check_experiment($affy_experiment, \%affy_filename, \%affy_chip, \%affy_organ, \%affy_stage);

        for my $name ( keys %affy_filename ){
            $all_filenames{lc($name)}++;
        }

        $affy_experiment = $tmp[1];
        %affy_filename   = ();
        $affy_filename{$tmp[0]}->{'count'}++;
        %affy_chip       = ();
        $affy_chip{$tmp[2]}++;
        if ( $tmp[3] eq 1 && $tmp[4] eq 1 ){
            $affy_filename{$tmp[0]}->{'path'}     = $processed_mas5;
            $affy_filename{$tmp[0]}->{'out'}      = '';
        }
        if ( $tmp[3] eq 2 && $tmp[4] eq 2 ){
            $affy_filename{$tmp[0]}->{'path'}     = $processed_schuster;
            $affy_filename{$tmp[0]}->{'out'}      = '.out';
            # raw data
            $affy_filename{$tmp[0]}->{'raw_path'} = $cel_data;
        }
        %affy_organ = ();
        $affy_organ{$tmp[7]}++;
        %affy_stage = ();
        $affy_stage{$tmp[9]}++;
    }
    else {
        $affy_experiment = $tmp[1];
        $affy_filename{$tmp[0]}->{'count'}++;
        $affy_chip{$tmp[2]}++;

        if ( $tmp[3] eq 1 && $tmp[4] eq 1 ){
            $affy_filename{$tmp[0]}->{'path'}     = $processed_mas5;
            $affy_filename{$tmp[0]}->{'out'}      = '';
        }
        if ( $tmp[3] eq 2 && $tmp[4] eq 2 ){
            $affy_filename{$tmp[0]}->{'path'}     = $processed_schuster;
            $affy_filename{$tmp[0]}->{'out'}      = '.out';
            # raw data
            $affy_filename{$tmp[0]}->{'raw_path'} = $cel_data;
        }
        $affy_organ{$tmp[7]}++;
        $affy_stage{$tmp[9]}++;
    }
    $start++;
}
close $IN5;

# Last experiment
check_experiment($affy_experiment, \%affy_filename, \%affy_chip, \%affy_organ, \%affy_stage);

# Check that experiments with incompatible or commented chips had actually
# ONLY incompatible or commented chips.
for my $expCompatible ( keys %nonCommentedAndCompatibleExperiments ){
    delete $commentedOrIncompatibleExperiments{$expCompatible};
}
# Then we remove the experiences with only incompatible or commented chips
# from the list of experiments to check.
for my $expIncompatible ( keys %commentedOrIncompatibleExperiments ){
    delete $experiments{$expIncompatible};
}
if ( (keys %experiments) > 0 ){
    print "Problem! These experiments are in microarrayExperiment and not in affymetrixChip:\n";
    for my $exp ( keys %experiments ){
        print "$exp\n";
    }
}

$dbh->disconnect;
exit 0;


sub check_experiment {
    my ($exp, $files, $chips, $organs, $stages) = @_;
    print "\t$exp\n";

    if ( !exists $experiments{$exp} ){
        print "Problem! Experiment not inserted in microarrayExperiment or experiment annotated twice in affymetrixChip [$exp]\n";
    }

    delete $experiments{$exp};

    # not needed anymore, gcRMA is now used even when there is only one cel file
    #if ( scalar(keys %{$files}) eq 1 ){
    #   for my $file ( keys %{$files} ){
    #       if ( $file =~ /\.cel$/i ){
    #           print "\tBe careful, this cannot be normalized by gcRMA, mas5 should be used... and affymetrixChip should be changed to 1 1\n";
    #       }
    #   }
    #}

    for my $file ( keys %{$files} ){
        if ( length($file) >= 255 ){
            print "Problem! File name is too long!\n";
        }
        if ( $$files{$file}->{'count'} > 1 ){
            die "Problem! two or more filenames are identical in this experiment\n";
        }
        # We remove this check, cel files can actually have the same name in different experiments
        #if ( exists $all_filenames{lc($file)} ){
        #   print "Problem! two or more filenames are identical between different experiments\n";
        #}
        my $temp = 0;
        for my $file2 ( keys %{$files} ){
            if ( lc $file eq lc $file2 ){
                $temp++;
            }
        }
        if ( $temp > 1 ){
            die "Problem! two or more filenames are identical in this experiment.\n";
        }
    }

    # Check normalized data
    if ( $ARGV[0] eq 'after' ){
        for my $file ( keys %{$files} ){
            my $flag = 0;
            #print $$files{$file}->{'path'}.$exp.'/', "\n";
            opendir(my $DIRA, $$files{$file}->{'path'}.$exp.'/')  or do {print "Problem! No directory for ", $$files{$file}->{'path'}.$exp.'/', "\n"; next;};
            while ( defined(my $file_dir = readdir($DIRA)) ){
                if ( $file_dir eq $file.$$files{$file}->{'out'} ){
                    $flag = 1;
                }
                #print $$files{$file}->{'path'}.$exp.'/'.$file_dir, "\n"
                if ( -z $$files{$file}->{'path'}.$exp.'/'.$file_dir && $file_dir !~ /^\./ ){
                    print "\t$file_dir seems to be empty !\n";
                }
            }
            closedir($DIRA);

            if ( $flag eq 0 ){
                print "\tProblem! Missing data file or problem in file name: ", $file.$$files{$file}->{'out'}, "\n";
            }
        }
    }

    # Check raw data (do nothing for MAS5 data)
    if ( $ARGV[0] eq 'before' ){
        for my $file ( keys %{$files} ){
            if ( defined $$files{$file}->{'raw_path'} ){
                my $flag = 0;
                #print $$files{$file}->{'raw_path'}.$exp.'/', "\n";
                opendir(my $DIRB, $$files{$file}->{'raw_path'}.$exp.'/')  or print "Problem! No directory for data\n";
                while ( defined(my $file_dir = readdir($DIRB)) ){
                    if ( $file_dir eq $file || $file_dir eq $file.'.gz' ){
                        $flag = 1;
                    }
                }
                closedir($DIRB);

                if ( $flag eq 0 ){
                    print $file, "\n";
                    print "Problem! Missing raw data file or problem in file name\n";
                }
            }
        }
    }

    for my $chip ( keys %{$chips} ){
        print "Problem! A new Id has to be inserted in chipType ($chip)\n"  if ( !exists $chip_id{$chip} );
    }

#    for my $organ ( keys %{$organs} ){
#        print "Problem! An organ annotated is not in Bgee ($organ)\n"       if ( !exists $organs{$organ} );
#    }
#
#    for my $stage ( keys %{$stages} ){
#        print "Problem! A stage annotated is not in Bgee ($stage)\n"        if ( !exists $stages{$stage} );
#    }

    return;
}

