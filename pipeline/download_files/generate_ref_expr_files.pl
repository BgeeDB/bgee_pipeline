#!/usr/bin/env perl

# Frederic Bastian, May 2015
# Julien Wollbrett, April 2024

# This script is reponsible for generating the download files
# containing our processed expression data (not the calls,
# but the CPM, TPM or FPKM values, read counts, log signal intensities, etc).
# See opt debug message below for information about parameters.
#
# Since Bgee 15.1 all RNASeq datatypes (bulk, full length and droplet based)
# are stored in the same tables in the database. The code has been refactored
# to be able to use the same function for each of these datatypes

#############################################################

use strict;
use warnings;
use diagnostics;

use Archive::Tar;
use Getopt::Long;
use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Copy qw(move);
use File::Basename qw(basename);
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;


$| = 1;

# Define arguments & their default value
my ($bgee_connector)                = ('');
my ($speciesArg)                    = ('');
my ($affyDir)                       = ('');
my ($estDir)                        = ('');
my ($inSituDir)                     = ('');
my ($rnaSeqDir)                     = ('');
my ($flScRnaSeqDir)                 = ('');
my ($dbScRnaSeqDir)                 = ('');
my ($bgeeVersion)                   = ('');
my ($debug)                         = (0);
my %opts = ('bgee=s'                => \$bgee_connector,     # Bgee connector string
            'speciesArg=s'         => \$speciesArg,
            'affyDir=s'            => \$affyDir,
            'estDir=s'             => \$estDir,
            'inSituDir=s'          => \$inSituDir,
            'rnaSeqDir=s'          => \$rnaSeqDir,
            'flScRnaSeqDir=s'      => \$flScRnaSeqDir,
            'dbScRnaSeqDir=s'      => \$dbScRnaSeqDir,
            'bgeeVersion=s'        => \$bgeeVersion,
            'debug'                => \$debug
           );

# Check arguments
my $emptyArg = '-';
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $speciesArg eq '' ||
     $affyDir eq '' || $estDir eq '' || $inSituDir eq '' || $rnaSeqDir eq '' || $flScRnaSeqDir eq '' || $dbScRnaSeqDir eq '' || $bgeeVersion eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEE_CMD) -speciesArg=\$(SPECIES_ARG) -affyDir=\$(AFFY_DIR) -estDir=\$(EST_DIR) -inSituDir=\$(IN_SITU_DIR) -rnaSeqDir=\$(RNA_SEQ_DIR) -flScRnaSeqDir=\$(SC_RNA_SEQ_FL_DIR)
\t-bgee                 Bgee    connector string
\t-speciesArg           A comma-separated list of species IDs to generate files for, or '-' to generate files for all species
\t-affyDir              The path to generate Affymetrix related files, or '-' for not generating data for Affymetrix
\t-estDir               The path to generate EST related files, or '-' for not generating data for EST
\t-inSituDir            The path to generate in situ hybridization related files, or '-' for not generating data for in situ hybridization
\t-rnaSeqDir            The path to generate RNA-Seq related files, or '-' for not generating data for RNA-Seq
\t flScRnaSeqDir        The path to generate full length single cell RNA-Seq related files, or '-' for not generating data for full length single cell RNA-Seq
\t dbScRnaSeqDir        The path to generate droplet based single cell RNA-Seq related files, or '-' for not generating data for droplet based single cell RNA-Seq
\t-bgeeVersion          The Bgee release for which the files are being generated, to generate FTP links to correct version, e.g. 'bgee_v15'.
\t-debug                more verbose output
\n";
    exit 1;
}

#############################################################
# Retrieve information needed for all species / all data types


# -----------------------
# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);
# species list
my @speciesList = ();
if ($speciesArg ne $emptyArg) {
    @speciesList = split(',', $speciesArg);
}

# -----------------------
# retrieve species information (notably to generate file names)
# $species{speciesId}{'name'} = genus . species
my %species;

my $speciesSql = 'SELECT speciesId, genus, species FROM species';
if (@speciesList) {
    $speciesSql .= ' WHERE speciesId IN (';
    for my $i (0 .. $#speciesList) {
        if ($i > 0) {
            $speciesSql .= ', ';
        }
        $speciesSql .= $speciesList[$i];
    }
    $speciesSql .= ')';
}
my $speciesStmt = $dbh->prepare($speciesSql);
$speciesStmt->execute()  or die $speciesStmt->errstr;
while ( my @data = $speciesStmt->fetchrow_array ){
    $species{$data[0]}{'name'} = $data[1].' '.$data[2];
}


# -----------------------
# Retrieve information about data sources
# $dataSources{dataSourceId}{'name'} = data source name
# $dataSources{dataSourceId}{'experimentBaseUrl'} = experiment base URL (replace '[experiment_id]' with exp ID)
# $dataSources{dataSourceId}{'evidenceBaseUrl'} = evidence base URL (replace '[evidence_id]' with evidence ID, such as chip ID or library ID)
my %dataSources;
my $expIdTag      = quotemeta '[experiment_id]';
my $evidenceIdTag = quotemeta '[evidence_id]';

my $sourceStmt = $dbh->prepare('SELECT dataSourceId, dataSourceName, experimentUrl, evidenceUrl FROM dataSource');
$sourceStmt->execute()  or die $sourceStmt->errstr;
while ( my @data = $sourceStmt->fetchrow_array ){
    $dataSources{$data[0]}{'name'}              = $data[1];
    $dataSources{$data[0]}{'experimentBaseUrl'} = $data[2];
    $dataSources{$data[0]}{'evidenceBaseUrl'}   = $data[3];
}


#############################################################
# Now, for each species, generate the files

# Link to FTP storing our files
my $ftpFilePath = 'https://www.bgee.org/ftp/';

for my $speId ( keys %species ){
    print "Generating files for species $speId...\n";
    # -----------------------
    # Affy data
    if ($affyDir ne $emptyArg) {
        print "Generating Affymetrix files...\n";
        my $absDir = $affyDir;
        if ( !File::Spec->file_name_is_absolute($affyDir) ){
            $absDir = File::Spec->rel2abs($affyDir) ;
        }
        generateAffyFiles($speId, $species{$speId}{'name'}, \%dataSources, $absDir);
        print "Done generating Affymetrix files...\n";
    }
    # -----------------------
    # bulk RNA-Seq data
    if ($rnaSeqDir ne $emptyArg) {
        print "Generating RNA-Seq files...\n";
        my $isSingleCell = 0;
        my $isSampleMultiplexing = 0;
        my $absDir = $rnaSeqDir;
        if ( !File::Spec->file_name_is_absolute($rnaSeqDir) ){
            $absDir = File::Spec->rel2abs($rnaSeqDir) ;
        }
        generateRnaSeqFiles($speId, $species{$speId}{'name'}, \%dataSources, $absDir, $isSingleCell, $isSampleMultiplexing);
        print "Done generating RNA-Seq files...\n";
    }
    # -----------------------
    # Full Length single cell RNA-Seq data
    if ($flScRnaSeqDir ne $emptyArg) {
        print "Generating Full Length single cell RNA-Seq files...\n";
        my $isSingleCell = 1;
        my $isSampleMultiplexing = 0;
        my $absDir = $flScRnaSeqDir;
        if ( !File::Spec->file_name_is_absolute($flScRnaSeqDir) ){
            $absDir = File::Spec->rel2abs($flScRnaSeqDir) ;
        }
        generateRnaSeqFiles($speId, $species{$speId}{'name'}, \%dataSources, $absDir, $isSingleCell, $isSampleMultiplexing);
        print "Done generating Full Length single cell RNA-Seq files...\n";
    }
    # -----------------------
    # Droplet based single cell RNA-Seq data
    if ($flScRnaSeqDir ne $emptyArg) {
        print "Generating Full Length single cell RNA-Seq files...\n";
        my $isSingleCell = 1;
        my $isSampleMultiplexing = 1;
        my $absDir = $dbScRnaSeqDir;
        if ( !File::Spec->file_name_is_absolute($dbScRnaSeqDir) ){
            $absDir = File::Spec->rel2abs($dbScRnaSeqDir) ;
        }
        generateRnaSeqFiles($speId, $species{$speId}{'name'}, \%dataSources, $absDir, $isSingleCell, $isSampleMultiplexing);
        print "Done generating Droplet based single cell RNA-Seq files...\n";
    }
    print "Done generating files for species $speId.\n";
}

#############################################################
# SUBS

$dbh->disconnect;
exit 0;

sub generateAffyFiles {
    my @args = @_;
    my $speciesId      = $args[0];
    my $speciesName    = $args[1];
    my $dataSourcesRef = $args[2];
    my $filesDir       = $args[3];

    my $speciesNameForFile = $speciesName;
    $speciesNameForFile =~ s/ /_/g;
    my $tmpDirName = $speciesNameForFile.'_tmp';
    my $speciesDirName = $speciesNameForFile;
    my $expFileName = $speciesNameForFile.'_Affymetrix_experiments.tsv';
    my $chipFileName = $speciesNameForFile.'_Affymetrix_chips.tsv';
	# recreate temp directory for experiment and chip information
	remove_tree(File::Spec->catdir( $filesDir, $tmpDirName));
    make_path(File::Spec->catdir( $filesDir, $tmpDirName));
    # path to cel files directory relative to FTP dir
    my $celFilePath = 'affymetrix_data/cel_files/';
    # path to MAS5 files directory relative to FTP dir
    my $mas5FilePath = 'affymetrix_data/mas5_files/';
    # path to processed expression data on FTP
    my $exprFilePath = $bgeeVersion.'/download/processed_expr_values/affymetrix/'
        .$speciesNameForFile.'/';

    # -------------------------------------------------
    # First, we retrieve information about experiments
    # $exp{'expId'}           = exp ID
    # $exp{'name'}            = exp name
    # $exp{'desc'}            = exp description
    # $exp{'sourceId'}        = data source ID
    # $exp{'chipCount'}       = chip count
    # $exp{'conditionCount'}  = condition count
    # $exp{'anatEntityCount'} = anatomical entity count
    # $exp{'stageCount'}      = stage count
    # $exp{'hasMAS5Chip'}     = 1 if the experiment includes MAS5 files
    # $exp{'hasCELChip'}      = 1 if the experiment includes CEL files
    my @experiments = ();

    my $sqlExpPart =
        'SELECT t1.microarrayExperimentId, t1.microarrayExperimentName, t1.microarrayExperimentDescription, '
              .'t1.dataSourceId, COUNT(DISTINCT bgeeAffymetrixChipId) AS chipCount, '
              # to count the number of conditions, we need to used the expression-mapped condition IDs, 
              # which merge non-informative annotations such as 'not annotated' and 'NA'
              .'COUNT(DISTINCT exprMappedConditionId) AS conditionCount, '
              .'COUNT(DISTINCT anatEntityId, stageId) AS anatEntityStageCount, '
              .'COUNT(DISTINCT anatEntityId) AS anatEntityCount, '
              .'COUNT(DISTINCT stageId) AS stageCount, '
              .'(SELECT COUNT(DISTINCT t10.sex) FROM cond AS t10 '
                  .'INNER JOIN affymetrixChip AS t11 ON t10.conditionId = t11.conditionId '
                  .'WHERE t10.sex NOT IN ('.join(', ', map { "'$_'" } @Utils::NO_ANNOT_SEX_INFO).') '
                  .'AND t11.microarrayExperimentId = t1.microarrayExperimentId) AS sexCount, '
              .'(SELECT COUNT(DISTINCT t10.strain) FROM cond AS t10 '
                  .'INNER JOIN affymetrixChip AS t11 ON t10.conditionId = t11.conditionId '
                  .'WHERE t10.strain NOT IN ('.join(', ', map { "'$_'" } @Utils::NO_ANNOT_STRAIN_INFO).') '
                  .'AND t11.microarrayExperimentId = t1.microarrayExperimentId) AS strainCount, '
              .'(SELECT 1 FROM affymetrixChip AS t3 '
                  .'WHERE t3.microarrayExperimentId = t1.microarrayExperimentId AND '
                  .'t3.normalizationType = "MAS5" LIMIT 1) AS hasMAS5Chip, '
              .'(SELECT 1 FROM affymetrixChip AS t3 '
                  .'WHERE t3.microarrayExperimentId = t1.microarrayExperimentId AND '
                  .'t3.normalizationType != "MAS5" LIMIT 1) AS hasCELChip '
              .'FROM microarrayExperiment AS t1 '
              .'INNER JOIN affymetrixChip AS t2 ON t1.microarrayExperimentId = t2.microarrayExperimentId '
              .'INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId '
              .'WHERE speciesId = ? '
              .'GROUP BY t1.microarrayExperimentId ';
    my $sql = $sqlExpPart
              .'ORDER BY chipCount DESC, conditionCount DESC, anatEntityStageCount DESC, '
              .'anatEntityCount DESC, stageCount DESC, sexCount DESC, strainCount DESC, '
              .'t1.microarrayExperimentId DESC';

    my $stmt = $dbh->prepare($sql);
    $stmt->execute($speciesId)  or die $stmt->errstr;

    while ( my @data = $stmt->fetchrow_array ){
        my %exp;
        $exp{'expId'}                = $data[0];
        $exp{'name'}                 = $data[1];
        $exp{'desc'}                 = $data[2];
        $exp{'sourceId'}             = $data[3];
        $exp{'chipCount'}            = $data[4];
        $exp{'conditionCount'}       = $data[5];
        $exp{'anatEntityStageCount'} = $data[6];
        $exp{'anatEntityCount'}      = $data[7];
        $exp{'stageCount'}           = $data[8];
        $exp{'sexCount'}             = $data[9];
        $exp{'strainCount'}          = $data[10];
        $exp{'hasMAS5Chip'} = 0;
        if (defined $data[11] && $data[11] == 1) {
            $exp{'hasMAS5Chip'} = 1;
        }
        $exp{'hasCELChip'} = 0;
        if (defined $data[12] && $data[12] == 1) {
            $exp{'hasCELChip'} = 1;
        }

        push @experiments, \%exp;
    }
    
    if( @experiments > 0 ){
	    # recreate directory for experiment and chip information
    	remove_tree(File::Spec->catdir( $filesDir, $speciesDirName));
    	make_path(File::Spec->catdir( $filesDir, $speciesDirName));
    }

    # Print experiment information into file
    my $expFile = File::Spec->catfile($filesDir, $tmpDirName, $expFileName);
    open(my $fh, '>', $expFile) or die "Could not open file '$expFile' $!";
    print $fh "Experiment ID\tExperiment name\t"
              ."Chip count\tCondition count\tOrgan-stage count\t"
              ."Organ count\tStage count\tSex count\tStrain count\t"
              ."Data source\tData source URL\tBgee normalized data URL\tBgee raw files URL\t"
              ."Experiment description\n";
    for my $exp ( @experiments ){
        print $fh $exp->{'expId'}."\t";
        
        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        my $name = $exp->{'name'};
        $name =~ s/"/'/g;
        print $fh '"'.$name.'"'."\t";
        
        print $fh $exp->{'chipCount'}."\t"
                  .$exp->{'conditionCount'}."\t"
                  .$exp->{'anatEntityStageCount'}."\t"
                  .$exp->{'anatEntityCount'}."\t"
                  .$exp->{'stageCount'}."\t"
                  .$exp->{'sexCount'}."\t"
                  .$exp->{'strainCount'}."\t";
        my $sourceUrl  = 'NA';
        my $sourceName = '';
        if (defined $dataSourcesRef->{$exp->{'sourceId'}}) {
            if ($dataSourcesRef->{$exp->{'sourceId'}}->{'experimentBaseUrl'}) {
                $sourceUrl = $dataSourcesRef->{$exp->{'sourceId'}}->{'experimentBaseUrl'};
                $sourceUrl =~ s/$expIdTag/$exp->{'expId'}/g;
            }
            $sourceName = $dataSourcesRef->{$exp->{'sourceId'}}->{'name'};
        }
        print $fh $sourceName."\t".$sourceUrl."\t";

        my $probesetFileName = $speciesNameForFile.'_Affymetrix_probesets_'.$exp->{'expId'}.'.tar.gz';
        $probesetFileName =~ s/ /_/g;
        print $fh $ftpFilePath.$exprFilePath.$probesetFileName."\t";

        if ($exp->{'hasCELChip'}) {
            print $fh $ftpFilePath.$celFilePath.$exp->{'expId'}.'/';
        }
        if ($exp->{'hasMAS5Chip'}) {
            if ($exp->{'hasCELChip'}) {
                print $fh ' ';
            }
            print $fh $ftpFilePath.$mas5FilePath.$exp->{'expId'}.'/';
        }
        print $fh "\t";

        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        my $desc = $exp->{'desc'};
        $desc =~ s/"/'/g;
        print $fh '"'.$desc.'"';

        print $fh "\n";

    }
    close $fh;

    # -------------------------------------------------
    # Now, information about chips.
    # $chip{'expId'}                   = experiment ID
    # $chip{'chipId'}                  = chip ID
    # $chip{'chipTypeId'}              = chip type ID
    # $chip{'chipTypeName'}            = chip type name
    # $chip{'chipTypeCdfName'}         = CDF name of chip type
    # $chip{'qualityScoreThreshold'}   = IQRray score threshold for the chip type
    # $chip{'percentPresentThreshold'} = percent present score threshold for the chip type
    # $chip{'scanDate'}                = scan date
    # $chip{'normalizationType'}       = normalization type
    # $lib{'anatEntityId'}             = anat. entity ID annotated for this chip
    # $lib{'anatEntityName'}           = anat. entity name annotated for this chip
    # $lib{'stageId'}                  = stage ID annotated for this chip
    # $lib{'stageName'}                = stage name annotated for this chip
    # $lib{'sex'}                      = annotated sex info (not mapped for expression table)
    # $lib{'strain'}                   = annotated strain info (not mapped for expression table)
    # $chip{'qualityScore'}            = quality score
    # $chip{'percentPresent'}          = percent present
    # $chip{'sourceId'}                = data source ID
    my @chips = ();

    $sql = 'SELECT t1.microarrayExperimentId, t1.affymetrixChipId, t1.chipTypeId, '
              .'t2.chipTypeName, t2.cdfName, t2.qualityScoreThreshold, t2.percentPresentThreshold, '
              .'t1.scanDate, t1.normalizationType, t1.qualityScore, t1.percentPresent, '
              .'t6.anatEntityId, t3.anatEntityName, t6.stageId, t4.stageName, t6.sex, t6.strain, '
              .'t5.dataSourceId '
              .'FROM affymetrixChip AS t1 '
              .'INNER JOIN cond AS t6 ON t1.conditionId = t6.conditionId '
              .'INNER JOIN chipType AS t2 ON t1.chipTypeId = t2.chipTypeId '
              .'INNER JOIN anatEntity AS t3 ON t6.anatEntityId = t3.anatEntityId '
              .'INNER JOIN stage AS t4 ON t6.stageId = t4.stageId '
              .'INNER JOIN ('.$sqlExpPart.') AS t5 ON t1.microarrayExperimentId = t5.microarrayExperimentId '
              .'WHERE t6.speciesId = ? '
              .'ORDER BY chipCount DESC, conditionCount DESC, anatEntityStageCount DESC, '
              .'anatEntityCount DESC, stageCount DESC, sexCount DESC, strainCount DESC, '
              .'t1.microarrayExperimentId, t1.chipTypeId, t6.anatEntityId, t6.stageId, t6.sex, t6.strain, '
              .'t1.affymetrixChipId';
    print $sql."\n";
    $stmt = $dbh->prepare($sql);
    $stmt->execute($speciesId, $speciesId)  or die $stmt->errstr;

    while ( my @data = $stmt->fetchrow_array ){
        my %chip;
        $chip{'expId'}                   = $data[0];
        $chip{'chipId'}                  = $data[1];
        $chip{'chipTypeId'}              = $data[2];
        $chip{'chipTypeName'}            = $data[3];
        $chip{'cdfName'}                 = $data[4];
        $chip{'qualityScoreThreshold'}   = $data[5];
        $chip{'percentPresentThreshold'} = $data[6];
        $chip{'scanDate'}                = $data[7];
        $chip{'normalizationType'}       = $data[8];
        $chip{'qualityScore'}            = $data[9];
        $chip{'percentPresent'}          = $data[10];
        $chip{'anatEntityId'}            = $data[11];
        $chip{'anatEntityName'}          = $data[12];
        $chip{'stageId'}                 = $data[13];
        $chip{'stageName'}               = $data[14];
        $chip{'sex'}                     = $data[15];
        $chip{'strain'}                  = $data[16];
        $chip{'sourceId'}                = $data[17];
        push @chips, \%chip;
    }

    # Print chip information into file
    my $chipFile = File::Spec->catfile($filesDir, $tmpDirName, $chipFileName);
    open($fh, '>', $chipFile) or die "Could not open file '$chipFile' $!";
    print $fh "Experiment ID\tChip ID\tAnatomical entity ID\tAnatomical entity name\t"
              ."Stage ID\tStage name\tSex\tStrain\tIQRray score\tMAS5 percent present\t"
              ."Normalization type\tScan date\tChip type ID\tCDF name\tChip type name\t"
              ."IQRray score threshold for the chip type\t"
              ."MAS5 percent present threshold for the chip type\t"
              ."Data source\tData source URL\tBgee normalized data URL\tBgee normalized data file\t"
              ."Bgee raw file URL\n";
    for my $chip ( @chips ){
        print $fh $chip->{'expId'}."\t"
                  .$chip->{'chipId'}."\t"
                  .$chip->{'anatEntityId'}."\t";
        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        my $toPrint = $chip->{'anatEntityName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        print $fh $chip->{'stageId'}."\t";
        $toPrint = $chip->{'stageName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";
        
        print $fh $chip->{'sex'}."\t";
        
        $toPrint = $chip->{'strain'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        if ($chip->{'qualityScore'} > 0) {
            print $fh $chip->{'qualityScore'};
        } else {
            print $fh 'NA';
        }
        print $fh "\t";
        if ($chip->{'percentPresent'} > 0) {
            print $fh $chip->{'percentPresent'};
        } else {
            print $fh 'NA';
        }
        print $fh "\t";

        print $fh $chip->{'normalizationType'}."\t";

        if (defined $chip->{'scanDate'} && $chip->{'scanDate'}) {
            print $fh $chip->{'scanDate'};
        } else {
            print $fh 'NA';
        }
        print $fh "\t";

        print $fh $chip->{'chipTypeId'}."\t";

        $toPrint = $chip->{'cdfName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";
        $toPrint = $chip->{'chipTypeName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        if ($chip->{'qualityScoreThreshold'} > 0) {
            print $fh $chip->{'qualityScoreThreshold'};
        } else {
            print $fh 'NA';
        }
        print $fh "\t";
        if ($chip->{'percentPresentThreshold'} > 0) {
            print $fh $chip->{'percentPresentThreshold'};
        } else {
            print $fh 'NA';
        }
        print $fh "\t";

        my $sourceUrl  = 'NA';
        my $sourceName = '';
        if (defined $dataSourcesRef->{$chip->{'sourceId'}}) {
            if ($dataSourcesRef->{$chip->{'sourceId'}}->{'evidenceBaseUrl'}) {
                $sourceUrl = $dataSourcesRef->{$chip->{'sourceId'}}->{'evidenceBaseUrl'};
                $sourceUrl =~ s/$expIdTag/$chip->{'expId'}/g;
                $sourceUrl =~ s/$evidenceIdTag/$chip->{'chipId'}/g;
            }
            $sourceName = $dataSourcesRef->{$chip->{'sourceId'}}->{'name'};
        }
        print $fh $sourceName."\t".$sourceUrl."\t";

        my $expFileName = $speciesNameForFile.'_Affymetrix_probesets_'.$chip->{'expId'}.'.tar.gz';
        $expFileName =~ s/ /_/g;
        print $fh $ftpFilePath.$exprFilePath.$expFileName."\t";

        my $probesetFileName = $speciesNameForFile.'_probesets_'
                    .$chip->{'expId'}.'_'.$chip->{'chipTypeId'}.'_'.$chip->{'normalizationType'}.'.tsv';
        $probesetFileName =~ s/ /_/g;
        print $fh $probesetFileName."\t";

        my $downloadUrl = $ftpFilePath;
        if ($chip->{'normalizationType'} ne 'MAS5') {
            $downloadUrl .= $celFilePath;
        } else {
            $downloadUrl .= $mas5FilePath;
        }
        $downloadUrl .= $chip->{'expId'}.'/'.$chip->{'chipId'};
        if ($chip->{'normalizationType'} ne 'MAS5') {
            $downloadUrl .= '.CEL.gz';
        }
        print $fh $downloadUrl."\n";
    }
    close $fh;

    # -------------------------------------------------
    # Now, we generate the probeset files,
    # one file, per experiment/chip type/normalization method.
    # First, we get the list of exp/chip type/normalization
    # $fileParams{expId}{chipTypeId}{normalization} = 1;
    my %fileParams;
    for my $chip ( @chips ){
        $fileParams{$chip->{'expId'}}{$chip->{'chipTypeId'}}{$chip->{'normalizationType'}} = 1;
    }

    # Now, generate the files. We'll store the path to all tar.gz files generated,
    # to make one giant tar.gz at the end.
    my @tarFileNames = ();
    for my $expId ( keys %fileParams ){
        # we will store all names of files generated for an experiment, to pack all of them together
        my @probesetFileNames = ();

        for my $chipTypeId ( sort keys %{$fileParams{$expId}} ){
            for my $normalizationType ( sort keys %{$fileParams{$expId}{$chipTypeId}} ){
                # Create file
                my $fileName = $speciesNameForFile.'_probesets_'
                    .$expId.'_'.$chipTypeId.'_'.$normalizationType.'.tsv';
                $fileName =~ s/ /_/g;
                push @probesetFileNames, $fileName;
                my $probesetFile = File::Spec->catfile($filesDir, $tmpDirName,
                    $fileName);
                open(my $fh, '>', $probesetFile) or die "Could not open file '$probesetFile' $!";
                print $fh "Experiment ID\tChip ID\tProbeset ID\tGene ID\t"
                          ."Anatomical entity ID\tAnatomical entity name\t"
                          ."Stage ID\tStage name\tSex\tStrain\tLog of normalized signal intensity\t"
                          ."Detection flag\tpValue\tState in Bgee\n";

                # Retrieve data from database.
                $sql = 'SELECT STRAIGHT_JOIN t1.microarrayExperimentId, t1.affymetrixChipId, '
                       .'t2.affymetrixProbesetId, t3.geneId, '
                       .'t4.anatEntityId, t5.anatEntityName, t4.stageId, t6.stageName, t4.sex, t4.strain, '
                       .'t2.normalizedSignalIntensity, t2.rawDetectionFlag, t2.pValue, '
                       .'t2.expressionId, t2.reasonForExclusion, '
                        # FIXME retrieve call type
                       .'IF(t2.expressionId IS NOT NULL, "data", "no data") AS globalAffymetrixData, '
                       .'t1.normalizationType '
                       .'FROM affymetrixChip AS t1 '
                       .'INNER JOIN cond AS t4 ON t1.conditionId = t4.conditionId '
                       .'INNER JOIN affymetrixProbeset AS t2 '
                           .'ON t1.bgeeAffymetrixChipId = t2.bgeeAffymetrixChipId '
                       .'INNER JOIN gene AS t3 ON t2.bgeeGeneId = t3.bgeeGeneId '
                       .'INNER JOIN anatEntity AS t5 ON t5.anatEntityId = t4.anatEntityId '
                       .'INNER JOIN stage AS t6 ON t6.stageId = t4.stageId '
                       .'LEFT OUTER JOIN expression AS t7 ON t2.expressionId = t7.expressionId '
                       .'WHERE t4.speciesId = ? AND t1.microarrayExperimentId = ? '
                           .'AND t1.chipTypeId = ? AND t1.normalizationType = ? '
                      .'ORDER BY t4.anatEntityId, t4.stageId, t4.sex, t4.strain, '
                          .'t1.affymetrixChipId, t3.geneId';

                $stmt = $dbh->prepare($sql);
                $stmt->execute($speciesId, $expId, $chipTypeId, $normalizationType)  or die $stmt->errstr;

                while ( my @data = $stmt->fetchrow_array ){
                    # we write the data directly to not store them in memory
                    print $fh $data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3]."\t"
                              .$data[4]."\t";
                    # we replace double quotes with simple quotes, and we surround with double quotes
                    # the values to escape potential special characters
                    my $toPrint = $data[5];
                    $toPrint =~ s/"/'/g;
                    print $fh '"'.$toPrint.'"'."\t";

                    print $fh $data[6]."\t";
                    $toPrint = $data[7];
                    $toPrint =~ s/"/'/g;
                    print $fh '"'.$toPrint.'"'."\t";

                    print $fh $data[8]."\t";
                    $toPrint = $data[9];
                    $toPrint =~ s/"/'/g;
                    print $fh '"'.$toPrint.'"'."\t";

                    # for signal intensities, if the normalization used was MAS5,
                    # the data are not log, while data normalized with gcRMA are.
                    if ($data[16] eq 'MAS5' && $data[10] > 0) {
                        print $fh sprintf("%.2f", log($data[10])/log(2));
                    } else {
                        print $fh $data[10];
                    }
                    print $fh "\t";

                    print $fh $data[11]."\t".$data[12]."\t";

                    if ($data[14] eq $Utils::CALL_NOT_EXCLUDED) {
                        print $fh 'Part of a call';
                        # TODO manage according to retrieved call type
                        #if (defined $data[12] && $data[12]) {
                        #    print $fh 'present ';
                        #} else {
                        #    print $fh 'absent ';
                        #}
                        #print $fh $data[15];
                    } else {
                        print $fh 'Result excluded, reason: '.$data[14];
                    }
                    print $fh "\n";
                }
                close $fh;
            }
        }
        # we compress all files for an experiment together, and delete the uncompressed files
        # to save disk space
        
        # We close the connection before compressing the files, it can take quite some time
        $dbh->disconnect;
        
        my $fileName = $speciesNameForFile.'_Affymetrix_probesets_'.$expId.'.tar.gz';
        $fileName =~ s/ /_/g;
        # remove old archive
        unlink File::Spec->catfile($filesDir, $speciesDirName, $fileName);
        # create new archive 
        my $experimentTar = Archive::Tar->new();
        my @probesetFilePaths = ();
        for my $probesetFileName ( @probesetFileNames ) {
        	push @probesetFilePaths, File::Spec->catfile($filesDir, $tmpDirName, $probesetFileName);
        }
        $experimentTar->add_files(@probesetFilePaths);
        my @file_objs = $experimentTar->get_files;
        for my $file_objs (@file_objs){
        	$experimentTar->rename( $file_objs, basename($file_objs->full_path)) 
        }
        #write the tar file as a compressed tar.gz file
        $experimentTar->write(File::Spec->catfile($filesDir, $speciesDirName, $fileName), COMPRESS_GZIP);
        
        # Remove uncompressed files to store disk space
        for my $probesetFileName ( @probesetFileNames ) {
            unlink File::Spec->catfile($filesDir, $tmpDirName, $probesetFileName);
        }
        # Store the name of this experiment tar.gz file, to make one giant tar.gz of all experiments
        # at the end.
        push @tarFileNames, $fileName;
        
        # reopen the connection for further use 
        $dbh = Utils::connect_bgee_db($bgee_connector);
    }

    
    # We close the connection before compressing the files, it can take quite some time
    $dbh->disconnect;
	if (scalar(keys %fileParams) > 0){
    	unlink File::Spec->catfile($filesDir, $speciesDirName, $speciesNameForFile.'_Affymetrix_experiments_chips.tar.gz');
    	my $tar = Archive::Tar->new();
    	$tar->add_files(File::Spec->catfile($filesDir, $tmpDirName, $expFileName),
    		File::Spec->catfile($filesDir, $tmpDirName, $chipFileName));

    	for my $file_objs ($tar->get_files){
            $tar->rename( $file_objs, basename($file_objs->full_path));
    	}
   		$tar->write(File::Spec->catfile($filesDir, $speciesDirName, $speciesNameForFile.'_Affymetrix_experiments_chips.tar.gz'), 
    		COMPRESS_GZIP);

	    # Compress each file with probesets independently, and add them to a global tar.gz file
    	unlink File::Spec->catfile($filesDir, $speciesDirName, $speciesNameForFile.'_Affymetrix_probesets.tar.gz');
    	$tar = Archive::Tar->new();
    	my @tarFilePaths = ();
    	for my $tarFileName ( @tarFileNames ) {
        	push @tarFilePaths, File::Spec->catfile($filesDir, $speciesDirName, $tarFileName);
    	}
    	$tar->add_files(@tarFilePaths);
    	for my $file_objs ($tar->get_files){
            $tar->rename( $file_objs, basename($file_objs->full_path));
    	}
    	$tar->write(File::Spec->catfile($filesDir, $speciesDirName, $speciesNameForFile.'_Affymetrix_probesets.tar.gz'), 
    		COMPRESS_GZIP);
	}
    remove_tree(File::Spec->catdir( $filesDir, $tmpDirName));
    
    # reopen the connection for further use 
    $dbh = Utils::connect_bgee_db($bgee_connector);
}

sub generateRnaSeqFiles {
    my @args = @_;
    my $speciesId             = $args[0];
    my $speciesName           = $args[1];
    my $dataSourcesRef        = $args[2];
    my $filesDir              = $args[3];
    my $isSingleCell          = $args[4];
    my $isSampleMultiplexing = $args[5];

    my $speciesNameForFile = $speciesName;
    $speciesNameForFile =~ s/ /_/g;
    my $tmpDirName = $speciesNameForFile.'_tmp';
    my $speciesDirName = $speciesNameForFile;
    my $expFileName = '';
    my $libFileName = '';
    my $fileNamePattern = '';
    my $exprFilePath = $bgeeVersion.'/download/processed_expr_values/';
    if (!$isSingleCell) {
        $expFileName = $speciesNameForFile.'_RNA-Seq_libraries.tsv';
        $libFileName = $speciesNameForFile.'_RNA-Seq_experiments.tsv';
        $fileNamePattern = '_RNA-Seq_read_counts_TPM_';
        $exprFilePath .= 'rna_seq/'.$speciesNameForFile.'/';
    } elsif ($isSingleCell) {
        if ($isSampleMultiplexing) {
            $expFileName = $speciesNameForFile.'_droplet_based_RNA-Seq_libraries.tsv';
            $libFileName = $speciesNameForFile.'_droplet_based_RNA-Seq_experiments.tsv';
            $fileNamePattern = '_Full-Length_SC_RNA-Seq_read_counts_TPM_';
            $exprFilePath .= 'full_length_sc_rna_seq/'.$speciesNameForFile.'/';
        } elsif (!$isSampleMultiplexing) {
            $expFileName = $speciesNameForFile.'_full_length_RNA-Seq_libraries.tsv';
            $libFileName = $speciesNameForFile.'_droplet_based_RNA-Seq_experiments.tsv';
            $fileNamePattern = '_Droplet-Based_SC_RNA-Seq_read_counts_CPM_';
            $exprFilePath .= 'droplet_based_sc_rna_seq/'.$speciesNameForFile.'/';
        }
    }
    if ($expFileName eq '' || $libFileName eq '') {
        die "unrecognized RNASeq datatype";
    }
    # recreate temp directory for experiment and library information
    remove_tree(File::Spec->catdir( $filesDir, $tmpDirName));
    make_path(File::Spec->catdir( $filesDir, $tmpDirName));
    
    # -------------------------------------------------
    # First, we retrieve information about experiments
    # $exp{'expId'}                = exp ID
    # $exp{'name'}                 = exp name
    # $exp{'desc'}                 = exp description
    # $exp{'sourceId'}             = data source ID
    # $exp{'libCount'}             = library count
    # $exp{'conditionCount'}       = informative condition count (expression-mapped condition IDs)
    # $exp{'anatEntityStageCount'} = count of couples organ/stage
    # $exp{'anatEntityCount'}      = anatomical entity count
    # $exp{'stageCount'}           = stage count
    # sex types that are not 'not annotated' and 'NA'
    # $exp{'sexCount'}             = informative sex information count
    # strain count that are not 'not annotated', 'NA', 'confidential_restricted_data'
    # $exp{'strainCount'}          = informative strain count
    my @experiments = ();

    my $sqlExpPart =
    'SELECT t1.rnaSeqExperimentId, t1.rnaSeqExperimentName, t1.rnaSeqExperimentDescription, '
              .'t1.dataSourceId, COUNT(DISTINCT rnaSeqLibraryId) AS libCount, '
              # to count the number of conditions, we need to used the expression-mapped condition IDs, 
              # which merge non-informative annotations such as 'not annotated' and 'NA'
              .'COUNT(DISTINCT exprMappedConditionId) AS conditionCount, '
              .'COUNT(DISTINCT anatEntityId, stageId) AS anatEntityStageCount, '
              .'COUNT(DISTINCT anatEntityId) AS anatEntityCount, '
              # distinct cellTypeIds are also exported to the bulk RNASeq files even if the
              # number should always be 1 (root of the ontology)
              .'COUNT(DISTINCT cellTypeId) AS cellTypeCount '
              .'COUNT(DISTINCT stageId) AS stageCount, '
              .'COUNT(DISTINCT sex)  AS sexCount, '
              .'COUNT(DISTINCT strain) AS strainCount '
              .'FROM rnaSeqExperiment AS t1 '
              .'INNER JOIN rnaSeqLibrary AS t2 ON t1.rnaSeqExperimentId = t2.rnaSeqExperimentId '
              .'INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId '
              .'INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId '
              .'WHERE speciesId = ? ';
    # filter data corresponding to the expected datatype
    if ($isSingleCell) {
        $sqlExpPart .=
               'AND rnaSeqTechnologyIsSingleCell = 1 ';
        if ($isSampleMultiplexing) {
            $sqlExpPart .=
               'AND libraryMultiplexing = 1 ';
        } else {
            $sqlExpPart .=
               'AND libraryMultiplexing = 0 ';
        }
    } else {
        $sqlExpPart .=
               'AND rnaSeqTechnologyIsSingleCell = 0 ';
    }
    $sqlExpPart .= 
               'GROUP BY t1.rnaSeqExperimentId ';
    my $sql = $sqlExpPart
              .'ORDER BY libCount DESC, conditionCount DESC, anatEntityStageCount DESC, '
              .'anatEntityCount DESC, cellTypeCount DESC, stageCount DESC, sexCount DESC, '
              .'strainCount DESC, t1.rnaSeqExperimentId';

    my $stmt = $dbh->prepare($sql);
    $stmt->execute($speciesId)  or die $stmt->errstr;

    while ( my @data = $stmt->fetchrow_array ){
        my %exp;
        $exp{'expId'}                = $data[0];
        $exp{'name'}                 = $data[1];
        $exp{'desc'}                 = $data[2];
        $exp{'sourceId'}             = $data[3];
        $exp{'libCount'}             = $data[4];
        $exp{'conditionCount'}       = $data[5];
        $exp{'anatEntityStageCount'} = $data[6];
        $exp{'anatEntityCount'}      = $data[7];
        $exp{'cellTypeCount'}        = $data[8];
        $exp{'stageCount'}           = $data[9];
        $exp{'sexCount'}             = $data[10];
        $exp{'strainCount'}          = $data[11];

        push @experiments, \%exp;
    }
    
    # recreate species directory for experiment and library information
    if( @experiments > 0){
    	remove_tree(File::Spec->catdir( $filesDir, $speciesDirName));
    	make_path(File::Spec->catdir( $filesDir, $speciesDirName));
    }

    # Print experiment information into file
    my $expFile = File::Spec->catfile($filesDir, $tmpDirName, $expFileName);
    open(my $fh, '>', $expFile) or die "Could not open file '$expFile' $!";
    print $fh "Experiment ID\tExperiment name\t"
              ."Library count\tCondition count\tOrgan-stage count\t"
              ."Organ count\tCelltype Count\tStage count\tSex count\tStrain count\t"
              ."Data source\tData source URL\tBgee normalized data URL\t"
              ."Experiment description\n";
    for my $exp ( @experiments ){
        print $fh $exp->{'expId'}."\t";
        
        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        my $name = $exp->{'name'};
        $name =~ s/"/'/g;
        print $fh '"'.$name.'"'."\t";
        
        print $fh $exp->{'libCount'}."\t"
                  .$exp->{'conditionCount'}."\t"
                  .$exp->{'anatEntityStageCount'}."\t"
                  .$exp->{'anatEntityCount'}."\t"
                  .$exp->{'cellTypeCount'}."\t"
                  .$exp->{'stageCount'}."\t"
                  .$exp->{'sexCount'}."\t"
                  .$exp->{'strainCount'}."\t";
        my $sourceUrl  = 'NA';
        my $sourceName = '';
        if (defined $dataSourcesRef->{$exp->{'sourceId'}}) {
            if ($dataSourcesRef->{$exp->{'sourceId'}}->{'experimentBaseUrl'}) {
                $sourceUrl = $dataSourcesRef->{$exp->{'sourceId'}}->{'experimentBaseUrl'};
                $sourceUrl =~ s/$expIdTag/$exp->{'expId'}/g;
            }
            $sourceName = $dataSourcesRef->{$exp->{'sourceId'}}->{'name'};
        }
        print $fh $sourceName."\t".$sourceUrl."\t";

        my $resultsFileName = $speciesNameForFile.$fileNamePattern.$exp->{'expId'}.'.tsv.gz';
        $resultsFileName =~ s/ /_/g;
        print $fh $ftpFilePath.$exprFilePath.$resultsFileName."\t";

        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        my $desc = $exp->{'desc'};
        $desc =~ s/"/'/g;
        print $fh '"'.$desc.'"';

        print $fh "\n";

    }
    close $fh;

    # -------------------------------------------------
    # Now, information about libraries.
    # TODO: add field allMappedReadsCount when database is updated
    # $lib{'expId'}                            = experiment ID
    # $lib{'libId'}                            = library ID
    # $lib{'platformId'}                       = RNA-Seq  platform ID
    # $lib{'tmmFactor'}                        = TMM normalization factor
    # $lib{'tpmThreshold'}                     = Threshold of TPM value to consider a gene as expressed
    # $lib{'fpkmThreshold'}                    = Threshold of FPKM value to consider a gene as expressed
    # $lib{'allGenesPercentPresent'}           = percentage of genes called present
    # $lib{'proteinCodingGenesPercentPresent'} = percentage of protein coding genes called present
    # $lib{'intergenicRegionsPercentPresent'}  = percentage of intergenic regions called present
    # $lib{'allReadsCount'}                    = total number of reads in library, including those not mapped.
    # $lib{'mappedReadsCount'}                 = number of reads mapped by pseudo alignement.
    # $lib{'minReadLength'}                    = min. length of the reads
    # $lib{'maxReadLength'}                    = max. length of the reads
    # $lib{'libraryType'}                      = enum('single','paired')
    # $lib{'libraryOrientation'}               = enum('forward','reverse','unstranded')
    # $lib{'anatEntityId'}                     = anat. entity ID annotated for this library
    # $lib{'anatEntityName'}                   = anat. entity name annotated for this library
    # $lib{'stageId'}                          = stage ID annotated for this library
    # $lib{'stageName'}                        = stage name annotated for this library
    # $lib{'sex'}                              = annotated sex info (not mapped for expression table)
    # $lib{'strain'}                           = annotated strain info (not mapped for expression table)
    # $lib{'exprMappedAnatEntityId'}           = anat. entity ID remapped for expression table for this library
    # $lib{'exprMappedAnatEntityName'}         = anat. entity name remapped for expression table for this library
    # $lib{'exprMappedStageId'}                = stage ID remapped for expression table for this library
    # $lib{'exprMappedStageName'}              = stage name remapped for expression table for this library
    # $lib{'exprMappedSex'}                    = sex info remapped for expression table
    # $lib{'exprMappedStrain'}                 = strain info remapped for expression table
    # $lib{'rnaSeqLibraryAnnotatedSampleDistinctRankCount'}         = count of distinct ranks in the library
    # $lib{'maxRank'}                          = maximum rank in the corresponding global condition
    # $lib{'sourceId'}                         = data source ID
    # $lib{'runIds'}                           = IDs of runs used, separated by '|'
    my @libs = ();

    $sql = 'SELECT t1.rnaSeqExperimentId, t1.rnaSeqLibraryId, t3.anatEntityId, t4.anatEntityName, ';
    if ($isSingleCell) {
        $sql .= 
               't3.cellTypeId, t5.anatEntityName, t2.cellTypeAuthorAnnotation, '
              .'t20.cellTypeId AS exprMappedCellTypeId, t40.anatEntityName AS exprMappedCellTypeName, ';
    }
    $sql .=    't3.stageId, t6.stageName, t3.sex, '
              .'t3.strain, t2.anatEntityAuthorAnnotation, t2.stageAuthorAnnotation, '
              .'t1.rnaSeqSequencerName, t1.strandSelection, t1.cellCompartment, '
              .'t1.sequencedTranscriptPart, t1.rnaSeqPopulationCapture, t1.genotype, t1.allReadsCount, '
              .'t1.mappedReadsCound, t1.minReadLength, t1.maxReadLength, t1.libraryType, t2.abundanceThreshold, '
              .'t2.abundanceUnit, t2.tmmFactor, t2.allGenesPercentPresent, t1.proteinCodingGenesPercentPresent, '
              .'t2.intergenicRegionsPercentPresent, '
              .'t20.anatEntityId AS exprMappedAnatEntityId, t30.anatEntityName AS exprMappedAnatEntityName, '
              .'t20.stageId AS exprMappedStageId, t50.stageName AS exprMappedStageName, '
              .'t20.sex AS exprMappedSex, t20.strain AS exprMappedStrain, '
              .'t2.rnaSeqLibraryAnnotatedSampleDistinctRankCount, '
              # TODO to change if we ever use globalMaxRank instead of maxRank?
              # But then we would have ranks not only for conditions with data,
              # so I guess it would not be present in this file. To rethink in this case.
              .'t2.rnaSeqLibraryAnnotatedSampleMaxRank, '
              .'t5.dataSourceId, '
              .'GROUP_CONCAT(DISTINCT t6.rnaSeqRunId ORDER BY t6.rnaSeqRunId SEPARATOR "|") AS runIds '
              .'FROM rnaSeqLibrary AS t1 '
              .'INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId '
              .'INNER JOIN cond AS t3 ON t3.conditionId = t2.conditionId '
              .'INNER JOIN anatEntity AS t4 ON t4.anatEntityId = t3.anatEntityId ';
    if ($isSingleCell) {
        $sql .= 
               'INNER JOIN anatEntity AS t5 ON t5.cellTypeId = t3.anatEntityId '
              .'INNER JOIN anatEntity AS t40 ON t20.cellTypeId = t40.anatEntityId ';
    }
    $sql .=    'INNER JOIN stage AS t6 ON t3.stageId = t6.stageId '
              .'INNER JOIN cond AS t20 ON t3.exprMappedConditionId = t20.conditionId '
              .'INNER JOIN anatEntity AS t30 ON t20.anatEntityId = t30.anatEntityId '
              .'INNER JOIN stage AS t50 ON t20.stageId = t40.stageId '
              .'INNER JOIN ('.$sqlExpPart.') AS t7 ON t1.rnaSeqExperimentId =   t7.rnaSeqExperimentId '
              .'LEFT OUTER JOIN rnaSeqRun AS t8 ON t1.rnaSeqLibraryId = t8.rnaSeqLibraryId '
              .'WHERE t3.speciesId = ? ';
    if ($isSingleCell) {
        $sql .=
              'AND rnaSeqTechnologyIsSingleCell = 1 ';
        if ($isSampleMultiplexing) {
            $sql .=
              'AND libraryMultiplexing = 1 ';
        } else {
            $sql .=
              'AND libraryMultiplexing = 0 ';
        }
    } else {
        $sql .=
              'AND rnaSeqTechnologyIsSingleCell = 0 ';
    }
    $sql .=     'GROUP BY t1.rnaSeqLibraryId '
              .'ORDER BY libCount DESC, conditionCount DESC, anatEntityStageCount DESC, '
              .'anatEntityCount DESC, cellTypeCount DESC, '
              .'stageCount DESC, sexCount DESC, strainCount DESC, '
              .'t1.rnaSeqExperimentId, t3.anatEntityId, t3.stageId, t3.cellTypeId, '
              .'t3.sex, t3.strain, t1.rnaSeqLibraryId';

    $stmt = $dbh->prepare($sql);
    $stmt->execute($speciesId, $speciesId)  or die $stmt->errstr;
    my $i = 0;

    while ( my @data = $stmt->fetchrow_array ){
        my %lib;

        $lib{'expId'}                            = $data[$i++];
        $lib{'libId'}                            = $data[$i++];
        $lib{'anatEntityId'}                     = $data[$i++];
        $lib{'anatEntityName'}                   = $data[$i++];
        if ($isSingleCell) {
            $lib{'cellTypeId'}                   = $data[$i++];
            $lib{'cellTypeName'}                 = $data[$i++];
            $lib{'cellTypeAuthorAnnotation'}     = $data[$i++];
            $lib{'exprMappedCellTypeId'}         = $data[$i++];
            $lib{'exprMappedCellTypeName'}       = $data[$i++];
        }

        $lib{'stageId'}                          = $data[$i++];
        $lib{'stageName'}                        = $data[$i++];
        $lib{'sex'}                              = $data[$i++];
        $lib{'strain'}                           = $data[$i++];
        $lib{'anatEntityAuthorAnnotation'}       = $data[$i++];
        $lib{'stageAuthorAnnotation'}            = $data[$i++];
        $lib{'sequencerName'}                    = $data[$i++];
        $lib{'strandSelection'}                  = $data[$i++];
        $lib{'cellCompartment'}                  = $data[$i++];
        $lib{'sequencedTranscriptPart'}          = $data[$i++];
        $lib{'rnaSeqPopulationCapture'}          = $data[$i++];
        $lib{'genotype'}                         = $data[$i++];
        $lib{'allReadsCount'}                    = $data[$i++];
        $lib{'mappedReadsCound'}                 = $data[$i++];
        $lib{'minReadLength'}                    = $data[$i++];
        $lib{'maxReadLength'}                    = $data[$i++];
        $lib{'libraryType'}                      = $data[$i++];
        $lib{'abundanceThreshold'}               = $data[$i++];
        $lib{'abundanceUnit'}                    = $data[$i++];
        $lib{'tmmFactor'}                        = $data[$i++];
        $lib{'allGenesPercentPresent'}           = $data[$i++];
        $lib{'proteinCodingGenesPercentPresent'} = $data[$i++];
        $lib{'intergenicRegionsPercentPresent'}  = $data[$i++];
        $lib{'exprMappedAnatEntityId'}           = $data[$i++];
        $lib{'exprMappedAnatEntityName'}         = $data[$i++];
        $lib{'exprMappedStageId'}                = $data[$i++];
        $lib{'exprMappedStageName'}              = $data[$i++];
        $lib{'exprMappedSex'}                    = $data[$i++];
        $lib{'exprMappedStrain'}                 = $data[$i++];
        $lib{'rnaSeqLibraryAnnotatedSampleDistinctRankCount'}  = $data[$i++];
        $lib{'rnaSeqLibraryAnnotatedSampleMaxRank'}            = $data[$i++];
        $lib{'sourceId'}                         = $data[$i++];
        $lib{'runIds'}                           = $data[$i++];
        push @libs, \%lib;
    }
    # Print library information into file
    my $libFile = File::Spec->catfile($filesDir, $tmpDirName, $libFileName);
    open($fh, '>', $libFile) or die "Could not open file '$libFile' $!";
    print $fh  "Experiment ID\tLibrary ID\tAnatomical entity ID\tAnatomical entity name\t"
              ."Anatomical entity author annotation\t";
    if ($isSingleCell) {
        print $fh "Celltype ID\tCelltype name\tCelltype author annotation\t";
    }
    print $fh  "Stage ID\tStage name\tStage author annotation\tSex\tStrain\t"
              ."Expression mapped anatomical entity ID\tExpression mapped anatomical entity name\t";
    if ($isSingleCell) {
        print $fh "Expression mapped celltype ID\tExpression mapped celltype name\t";
    }
    print $fh  "Expression mapped stage ID\tExpression mapped stage name\t"
              ."Expression mapped sex\tExpression mapped strain\t"
              ."Platform ID\tProtocol\tLibrary type\tLibrary orientation\t"
              ."TMM normalization factor\tExpression threshold\tExpression unit\t"
              ."Cell compartment\tSequenced transcript part\tGenotype\t"
              ."Read count\tMapped read count\t"
              ."Min. read length\tMax. read length\tAll genes percent present\t"
              ."Protein coding genes percent present\tIntergenic regions percent present\t"
              ."Distinct rank count\tMax rank in the expression mapped condition\tRun IDs\t"
              ."Data source\tData source URL\tBgee normalized data URL\tRaw file URL\n";

    for my $lib ( @libs ){
        print $fh $lib->{'expId'}."\t"
                  .$lib->{'libId'}."\t";

        print $fh $lib->{'anatEntityId'}."\t";
        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        my $toPrint = $lib->{'anatEntityName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";
        if ($isSingleCell) {
            print $fh $lib->{'cellTypeId'}."\t";
            $toPrint = $lib->{'cellTypeName'};
            $toPrint =~ s/"/'/g;
            print $fh '"'.$toPrint.'"'."\t";
            $toPrint = $lib->{'cellTypeAuthorAnnotation'};
            $toPrint =~ s/"/'/g;
            print $fh '"'.$toPrint.'"'."\t";
        }

        print $fh $lib->{'stageId'}."\t";
        $toPrint = $lib->{'stageName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";
        
        print $fh $lib->{'sex'}."\t";
        
        $toPrint = $lib->{'strain'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        print $fh $lib->{'exprMappedAnatEntityId'}."\t";
        # we replace double quotes with simple quotes, and we surround with double quotes
        # the values to escape potential special characters
        $toPrint = $lib->{'exprMappedAnatEntityName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        print $fh $lib->{'exprMappedStageId'}."\t";
        $toPrint = $lib->{'exprMappedStageName'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        print $fh $lib->{'exprMappedSex'}."\t";

        $toPrint = $lib->{'exprMappedStrain'};
        $toPrint =~ s/"/'/g;
        print $fh '"'.$toPrint.'"'."\t";

        print $fh $lib->{'sequencerName'}."\t".$lib->{'rnaSeqPopulationCapture'}."\t".$lib->{'libraryType'}."\t".$lib->{'strandSelection'}."\t"
            .$lib->{'tmmNormalisation'}."\t".$lib->{'abundanceThreshold'}."\t".$lib->{'abundanceUnit'}."\t"
            .$lib->{'cellCompartment'}."\t".$lib->{'sequencedTranscriptPart'}."\t".$lib->{'genotype'}."\t"
            .$lib->{'allReadsCount'}."\t".$lib->{'mappedReadsCount'}."\t";

        print $fh $lib->{'minReadLength'}."\t".$lib->{'maxReadLength'}."\t"
            .$lib->{'allGenesPercentPresent'}."\t".$lib->{'proteinCodingGenesPercentPresent'}."\t"
            .$lib->{'intergenicRegionsPercentPresent'}."\t";

        my $maxRank = "NA";
        if(defined $lib->{'maxRank'}) {
            $maxRank = $lib->{'maxRank'};
        }
        print $fh $lib->{'rnaSeqLibraryAnnotatedSampleDistinctRankCount'}."\t".$maxRank."\t";

        if (defined $lib->{'runIds'} && $lib->{'runIds'}) {
            print $fh $lib->{'runIds'};
        } else {
            print $fh 'NA';
        }
        print $fh "\t";

        my $sourceUrl  = 'NA';
        my $sourceName = '';
        if (defined $dataSourcesRef->{$lib->{'sourceId'}}) {
            if ($dataSourcesRef->{$lib->{'sourceId'}}->{'evidenceBaseUrl'}) {
                $sourceUrl = $dataSourcesRef->{$lib->{'sourceId'}}->{'evidenceBaseUrl'};
                $sourceUrl =~ s/$expIdTag/$lib->{'expId'}/g;
                $sourceUrl =~ s/$evidenceIdTag/$lib->{'libId'}/g;
            }
            $sourceName = $dataSourcesRef->{$lib->{'sourceId'}}->{'name'};
        }
        print $fh $sourceName."\t".$sourceUrl."\t";

        my $resultsFileName = $speciesNameForFile.$fileNamePattern.$lib->{'expId'}.'.tsv.gz';
        $resultsFileName =~ s/ /_/g;
        print $fh $ftpFilePath.$exprFilePath.$resultsFileName."\t";

        print $fh 'https://trace.ncbi.nlm.nih.gov/Traces/study/?acc='.$lib->{'libId'}."\n";
    }
    close $fh;

    # -------------------------------------------------
    # Now, we generate the RNA-Seq results files,
    # one file, per experiment.
    # But with GTEx, we can't load all rnaSeqResults at once in memory, 
    # so we query one library at a time.
    # XXX: Do we normalize libraries of different types independently?
    # First, we get the list of exp (and not the library types for now)
    # $fileParams{expId} = 1;
    my %fileParams;
    for my $lib ( @libs ){
    	if (!exists $fileParams{$lib->{'expId'}}) {
    		$fileParams{$lib->{'expId'}} = ();
    	}
    	push @{$fileParams{$lib->{'expId'}}}, $lib->{'libId'};
    }

    # Now, generate the files. We'll store the path to all gz files generated,
    # to make one giant tar.gz at the end.
    my @tarFileNames = ();
    for my $expId ( keys %fileParams ){
        # we will store all names of files generated for an experiment, to pack all of them together
        my @resultsFileNames = ();
    
        # UPDATES Bgee 15.2
        # file header to update
        # - read count -> can be "read count" or "UMI count"
        # - TPM -> can be "TPM" or "CPM"
        # - FPKM -> removed
        # - cellTypeId -> new column
        # - cellTypeName -> new column

        # Because we reopen the connection in the expId loop we cannot prepare the query outside the loop.
        # XXX: left outer join to expression to retrieve the global call quality?
        $sql =           'SELECT t4.rnaSeqExperimentId, t3.rnaSeqLibraryId, t4.libraryType, t2.geneId, '
                        .'t5.anatEntityId, t6.anatEntityName, ';
        if ($isSingleCell) {
            $sql .=      't5.cellTypeId, t7.anatEntityName, ';
        }
        $sql .=          't5.stageId, t8.stageName, t5.sex, t5.strain, ';
        if ($isSingleCell && $isSampleMultiplexing) {
            $sql .=      't1.UMIsCount, ';
        } else {
            $sql .=      't1.ReadsCount, ';
        }
        $sql .=          't1.abundance, t1.rawRank, t1.pValue, t1.reasonForExclusion, '
                        .'FROM rnaSeqLibraryAnnotatedSampleGeneResult AS t1 '
                        .'INNER JOIN gene AS t2 ON t2.bgeeGeneId = t1.bgeeGeneId '
                        .'INNER JOIN rnaSeqLibraryAnnotatedSample AS t3 ON t3.rnaSeqLibraryAnnotatedSampleId = t1.rnaSeqLibraryAnnotatedSampleId '
                        .'INNER JOIN rnaSeqLibrary AS t4 ON t4.rnaSeqLibraryId = t3.rnaSeqLibraryId '
                        .'INNER JOIN cond AS t5 ON t5.conditionId = t3.conditionId '
                        .'INNER JOIN anatEntity AS t6 ON t6.anatEntityId = t5.anatEntityId ';
        if ($isSingleCell) {
            $sql .=      'INNER JOIN anatEntity AS t7 ON t7.anatEntityId = t5.cellTypeId ';
        }
        $sql .=          'INNER JOIN stage AS t8 ON t8.stageId = t5.stageId '
                        # removed that useless join as part of Bgee 15.2
                        #.'LEFT OUTER JOIN expression AS t7 ON t1.expressionId = t7.expressionId '
                        .'WHERE t1.rnaSeqLibraryId = ? '
                        .'AND t1.expressionId IS NOT NULL '
                        .'AND t5.speciesId = ? '
                        .'ORDER BY t2.geneId';
        $stmt = $dbh->prepare($sql);
                
        # Note: actually, it is not enough to simply sort the libraries based on their ID, 
        # we need to apply the same sorting as when generating the library info file. 
        # A simple solution is to redo a SQL query with the same ORDER BY clause. 
        my $selectOrderedLibraryIds =   'SELECT t1.rnaSeqLibraryId FROM rnaSeqLibrary AS t1 '
                                       .'INNER JOIN rnaSeqLibraryAnnotatedSample AS t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId '
                                       .'INNER JOIN cond AS t3 ON t2.conditionId = t3.conditionId '
                                       .'WHERE t1.rnaSeqExperimentId = ? '
                                       .'AND t3.speciesId = ? '
                                       .'ORDER BY t3.anatEntityId, ';
        if ($isSingleCell) {
            $selectOrderedLibraryIds .= 't3.cellTypeId, ';
        }
        $selectOrderedLibraryIds .=  't3.stageId, t3.sex, t3.strain, t1.rnaSeqLibraryId';
        
        my $getExpLibs = $dbh->prepare($selectOrderedLibraryIds);
        my $resultsFileName = $speciesNameForFile.$fileNamePattern
            .$expId.'.tsv';
        $resultsFileName =~ s/ /_/g;
        push @resultsFileNames, $resultsFileName;
        my $resultsFile = File::Spec->catfile($filesDir, $tmpDirName,
            $resultsFileName);
        open(my $fh, '>', $resultsFile) or die "Could not open file '$resultsFile' $!";
        print $fh "Experiment ID\tLibrary ID\tLibrary type\tGene ID\t"
                  ."Anatomical entity ID\tAnatomical entity name\t";
        if ($isSingleCell) {
            print $fh "Celltype ID\tCelltype name\t";
        }
        print $fh "Stage ID\tStage name\tSex\tStrain\t";
        if ($isSampleMultiplexing) {
            print $fh "UMI count\tCPM\t";
        } else {
            print $fh "Read count\tTPM\t";
        }
        print $fh "Rank\tDetection flag\tpValue\tState in Bgee\n";
        
        $getExpLibs->execute($expId, $speciesId) or die $getExpLibs->errstr;

        while ( my @libs = $getExpLibs->fetchrow_array ){
            my $libId = $libs[0];
            # Retrieve data from database.
            $stmt->execute($libId, $speciesId) or die $stmt->errstr;
            my $i = 0;
            while ( my @data = $stmt->fetchrow_array ){
                # we write the data directly to not store them in memory
                print $fh $data[$i++]."\t".$data[$i++]."\t".$data[$i++]."\t".$data[$i++]."\t"
                    .$data[$i++]."\t";
                # we replace double quotes with simple quotes, and we surround with double quotes
                # the values to escape potential special characters
                my $toPrint = $data[$i++];
                $toPrint =~ s/"/'/g;
                print $fh '"'.$toPrint.'"'."\t";
                #print celltypeId and celltypeName
                if ($isSingleCell) {
                    print $fh $data[$i++]."\t";
                    my $toPrint = $data[$i++];
                    $toPrint =~ s/"/'/g;
                    print $fh '"'.$toPrint.'"'."\t";
                }
                print $fh $data[$i++]."\t";
                $toPrint = $data[$i++];
                $toPrint =~ s/"/'/g;
                print $fh '"'.$toPrint.'"'."\t";

                print $fh $data[8]."\t";
                # print strain
                $toPrint = $data[9];
                $toPrint =~ s/"/'/g;
                print $fh '"'.$toPrint.'"'."\t";

                # Read count and abundance (TPM or CPM), rank, detection
                print $fh $data[$i++]."\t".$data[$i++]."\t";
                # rank, detection flag and pValue
                my $rank = "NA";
                if (defined $data[$i]) {
                    $rank = $data[$i];
                }
                $i++;
                print $fh $rank."\t".$data[$i++]."\t";
                my $pValue = "NA";
                if(defined $data[$i]) {
                    $pValue = $data[$i];
                }
                $i++;
                print $fh $pValue."\t";
                if ($pValue eq "NA" || $pValue > 0.05) {
                    print $fh "absent\t";
                } else {
                    print $fh "present\t";
                }
                if ($data[$i] eq $Utils::CALL_NOT_EXCLUDED) {
                	print $fh 'Part of a call';
                } else {
                    print $fh 'Result excluded, reason: '.$data[$i];
                }
                $i++;
                print $fh "\n";
            }
        }
        close $fh;
        
        # We close the connection before compressing the files, it can take quite some time
        $dbh->disconnect;
        # we compress all files for an experiment together, and delete the uncompressed files
        # to save disk space
        my $fileName = $speciesNameForFile.'_RNA-Seq_read_counts_TPM_FPKM_'.$expId.'.tsv.gz';
        $fileName =~ s/ /_/g;
        unlink File::Spec->catfile($filesDir, $speciesDirName, $fileName);
        # use shell command to compress files because perl compression modules load data in memory
        # to compress. This is not possible because GTeX experiment is too big.
        # https://metacpan.org/pod/distribution/Archive-Tar/lib/Archive/Tar.pm#FAQ
        # https://www.perlmonks.org/?node_id=201013
        system("gzip -cvf ".File::Spec->catfile($filesDir, $tmpDirName, $resultsFileName)
            ." > ".File::Spec->catfile($filesDir, $speciesDirName, $fileName));
        #it is now safe to remove uncompressed file to save disk space
        unlink File::Spec->catfile($filesDir, $tmpDirName, $resultsFileName); 
        # Store the name of this experiment gz file, to make one giant tar.gz of all experiments
        # at the end.
        push @tarFileNames, $fileName;

        # reopen the connection for further use 
        $dbh = Utils::connect_bgee_db($bgee_connector);
    }


    # -------------------------------------------------
    #everything went fine, we move and tar.gz the tmp files
    
    # We close the connection before compressing the files, it can take quite some time
    $dbh->disconnect;

	if (scalar(keys %fileParams) > 0){
    	unlink File::Spec->catfile($filesDir, $speciesDirName, $speciesNameForFile.'_RNA-Seq_experiments_libraries.tar.gz');
	    my $tar = Archive::Tar->new();
	    $tar->add_files(File::Spec->catfile($filesDir, $tmpDirName, $expFileName),
	    	File::Spec->catfile($filesDir, $tmpDirName, $libFileName));
	    for my $file_objs ($tar->get_files){
            $tar->rename( $file_objs, basename($file_objs->full_path));
	    }
	    $tar->write(
	        File::Spec->catfile($filesDir, $speciesDirName, $speciesNameForFile.'_RNA-Seq_experiments_libraries.tar.gz'), 
		        COMPRESS_GZIP);

 	   	# Compress each file with RNA-Seq results independently, and add them to a global tar.gz file
 	   	unlink File::Spec->catfile($filesDir, $speciesDirName, $speciesNameForFile.'_RNA-Seq_read_counts_TPM_FPKM.tar.gz');
   		$tar = Archive::Tar->new();
    
    	my @tarFilePaths = ();
	    for my $tarFileName ( @tarFileNames ) {
	       	push @tarFilePaths, File::Spec->catfile($filesDir, $speciesDirName, $tarFileName);
	    }
	    $tar->add_files(@tarFilePaths);
	    for my $file_objs ($tar->get_files){
            $tar->rename( $file_objs, basename($file_objs->full_path));
	    }
	    $tar->write(
	        File::Spec->catfile($filesDir, $speciesDirName, $speciesNameForFile.'_RNA-Seq_read_counts_TPM_FPKM.tar.gz'),
	        COMPRESS_GZIP);
	}

    remove_tree(File::Spec->catdir( $filesDir, $tmpDirName));
    
    # reopen the connection for further use 
    $dbh = Utils::connect_bgee_db($bgee_connector);
}