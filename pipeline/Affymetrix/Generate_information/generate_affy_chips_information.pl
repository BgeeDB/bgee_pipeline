#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Frederic Bastian, created July 2012
#
# This script generates various information for affymetrix chips:
#   compute sha512 checksums and "unique strings" to check for duplicated files
#   compute quality scores and percent present, extract cdf_name
# This will be done only for chips for which this information is not already generated (not already in the information file)
#
# This script will then inform about duplicated chips and low quality chips
# for which this information was not already computed.
# And then, will inform about all duplicated chips and low quality chips,
# even for which this information was already computed.
#############################################################

$| = 1;

use Digest::SHA;
use File::Spec;
use FindBin;
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

require 'affy_utils.pl';
require 'bgee_utils.pl';


# Define arguments & their default value
my ($affyChipFilesDir, $affymetrixChip, $affymetrixChipInformation, $chipTypeQual) = ('', '', '', '');
my ($debug) = (0);
my %opts = ('affyChipFilesDir=s'            => \$affyChipFilesDir,
            'affymetrixChip=s'              => \$affymetrixChip,
            'affymetrixChipInformation=s'   => \$affymetrixChipInformation,
            'chipTypeQual=s'                => \$chipTypeQual,
            'debug'                         => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $affyChipFilesDir eq '' || $affymetrixChip eq '' || $affymetrixChipInformation eq '' || $chipTypeQual eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -affyChipFilesDir=\$(AFFYDATAPATH)  -affymetrixChip=\$(AFFY_CHIP_FILEPATH)  -affymetrixChipInformation=\$(AFFY_CHIPINFO_FILEPATH)  -chipTypeQual=\$(AFFY_CHIPTYPEQUAL_FILEPATH)
\t-affyChipFilesDir             Main Affymetrix directory
\t-affymetrixChip               affymetrixChip                              annotation file
\t-affymetrixChipInformation    affymetrixChipInformation                   pipeline   file
\t-chipTypeQual                 chipTypeCorrespondencesAndQualityThresholds pipeline   file
\n";
    exit 1;
}

if ( !-e $affyChipFilesDir ){
    die "Could not find [$affyChipFilesDir] (should be '\$(AFFYDATAPATH)')\n";
}


##########################################
##   START OF THE SCRIPT
##########################################

# First, we get all chips from the affymetrixChip file
# $affyChips{expId}{chipId}
my %affyChips = getAllChips( $affymetrixChip );

# then we get all chips info already computed (checksums, quality scores, ...)
# $affyChipsInfo{expId}{chipId}
my %affyChipsInfo = getAllChipsInfo( $affymetrixChipInformation );
my $countExistingInfo = 0;
for my $expId ( keys %affyChipsInfo ){
    for my $chipId ( keys %{$affyChipsInfo{$expId}} ){
        $countExistingInfo++;
    }
}

# Now we generate information for chips for which it hasn't been computed yet
# (checksums, quality score, ...)
my %newAffyChipsInfo = generateNewInfo(\%affyChips, \%affyChipsInfo);
if ( scalar keys %newAffyChipsInfo != 0 ){
    print "Some new files were detected, information was computed for them\n";
    print "Write new info in file...\n";
    my $count = writeNewInfoInFile(\%newAffyChipsInfo);
    print "$count lines written\n";
    print 'This should make, in affymetrixChipInformation: 1 header + ',
          "$countExistingInfo already existing lines + $count new lines = ",
          (1 + $countExistingInfo + $count), '. You can check it by using: ',
          'wc ', $affymetrixChipInformation,
          ". Also, check the new entries carefully.\n";
}
else {
    print "No new info generated\n"  if ( $debug );
}


# Get quality score thresholds for cdf names
my %chipTypeInfo = getChipTypesInformation( $chipTypeQual );

# Display new duplicated chips
if ( scalar keys %newAffyChipsInfo != 0 ){
    print "-------------------\n";
    print "Try to find duplicates amongst new chips:\n";
    displayDuplicatedChips(\%newAffyChipsInfo, \%affyChips, \%newAffyChipsInfo, 0);
    print "\n-------------------\n";
    print "Try to find duplicates between new chips and already annotated chips:\n";
    displayDuplicatedChips(\%newAffyChipsInfo, \%affyChips, \%affyChipsInfo, 0);

    # Display new chips with quality score below thresholds or incompatible chipType
    print "-------------------\n";
    print "Try to find new chips with quality score below thresholds, or incompatible chip type:\n";
    displayLowQualityAndIncompatibleChips(\%newAffyChipsInfo, \%affyChips, \%chipTypeInfo, 0, 1);
}


# Display already annotated duplicated chips
print "\n-------------------\n";
print "Try to find duplicates amongst already annotated chips:\n";
displayDuplicatedChips(\%affyChipsInfo, \%affyChips, \%affyChipsInfo, 0);


# Display already annotated with quality score below thresholds
print "\n-------------------\n";
print "Try to find already annotated chips with quality score below thresholds:\n";
displayLowQualityChips(\%affyChipsInfo, \%affyChips, \%chipTypeInfo, 0);


print "\n-------------------\n";
print "Try to find all incompatible annotated chips:\n";
displayIncompatibleChips(\%affyChips, \%chipTypeInfo, 0);


if ( scalar keys %newAffyChipsInfo != 0 ){
    print "\n-------------------\n";
    print 'Now, for information, try to find duplicates amongst new chips that are used, ',
          "compared to new chips that are commented (should not even exist...), if we have info for them:\n";
    displayDuplicatedChips(\%newAffyChipsInfo, \%affyChips, \%newAffyChipsInfo, 1);
    print "\n-------------------\n";
    print 'Now, for information, try to find duplicates amongst new chips that are used, ',
          "compared to already annotated chips that are commented, if we have info for them:\n";
    displayDuplicatedChips(\%newAffyChipsInfo, \%affyChips, \%affyChipsInfo, 1);

    # Display new chips with quality score below thresholds
    print "\n-------------------\n";
    print 'Now, for information, try to find new chips that are commented (should not even exist...), incompatible or ',
          "with quality score below thresholds, if we have info for them:\n";
    displayLowQualityAndIncompatibleChips(\%newAffyChipsInfo, \%affyChips, \%chipTypeInfo, 1, 1);
}

# Next things are too verbose
exit 0  if ( !$debug );
print "\n\n============================================================\n";
print "Now all this part below is just for your information, comparing used chips to commented chips, not really important\n";
print "============================================================\n";
print "\n\n";

print "\n-------------------\n";
print 'Now, for information, try to find duplicates amongst already annotated chips that are used, ',
      "compared to already annotated chips that are commented, if we have info for them:\n";
displayDuplicatedChips(\%affyChipsInfo, \%affyChips, \%affyChipsInfo, 1);

# Display already annotated with quality score below thresholds
print "\n-------------------\n";
print 'Now, for information, try to find already annotated chips that are commented, ',
      "with quality score below thresholds, if we have info for them:\n";
displayLowQualityChips(\%affyChipsInfo, \%affyChips, \%chipTypeInfo, 1);


print "\n-------------------\n";
print "Now, for information, try to find all incompatible commented chips, if we have info for them:\n";
displayIncompatibleChips(\%affyChips, \%chipTypeInfo, 1);

exit 0;



# Sub to get the affymetrixChip file, or the affymetrixChipInformation file,
# depending on the argument provided to the sub (either 'affymetrixChip', or 'affymetrixChipInformation',
# or 'qualityThresholds')
sub get_output_file {
    my ($fileName) = @_;

    my ($volume, $directories, $file) = File::Spec->splitpath($affyChipFilesDir, 1);
    my @dirs = File::Spec->splitdir($directories);
    push(@dirs, 'chip_information');
    push(@dirs, 'logs');

    my $filePath = File::Spec->catfile(@dirs, $fileName);
    if ( !File::Spec->file_name_is_absolute($filePath) ){
        $filePath = File::Spec->rel2abs($filePath);
    }
    return $filePath;
}

sub get_result_file {
    my ($fileName) = @_;

    my ($volume, $directories, $file) = File::Spec->splitpath($affyChipFilesDir, 1);
    my @dirs = File::Spec->splitdir($directories);
    push(@dirs, 'chip_information');
    push(@dirs, 'results');

    my $filePath = File::Spec->catfile(@dirs, $fileName);
    if ( !File::Spec->file_name_is_absolute($filePath) ){
        $filePath = File::Spec->rel2abs($filePath);
    }
    return $filePath;
}

# Sub to get an affymetrix file (like a CEL file, or a processed_mas5 file)
# arguments are: experimentId, affymetrixChipId, fileType (either 'cel_data', 'processed_mas5_original_files', or 'processed_mas5')
sub get_affy_chip_file {
    my ($expId, $chipId, $fileDir) = @_;

    my ($volume, $directories, $file) = File::Spec->splitpath($affyChipFilesDir, 1);
    my @dirs = File::Spec->splitdir($directories);
    push(@dirs, $fileDir);
    push(@dirs, $expId);

    my $filePath = File::Spec->catfile(@dirs, $chipId);
    if ( !File::Spec->file_name_is_absolute($filePath) ){
        $filePath = File::Spec->rel2abs($filePath);
    }
    return $filePath;
}

sub generateCheckSum {
    my ($file, $removeFirstLine) = @_;
    # To define whether the first line of the file should be removed before generating the checksum
    # It is the case for filtered processed mas5, where lines are ordered by probeset IDs, but header is not standardized
    # Better chance to detect duplicates by removing it

    my $sha = Digest::SHA->new(512);

    if ( $removeFirstLine ){
        open(my $TEMPOUT, '>', 'firtLineRemovedTempFile')  or die "Could not open temp file\n";
        open(my $IN,      '<', "$file")                    or die "Could not read [$file]\n";
        my $line = <$IN>;#first line removed
        while ( defined ($line = <$IN>) ){
            print ${TEMPOUT} $line;
        }
        close $IN;
        close $TEMPOUT;

        $sha->addfile('firtLineRemovedTempFile', 'p');
        unlink('firtLineRemovedTempFile');
    }
    else {
        $sha->addfile($file, 'p');
    }

    return $sha->hexdigest();
}

sub generateChipInformation {
    # Path to the cel file or mas5 filtered file
    my ($file, $expId, $chipId, $normalizationTypeId) = @_;

    my $cdfName;
    my $scanDate;
    my $qualityScore;
    my $percentPresent;
    my $expressionBasedUniqueString;
    my $errorCode;
    my $errorMessage;

    # Marta R script here, read results from a temp file
    # input:
    my $celDataRScriptFile = $FindBin::Bin.'/get_celfile_info.R'; # Same directory than this script
    my $mas5RScriptFile    = $FindBin::Bin.'/get_MAS5_info.R';

    # output:
    my $logOutputFile      = get_output_file($expId.'_'.$chipId.'.out');
    # tmp file to retrieve results generated by the R scripts
    my $resultOutput       = get_result_file($expId.'_'.$chipId.'.tsv');

    my $problem = 0;
    # cel file
    if ( $normalizationTypeId eq '2' ){
        my $chipTypeInfoFile = $chipTypeQual;
        # Launch R script
        ## R CMD BATCH --no-save --no-restore '--args celfile_path=$celfile_path chipType_info_path=$chipType_info_path output_file=$output_file' get_celfile_info.R $Rout_path
        my @args = ("R CMD BATCH --no-save --no-restore '--args celfile_path=\"".$file."\" chipType_info_path=\"".$chipTypeInfoFile."\" output_file=\"".$resultOutput."\"' ".$celDataRScriptFile.' '.$logOutputFile);
        system(@args)==0  or do{warn "\tsystem @args failed: $? "; system("mv $logOutputFile $logOutputFile.PROB"); $problem = 1;};

        if ( !$problem ){
            ($cdfName, $scanDate, $qualityScore, $percentPresent, $expressionBasedUniqueString, $errorCode, $errorMessage) = readCelFileInformationResult($resultOutput);
            if ( $errorCode ne '' ){
                if ( $errorCode eq 'reading error' || $errorCode eq 'processing error' ){
                    print "Error while analyzing $file: $errorCode - $errorMessage\n";
                }
                elsif ( $errorCode eq 'unknown cdfname' ){
                    print "Error while analyzing $file: new cdf called '$cdfName'. You should add the corresponding chip type Id to chipType in annotation.xls, ".
                        "and the required information in chipTypeCorrespondencesAndQualityThresholds\n";
                }
                elsif ( $errorCode eq 'incompatible chipType' ){
                    print "Warning while analyzing $file: the chip type with cdf name '$cdfName' is not compatible with the current pipeline.\n";
                }
            }
        }
    }
    else { # MAS5
        # Launch R script
        ## R CMD BATCH --no-save --no-restore '--args celfile_path=$celfile_path chipType_info_path=$chipType_info_path output_file=$output_file' get_celfile_info.R $Rout_path
        my @args = ("R CMD BATCH --no-save --no-restore '--args MAS5_path=\"".$file."\" output_file=\"".$resultOutput."\"' ".$mas5RScriptFile." ".$logOutputFile);
        system(@args)==0  or do{warn "\tsystem @args failed: $? "; system("mv $logOutputFile $logOutputFile.PROB"); $problem = 1;};

        if ( !$problem ){
            ($percentPresent, $expressionBasedUniqueString, $errorCode, $errorMessage) = readMAS5FileInformationResult($resultOutput);
            if ( $errorCode ne '' ){
                print "Error while analyzing $file: $errorCode - $errorMessage\n";
            }
        }
    }

    $cdfName                     = ''
        if ( !defined $cdfName                     || $cdfName                     =~ /^false$/i || $cdfName                     =~ /^n\/?a$/i );
    $scanDate                    = ''
        if ( !defined $scanDate                    || $scanDate                    =~ /^false$/i || $scanDate                    =~ /^n\/?a$/i );
    $qualityScore                = ''
        if ( !defined $qualityScore                || $qualityScore                =~ /^false$/i || $qualityScore                =~ /^n\/?a$/i );
    $percentPresent              = ''
        if ( !defined $percentPresent              || $percentPresent              =~ /^false$/i || $percentPresent              =~ /^n\/?a$/i );
    $expressionBasedUniqueString = ''
        if ( !defined $expressionBasedUniqueString || $expressionBasedUniqueString =~ /^false$/i || $expressionBasedUniqueString =~ /^n\/?a$/i );
    $expressionBasedUniqueString =~ s{ }{}g;

    return ($cdfName, $scanDate, $qualityScore, $percentPresent, $expressionBasedUniqueString);
}

sub readMAS5FileInformationResult {
    # Path to the output result file
    my ($file) = @_;

    my $percentPresent;
    my $expressionBasedUniqueString;
    my $errorCode;
    my $errorMessage;

    open(my $IN, '<', "$file")  or die "Could not read result file [$file]\n";
    my $line = <$IN>; #header
    # just read the first line, supposed to be the only one
    if ( defined ($line = <$IN>) ){
        my @tmp = map { bgeeTrim($_) }
                  split(/\t/, $line);
        $percentPresent               = $tmp[3];
        $expressionBasedUniqueString  = $tmp[4];
        $errorCode                    = $tmp[5];
        $errorMessage                 = $tmp[6];
    }
    close $IN;

    $percentPresent              = ''
        if ( !defined $percentPresent              || $percentPresent              =~ /^false$/i || $percentPresent              =~ /^n\/?a$/i );
    $expressionBasedUniqueString = ''
        if ( !defined $expressionBasedUniqueString || $expressionBasedUniqueString =~ /^false$/i || $expressionBasedUniqueString =~ /^n\/?a$/i );
    $errorCode                   = ''
        if ( !defined $errorCode                   || $errorCode                   =~ /^false$/i || $errorCode                   =~ /^n\/?a$/i );
    $errorMessage                = ''
        if ( !defined $errorMessage                || $errorMessage                =~ /^false$/i || $errorMessage                =~ /^n\/?a$/i );

    return ($percentPresent, $expressionBasedUniqueString, $errorCode, $errorMessage);
}

sub readCelFileInformationResult {
    # Path to the output result file
    my ($file) = @_;

    my $cdfName;
    my $scanDate;
    my $qualityScore;
    my $percentPresent;
    my $expressionBasedUniqueString;
    my $errorCode;
    my $errorMessage;

    open(my $IN, '<', "$file")  or die "Could not read result file [$file]\n";
    my $line = <$IN>; #header
    # just read the first line, supposed to be the only one
    if ( defined ($line = <$IN>) ){
        my @tmp  = map { bgeeTrim($_) }
                   split(/\t/, $line);
        $cdfName                     = $tmp[3];
        $scanDate                    = $tmp[4];
        $percentPresent              = $tmp[6];
        $qualityScore                = $tmp[7];
        $expressionBasedUniqueString = $tmp[8];
        $errorCode                   = $tmp[9];
        $errorMessage                = $tmp[10];
    }
    close $IN;

    $cdfName                     = ''
        if ( !defined $cdfName                     || $cdfName                     =~ /^false$/i || $cdfName                     =~ /^n\/?a$/i );
    $scanDate                    = ''
        if ( !defined $scanDate                    || $scanDate                    =~ /^false$/i || $scanDate                    =~ /^n\/?a$/i );
    $percentPresent              = ''
        if ( !defined $percentPresent              || $percentPresent              =~ /^false$/i || $percentPresent              =~ /^n\/?a$/i );
    $qualityScore                = ''
        if ( !defined $qualityScore                || $qualityScore                =~ /^false$/i || $qualityScore                =~ /^n\/?a$/i );
    $expressionBasedUniqueString = ''
        if ( !defined $expressionBasedUniqueString || $expressionBasedUniqueString =~ /^false$/i || $expressionBasedUniqueString =~ /^n\/?a$/i );
    $errorCode                   = ''
        if ( !defined $errorCode                   || $errorCode                   =~ /^false$/i || $errorCode                   =~ /^n\/?a$/i );
    $errorMessage                = ''
        if ( !defined $errorMessage                || $errorMessage                =~ /^false$/i || $errorMessage                =~ /^n\/?a$/i );

    return ($cdfName, $scanDate, $qualityScore, $percentPresent, $expressionBasedUniqueString, $errorCode, $errorMessage);
}

sub generateNewInfo {
    my ($affyChipsRef, $affyChipsInfoRef) = @_;

    # Store chips for which we generate new info
    my %newChipsInfo;
    # Get the chip type compatibility
    my %chipTypeInfo = getChipTypesInformation( $chipTypeQual );

    # Search for chips with no info generated
    for my $expId ( sort keys %$affyChipsRef ){
        for my $chipId ( sort keys %{$affyChipsRef->{$expId}} ){
            # Do we already have info for this chip?
            if ( !defined $affyChipsInfoRef->{$expId}->{$chipId} ){
                # Experiment is commented out
                next  if ( $affyChipsRef->{$expId}->{$chipId}->{'commented'} eq '1');
                next  if ( !$affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'} && $affyChipsRef->{$expId}->{$chipId}->{'comments'} =~ /Supplementary data files not provided/i );

                my $expressionBasedUniqueString = '';
                my $celDataChecksum             = '';
                my $mas5OriginalFileChecksum    = '';
                my $mas5Checksum                = '';
                my $qualityScore                = '';
                my $percentPresent              = '';
                my $cdfName                     = '';
                my $scanDate                    = '';

                my $infoToWrite = 0;
                printf('Generating info for chipId %-17s - expId %-10s ... ', $chipId, $expId);

                if ( !defined $chipTypeInfo{$affyChipsRef->{$expId}->{$chipId}->{'chipTypeId'}}{'correspondence'} || $chipTypeInfo{$affyChipsRef->{$expId}->{$chipId}->{'chipTypeId'}}{'correspondence'} eq '' ){
                    print "Unknown chipTypeId, you should edit chipTypeCorrespondencesAndQualityThresholds\n";
                    next;
                }

                # Test if the chip type is supposed to be compatible
                if ( isChipIncompatibleOrLowQuality($affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'}, undef, \%chipTypeInfo, $affyChipsRef->{$expId}->{$chipId}->{'chipTypeId'}) ){
                    print "Incompatible chip type, skipped.\n";
                    next;
                }

                # Now we try to get cel file, mas5 original file, and mas5 filtered file,
                # whatever the normalization type, just in case a cel file has been normalized using MAS5, etc.

                # Try to get a cel file
                my $file = get_affy_chip_file($expId, $chipId, 'cel_data');
                # if the file was gziped
                my $wasGzipped = 0;
                if ( -e $file.'.gz' ){
                    print "gunzip -d ...";
                    gunzip $file.'.gz' => $file  or do{
                        print "Decompression failed, file [$file] skipped: $GunzipError\n";
                        next;
                    };
                    print 'Done';
                    $wasGzipped = 1;
                }

                print 'Generating checksum...';
                # if the cel file doesn't exist
                if ( !-e $file ){
                    # It's a problem only if the file is not commented,
                    # and if it is declared as having a gcRMA normalization type
                    if ( !$affyChipsRef->{$expId}->{$chipId}->{'commented'} && $affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'} == 2 ){
                        print "Error, could not find CEL file [$file]\n";
                        # This needs to be fixed, do not generate chip info
                        next;
                    }
                }
                else {
                    $celDataChecksum = generateCheckSum($file, 0);
                    if ( $celDataChecksum eq '' ){
                        print "Error, problem while computing checksum for [$file]\n";
                        next;
                    }
                    $infoToWrite = 1;
                }

                # Try to get a MAS5 original file
                $file = get_affy_chip_file($expId, $chipId, 'processed_mas5_original_files');
                # if the file doesn't exist
                if ( !-e $file ){
                    # It's a problem only if the file is not commented,
                    # and if it is declared as having a MAS5 normalization type
                    if ( !$affyChipsRef->{$expId}->{$chipId}->{'commented'} && $affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'} == 1 ){
                        print "Error, could not find MAS5 original file [$file]\n";
                        # This needs to be fixed, do not generate chip info
                        next;
                    }
                }
                else {
                    $mas5OriginalFileChecksum = generateCheckSum($file, 0);
                    if ( $mas5OriginalFileChecksum eq '' ){
                        print "Error, problem while computing checksum for [$file]\n";
                        next;
                    }
                    $infoToWrite = 1;
                }

                # Try to get a MAS5 filtered file
                $file = get_affy_chip_file($expId, $chipId, 'processed_mas5');
                # if the file doesn't exist
                if ( !-e $file ){
                    # It's a problem only if the file is not commented,
                    # and if it is declared as having a MAS5 normalization type
                    if ( !$affyChipsRef->{$expId}->{$chipId}->{'commented'} && $affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'} == 1 ){
                        print "Error, could not find MAS5 filtered file [$file]\n";
                        # This needs to be fixed, do not generate chip info
                        next;
                    }
                }
                else {
                    $mas5Checksum = generateCheckSum($file, 1);
                    if ( $mas5Checksum eq '' ){
                        print "Error, problem while computing checksum for [$file]\n";
                        next;
                    }
                    $infoToWrite = 1;
                }
                print 'Done   ';

                print 'Generating chip info...';
                # Now we generate quality score, etc, only for the file actually supposed to be used
                if ( $affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'} == 1 ){
                    # MAS5 filtered file
                    $file = get_affy_chip_file($expId, $chipId, 'processed_mas5');
                    if ( -e $file ){
                        ($cdfName, $scanDate, $qualityScore, $percentPresent, $expressionBasedUniqueString) =
                            generateChipInformation($file, $expId, $chipId, $affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'});
                        if ( $expressionBasedUniqueString eq '' || $percentPresent eq '' ){
                            print "Error, could not get chip info for file [$file]\n";
                            next;
                        }
                        $infoToWrite = 1;
                    }
                }
                elsif ( $affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'} == 2 ){
                    # CEL file
                    $file = get_affy_chip_file($expId, $chipId, 'cel_data');
                    if ( -e $file ){
                        ($cdfName, $scanDate, $qualityScore, $percentPresent, $expressionBasedUniqueString) =
                            generateChipInformation($file, $expId, $chipId, $affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'});
                        if ( $expressionBasedUniqueString eq '' || $percentPresent eq '' || $qualityScore eq '' || $cdfName eq '' || $scanDate eq '' ){
                            print "Error, could not get chip info for file [$file]\n";
                            next;
                        }
                        if ( $chipTypeInfo{$affyChipsRef->{$expId}->{$chipId}->{'chipTypeId'}}{'correspondence'} ne $cdfName ){
                            print "Error, wrong chipTypeId provided\n";
                            next;
                        }
                        $infoToWrite = 1;
                    }
                }
                print 'Done ';

                if ( $wasGzipped == 1 && -e $file.'.gz' ){
                    print 'Removing uncompressed file...';
                    if ( !unlink($file) ){
                        print "Error: could not remove [$file]\n";
                    }
                    else {
                        print 'Done ';
                    }
                }

                if ( $infoToWrite ){
                    $newChipsInfo{$expId}{$chipId}{'cdfName'}                     = $cdfName;
                    $newChipsInfo{$expId}{$chipId}{'scanDate'}                    = $scanDate;
                    $newChipsInfo{$expId}{$chipId}{'qualityScore'}                = $qualityScore;
                    $newChipsInfo{$expId}{$chipId}{'percentPresent'}              = $percentPresent;
                    $newChipsInfo{$expId}{$chipId}{'expressionBasedUniqueString'} = $expressionBasedUniqueString;
                    $newChipsInfo{$expId}{$chipId}{'celDataChecksum'}             = $celDataChecksum;
                    $newChipsInfo{$expId}{$chipId}{'mas5OriginalFileChecksum'}    = $mas5OriginalFileChecksum;
                    $newChipsInfo{$expId}{$chipId}{'mas5Checksum'}                = $mas5Checksum;
                    print "\n";
                }
                else {
                    print "Problem, nothing computed for this chip\n";
                }
            }
        }
    }

    return %newChipsInfo;
}

sub writeNewInfoInFile {
    # a ref to the hash containing chips and their info (generated by generateNewInfo)
    my ($affyChipsInfoRef) = @_;

    open(my $OUT, '>>', "$affymetrixChipInformation")  or die "Could not open file [$affymetrixChipInformation]\n";
    my $count = 0;
    for my $expId ( keys %$affyChipsInfoRef ){
        for my $chipId ( keys %{$affyChipsInfoRef->{$expId}} ){
            print ${OUT} $chipId, "\t", $expId, "\t",
                         $affyChipsInfoRef->{$expId}->{$chipId}{'cdfName'}, "\t",
                         $affyChipsInfoRef->{$expId}->{$chipId}{'scanDate'}, "\t",
                         $affyChipsInfoRef->{$expId}->{$chipId}{'qualityScore'}, "\t",
                         $affyChipsInfoRef->{$expId}->{$chipId}{'percentPresent'}, "\t",
                         $affyChipsInfoRef->{$expId}->{$chipId}{'expressionBasedUniqueString'}, "\t",
                         $affyChipsInfoRef->{$expId}->{$chipId}{'celDataChecksum'}, "\t",
                         $affyChipsInfoRef->{$expId}->{$chipId}{'mas5OriginalFileChecksum'}, "\t",
                         $affyChipsInfoRef->{$expId}->{$chipId}{'mas5Checksum'}, "\t\n";
            $count++;
        }
    }
    close $OUT;

    return $count;
}

