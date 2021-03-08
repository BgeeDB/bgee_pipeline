#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Frederic Bastian, created July 2012

# USAGE: perl check_affy_info.pl
# a script to check that the required info were loaded for used chips
# (quality scores, percent present, ...)
##############################################

use Getopt::Long;
use lib '.';
require 'affy_utils.pl';

# Define arguments & their default value
my ($affymetrixChip) = ('');
my ($affyChipInformation, $chipTypeQual) = ('', '');
my %opts = (
            'affyChipInformation=s'  => \$affyChipInformation,
            'chipTypeQual=s'         => \$chipTypeQual,
            'affymetrixChip=s'       => \$affymetrixChip,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $affymetrixChip eq '' || $affyChipInformation eq '' || $chipTypeQual eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -affyChipInformation=\$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=\$(AFFY_CHIPTYPEQUAL_FILEPATH) -affymetrixChip=\$(AFFY_CHIP_FILEPATH)
\t-affyChipInformation        affymetrixChipInformation                     pipeline   file
\t-chipTypeQual               chipTypeCorrespondencesAndQualityThresholds   pipeline   file
\t-affymetrixChip             affymetrixChip                                annotation file
\n";
    exit 1;
}

$| = 1;


# Get chip annotations
my %affyChips     = getAllChips($affymetrixChip);
# Get chip info
my %affyChipsInfo = getAllChipsInfo($affyChipInformation);
# Get the incompatible chip type
my %chipTypeInfo  = getChipTypesInformation($chipTypeQual);

# Now, check that we have all required info
for my $expId ( %affyChips ){
    for my $chipId ( keys %{$affyChips{$expId}} ){
        # if the chip is commented, skip
        next  if ( $affyChips{$expId}{$chipId}{'commented'} );
        # incompatible chip?
        # We do not check for low quality here (undef value used)
        next  if ( isChipIncompatibleOrLowQuality($affyChips{$expId}{$chipId}{'normalizationTypeId'}, undef, \%chipTypeInfo, $affyChips{$expId}{$chipId}{'chipTypeId'}) );
        # No info for this chip?
        if ( !defined $affyChipsInfo{$expId}{$chipId} ){
            print "Warning, information was not generated for chipId: $chipId - expId: $expId\n";
            next;
        }
        # normalization method defined?
        if ( !defined $affyChips{$expId}{$chipId}{'normalizationTypeId'} || ($affyChips{$expId}{$chipId}{'normalizationTypeId'} ne '1' &&  $affyChips{$expId}{$chipId}{'normalizationTypeId'} ne '2') ){
            print "Warning, normalization type not defined for chipId: $chipId - expId: $expId\n";
            next;
        }

        # for CEL file
        if ( $affyChips{$expId}{$chipId}{'normalizationTypeId'} eq '2' ){
            # Do we have the quality score?
            if ( !defined $affyChipsInfo{$expId}{$chipId}{'qualityScore'}    || $affyChipsInfo{$expId}{$chipId}{'qualityScore'} eq ''    || !$affyChipsInfo{$expId}{$chipId}{'qualityScore'} ){
                print "Warning, quality score arIQR not defined for chipId: $chipId - expId: $expId\n";
                next;
            }
            # Do we have the cel file checksum?
            if ( !defined $affyChipsInfo{$expId}{$chipId}{'celDataChecksum'} || $affyChipsInfo{$expId}{$chipId}{'celDataChecksum'} eq '' || !$affyChipsInfo{$expId}{$chipId}{'celDataChecksum'} ){
                print "Warning, celDataChecksum not defined for chipId: $chipId - expId: $expId\n";
                next;
            }
            # Do we have the scan date?
            if ( !defined $affyChipsInfo{$expId}{$chipId}{'scanDate'}        || $affyChipsInfo{$expId}{$chipId}{'scanDate'} eq ''        || !$affyChipsInfo{$expId}{$chipId}{'scanDate'} ){
                print "Warning, scan date not defined for chipId: $chipId - expId: $expId (it happens sometimes, so it might be harmless)\n";
                #next;
            }
            # Do we have the cdf name?
            if ( !defined $affyChipsInfo{$expId}{$chipId}{'cdfName'}         || $affyChipsInfo{$expId}{$chipId}{'cdfName'} eq ''         || !$affyChipsInfo{$expId}{$chipId}{'cdfName'} ){
                print "Warning, cdfName not defined for chipId: $chipId - expId: $expId\n";
                next;
            }
            # Check that we have the mapping cdf name -> chipTypeId
            if ( !defined $chipTypeInfo{$affyChipsInfo{$expId}{$chipId}{'cdfName'}}{'correspondence'} || $chipTypeInfo{$affyChipsInfo{$expId}{$chipId}{'cdfName'}}{'correspondence'} eq '' ){
                print "Warning, no correspondence provided between cdf name and chip type ID for chipId: $chipId - expId: $expId\n";
                next;
            }
            # Check that the cdf name corresponds to the chipTypeId provided by annotators
            if ( $chipTypeInfo{$affyChipsInfo{$expId}{$chipId}{'cdfName'}}{'correspondence'} ne $affyChips{$expId}{$chipId}{'chipTypeId'} ){
                print 'Warning, the cdf name ', $affyChipsInfo{$expId}{$chipId}{'cdfName'},
                      'does not correspond to the chipTypeId ', $affyChips{$expId}{$chipId}{'chipTypeId'},
                      "provided by curators for chipId: $chipId - expId: $expId\n";
                next;
            }
        }
        # For MAS5 files
        else {
            # Do we have the MAS5 original file checksum?
            if ( !defined $affyChipsInfo{$expId}{$chipId}{'mas5OriginalFileChecksum'} || $affyChipsInfo{$expId}{$chipId}{'mas5OriginalFileChecksum'} eq '' || !$affyChipsInfo{$expId}{$chipId}{'mas5OriginalFileChecksum'} ){
                print "Warning, mas5OriginalFileChecksum not defined for chipId: $chipId - expId: $expId\n";
                next;
            }
            # Do we have the MAS5 filtered file checksum?
            if ( !defined $affyChipsInfo{$expId}{$chipId}{'mas5Checksum'}             || $affyChipsInfo{$expId}{$chipId}{'mas5Checksum'} eq ''             || !$affyChipsInfo{$expId}{$chipId}{'mas5Checksum'} ){
                print "Warning, mas5Checksum not defined for chipId: $chipId - expId: $expId\n";
                next;
            }
            # Check that we  have the chipTypeId in chipTypeInfo
            if ( !defined $chipTypeInfo{$affyChips{$expId}{$chipId}{'chipTypeId'}} ){
                print 'Warning, no information in chipTypeCorrespondencesAndQualityThresholds ',
                      'for the chipTypeId ', $affyChips{$expId}{$chipId}{'chipTypeId'},
                      " for chipId: $chipId - expId: $expId\n";
                next;
            }
        }

        # Do we have percent present score?
        if ( !defined $affyChipsInfo{$expId}{$chipId}{'percentPresent'} || $affyChipsInfo{$expId}{$chipId}{'percentPresent'} eq '' || !$affyChipsInfo{$expId}{$chipId}{'percentPresent'} ){
            print "Warning, percent present not defined for chipId: $chipId - expId: $expId\n";
            next;
        }
        # Do we have expression based unique ID?
        if ( !defined $affyChipsInfo{$expId}{$chipId}{'expressionBasedUniqueString'} || $affyChipsInfo{$expId}{$chipId}{'expressionBasedUniqueString'} eq '' || !$affyChipsInfo{$expId}{$chipId}{'expressionBasedUniqueString'} ){
            print "Warning, expressionBasedUniqueString not defined for chipId: $chipId - expId: $expId\n";
            next;
        }
    }
}

# Now we check for duplicated files
print "Checking for duplicated files (it needs to be fixed if you find some)... \n";
displayDuplicatedChips(\%affyChipsInfo, \%affyChips, \%affyChipsInfo, 0);

# and for incompatible files
#print "\nChecking for incompatible files (just for your information, the pipeline will reject them anyway). Incompatible chips do not appear here. \n";
displayIncompatibleChips(\%affyChips, \%chipTypeInfo, 0, \%affyChipsInfo);

# and for low quality files
print "\nChecking for low quality files (just for your information, the pipeline will reject them anyway). Incompatible chips do not appear here. \n";
displayLowQualityChips(\%affyChipsInfo, \%affyChips, \%chipTypeInfo, 0, 1);

exit 0;

