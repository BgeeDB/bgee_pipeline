#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Frederic Bastian, created July 2012
#############################################################

$| = 1;
require 'bgee_utils.pl';
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;


sub getAllChipsInfo {
    my ($affyFile) = @_;
    # $affyChips{expId}{chipId}
    my %affyChips;

    open(my $IN, '<', "$affyFile")  or die "Could not read file $affyFile\n";
    my $line = <$IN>;    #header
    while ( defined ($line = <$IN>) ){

        my @tmp = map { bgeeTrim($_) }
                  split(/\t/, $line);
        my $chipId = $tmp[0];
        my $expId  = $tmp[1];

        $affyChips{$expId}{$chipId} = ();
        $affyChips{$expId}{$chipId}{'cdfName'}                     = $tmp[2]   if ( $tmp[2]  ne '' );
        $affyChips{$expId}{$chipId}{'scanDate'}                    = $tmp[3]   if ( $tmp[3]  ne '' );
        $affyChips{$expId}{$chipId}{'qualityScore'}                = $tmp[4]   if ( $tmp[4]  ne '' );
        $affyChips{$expId}{$chipId}{'percentPresent'}              = $tmp[5]   if ( $tmp[5]  ne '' );
        $affyChips{$expId}{$chipId}{'expressionBasedUniqueString'} = $tmp[6]   if ( $tmp[6]  ne '' );
        $affyChips{$expId}{$chipId}{'celDataChecksum'}             = $tmp[7]   if ( $tmp[7]  ne '' );
        $affyChips{$expId}{$chipId}{'mas5OriginalFileChecksum'}    = $tmp[8]   if ( $tmp[8]  ne '' );
        $affyChips{$expId}{$chipId}{'mas5Checksum'}                = $tmp[9]   if ( $tmp[9]  ne '' );
        $affyChips{$expId}{$chipId}{'status'}                      = $tmp[10]  if ( $tmp[10] ne '' );
    }
    close $IN;

    return %affyChips;
}

sub getAllChips {
    my ($affyFile) = @_;

    # $affyChips{expId}{chipId}
    my %affyChips;
    open(my $IN2, '<', "$affyFile")  or die "Could not read file $affyFile\n";
    my $line = <$IN2>;  # header
    while ( defined ($line = <$IN2>) ){
        # Detect commented lines,
        # but we still want the information
        my $commented = 0;
        #FIXME if Utils::read_spreadsheet is used we lose the line information per se because the columns are in hash

        # should work if the first element is quoted or not
        if (( $line =~ /^#(.+)/ ) or ($line =~ /^\"#(.+)/ )){
            $commented = 1;
            $line      = $1;
        }

        my @tmp = map { bgeeTrim($_) }
                  split(/\t/, $line);
        $affyChips{$tmp[1]}{$tmp[0]}{'commented'}           = $commented;
        $affyChips{$tmp[1]}{$tmp[0]}{'chipTypeId'}          = $tmp[2];
        $affyChips{$tmp[1]}{$tmp[0]}{'normalizationTypeId'} = $tmp[3];
        $affyChips{$tmp[1]}{$tmp[0]}{'comments'}            = $tmp[24];
    }
    close $IN2;

    return %affyChips;
}

# Sub to extract probesets from a MAS5 of Schuster Affymetrix file.
# This sub also logs some warnings if incorrectly formatted files are seen.
#
# Arguments:
#   * $experimentId: ID of the experiment
#   * $chipId: ID of the chip
#   * $normType: normalization type
#   * $processedMas5Dir: directory where processed MAS5 files are stored
#   * $processedSchusterDir: directory where processed Schuster files are stored
#   * $logWarn: whether warnings should be logged.
#
# Returns a reference to a hash defined as:
# $hashRef->{$probesetId}->{'call'}
# $hashRef->{$probesetId}->{'signal'}
# $hashRef->{$probesetId}->{'quality'}
sub extractProbesetsFromFile {
    my ($experimentId, $chipId, $normType, $processedMas5Dir, $processedSchusterDir, $logWarn) = @_;

    my %pbsets;
    my %correspCall = get_corresp_call(); ## mas5_utils.pl
    my $problems = 0; #count the number of problems for harmonizing the calls
    my $nbr_line = 0;
    my $nbr_null = 0;

    # mas5 normalization & detection
    if ( $normType eq 'MAS5' ) {
        my $file_name = $processedMas5Dir.'/'.$experimentId.'/'.$chipId;
        my %column_indexes = get_mas5_columns($file_name);
        if ( $column_indexes{'call'} != -1 && $column_indexes{'signal'} != -1 && $column_indexes{'probeset_id'} != -1 ) {
            open(my $INMAS52, '<', "$file_name")  or warn "\tWarning! Can't read file [$file_name] for [$chipId][$experimentId]\n";
            my $line = <$INMAS52>;   #header
            while ( defined ($line = <$INMAS52>) ) {
                if ( !is_valid_processed_mas5_line($line, \%column_indexes) ) {
                  print "\tWarning, unrecognized line: $line for [$chipId][$experimentId]\n" if ( $logWarn );
                  next;
                }
                chomp $line;
                $nbr_line++;
                my @tmp = split(/\t/, $line);
                my $call          = bgeeTrim($tmp[$column_indexes{'call'}]);
                my $signal        = bgeeTrim($tmp[$column_indexes{'signal'}]);
                my $probesetId    = bgeeTrim($tmp[$column_indexes{'probeset_id'}]);
                my $adjusted_call = bgeeTrim($tmp[$column_indexes{'adjusted_call'}]);
                my $q_value       = bgeeTrim($tmp[$column_indexes{'q_value'}]);

                # check that the call  and the signal intensity are defined
                if ( exists $correspCall{$call} && defined $signal && $signal ne 'null' ) {
                    #pbset -> quality, call and signal
                    $pbsets{$probesetId}->{'call'}          = $correspCall{$call};
                    $pbsets{$probesetId}->{'p_value'}       = 1.0; # Because no p-value for MAS5
                    $pbsets{$probesetId}->{'q_value'}       = $q_value;
                    $pbsets{$probesetId}->{'adjusted_call'} = $adjusted_call;
                    # a way to convert signal to numeric
                    if ( $signal + 0 >= 0 ) {
                        $pbsets{$probesetId}->{'signal'} = $signal;
                    }
                    else {
                        $pbsets{$probesetId}->{'signal'} = 0;
                    }
                    # mas5 quality is always set to low
                    $pbsets{$probesetId}->{'quality'} = $Utils::LOW_QUAL;
                }
                elsif ( exists $correspCall{$call} && defined $signal && $signal eq 'null' ) {
                    $nbr_null++;
                }
                else {
                    $problems++;
                }
            }
            close $INMAS52;
        }
        else {
            #print "\tWarning! Incorrectly formated file: $file_name\n"; #already printed earlier
        }
    }

    # gcRMA normalization / Schuster et al. detection
    # /!\ the order of the columns is inverted !
    elsif ( $normType eq 'gcRMA' ) {
        (open(INGCRMA2, '<', $processedSchusterDir.'/'.$experimentId.'/'.$chipId.'.cel.out') or
         open(INGCRMA2, '<', $processedSchusterDir.'/'.$experimentId.'/'.$chipId.'.CEL.out') or
         open(INGCRMA2, '<', $processedSchusterDir.'/'.$experimentId.'/'.$chipId.'.out')) or
         warn "\tWarning! Can't read data file for [$chipId][$experimentId]!\n";
        while ( defined (my $line = <INGCRMA2>) ) {
            chomp $line;
            $nbr_line++;
            my @tmp = split(/\t/, $line);
            my $probesetId    = bgeeTrim($tmp[0]);
            my $signal        = bgeeTrim($tmp[1]);
            my $p_value       = bgeeTrim($tmp[2]);
            my $call          = bgeeTrim($tmp[3]);
            my $q_value       = bgeeTrim($tmp[4]);
            my $adjusted_call = bgeeTrim($tmp[5]);

            if ( exists $correspCall{$call} && defined $signal && $signal ne 'null' ) {
                #pbset -> quality, call and signal
                $pbsets{$probesetId}->{'call'} = $correspCall{$call};
                if ( $signal + 0 >= 0 ) {
                    $pbsets{$probesetId}->{'signal'} = $signal;
                }
                else {
                    $pbsets{$probesetId}->{'signal'} = 0;
                }

                if ( $correspCall{$call} eq $Utils::PRESENT_CALL || $correspCall{$call} eq $Utils::ABSENT_CALL ) {
                    $pbsets{$probesetId}->{'quality'} = $Utils::HIGH_QUAL;
                }
                elsif ( $correspCall{$call} eq $Utils::MARGINAL_CALL ) {
                    $pbsets{$probesetId}->{'quality'} = $Utils::LOW_QUAL;
                }

                $pbsets{$probesetId}->{'p_value'}       = $p_value;
                $pbsets{$probesetId}->{'q_value'}       = $q_value;
                $pbsets{$probesetId}->{'adjusted_call'} = $adjusted_call;
            }
            elsif ( exists $correspCall{$call} && defined $signal && $signal eq 'null' ) {
                $nbr_null++;
            }
            else {
                $problems++;
            }
        }
        close INGCRMA2;
    }
    else {
        die "Unrecognized normalization type: $normType\n";
    }

    # Data are cleaner now, so should rarely meet this warning!
    if ( $logWarn ) {
        print "\tWarning! Some expression calls are not standard ($problems) for [$chipId][$experimentId]\n" if ( $problems > 0 );
        print "\tWarning! No line was apparently read for [$chipId][$experimentId]\n" if ( $nbr_line eq 0 );
        # Threshold at 2% for warning on proportion of null calls
        if ($nbr_line > 0){
            print "\tWarning! Some expression calls have too many null values (".(100*$nbr_null/$nbr_line)." %) for [$chipId][$experimentId]\n"  if ( (100*$nbr_null/$nbr_line) > 2);
        }
    }

    return \%pbsets;
}

sub displayIncompatibleChips {
    # The main hash containing normalization method ID and chipTypeId
    my ($affyChipsRef, $chipTypeInfoRef, $commentedChips, $affyChipsInfoRef) = @_;
    # A hash containing quality score and percent present thresholds for cdf names
    #   {cdfName}{'qualityScore} = ...
    #   {cdfName}{'percentPresent} = ...
    # A boolean to define whether we want to study commented chips
    #   must be equal to 1 or 0, only.
    # Optional, to compare to Marta's original status

    # To know if we have found at least one chip with quality below thresholds
    my $found      = 0;
    my $totalCount = 0;

    for my $expId ( sort keys %$affyChipsRef ){
        for my $chipId ( sort keys %{$affyChipsRef->{$expId}} ){
            next  if ( $chipId =~ /^#/ ); # Commented chips

            my $status = 'compatible';
            if ( (!$commentedChips && !defined $affyChipsRef->{$expId}->{$chipId}->{'commented'}) ||
                 (defined $affyChipsRef->{$expId}->{$chipId}->{'commented'} && $affyChipsRef->{$expId}->{$chipId}->{'commented'} != $commentedChips) ){
                next;
            }
            if ( isChipIncompatibleOrLowQuality($affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'}, undef, $chipTypeInfoRef, $affyChipsRef->{$expId}->{$chipId}->{'chipTypeId'}) ){
                if ( !$found ){
                    print "====================\nIncompatible chip found:\n====================\n";
                }
                print "ChipId: $chipId - expId: $expId\n";

                $found  = 1;
                $totalCount++;
                $status = 'incompatible';
            }

            if ( defined $affyChipsInfoRef && defined $affyChipsInfoRef->{$expId}->{$chipId}->{'status'} && $affyChipsInfoRef->{$expId}->{$chipId}->{'status'} ne '' ){
                if ( ($status eq 'incompatible' && ($affyChipsInfoRef->{$expId}->{$chipId}->{'status'} eq 'ok' ||
                                                    $affyChipsInfoRef->{$expId}->{$chipId}->{'status'} eq 'quality_unknown' ||
                                                    $affyChipsInfoRef->{$expId}->{$chipId}->{'status'} eq 'low_quality'))
                     || ($status eq 'compatible' && $affyChipsInfoRef->{$expId}->{$chipId}->{'status'} eq 'incompatible') ){
                    print "\tChip type compatibility incoherent to Marta's status: chipId: $chipId - expId: $expId\n";
                }
            }
        }
    }

    if ( !$found ){
        print "None found\n";
    }
    else {
        print "$totalCount chips found\n";
    }
}

sub displayLowQualityChips {
    # A ref to the hash containing chips and their info (generated by getAllChipsInfo or generateNewInfo)
    #   for which we want to find chips with quality scores below thresholds
    my ($affyChipsInfoRef, $affyChipsRef, $chipTypeInfoRef, $commentedChips, $incoherentStatus) = @_;
    # The main hash containing normalization method ID and chipTypeId
    # A hash containing quality score and percent present thresholds for cdf names
    #   {cdfName}{'qualityScore} = ...
    #   {cdfName}{'percentPresent} = ...
    # A boolean to define whether we want to study commented chips
    #   must be equal to 1 or 0, only.
    # Optional, to compare to Marta's status the first time

    displayLowQualityAndIncompatibleChips($affyChipsInfoRef, $affyChipsRef, $chipTypeInfoRef, $commentedChips, 0, $incoherentStatus);
}

sub displayLowQualityAndIncompatibleChips {
    # A ref to the hash containing chips and their info (generated by getAllChipsInfo or generateNewInfo)
    #   for which we want to find chips with quality scores below thresholds
    my ($affyChipsInfoRef, $affyChipsRef, $chipTypeInfoRef, $commentedChips, $incompatibleChips, $incoherentStatus) = @_;
    # The main hash containing normalization method ID and chipTypeId
    # A hash containing quality score and percent present thresholds for cdf names
    #   {cdfName}{'qualityScore} = ...
    #   {cdfName}{'percentPresent} = ...
    # A boolean to define whether we want to study commented chips
    #   must be equal to 1 or 0, only.
    # To know if we should display only low quality chips
    #   if 1, display also incompatible chips
    # Optional, to compare to Marta's status the first time

    my $totalCount = 0;

    # To know if we have found at least one chip with quality below thresholds
    my $found = 0;

    for my $expId ( sort keys %$affyChipsInfoRef ){
        for my $chipId ( sort keys %{$affyChipsInfoRef->{$expId}} ){
            next  if ( $chipId =~ /^#/ ); # Commented chips

            if ( (!$commentedChips && !defined $affyChipsRef->{$expId}->{$chipId}->{'commented'}) ||
                 (defined $affyChipsRef->{$expId}->{$chipId}->{'commented'} && $affyChipsRef->{$expId}->{$chipId}->{'commented'} != $commentedChips) ){
                next;
            }

            # Should we check only for low quality chips?
            # if $incompatibleChips is false, and chipType is incompatible, next!
            if ( !$incompatibleChips &&
                 isChipIncompatibleOrLowQuality($affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'}, undef, $chipTypeInfoRef, $affyChipsRef->{$expId}->{$chipId}->{'chipTypeId'}) ){
                next;
            }

            my $lowQualOrIncompatible = 0;
            if ( isChipIncompatibleOrLowQuality($affyChipsRef->{$expId}->{$chipId}->{'normalizationTypeId'}, $affyChipsInfoRef->{$expId}->{$chipId}, $chipTypeInfoRef, $affyChipsRef->{$expId}->{$chipId}->{'chipTypeId'}) ){
                if ( !$found ){
                    print "====================\nLow quality ";
                    if ( $incompatibleChips ){
                        print 'or incompatible ';
                    }
                    print "chip found: \n====================\n";
                }
                print "ChipId: $chipId - expId: $expId\n";
                $found = 1;
                $totalCount++;
                $lowQualOrIncompatible = 1;
            }

            if ( defined $incoherentStatus && $incoherentStatus && defined $affyChipsInfoRef->{$expId}->{$chipId}->{'status'} && $affyChipsInfoRef->{$expId}->{$chipId}->{'status'} ne '' ){
                if ( ($lowQualOrIncompatible && ($affyChipsInfoRef->{$expId}->{$chipId}->{'status'} eq 'ok' || $affyChipsInfoRef->{$expId}->{$chipId}->{'status'} eq 'quality_unknown')) ||
                     (!$lowQualOrIncompatible && ($affyChipsInfoRef->{$expId}->{$chipId}->{'status'} eq 'low_quality' || $affyChipsInfoRef->{$expId}->{$chipId}->{'status'} eq 'incompatible')) ){
                    print "\tQuality incoherent to Marta's status: chipId: $chipId - expId: $expId\n";
                }
            }
        }
    }

    if ( !$found ){
        print "None found\n";
    }
    else {
        print "$totalCount chips found\n";
    }
}


sub displayDuplicatedChips {
    # A ref to the hash containing chips and their info (generated by getAllChipsInfo or generateNewInfo)
    #   for which we want to find duplicates in the $targetAffyChipsInfoRef.
    #   $sourceAffyChipsInfoRef and $targetAffyChipsInfoRef can actually be a ref to the same hash,
    #   when we want to compare new chips to new chips, and old chips to old chips
    my ($sourceAffyChipsInfoRef, $affyChipsRef, $targetAffyChipsInfoRef, $compareToCommentedChips) = @_;
    # The main hash containing 'commented' status
    # A ref to the hash containing chips and their info (generated by getAllChipsInfo or generateNewInfo)
    #   in which we want to retrieve duplicated chips as compared to chips in $sourceAffyChipsInfoRef.
    #   $sourceAffyChipsInfoRef and $targetAffyChipsInfoRef can actually be a ref to the same hash,
    #   when we want to compare new chips to new chips, and old chips to old chips
    # A boolean to define whether the comparison should be done to commented chips
    #   must be equal to 1 or 0, only.

    # A hash to store duplicated chips, to not display them several times
    # {expId1}{chipId1}{expId2}{chipId2}
    my %identifiedDuplicatedChips;
    my @infos = ('qualityScore', 'percentPresent', 'expressionBasedUniqueString', 'celDataChecksum', 'mas5OriginalFileChecksum', 'mas5Checksum', 'scanDate');

    my $duplicatesFound = 0;
    # Iterate the source chips
    for my $sourceExpId ( sort keys %$sourceAffyChipsInfoRef ){
        for my $sourceChipId ( sort keys %{$sourceAffyChipsInfoRef->{$sourceExpId}} ){
            # Is the source chip commented? we never want them to be source chips
            next  if ( !defined $affyChipsRef->{$sourceExpId}->{$sourceChipId}->{'commented'} || $affyChipsRef->{$sourceExpId}->{$sourceChipId}->{'commented'} );

            my $duplicateRank = 0;
            # A hash to be able to compare the last target duplicate found for this source chip
            # to all previously identified target duplicates during this loop
            # {duplicateRank}{expId}{chipId}
            my %duplicatesToCurrentSourceChip = ();

            # Now, let's iterate the target chips to search for duplicates
            for my $targetExpId ( sort keys %$targetAffyChipsInfoRef ){
                for my $targetChipId ( sort keys %{$targetAffyChipsInfoRef->{$targetExpId}} ){
                    # The target and source chips can actually be the same,
                    # when we want to compare new chips to new chips, or old chips to old chips
                    next  if ( $sourceExpId eq $targetExpId && $sourceChipId eq $targetChipId );

                    # Do we want to compare to commented chips, or actually used chips?
                    if ( (!$compareToCommentedChips && !defined $affyChipsRef->{$targetExpId}->{$targetChipId}->{'commented'}) ||
                         (defined $affyChipsRef->{$targetExpId}->{$targetChipId}->{'commented'} && $affyChipsRef->{$targetExpId}->{$targetChipId}->{'commented'} != $compareToCommentedChips) ){
                        next;
                    }

                    # Do we already have displayed that these chips were duplicated?
                    if ( defined $identifiedDuplicatedChips{$sourceExpId}{$sourceChipId}{$targetExpId}{$targetChipId} || defined $identifiedDuplicatedChips{$targetExpId}{$targetChipId}{$sourceExpId}{$sourceChipId} ){
                        next;
                    }

                    # Now, try to find some equality
                    # some conditions are good enough to believe it's a duplicated chip,
                    # even if it is the only evidence
                    my $enoughConditions = 0;

                    # Other conditions are good only if not alone (quality_score, percent_present)
                    my $conditionsCount = 0;
                    for my $info ( @infos ){
                        if ( defined $sourceAffyChipsInfoRef->{$sourceExpId}->{$sourceChipId}->{$info} &&
                             defined $targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} &&
                             $sourceAffyChipsInfoRef->{$sourceExpId}->{$sourceChipId}->{$info} ne '' &&
                             $targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} ne '' ){
                            if ( $info eq 'qualityScore' || $info eq 'percentPresent' ){
                                if ( ($sourceAffyChipsInfoRef->{$sourceExpId}->{$sourceChipId}->{$info} + 0) == ($targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} + 0) ){
                                    $conditionsCount++;
                                }
                            }
                            elsif ( $sourceAffyChipsInfoRef->{$sourceExpId}->{$sourceChipId}->{$info} eq $targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} ){
                                if ( $info eq 'expressionBasedUniqueString' || $info eq 'celDataChecksum' || $info eq 'mas5OriginalFileChecksum' || $info eq 'mas5Checksum' || $info eq 'scanDate' ){
                                    $enoughConditions = 1;
                                }
                                $conditionsCount++;
                            }
                        }
                    }
                    if ( $conditionsCount < 2 && !$enoughConditions ){ # not duplicated
                        next;
                    }
                    $duplicatesFound = 1;

                    if ( $duplicateRank == 0 ){
                        print "====================\nPotential duplicates found:\n====================\n";
                        print "chip 1 - chip ID: $sourceChipId - experiment ID: $sourceExpId\n";
                    }

                    if ( $compareToCommentedChips ){
                        print 'Commented ';
                    }
                    print 'chip ', ($duplicateRank + 2),
                          " - chip ID: $targetChipId - experiment ID: $targetExpId\n";

                    print "\tEquality to chip 1 based on: ";
                    for my $info ( @infos ){
                        if ( defined $sourceAffyChipsInfoRef->{$sourceExpId}->{$sourceChipId}->{$info} &&
                             defined $targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} &&
                             $sourceAffyChipsInfoRef->{$sourceExpId}->{$sourceChipId}->{$info} ne '' &&
                             $targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} ne '' ){
                            if ( $info eq 'qualityScore' || $info eq 'percentPresent' ){
                                if ( ($sourceAffyChipsInfoRef->{$sourceExpId}->{$sourceChipId}->{$info} + 0) == ($targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} + 0) ){
                                    print "$info ";
                                }
                            }
                            elsif ( $sourceAffyChipsInfoRef->{$sourceExpId}->{$sourceChipId}->{$info} eq $targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} ){
                                print "$info ";
                            }
                        }
                    }
                    print "\n";
                    $identifiedDuplicatedChips{$sourceExpId}{$sourceChipId}{$targetExpId}{$targetChipId} = 1;
                    $identifiedDuplicatedChips{$targetExpId}{$targetChipId}{$sourceExpId}{$sourceChipId} = 1;

                    # Now we want to display the equality
                    # to all other target chips previously identified
                    # as duplicates of the current source chip
                    for ( my $i = 0; $i < $duplicateRank; $i++ ){
                        for my $iterateExpId ( sort keys %{$duplicatesToCurrentSourceChip{$i}} ){
                            for my $iterateChipId ( sort keys %{$duplicatesToCurrentSourceChip{$i}{$iterateExpId}} ){
                                if ( defined $identifiedDuplicatedChips{$iterateExpId}{$iterateChipId}{$targetExpId}{$targetChipId} ||
                                     defined $identifiedDuplicatedChips{$targetExpId}{$targetChipId}{$iterateExpId}{$iterateChipId} ){
                                    next;
                                }

                                print "\tEquality to ";
                                if ( $compareToCommentedChips ){
                                    print 'commented ';
                                }
                                print 'chip ', ($i + 2), ' based on: ';
                                my $infoFound = 0;
                                for my $info ( @infos ){
                                    if ( defined $targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} &&
                                         defined $targetAffyChipsInfoRef->{$iterateExpId}->{$iterateChipId}->{$info} &&
                                         $targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} ne '' &&
                                         $targetAffyChipsInfoRef->{$iterateExpId}->{$iterateChipId}->{$info} ne '' ){
                                        if ( $info eq 'qualityScore' || $info eq 'percentPresent' ){
                                            if ( ($targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} + 0) == ($targetAffyChipsInfoRef->{$iterateExpId}->{$iterateChipId}->{$info} + 0) ){
                                                print "$info ";
                                                $infoFound = 1;
                                            }
                                        }
                                        elsif ( $targetAffyChipsInfoRef->{$targetExpId}->{$targetChipId}->{$info} eq $targetAffyChipsInfoRef->{$iterateExpId}->{$iterateChipId}->{$info} ){
                                            print "$info ";
                                            $infoFound = 1;
                                        }
                                    }
                                }
                                $identifiedDuplicatedChips{$iterateExpId}{$iterateChipId}{$targetExpId}{$targetChipId} = 1;
                                $identifiedDuplicatedChips{$targetExpId}{$targetChipId}{$iterateExpId}{$iterateChipId} = 1;
                                if ( !$infoFound ){
                                    print 'none!';
                                }
                                print "\n";
                            }
                        }
                    }

                    $duplicatesToCurrentSourceChip{$duplicateRank}{$targetExpId}{$targetChipId} = 1;
                    $duplicateRank++;
                }
            }
        }
    }
    if ( !$duplicatesFound ){
        print "None found\n";
    }
}

sub getChipTypesInformation {
    #TODO use spreadsheet_reader()
    my ($qualityFile) = @_;
    # {cdfName}{'qualityScore'}   = ...
    # {cdfName}{'percentPresent'} = ...
    my %chipTypeInfo;

    open(my $IN3, '<', "$qualityFile")  or die "Could not read file [$qualityFile]\n";
    my $line = <$IN3>;  # header
    # E.g. cdfName    ArrayExpressChipTypeId    GEOChipTypeId    ChipTypeStatus    arIQR    PPresent
    while ( defined ($line = <$IN3>) ){
        my @tmp = map { $_ = bgeeTrim($_); $_ || '' }
                  split(/\t/, $line);
        # Use GEOChipTypeId as ChipTypeId if ArrayExpressChipTypeId is empty
        $tmp[1] = $tmp[1] ne '' ? $tmp[1]
                :                 $tmp[2];

        next  if ( $tmp[0] eq '' || $tmp[1] eq '' );

        $tmp[3] = ''  if ( $tmp[3] =~ /^incompatible$/i || $tmp[3] =~ /^false$/i || $tmp[3] =~ /^n\/?a$/i );
        $tmp[4] = ''  if ( $tmp[4] =~ /^false$/i || $tmp[4] =~ /^n\/?a$/i );
        $tmp[5] = ''  if ( $tmp[5] =~ /^false$/i || $tmp[5] =~ /^n\/?a$/i );

        # qualitythresholds based on cdf name
        $chipTypeInfo{$tmp[0]}{'status'}         = $tmp[3];
        $chipTypeInfo{$tmp[0]}{'qualityScore'}   = $tmp[4];
        $chipTypeInfo{$tmp[0]}{'percentPresent'} = $tmp[5];

        # qualitythresholds based on chipTypeId
        $chipTypeInfo{$tmp[1]}{'status'}         = $tmp[3];
        $chipTypeInfo{$tmp[1]}{'qualityScore'}   = $tmp[4];
        $chipTypeInfo{$tmp[1]}{'percentPresent'} = $tmp[5];

        # Correspondences
        $chipTypeInfo{$tmp[1]}{'correspondence'} = $tmp[0];
        $chipTypeInfo{$tmp[0]}{'correspondence'} = $tmp[1];
    }
    close $IN3;

    return %chipTypeInfo;
}

sub isChipIncompatibleOrLowQuality {
    my ($normalizationTypeId, $affyChipInfoRef, $chipTypeInfoRef, $chipTypeId) = @_;
    # Optional: chipTypeId (for MAS5 file we often don't have the cdf name.
    # cdf name and chipTypeId can be both used to get quality thresholds anyway)

    my $cdfName;
    if ( defined $chipTypeId && $chipTypeId ne '' ){
        $cdfName = $chipTypeId;
    }
    elsif ( defined $affyChipInfoRef->{'cdfName'} && $affyChipInfoRef->{'cdfName'} ne '' ){
        $cdfName = $affyChipInfoRef->{'cdfName'};
    }
    else { # commented file?
        return 0;
    }

    if ( !defined $cdfName || !defined $chipTypeInfoRef->{$cdfName} || $chipTypeInfoRef->{$cdfName}->{'status'} ne 'compatible' ){
        return 1;
    }

    # In case we just want to test the chip type (previous condition)
    if ( !defined $affyChipInfoRef || !$affyChipInfoRef ){
        return 0;
    }

    # commented file?
    if ( !defined $normalizationTypeId || ($normalizationTypeId ne '1' && $normalizationTypeId ne '2') ){
        if ( defined $affyChipInfoRef->{'qualityScore'} && $affyChipInfoRef->{'qualityScore'} ne '' && $affyChipInfoRef->{'qualityScore'} ne '0' &&
             defined $chipTypeInfoRef->{$cdfName}->{'qualityScore'} && $chipTypeInfoRef->{$cdfName}->{'qualityScore'} ne '' ){
            if ( ($affyChipInfoRef->{'qualityScore'} + 0) < ($chipTypeInfoRef->{$cdfName}->{'qualityScore'} + 0) ){
                return 1;
            }
        }

        if ( defined $affyChipInfoRef->{'percentPresent'} && $affyChipInfoRef->{'percentPresent'} ne '' && $affyChipInfoRef->{'percentPresent'} ne '0' &&
             defined $chipTypeInfoRef->{$cdfName}->{'percentPresent'} && $chipTypeInfoRef->{$cdfName}->{'percentPresent'} ne '' ){
            if ( ($affyChipInfoRef->{'percentPresent'} + 0) < ($chipTypeInfoRef->{$cdfName}->{'percentPresent'} + 0) ){
                return 1;
            }
        }
    }
    else { # CEL file
        if ( $normalizationTypeId eq '2' &&
             # if we have a quality score threshold defined for this cdf name
             defined $chipTypeInfoRef->{$cdfName}->{'qualityScore'} && $chipTypeInfoRef->{$cdfName}->{'qualityScore'} ne '' &&
             # commented file with parameters set?
             defined $affyChipInfoRef->{'qualityScore'} && $affyChipInfoRef->{'qualityScore'} ne '' ){
            # is the quality score of the current chip below the threshold?
            if ( ($affyChipInfoRef->{'qualityScore'} + 0) < ($chipTypeInfoRef->{$cdfName}->{'qualityScore'} + 0) ){
                return 1;
            }
        }
        if ( defined $chipTypeInfoRef->{$cdfName}->{'percentPresent'} && $chipTypeInfoRef->{$cdfName}->{'percentPresent'} ne '' &&
             # commented file with parameters set?
             defined $affyChipInfoRef->{'percentPresent'} && $affyChipInfoRef->{'percentPresent'} ne '' ){
            # is the quality score of the current chip below the threshold?
            if ( ($affyChipInfoRef->{'percentPresent'} + 0) < ($chipTypeInfoRef->{$cdfName}->{'percentPresent'} + 0) ){
                return 1;
            }
        }
    }

    return 0;
}

1;

