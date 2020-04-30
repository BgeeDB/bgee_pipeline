#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Frederic Bastian, created June 2012
# USAGE: just a function file, called by other scripts
#############################################################

$| = 1;

require('bgee_utils.pl');
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

sub get_corresp_call {
    my %corresp_call = (
            'present'   => $Utils::PRESENT_CALL,
            'absent'    => $Utils::ABSENT_CALL,
            'marginal'  => $Utils::MARGINAL_CALL,
            'undefined' => $Utils::UNDEFINED_CALL,
            'p'         => $Utils::PRESENT_CALL,
            'a'         => $Utils::ABSENT_CALL,
            'm'         => $Utils::MARGINAL_CALL,
            'u'         => $Utils::UNDEFINED_CALL,
            'Present'   => $Utils::PRESENT_CALL,
            'Absent'    => $Utils::ABSENT_CALL,
            'Marginal'  => $Utils::MARGINAL_CALL,
            'Undefined' => $Utils::UNDEFINED_CALL,
            'P'         => $Utils::PRESENT_CALL,
            'A'         => $Utils::ABSENT_CALL,
            'M'         => $Utils::MARGINAL_CALL,
            'U'         => $Utils::UNDEFINED_CALL,
            'RP'        => $Utils::UNDEFINED_CALL,
            );

    return %corresp_call;
}

# function to find the column containing mas5 calls and signals,
# based on header names
sub get_mas5_columns_from_header {
    my ($header) = @_;
    $header = bgeeTrim($header);

    my %proper_columns;
    $proper_columns{'call'}        = -1;
    $proper_columns{'signal'}      = -1;
    $proper_columns{'probeset_id'} = -1;

    my @tmp = split(/\t/, $header);
    my $error = 0;

    # first, let's try to identify proper columns with real proper names
    for ( my $y = 0; $y <= $#tmp; $y++ ){
        $tmp[$y] = bgeeTrim($tmp[$y]);
        if ( $tmp[$y] =~ /GEO:AFFYMETRIX(_|\s)VALUE$/ || $tmp[$y] =~ /Affymetrix:CHPSignal$/){
            if ( $proper_columns{'signal'} != -1 && $proper_columns{'signal'} != $y ){
                $error = 1;
            }
            $proper_columns{'signal'} = $y;
        }
        elsif ( $tmp[$y] =~ /GEO:AFFYMETRIX(_ABS)?(_|\s)CALL$/ || $tmp[$y] =~ /Affymetrix:CHPDetection$/ || $tmp[$y] =~ /GEO:AFFYMETRIX_Detection$/){
            if ( $proper_columns{'call'} != -1 && $proper_columns{'call'} != $y ){
                $error = 1;
            }
            $proper_columns{'call'} = $y;
        }
        elsif ( $tmp[$y] eq 'CompositeSequence name' || $tmp[$y] eq 'CompositeSequence_name' || $tmp[$y] eq 'Affymetrix:CHPProbeSetName'){
            if ( $proper_columns{'probeset_id'} != -1 && $proper_columns{'probeset_id'} != $y ){
                $error = 1;
            }
            $proper_columns{'probeset_id'} = $y;
        }
    }
    if ( $error != 0 ){
        $proper_columns{'call'}        = -1;
        $proper_columns{'signal'}      = -1;
        $proper_columns{'probeset_id'} = -1;
    }

    # then, if we could not do it with real proper names...
    if ( $proper_columns{'call'} == -1 || $proper_columns{'signal'} == -1 || $proper_columns{'probeset_id'} == -1 ){
        # let's try with more relaxed constraints on names,
        # but more on columns order and total count
        if ( $tmp[0] eq 'ID_REF' && $tmp[1] eq 'VALUE' && ($tmp[2] eq 'ABS_CALL' || $tmp[2] eq 'ABS CALL' || $tmp[2] eq 'Detection') ){
            $proper_columns{'probeset_id'} = 0;
            $proper_columns{'call'}        = 2;
            $proper_columns{'signal'}      = 1;
        }
        if ( $tmp[0] eq 'ID_REF' && $tmp[2] eq 'VALUE' && ($tmp[1] eq 'ABS_CALL' || $tmp[1] eq 'ABS CALL' || $tmp[1] eq 'Detection') ){
            $proper_columns{'probeset_id'} = 0;
            $proper_columns{'call'}        = 1;
            $proper_columns{'signal'}      = 2;
        }
        elsif ( $tmp[0] eq 'Name' && $tmp[1] =~ /GEO:AFFYMETRIX_ABS_CALL$/ && $tmp[2] =~ /GEO:AFFYMETRIX_VALUE$/ ){
            $proper_columns{'probeset_id'} = 0;
            $proper_columns{'call'}        = 1;
            $proper_columns{'signal'}      = 2;
        }
        $proper_columns{'q_value'}       = 3;
        $proper_columns{'adjusted_call'} = 4;
    }

    return %proper_columns;
}

# function to find the column containing mas5 calls and signals.
# it's not always in the same columns
sub get_mas5_columns {
    my ($mas5_file) = @_;

    my %corresp_call = get_corresp_call();
    open(my $IN, '<', "$mas5_file")  or do {print "\tWarning! Can't read file to find the column containing mas5 calls!\n"; return;};
    my $i = 0;
    my %proper_columns;
    $proper_columns{'call'}        = -1;
    $proper_columns{'signal'}      = -1;
    $proper_columns{'probeset_id'} = -1;
    my $error = 0;

    my $firstLine = 1;
    while ( defined (my $line = <$IN>) && $i < 500 && $error == 0 ){
        $line =~ s/(\r\n|\r)$/\n/;
        chomp $line;

        if ( $firstLine == 1 ){
            if ( $line =~ /^#/ ){
                # it is an additional header line, no problem
                # but it is still not the header, we do not set $firstLine to 0
            }
            else { # it is the header
                $firstLine = 0;
                # try to get proper columns from the header
                %proper_columns = get_mas5_columns_from_header($line);
                # if we were successfull in getting the columns from the header
                if ( $proper_columns{'call'} != -1 && $proper_columns{'signal'} != -1 && $proper_columns{'probeset_id'} != -1 ){
                    # we can stop here
                    return %proper_columns;
                }
                else {
                    # print a warning because the header was weird,
                    # the file will then be inspected to try to auto-detect columns
                    warn "\t\tWarning, unrecognized header, in $mas5_file: $line\n";
                }
            }
            next;
        }

        # If we reach this part of the code, it means that the current line is not a header ($firstLine = 0),
        # and that we were not able to detect the proper columns from the header.
        # We try to auto-detect the columns here
        # We check the 500 first lines, just to be sure (see the while loop)
        my @tmp = split(/\t/, $line);

        #trying to find the proper columns
        for ( my $y = 0; $y <= $#tmp; $y++ ){
            $tmp[$y] = bgeeTrim($tmp[$y]);
            # Check that the column corresponds to a mas5 call
            if ( exists $corresp_call{$tmp[$y]} ){
                if ( $proper_columns{'call'} != -1 && $proper_columns{'call'} != $y ){
                    $error = 1;
                }
                $proper_columns{'call'} = $y;
                # In some files, the p-value is provided. We have to distinguish signals and p-values then
                # the p-value is always lesser than 1. Signals are most of the time greater. So it is fair to assume than in 500 lines,
                # we should see in the signal column a number greater than 1 at least once.
            }
            elsif ( $tmp[$y] =~ /^[\+-]*[0-9]*\.[0-9]+$/ && $tmp[$y] !~ /^[\. ]*$/ && ($tmp[$y] + 0) > 1 ){
                if ( $proper_columns{'signal'} != -1 && $proper_columns{'signal'} != $y ){
                    $error = 1;
                }
                $proper_columns{'signal'} = $y;
            }
            elsif ( $tmp[$y] =~ /^[\w\-]+?_(at|st)$/ ){
                if ( $proper_columns{'probeset_id'} != -1 && $proper_columns{'probeset_id'} != $y ){
                    $error = 1;
                }
                $proper_columns{'probeset_id'} = $y;
            }
        }
        $i++;
    }
    close $IN;

    if ( $proper_columns{'call'} == -1 || $proper_columns{'signal'} == -1 || $proper_columns{'probeset_id'} == -1 ){
        $error = 1;
    }
    if ( $error != 0 ){
        $proper_columns{'call'}        = -1;
        $proper_columns{'signal'}      = -1;
        $proper_columns{'probeset_id'} = -1;
    }
    return %proper_columns;
}

# Detect whether this line in a processed_mas5 file is valid,
# i.e. not a header, not a comment
sub is_valid_processed_mas5_line {
    my ($line, $column_index_ref) = @_;
    my %corresp_call     = get_corresp_call();
    $line =~ s/(\r\n|\r)$/\n/;
    chomp $line;

    # comment line
    if ( $line =~ /^#/ ){
        return 0;
    }
    my @tmp = split(/\t/, $line);
    $tmp[$column_index_ref->{'call'}]   = bgeeTrim($tmp[$column_index_ref->{'call'}]);
    $tmp[$column_index_ref->{'signal'}] = bgeeTrim($tmp[$column_index_ref->{'signal'}]);

    # Check for the presence of a mas5 call
    # if not present, most likely a header
    if ( !exists $corresp_call{$tmp[$column_index_ref->{'call'}]} ){
        return 0
    }

    # Check for the presence of a signal column
    # if not present, most likely a header, or a 'null' signal (usual)
    if ( !($tmp[$column_index_ref->{'signal'}] =~ /^[\+-]*\d*\.*\d*$/ && $tmp[$column_index_ref->{'signal'}] !~ /^[\. ]*$/) && $tmp[$column_index_ref->{'signal'}] ne 'null' ){
        return 0;
    }
    return 1;
}

sub is_a_valid_filtered_file {
    my ($file) = @_;

    open(my $IN2, '<', "$file")  or do {print "\nCan't read file $file!\n"; return;};
    my $line = <$IN2>;  # header
    my $valid = 1;
    while (defined ($line = <$IN2>)) {
        if ( !is_valid_filtered_processed_mas5_line($line) ){
            $valid = 0;
            last;
        }
    }
    close $IN2;
    return $valid;
}

sub is_valid_filtered_processed_mas5_line {
    my ($line) = @_;
    chomp $line;

    # comment line
    if ( $line =~ /^#/ ){
        return 0;
    }

    my @tmp = map { bgeeTrim($_) }
              split(/\t/, $line);
    if ( scalar(@tmp) != 3 ){
        return 0;
    }

    # Check for the presence of a mas5 call
    # we don't use global vars from Utils.pm, since these values are MAS5 specific, not Bgee specific.
    # XXX: but why not using sub get_corresp_call then?
    if ( $tmp[1] ne 'present' && $tmp[1] ne 'absent' && $tmp[1] ne 'marginal' && $tmp[1] ne 'undefined' ){
        return 0
    }

    # Check for the presence of a signal column
    if ( !($tmp[2] =~ /^[\+-]*\d*\.*\d*$/ && $tmp[2] !~ /^[\. ]*$/) && $tmp[2] ne 'null' ){
        return 0;
    }
    return 1;
}

