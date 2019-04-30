#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
require('rna_seq_utils.pl');
$| = 1; # no buffering of output
# Julien Roux, Nov 2016

# script aimed at parsing the .report files in each folder and collect the % reads aligned and read length infos
# Lines in .report files look like:
#  - "Kallisto pseudo-aligned 7063465 reads out of 8564610 (82.472709968141%)"
#  - "Minimum read length across runs: XXXX"
#  - "Maximum read length across runs: XXXX"

# Output in tab-delimited file reports_info_all_samples.txt
#  - columns libraryId, allReadsCount, mappedReadsCount, minReadLength, maxReadLength
#  - if no info: NULL
#####################################################################

# Define arguments & their default value
my ($all_results, $library_info, $excluded_libraries, $report_info) = ('', '', '', '');
my %opts = ('all_results=s'         => \$all_results,        # /var/bgee/extra/pipeline/rna_seq/all_results_bgee_v14/
            'library_info=s'        => \$library_info,       # rna_seq_sample_info.txt file
            'excluded_libraries=s'  => \$excluded_libraries, # rna_seq_sample_excluded.txt file
            'report_info=s'         => \$report_info,        # reports_info_all_samples.txt
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $library_info eq ''  || $excluded_libraries eq '' || $report_info eq '' || $all_results eq '' ){
    print "\n\tInvalid or missing argument:
\te.g., $0 -library_info=\$(RNASEQ_SAMPINFO_FILEPATH) -excluded_libraries=\$(RNASEQ_SAMPEXCLUDED_FILEPATH) -report_info=\$(RNASEQREPORTINFO) -all_results=\$(RNASEQALLRES) > $@.tmp 2>warnings.$@
\t-library_info        rna_seq_sample_info.txt file
\t-excluded_libraries  rna_seq_sample_excluded.txt file
\t-report_info         reports_info_all_samples.txt
\t-all_results         all_results directory
\n";
    exit 1;
}

#####################################################################

# Library info used to launch the pipeline
my %libraries         = getAllRnaSeqLibrariesInfo($library_info);
print "\t", scalar keys %libraries, " experiments with libraries mapped.\n";

# Excluded libraries (after mapping step)
my %excludedLibraries = getExcludedLibraries($excluded_libraries);
print "\t", scalar keys %excludedLibraries, " libraries excluded.\n";

my $count_libs = 0;
foreach my $expId ( sort keys %libraries ){
    foreach my $libraryId ( sort keys %{$libraries{$expId}} ){
        if ( exists($excludedLibraries{$libraryId}) ){
            delete $libraries{$expId}->{$libraryId};
        } else {
            $count_libs++;
        }
    }
}
print "\t", $count_libs, " libraries mapped and to be inserted.\n\n";


# Read the .report file for each library, extract infos and print them out
#NOTE $report_info is rewritten everytime this script runs!
open (my $OUT, '>', $report_info)  or die "Cannot write [$report_info]\n";
# header
print {$OUT} "#libraryId\tallReadsCount\tmappedReadsCount\tminReadLength\tmaxReadLength\n";

foreach my $expId ( sort keys %libraries ){
    LIBRARY:
    foreach my $libraryId ( sort keys %{$libraries{$expId}} ){
        print "\t$expId $libraryId\n";
        my ($allReadsCount, $mappedReadsCount, $minReadLength, $maxReadLength) = ('NULL', 'NULL', 'NULL', 'NULL');

        # test if .report file exists
        if ( -s $all_results.'/'.$libraryId.'/'.$libraryId.'.report' ){
            open(my $IN, '<', $all_results.'/'.$libraryId.'/'.$libraryId.'.report') or die "could not read .report file\n";
            while ( defined (my $line = <$IN>) ){
                if ( $line =~ m/Kallisto\spseudo\-aligned\s(\d+)\sreads\sout\sof\s(\d+)\s/ ){
                    $allReadsCount    = $2;
                    $mappedReadsCount = $1;
                }
                elsif ( $line =~ m/Minimum read length across runs\:\s(\d+)/ ){
                    $minReadLength = $1;
                }
                elsif ( $line =~ m/Maximum read length across runs\:\s(\d+)/ ){
                    $maxReadLength = $1;
                }
            }
            close $IN;

            # output infos
            print {$OUT} "$libraryId\t$allReadsCount\t$mappedReadsCount\t$minReadLength\t$maxReadLength\n";
        }
        else {
            print "Missing .report file for library $libraryId!\n";
        }
    }
}
close $OUT;
exit;
# Warning: if the info is several time in the report file (e.g., if sample was launched several time), only the last occurence will be considered

