#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use Sort::Naturally;
use FindBin;

require("$FindBin::Bin/rna_seq_utils.pl");
$| = 1; # no buffering of output

# Julien Roux, created Nov 2016
# This simple script exports the length of transcripts from all species
# into a single file to be sent to our servers for insertion into the
# database.
#####################################################################

my ($all_results, $library_info, $excluded_libraries, $length_info) = ('', '', '', '');
my %opts = ('library_info=s'        => \$library_info,       # rna_seq_sample_info.txt file
            'excluded_libraries=s'  => \$excluded_libraries, # rna_seq_sample_excluded.txt file
            'all_results=s'         => \$all_results,        # /var/bgee/extra/pipeline/rna_seq/all_results_bgee_v14/
            'length_info=s'         => \$length_info,        # rna_seq_length_info.txt file
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $library_info eq '' || $excluded_libraries eq '' || $all_results eq '' || $length_info eq '' ){
    print "\n\tInvalid or missing argument:
\te.g., $0 -library_info=\$(RNASEQ_SAMPINFO_FILEPATH) -excluded_libraries=\$(RNASEQ_SAMPEXCLUDED_FILEPATH) -all_results=\$(RNASEQ_CLUSTER_ALL_RES) -length_info=\$(RNASEQ_LENGTH_FILEPATH) > $@.tmp 2>warnings.$@
\t-library_info        rna_seq_sample_info.txt file
\t-excluded_libraries  rna_seq_sample_excluded.txt file
\t-all_results         all_results directory
\t-length_info         path to exported file
\n";
    exit 1;
}

# Library info used to launch the pipeline
my %libraries         = getAllRnaSeqLibrariesInfo($library_info);
# Excluded libraries (after mapping step)
my %excludedLibraries = getExcludedLibraries($excluded_libraries);


# Chose one randome library per species and record the lengths
my %all_species;
foreach my $expId ( sort keys %libraries ){
    foreach my $libraryId ( sort keys %{$libraries{$expId}} ){
        next  if ( exists($excludedLibraries{$libraryId}) );

        # use only for first library found for each species
        unless ( exists $all_species{$libraries{$expId}->{$libraryId}->{'speciesId'}} ){
            print 'Recording length for species ', $libraries{$expId}->{$libraryId}->{'speciesId'}, " [$libraries{$expId}->{$libraryId}->{'organism'}]\n";

            # check if data file exists
            my $length_file = $all_results.'/'.$libraryId.'/abundance+gene_id+new_genic_tpm+new_genic_fpkm.tsv';
            unless ( -s $length_file ){
                die "Missing or empty processed data file ($length_file) for library $libraryId! This library should maybe be added to the file of excluded libraries?\n";
            }

            open(my $IN, '<', $length_file)  or die "Could not read file [$length_file]\n";
            my $line = <$IN>;    #header
            while ( defined ($line = <$IN>) ){
                chomp $line;
                # file format: target_id    gene_id    length    eff_length    est_counts    tpm    fpkm    biotype
                my @tmp = map { bgeeTrim($_) } split(/\t/, $line);
                my $transcriptId = $tmp[0];
                my $geneId       = $tmp[1];
                my $length       = $tmp[2];

                # record the length
                $all_species{$libraries{$expId}->{$libraryId}->{'speciesId'}}->{$transcriptId}->{'geneId'} = $geneId;
                $all_species{$libraries{$expId}->{$libraryId}->{'speciesId'}}->{$transcriptId}->{'length'} = $length;
            }
            close $IN;
        }
    }
}
print "\n";

#################
# EXPORT LENGTH #
#################
open (my $OUT, '>', $length_info)  or die "Cannot write [$length_info]\n";
print {$OUT} join("\t", '#taxid', 'transcript', 'geneId', 'transcript_length'), "\n";
foreach my $species ( nsort keys %all_species ){
    print 'Exporting length for species ', $species, "\n";
    foreach my $transcript ( sort keys %{$all_species{$species}} ){
        print {$OUT} $species, "\t", $transcript, "\t", $all_species{$species}->{$transcript}->{'geneId'}, "\t", $all_species{$species}->{$transcript}->{'length'}, "\n";
    }
}
close $OUT;
exit 0;

