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

my ($all_results, $library_info, $excluded_libraries, $length_info, $tx2gene_dir) = ('', '', '', '', '');
my %opts = ('library_info=s'        => \$library_info,       # rna_seq_sample_info.txt file
            'excluded_libraries=s'  => \$excluded_libraries, # rna_seq_sample_excluded.txt file
            'tx2gene_dir=s'         => \$tx2gene_dir,        # directory containing transcript to gene mapping file for all species
            'all_results=s'         => \$all_results,        # /var/bgee/extra/pipeline/rna_seq/all_results_bgee_v15/
            'length_info=s'         => \$length_info,        # rna_seq_length_info.txt file
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $library_info eq '' || $excluded_libraries eq '' || $tx2gene_dir eq '' || $all_results eq '' || $length_info eq '' ){
    print "\n\tInvalid or missing argument:
\te.g., $0 -library_info=\$(RNASEQ_SAMPINFO_FILEPATH) -excluded_libraries=\$(RNASEQ_SAMPEXCLUDED_FILEPATH) -tx2gene_dir=\$(OUTPUT_DIR) -all_results=\$(RNASEQ_VITALIT_ALL_RES) -length_info=\$(RNASEQ_LENGTH_FILEPATH) > $@.tmp 2>warnings.$@
\t-library_info        rna_seq_sample_info.txt file
\t-excluded_libraries  rna_seq_sample_excluded.txt file
\t-tx2gene_dir         directory containing transcript to gene mapping file for all species
\t-all_results         all_results directory
\t-length_info         path to exported file
\n";
    exit 1;
}

# Library info used to launch the pipeline
my %libraries         = getAllRnaSeqLibrariesInfo($library_info);
# Excluded libraries (after mapping step)
my %excludedLibraries = getExcludedLibraries($excluded_libraries);


# Chose one random library per species and record the lengths
my %all_species;
foreach my $expId ( sort keys %libraries ){
    foreach my $libraryId ( sort keys %{$libraries{$expId}} ){
        next  if ( exists($excludedLibraries{$libraryId}) );

        # use only for first library found for each species
        unless ( exists $all_species{$libraries{$expId}->{$libraryId}->{'speciesId'}} ){
        	my $speciesId = $libraries{$expId}->{$libraryId}->{'speciesId'};
            print 'Recording length for species ', $speciesId, " [$libraries{$expId}->{$libraryId}->{'organism'}]\n";

            ## retrieve mapping between transcript ids and gene ids (without intergenic regions)
            my $genomeFilePath = $libraries{$expId}->{$libraryId}->{'genomeFilePath'};

            my %tx2gene_mapping;
            $genomeFilePath =~ m/.+\/(.+)/;
            open(my $IN_TX2GENE, '<', $tx2gene_dir.'/'.$1.'.tx2gene') or die "Could not read tx2gene file\n";

            my $line_tx2gene = <$IN_TX2GENE>;    #header
            while ( defined ($line_tx2gene = <$IN_TX2GENE>) ){
            	chomp $line_tx2gene;
                my @column = split(/\t/, $line_tx2gene);
                # intergenic regions have same transcript id and gene id
                if ($column[0] ne $column[1]) {
                	$tx2gene_mapping{$column[0]} = $column[1];
                }
            }
            close $IN_TX2GENE;

			## open kallisto abundance file to retrieve transcript id and corresponding read length
			my $kallisto_file = $all_results.'/'.$libraryId.'/abundance.tsv';
			unless ( -s $kallisto_file ){
                die "Missing or empty processed data file ($kallisto_file) for library $libraryId! This library should maybe be added to the file of excluded libraries?\n";
            }
            open(my $IN, '<', $kallisto_file)  or die "Could not read file [$kallisto_file]\n";
            my $line = <$IN>;    #header
            while ( defined ($line = <$IN>) ){
                chomp $line;
                # file format: target_id    gene_id    length    eff_length    est_counts    tpm    fpkm    biotype
                my @tmp = map { bgeeTrim($_) } split(/\t/, $line);

                my $transcriptId = $tmp[0];
                my $geneId       = $tx2gene_mapping{$transcriptId};
                my $length       = $tmp[2];

                # record the length
                if (!exists($tx2gene_mapping{$transcriptId})) {
                    $all_species{$speciesId}->{$transcriptId}->{'geneId'} = $geneId;
                    $all_species{$speciesId}->{$transcriptId}->{'length'} = $length;
                }
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

