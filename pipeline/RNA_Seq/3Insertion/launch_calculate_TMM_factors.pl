#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Roux, created Sep 2015, updated Nov 2016
# launch the calculation of TMM normalization factors
# USAGE: perl launch_calculate_TMM_factors.pl
#        perl launch_calculate_TMM_factors.pl <EXP_ID>
###################################################

use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;
$| = 1;


my $abundance_file = 'abundance_gene_level+new_tpm+new_fpkm+calls.tsv';

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug) = (0);
my ($path_generes, $path_target, $path_processed) = ('', '', '');
my %opts = ('debug'             => \$debug,            # more verbose
            'bgee=s'            => \$bgee_connector,   # Bgee connector string
            'path_generes=s'    => \$path_generes,     # rna_seq/all_results/
            'path_target=s'     => \$path_target,      # rna_seq/bioconductor/targets_TMM/
            'path_processed=s'  => \$path_processed,   # final result dir rna_seq/processed_TMM/
           );

#FIXME should mkdir $path_target have to be done?
# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $path_generes eq '' || $path_target eq '' || $path_processed eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -path_generes=\$(RNASEQALLRES) -path_target=\$(RNASEQTMMTARG) -path_processed=\$(RNASEQTMMPATH)  <expId>
\t-bgee             Bgee connector string
\t-path_generes     rna_seq/all_results/ directory path to read counts files for genic and intergenic features
\t-path_target      rna_seq/bioconductor/targets_TMM/ directory path for writing target files
\t-path_processed   rna_seq/processed_TMM/ directory path for writing output files
\t-debug            More verbose
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


my %experiments;
my $selExpr = $dbh->prepare('SELECT (SELECT t5.speciesId FROM gene AS t5 WHERE t5.bgeeGeneId = (SELECT t4.bgeeGeneId FROM rnaSeqResult AS t4 WHERE t4.rnaSeqLibraryId = t1.rnaSeqLibraryId LIMIT 1)) AS speciesId, rnaSeqExperimentId, rnaSeqLibraryId, rnaSeqPlatformId, libraryType, libraryOrientation FROM rnaSeqLibrary AS t1 ORDER BY speciesId, rnaSeqExperimentId, rnaSeqPlatformId, libraryType, libraryOrientation');
# Output looks like:
#speciesId rnaSeqExperimentId rnaSeqLibraryId rnaSeqPlatformId             libraryType libraryOrientation
#     6239 GSE16552           SRX006985       Illumina Genome Analyzer     single      NA
#     6239 GSE16552           SRX006988       Illumina Genome Analyzer     single      NA
#     6239 GSE16552           SRX006984       Illumina Genome Analyzer     single      NA
# ...

$selExpr->execute()  or die $selExpr->errstr;
while ( my @data = $selExpr->fetchrow_array ){
    # if one experiment ID specified as argument, run the script only for this experiment
    if ( !defined $ARGV[0] || (defined $ARGV[0] && $data[1] eq $ARGV[0]) ){

        # check that the read counts were generated for this sample in this experiment
        my $fileName;
        unless ( -e "$path_generes/$data[2]/$abundance_file" && -s "$path_generes/$data[2]/$abundance_file" ){
            die "Error, no processed file found for expId: [$data[1]] - libId: [$data[2]]\n";
            next;
        }

        # Store the samples
        # Experiment -> rnaSeqPlatformId -> libraryType -> libraryOrientation -> speciesId -> rnaSeqLibraryId = fileName
        $experiments{$data[1]}->{$data[3]}->{$data[4]}->{$data[5]}->{$data[0]}->{$data[2]}++;
    }
}
$selExpr->finish;
$dbh->disconnect;


for my $experimentId ( keys %experiments ){
    for my $rnaSeqPlatformId ( keys %{$experiments{$experimentId}} ){
        my $rnaSeqPlatformIdForFileOutput = $rnaSeqPlatformId;
        $rnaSeqPlatformIdForFileOutput =~ s/\s/\_/g; # for output file name, we need to replace spaces by _ in the platform IDs.

        for my $libraryType ( keys %{$experiments{$experimentId}->{$rnaSeqPlatformId}} ){
            for my $libraryOrientation ( keys %{$experiments{$experimentId}->{$rnaSeqPlatformId}->{$libraryType}} ){
                for my $speciesId ( keys %{$experiments{$experimentId}->{$rnaSeqPlatformId}->{$libraryType}->{$libraryOrientation}} ){
                    unless ( -e $path_processed ){
                        mkdir $path_processed;
                    }

                    # if there is only 1 sample in this experiment/platform/species/etc: TMM factor = 1 for this sample
                    if ( scalar keys %{$experiments{$experimentId}->{$rnaSeqPlatformId}->{$libraryType}->{$libraryOrientation}->{$speciesId}} eq 1 ){
                        for my $rnaSeqLibraryId ( keys %{$experiments{$experimentId}->{$rnaSeqPlatformId}->{$libraryType}->{$libraryOrientation}->{$speciesId}} ){
                            my $result = "$path_processed/${experimentId}_${rnaSeqPlatformIdForFileOutput}_${libraryType}_${libraryOrientation}_${speciesId}.tsv";
                            if ( $debug ){
                                print "Result file: $result\nrnaSeqExperimentId\trnaSeqLibraryId\ttmmFactor\n";
                                print "$experimentId\t$rnaSeqLibraryId\t1\n\n";
                            } else {
                                open(my $TMM, '>', $result)  or die 'Cannot open TMM file';
                                print {$TMM} "rnaSeqExperimentId\trnaSeqLibraryId\ttmmFactor\n";
                                print {$TMM} "$experimentId\t$rnaSeqLibraryId\t1\n";
                                close $TMM;
                            }
                        }
                    }
                    else {
                        print "Launching TMM factor calculation for $experimentId / $rnaSeqPlatformId / $libraryType / $libraryOrientation / $speciesId \n";
                        # Create file .targets
                        my $target = "$path_target/${experimentId}_${rnaSeqPlatformIdForFileOutput}_${libraryType}_${libraryOrientation}_${speciesId}.target";

                        if ( $debug ){
                            print "Target file: $target\nrnaSeqExperimentId\trnaSeqLibraryId\tfile\n";
                            for my $rnaSeqLibraryId ( keys %{$experiments{$experimentId}->{$rnaSeqPlatformId}->{$libraryType}->{$libraryOrientation}->{$speciesId}} ){
                                print "$experimentId\t$rnaSeqLibraryId\t", $path_generes.$rnaSeqLibraryId."/$abundance_file\n";
                            }
                            print "\n";
                        } else {
                            open(my $TARGET, '>', "$target")  or die 'Cannot open TARGET file';
                            print {$TARGET} "rnaSeqExperimentId\trnaSeqLibraryId\tfile\n";
                            for my $rnaSeqLibraryId ( keys %{$experiments{$experimentId}->{$rnaSeqPlatformId}->{$libraryType}->{$libraryOrientation}->{$speciesId}} ){
                                print {$TARGET} "$experimentId\t$rnaSeqLibraryId\t", $path_generes.$rnaSeqLibraryId."/$abundance_file\n";
                            }
                            close $TARGET;

                            # Launch R script
                            my $log = "$path_processed/${experimentId}_${rnaSeqPlatformIdForFileOutput}_${libraryType}_${libraryOrientation}_${speciesId}.log";
                            my $cmd = "R CMD BATCH --vanilla '--args target_file_path=\"$target\" output_file_path=\"$path_processed/${experimentId}_${rnaSeqPlatformIdForFileOutput}_${libraryType}_${libraryOrientation}_${speciesId}.tsv\"' $FindBin::Bin/calculate_TMM_factors.R $log";
                            system($cmd)==0  or do{warn "\tsystem [$cmd] failed: $?\n"; map { system("mv $_ $_.PROB") } glob("$path_processed/experimentId*.out");};
                        }
                    }
                }
            }
        }
    }
}

exit 0;

