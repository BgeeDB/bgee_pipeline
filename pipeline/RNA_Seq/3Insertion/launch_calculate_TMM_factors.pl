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
use File::Path qw(make_path);
use Parallel::ForkManager;
$| = 1;

#calls files from BgeeCall
my $abundance_file = 'gene_level_abundance+calls.tsv';

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug, $parallel_jobs) = (0, 0);
my ($path_generes, $path_target, $path_processed) = ('', '', '');
my %opts = ('debug'             => \$debug,            # more verbose
            'bgee=s'            => \$bgee_connector,   # Bgee connector string
            'path_generes=s'    => \$path_generes,     # rna_seq/all_results/
            'path_target=s'     => \$path_target,      # rna_seq/bioconductor/targets_TMM/
            'path_processed=s'  => \$path_processed,   # final result dir rna_seq/processed_TMM/
            'parallel_jobs=i'   => \$parallel_jobs     # number of jobs run in parallel
           );

#FIXME should mkdir $path_target have to be done?
# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $path_generes eq '' || $path_target eq '' || $path_processed eq '' ||
    $parallel_jobs == 0){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -path_generes=\$(RNASEQALLRES) -path_target=\$(RNASEQTMMTARG) -path_processed=\$(RNASEQTMMPATH) -parallel_jobs=10  <expId>
\t-bgee             Bgee connector string
\t-path_generes     rna_seq/all_results/ directory path to read counts files for genic and intergenic features
\t-path_target      rna_seq/bioconductor/targets_TMM/ directory path for writing target files
\t-path_processed   rna_seq/processed_TMM/ directory path for writing output files
\t-parallel_jobs    number of experiments for which TMM normalization is done in parallel
\t-debug            More verbose
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


my %experiments;
my $selExpr = $dbh->prepare('SELECT t3.speciesId, t2.rnaSeqExperimentId, t1.rnaSeqLibraryId, t2.rnaSeqSequencerName, '.
                            't2.libraryType, t2.rnaSeqTechnologyName, t2.strandSelection, t2.rnaSeqTechnologyIsSingleCell, t2.strandSelection '.
                            'FROM rnaSeqLibraryAnnotatedSample AS t1 '.
                            'INNER JOIN rnaSeqLibrary AS t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId '.
                            'INNER JOIN cond AS t3 ON t1.conditionId = t3.conditionId '.
                            'WHERE t2.rnaSeqTechnologyIsSingleCell = 0 '.
                            'AND NOT EXISTS('.
                            '  SELECT 1 from rnaSeqLibraryAnnotatedSample as t4 '.
                            '  WHERE t4.rnaSeqLibraryAnnotatedSampleId = t1.rnaSeqLibraryAnnotatedSampleId '.
                            '  AND t4.tmmFactor != 1'.
                            ') '.
                            'ORDER BY speciesId, rnaSeqExperimentId, rnaSeqSequencerName, libraryType, rnaSeqTechnologyName, strandSelection');
# Output looks like:
#speciesId rnaSeqExperimentId rnaSeqLibraryId rnaSeqSequencerName          libraryType rnaSeqTechnologyName  libraryOrientation
#     6239 GSE16552           SRX006985       Illumina Genome Analyzer     single      Illumina              NA
#     6239 GSE16552           SRX006988       Illumina Genome Analyzer     single      Illumina              NA
#     6239 GSE16552           SRX006984       Illumina Genome Analyzer     single      Illumina              NA
# ...

$selExpr->execute()  or die $selExpr->errstr;
while ( my @data = $selExpr->fetchrow_array ){
    # if one experiment ID specified as argument, run the script only for this experiment
    if ( !defined $ARGV[0] || (defined $ARGV[0] && $data[1] eq $ARGV[0]) ){

        # check that the read counts were generated for this sample in this experiment
        my $fileName;

        # Store the samples
        # Experiment -> rnaSeqTechnologyName -> rnaSeqSequencerName -> libraryType -> libraryOrientation -> speciesId -> rnaSeqLibraryId = fileName
        $experiments{$data[1]}->{$data[5]}->{$data[3]}->{$data[4]}->{$data[6]}->{$data[0]}->{$data[2]}++;
    }
}
$selExpr->finish;
# close the connection before using forkManager
 $dbh->disconnect;

my $pm = new Parallel::ForkManager($parallel_jobs);
for my $experimentId ( keys %experiments ){
    # Forks and returns the pid for the child
    my $pid = $pm->start and next;

    my $bgee_thread = Utils::connect_bgee_db($bgee_connector);

    my $libMissingCounts = $bgee_thread->prepare(
        'SELECT distinct t4.geneId, t3.readsCount '.
        'FROM rnaSeqLibrary as t1 '.
        'INNER JOIN rnaSeqLibraryAnnotatedSample as t2 ON t1.rnaSeqLibraryId = t2.rnaSeqLibraryId '.
        'INNER JOIN rnaSeqLibraryAnnotatedSampleGeneResult as t3 ON t2.rnaSeqLibraryAnnotatedSampleId = t3.rnaSeqLibraryAnnotatedSampleId '.
        'INNER JOIN gene AS t4 ON t3.bgeegeneId = t4.bgeeGeneId '.
        'WHERE t1.rnaSeqLibraryId = ? '.
        'ORDER BY t4.geneId');

    for my $rnaSeqTechnologyName ( keys %{$experiments{$experimentId}} ) {
        my $rnaSeqTechnologyNameForFileOutput = $rnaSeqTechnologyName;
        $rnaSeqTechnologyNameForFileOutput =~ s/\s/\-/g; # for output file name, we need to replace spaces by _ in the technology name.

        for my $rnaSeqSequencerName ( keys %{$experiments{$experimentId}->{$rnaSeqTechnologyName}} ){
            my $rnaSeqSequencerNameForFileOutput = $rnaSeqSequencerName;
            $rnaSeqSequencerNameForFileOutput =~ s/\s/\-/g; # for output file name, we need to replace spaces by _ in the platform IDs.

            for my $libraryType ( keys %{$experiments{$experimentId}->{$rnaSeqTechnologyName}->{$rnaSeqSequencerName}} ){
                for my $libraryOrientation ( keys %{$experiments{$experimentId}->{$rnaSeqTechnologyName}->{$rnaSeqSequencerName}->{$libraryType}} ){
                    for my $speciesId ( keys %{$experiments{$experimentId}->{$rnaSeqTechnologyName}->{$rnaSeqSequencerName}->{$libraryType}->{$libraryOrientation}} ){
                        unless ( -e $path_processed ){
                            mkdir $path_processed;
                        }

                        # if there is only 1 sample in this experiment/platform/species/etc: TMM factor = 1 for this sample
                        if ( scalar keys %{$experiments{$experimentId}->{$rnaSeqTechnologyName}->{$rnaSeqSequencerName}->{$libraryType}->{$libraryOrientation}->{$speciesId}} eq 1 ){
                            for my $rnaSeqLibraryId ( keys %{$experiments{$experimentId}->{$rnaSeqTechnologyName}->{$rnaSeqSequencerName}->{$libraryType}->{$libraryOrientation}->{$speciesId}} ){
                                my $result = "$path_processed/${experimentId}_${rnaSeqTechnologyNameForFileOutput}_${rnaSeqSequencerNameForFileOutput}_${libraryType}_${libraryOrientation}_${speciesId}.tsv";
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
                            print "Launching TMM factor calculation for $experimentId / $rnaSeqSequencerName / $libraryType / $libraryOrientation / $speciesId \n";
                            # Create file .targets
                            my $target = "$path_target/${experimentId}_${rnaSeqTechnologyNameForFileOutput}_${rnaSeqSequencerNameForFileOutput}_${libraryType}_${libraryOrientation}_${speciesId}.target";
                            if ( $debug ){
                                print "Target file: $target\nrnaSeqExperimentId\trnaSeqLibraryId\tfile\n";
                                for my $rnaSeqLibraryId ( keys %{$experiments{$experimentId}->{$rnaSeqTechnologyName}->{$rnaSeqSequencerName}->{$libraryType}->{$libraryOrientation}->{$speciesId}} ){
                                    print "$experimentId\t$rnaSeqLibraryId\t", $path_generes.$rnaSeqLibraryId."/$abundance_file\n";
                                }
                                print "\n";
                            } else {
                                open(my $TARGET, '>', "$target")  or die 'Cannot open TARGET file';
                                print {$TARGET} "rnaSeqExperimentId\trnaSeqLibraryId\tfile\n";
                                for my $rnaSeqLibraryId ( keys %{$experiments{$experimentId}->{$rnaSeqTechnologyName}->{$rnaSeqSequencerName}->{$libraryType}->{$libraryOrientation}->{$speciesId}} ){
                                    my $pathAbundanceFile = "$path_generes$rnaSeqLibraryId/$abundance_file";
                                    if (! -e $pathAbundanceFile) {
                                        # if file does not exist it means the library has been inserted in the database in a previous release but has to be
                                        # used again to calculate tmmfactor because an other library from the same experiment/species/etc. has been
                                        # inserted in the current release. 
                                        # retrieve missing counts from the database
                                        my $tempCountFolder = "$path_target/temp_count_dir/$rnaSeqLibraryId";
                                        make_path($tempCountFolder);
                                        $pathAbundanceFile = "$tempCountFolder/gene_read_counts.tsv";
                                        warn "abundance file generated by the pipeline missing for library $rnaSeqLibraryId. The abundance ".
                                          "file $pathAbundanceFile is created using counts from the database in order to calculate the TMM normalization.";

                                        open(my $COUNTS, '>', $pathAbundanceFile)  or die 'Cannot open COUNTS file';
                                        print {$COUNTS} "id\tcounts\n";
                                        $libMissingCounts->execute($rnaSeqLibraryId)  or die $libMissingCounts->errstr;
                                        while ( my @data = $libMissingCounts->fetchrow_array ){
                                            print {$COUNTS} "$data[0]\t$data[1]\n";
                                        }
                                        close $COUNTS;
                                    }
                                    print {$TARGET} "$experimentId\t$rnaSeqLibraryId\t$pathAbundanceFile\n";
                                }
                                close $TARGET;

                                # Launch R script
                                my $log = "$path_processed/${experimentId}_${rnaSeqTechnologyNameForFileOutput}_${rnaSeqSequencerNameForFileOutput}_${libraryType}_${libraryOrientation}_${speciesId}.log";
                                my $cmd = "R CMD BATCH --vanilla '--args target_file_path=\"$target\" output_file_path=\"$path_processed/${experimentId}_${rnaSeqTechnologyNameForFileOutput}_${rnaSeqSequencerNameForFileOutput}_${libraryType}_${libraryOrientation}_${speciesId}.tsv\"' $FindBin::Bin/calculate_TMM_factors.R $log";
                                system($cmd)==0  or do{warn "\tsystem [$cmd] failed: $?\n"; map { system("mv $_ $_.PROB") } glob("$path_processed/experimentId*.out");};
                            }
                        }
                    }
                }
            }
        }
    }
    $libMissingCounts->finish;
    $bgee_thread->disconnect;
    $pm->finish;
}
$pm->wait_all_children;

exit 0;

