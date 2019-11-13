#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use File::Slurp;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;


my $rnaseq_file_path = '/var/bgee/extra/pipeline/rna_seq/gene_results';
my $ensembl_version  = 75;

# Fake annotation files to produce
my $RNAseqExperiment = 'RNAseqExperiment';
my $RNAseqLibrary    = 'RNAseqLibrary';


# Read Marta' files + WormBase filtered TSV
my $tsv      = 'wormbase.rnaseq.tsv';
my $linfo    = 'library_info_file.txt';
my $rinfo    = 'run_info_file.txt';
my $platform = 'platforms.txt';


die "No wormbase.rnaseq.tsv file provided\n"    if ( !-e "$tsv"      || -z "$tsv" );
die "No library_info_file.txt file provided\n"  if ( !-e "$linfo"    || -z "$linfo" );
die "No run_info_file.txt file provided\n"      if ( !-e "$rinfo"    || -z "$rinfo" );
die "No platforms.txt file provided\n"          if ( !-e "$platform" || -z "$platform" );

# XXX: is the quote mode '"' when launching read_spreadsheet?

my %tsv      = %{ Utils::read_spreadsheet("$tsv",      "\t", 'csv', '', 1) };
#Paper  PMID  GSE  Platform  Experiment  Type  GSM  Tissue  Life_stage  Genotype  Species

my %linfo    = %{ Utils::read_spreadsheet("$linfo",    "\t", 'csv', '', 1) };
#exp_id  sample_id  library_id  cutoff  PP_all  PP_protein  PP_intronic  PP_intergenic  all_reads  mapped_left  mapped_right  pr_mapped_left  pr_mapped_right  min_read_length  max_read_length  library_type  library_orientation

my %rinfo    = %{ Utils::read_spreadsheet("$rinfo",    "\t", 'csv', '', 1) };
#exp_id  sample_id  library_id  run_id  checksum  size_in_bytes

my %platform = %{ Utils::read_spreadsheet("$platform", "\t", 'csv', '', 1) };
#exp_id  sample_id  library_id  platform


# Fake RNAseqLibrary file
my $to_write = "#Chip ID\tExperiment ID\tchipTypeId\torganId\torganName\tuberonId\tuberonName\tstageId\tstageName\tinfoOrgan\tinfoStage\tsampleTitle\tsampleSource\tSampleDescription\tSampleCharacteristics\torganAnnotationStatus\torganBiologicalStatus\tstageAnnotationStatus\tstageBiologicalStatus\tsex\tstrain\tTaxID\tcomment\tannotatorId\tlastModificationDate\t\n";
my %GSE;
for my $tsv_line ( 0..$#{$tsv{'GSM'}} ){
    # Find new RNA Seq file name based on this GSM
    for my $linfo_line ( 0..$#{$linfo{'library_id'}} ){
        if ( $tsv{'GSM'}[$tsv_line] eq $linfo{'library_id'}[$linfo_line] ){
            die "Problem with $linfo{'exp_id'}[$linfo_line]/$linfo{'sample_id'}[$linfo_line]\n"
                if ( !-e "$rnaseq_file_path/$linfo{'exp_id'}[$linfo_line]/$linfo{'sample_id'}[$linfo_line]_out_$ensembl_version" ||
                      -z "$rnaseq_file_path/$linfo{'exp_id'}[$linfo_line]/$linfo{'sample_id'}[$linfo_line]_out_$ensembl_version" );
            my $chipTypeId = '?????';
            for my $line ( 0..$#{$platform{'library_id'}} ){
                if ( $linfo{'sample_id'}[$linfo_line] eq $platform{'sample_id'}[$line] ){
                    $chipTypeId = $platform{'platform'}[$line];
                    last;
                }
            }
            $to_write .= "$linfo{'sample_id'}[$linfo_line]\t$linfo{'exp_id'}[$linfo_line]\t$chipTypeId\t\t\t$tsv{'Tissue'}[$tsv_line]\t\t$tsv{'Life_stage'}[$tsv_line]\t\t\t\t\t\t\t\t\t\t\t\t\t6239\tCaenorhabditis elegans\t\t\t\n";
            $GSE{ $linfo{'exp_id'}[$linfo_line] } = 1;
        }
    }
}
write_file("$RNAseqLibrary", $to_write);


# Fake RNAseqExperiment file
$to_write = "#EXPERIMENT_ID\tEXPERIMENT_NAME\tEXPERIMENT_DESCRIPTION\tEXPERIMENT_SOURCE\tEXPERIMENT_STATUS\tCOMMENT\n";
for my $gse ( sort keys %GSE ){
    my $source = $gse =~ /^GSE/ ? 'GEO'
               : $gse =~ /^SRP/ ? 'SRA'
               :                  '???';
    $to_write .= "$gse\t\t\t$source\ttotal\t\n";
}
write_file("$RNAseqExperiment", $to_write);

exit 0;

