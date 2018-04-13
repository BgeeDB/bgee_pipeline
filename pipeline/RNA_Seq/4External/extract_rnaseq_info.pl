#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use File::Slurp;


my $tsv_file = $ARGV[0]  or die "Cannot read tsv file.\n$0 tsv_file\n";

print "#Paper\tPMID\tGSE\tPlatform\tExperiment\tType\tGSM\tTissue\tLife_stage\tGenotype\tSpecies\n";
for my $line ( read_file("$tsv_file", chomp => 1) ){
    next  if ( $line =~ /^Paper\t/ || $line =~/^#/ ); # Header
    # Paper    PMID    GSE    Platform    Experiment    Type    GSM    Tissue    Life_stage    Genotype    Species

    my @field = split("\t", $line);
# grep 'SRA ID: SRX' Other_species/MrExpTable.csv | grep -P '\t(N2|C. elegans wild isolate\.)' | grep 'Caenorhabditis elegans' | grep 'RNASeq'
# grep 'GEO record: GSM' Other_species/MrExpTable.csv | grep 'microarray' | grep -P '\t(N2|C. elegans wild isolate\.)' | grep 'Caenorhabditis elegans' | grep 'GPL200'
    next  if ( $field[10] ne 'Caenorhabditis elegans' ); # ONLY "Caenorhabditis elegans" species
    next  if ( $field[9] ne 'N2' && $field[9] ne 'N2 Bristol' && $field[9] ne 'C. elegans wild isolate.' ); # ONLY wild types

    next  if ( $field[5] ne 'RNASeq' );
    next  if ( $field[6] !~ /^SRA ID: SRX\d+$/ ); # ONLY experiment with GEO id

    # Keep only single tissue and single stages
    if ( $field[7] =~ /,/ ){
        warn "\tSeveral tissues: [$field[7]]\n";
        next;
    }
    if ( $field[8] =~ /,/ ){
        warn "\tSeveral stages:  [$field[8]]\n";
        next;
    }

    # Replace organ N.A. by whole organism
    if ( $field[7] eq 'N.A.' ){
        $field[7] = 'UBERON:0001062';
    }
    # Replace stage N.A. by development
    if ( $field[8] eq 'N.A.' ){
        $field[8] = 'UBERON:0000104';
    }

    $field[7] =~ s{\([^\)]+?\)}{}g; # Remove text label, keep only ontology ids
    $field[8] =~ s{\([^\)]+?\)}{}g;
    $field[6] =~ s{^SRA ID: }{};

    print join("\t", @field), "\n";
}


exit 0;

