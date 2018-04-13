#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use File::Slurp;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;


my $cel_file_path = '/var/bgee/extra/pipeline/Affymetrix/cel_data';


# Fake annotation files to produce
my $affymetrixChip       = 'affymetrixChip';
my $chipType             = 'chipType';
my $microarrayExperiment = 'microarrayExperiment';


# Read Marta' files + WormBase filtered TSV
my $tsv     = 'wormbase.affy.tsv';
my $info    = 'wormbase_info.txt';
my $extra   = 'wormbase_info_more_cel.txt';
#my $comment = 'wormbase_exp_comments.txt';

die "No wormbase.affy.tsv file provided\n"          if ( !-e "$tsv"     || -z "$tsv" );
die "No wormbase_info.txt file provided\n"          if ( !-e "$info"    || -z "$info" );
die "No wormbase_info_more_cel.txt file provided\n" if ( !-e "$extra"   || -z "$extra" );
#die "No wormbase_exp_comments.txt file provided\n"  if ( !-e "$comment" || -z "$comment" );


# XXX: is the quote mode '"' when launching read_spreadsheet?
my %tsv     = %{ Utils::read_spreadsheet("$tsv",     "\t", 'csv', '', 1) };
#Experiment Tissue  Life_stage  Sex Treatment   GSM Platform

my %info    = %{ Utils::read_spreadsheet("$info",    "\t", 'csv', '', 1) };
#GSM  GSE  ftp_adress  old_cel_name  new_cel_name

my %extra   = %{ Utils::read_spreadsheet("$extra",   "\t", 'csv', '', 1) };
#GSM  GSE  ftp_adress  old_cel_name  new_cel_name

#my %comment = %{ Utils::read_spreadsheet("$comment", "\t", 'csv', '', 1) };
#sample_id  new_exp_id  WarmBase_exp  GEO_exp  comment



# Fake affymetrixChip file
my $to_write = "#Chip ID\tExperiment ID\tchipTypeId\tnormalizationTypeId\tdetectionTypeId\tanatEntityId\torganName\tUberonId\tUberonName\tstageId\n";
my $chipTypeId          = 'A-AFFY-60'; # https://www.ebi.ac.uk/arrayexpress/arrays/A-AFFY-60/?keywords=elegans
my $normalizationTypeId = 2; # gcRMA    only
my $detectionTypeId     = 2; # Schuster only
my %GSE;
for my $tsv_line ( 0..$#{$tsv{'GSM'}} ){
    # Find new cel file name based on this GSM
    for my $info_line ( 0..$#{$info{'GSM'}} ){
        if ( $tsv{'GSM'}[$tsv_line] eq $info{'GSM'}[$info_line] ){
            my $cel = $info{'new_cel_name'}[$info_line];
            $cel   =~ s{\.gz$}{};
            die "Problem with $info{'GSE'}[$info_line]/$cel\n"  if ( !-e "$cel_file_path/$info{'GSE'}[$info_line]/$cel" || -z "$cel_file_path/$info{'GSE'}[$info_line]/$cel" );
            $to_write .= "$cel\t$info{'GSE'}[$info_line]\t$chipTypeId\t$normalizationTypeId\t$detectionTypeId\t\t\t$tsv{'Tissue'}[$tsv_line]\t\t$tsv{'Life_stage'}[$tsv_line]\n";
            $GSE{ $info{'GSE'}[$info_line] } = 1;
        }
    }
    for my $extra_line ( 0..$#{$extra{'GSM'}} ){
        if ( $tsv{'GSM'}[$tsv_line] eq $extra{'GSM'}[$extra_line] ){
            my $cel = $extra{'new_cel_name'}[$extra_line];
            $cel   =~ s{\.gz$}{};
            die "Problem with $extra{'GSE'}[$extra_line]/$cel\n"  if ( !-e "$cel_file_path/$extra{'GSE'}[$extra_line]/$cel" || -z "$cel_file_path/$extra{'GSE'}[$extra_line]/$cel" );
            $to_write .= "$cel\t$extra{'GSE'}[$extra_line]\t$chipTypeId\t$normalizationTypeId\t$detectionTypeId\t\t\t$tsv{'Tissue'}[$tsv_line]\t\t$tsv{'Life_stage'}[$tsv_line]\n";
            $GSE{ $extra{'GSE'}[$extra_line] } = 1;
        }
    }
}
write_file("$affymetrixChip", $to_write);


# Fake chipType file
$to_write  = "#chipTypeId\tchipTypeName\tspeciesId\tEnsembl_Xref_table\tComments\t\n";
$to_write .= "$chipTypeId\tAffymetrix GeneChip C. elegans Genome Array [Celegans]\t6239\tC_elegans\t\t\n";
write_file("$chipType", $to_write);


# Fake microarrayExperiment file
$to_write = "#experimentId\texperimentName\texperimentDescription\texperimentSource\texperimentStatus\tcomment\n";
for my $gse ( sort keys %GSE ){
    $to_write .= "$gse\t\t\tGEO\tcomplete\t\n";
}
write_file("$microarrayExperiment", $to_write);

exit 0;

