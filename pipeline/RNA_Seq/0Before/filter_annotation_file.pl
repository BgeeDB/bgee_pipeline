#!/usr/bin/env perl
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use File::Slurp;
use Data::Dumper;


# Define arguments & their default value
my ($RNAlib)            = ('');
my ($RNAlibFiltered)    = ('');
my ($minConditionNbr)   = 6;
my ($debug)             = (0);
my %opts = ('RNAlib=s'            => \$RNAlib,
            'RNAlibFiltered=s'    => \$RNAlibFiltered,
            'minConditionNbr=i'   => \$minConditionNbr,
            'debug'               => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $RNAlib eq '' || $RNAlibFiltered eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -RNAlib=\$(RNASEQ_LIB_FILEPATH_FULL).ori -RNAlibFiltered=\$(RNASEQ_LIB_FILEPATH_FULL) -minConditionNbr=$minConditionNbr
\t-RNAlib              RNASeqLibrary_full.tsv before filtering
\t-RNAlibFiltered      RNASeqLibrary_full.tsv that will be used, after condition number filtering
\t-minConditionNbr     Minimal number of unique conditions per species
\t-debug               more verbose output
\n";
    exit 1;
}

# Read RNASeqLibrary_full.tsv before filtering
my $header = '';
my $store;
ANNOT:
for my $annot ( read_file("$RNAlib", chomp=>1) ){
    #"#libraryId"	"experimentId"	"platform"	"SRSId"	"anatId"	"anatName"	"stageId"	"stageName"	"url_GSM"	"infoOrgan"	"infoStage"	"anatAnnotationStatus"	"anatBiologicalStatus"	"stageAnnotationStatus"	"sex"	"strain"	"genotype"	"speciesId"	"protocol"	"protocolType"	"RNASelection"	"globin_reduction"	"replicate"	"sampleTitle"	"PATOid"	"PATOname"	"comment"	"condition"	"annotatorId"	"lastModificationDate"
    #"SRX5028084"	"SRP169832"	"Illumina HiSeq 2500"	"SRS4059375"	"UBERON:0001046"	"hindgut"	"SsalDv:0000065"	"parr stage"		"Hindgut"	"1 yr"	"perfect match"	"partial sampling"	"missing child term"	"NA"	"Aqua Gen strain"		8030			"polyA"			"Control hindgut 2"			"PRJNA506138,1 year, parr,Ctrl_Hind_2,Fragments are taken to do histological analyses"		"WAH"	"29/03/2023"

    if ( $annot =~ /^"#/ && $header eq '' ){
        $header = $annot;
        next ANNOT;
    }
    next ANNOT  if ( $annot =~ /^#/ || $annot =~ /^"#/ );
    my ($libraryId, $experimentId, $platform, undef, $uberonId, undef, $stageId, undef, undef,
        undef, undef, undef, undef, undef, $sex, $strain, undef, $speciesId)
        = split(/\t/, $annot, 30);

    push @{ $store->{$speciesId}->{'lib'} }, $annot;
    $store->{$speciesId}->{'cond'}->{"$uberonId--$stageId--$sex--$strain"}++;
    $store->{$speciesId}->{'organ'}->{$uberonId}++;
    $store->{$speciesId}->{'stage'}->{$stageId}++;
    $store->{$speciesId}->{'sex'}->{$sex}++;
    $store->{$speciesId}->{'strain'}->{$strain}++;
    $store->{$speciesId}->{'plateform'}->{$platform}++;
}


# Stats
print join("\t", '#speciesId', 'conditions', 'organs', 'stages', 'sexes', 'strains', 'plateforms'), "\n";
for my $sp ( sort {scalar keys %{ $store->{$a}->{'cond'} } <=> scalar keys %{ $store->{$b}->{'cond'} }} keys %$store ){
    print join("\t", $sp,
                     scalar keys %{ $store->{$sp}->{'cond'} },
                     scalar keys %{ $store->{$sp}->{'organ'} },
                     scalar keys %{ $store->{$sp}->{'stage'} },
                     scalar keys %{ $store->{$sp}->{'sex'} },
                     scalar keys %{ $store->{$sp}->{'strain'} },
                     scalar keys %{ $store->{$sp}->{'plateform'} },
              ), "\n";
}


# Print filtered RNASeqLibrary_full.tsv
write_file("$RNAlibFiltered", $header."\n");
for my $sp ( grep { scalar keys %{ $store->{$_}->{'cond'} } >= $minConditionNbr } sort {$a <=> $b} keys %$store ){
    for my $line ( sort @{ $store->{$sp}->{'lib'} } ){
        append_file("$RNAlibFiltered", $line, "\n");
    }
}

exit 0;

