#!/usr/bin/env perl
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use File::Slurp;


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
    #"#libraryId"    "experimentId"  "platform"  "organId"   "organName" "uberonId"  "uberonName"    "stageId"   "stageName" "infoOrgan" "infoStage" "sampleTitle"   "sampleSource"  "sampleDescription" "sampleCharacteristics" "organAnnotationStatus" "organBiologicalStatus" "stageAnnotationStatus" "stageBiologicalStatus" "sex"   "strain"    "speciesId" "comment"   "annotatorId"   "lastModificationDate"  "replicate" "infoReplicate" "SRSId" "tags"RNASeqProtocol"   "physiological status"  "globin_reduction"  "PATOid"    "PATOname"
    #"SRX843135" "SRP051959" "Illumina HiSeq 2000"           "UBERON:0001871"    "temporal lobe" "UBERON:0000113"    "post-juvenile adult stage" "Brain Temporal Lobe"   2.91                    "perfect match" "not documented"    "missing child term"    "not documented"    "F" "NA"    9555    "library selection other but ok see PMID:25392405"  "AUC"   "2018-06-25"            "SRS819312"     "polyA"             

    if ( $annot =~ /^"#/ && $header eq '' ){
        $header = $annot;
        next ANNOT;
    }
    next ANNOT  if ( $annot =~ /^#/ || $annot =~ /^"#/ );

    my ($libraryId, $experimentId, $platform, undef, undef, $uberonId, undef, $stageId, undef,
        undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, $sex, $strain, $speciesId, undef)
        = split(/\t/, $annot, 23);
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

