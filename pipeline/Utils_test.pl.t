#!/usr/bin/env perl

# This files contains the tests for the fractionnal ranking implementation
# Author: Philippe Moret

use strict;
use warnings;
use Test::More;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

# Test fractional ranking

# Test 1
my @data = ({"id" =>"A", "val" => 1.0},
    {"id" => "B", "val" => 42.1},
    {"id" => "C", "val" => 42.1}
);

my %expected = ("A", 1.0, "B",2.5,"C",2.5);
my %r = Utils::fractionnal_ranking(@data);

for my $k(keys(%expected)) {
    ok (  $r{$k} eq $expected{$k} , "Test1: Expecting $expected{$k}  for key $k");
}

# Test 2
@data = (
    {"id" => "A", "val" =>1.0},
    {"id" => "B", "val" =>1.0},
    {"id" => "C", "val" => 1.0},
    {"id" => "D", "val" =>1.0}
);
%expected = ("A",2.5,"B",2.5,"C",2.5,"D",2.5);
%r = Utils::fractionnal_ranking(@data);

for my $k(keys(%expected)) {
    ok (  $r{$k} eq $expected{$k} , "Test1: Expecting $expected{$k}  for key $k");
}

# Test 3
@data = ({"id" => "A", "val" => 1.0},
    {"id" => "B", "val" => 1.0},
    {"id" => "C" ,"val" => 1.0},
    {"id" => "D", "val" => 2.5},
    {"id" => "E","val" =>3.5});

%expected = ("A",2.0,"B",2.0,"C",2.0,"D",4.0,"E", 5.0);
%r = Utils::fractionnal_ranking(@data);

for my $k(keys(%expected)) {
    ok ($r{$k} eq $expected{$k} , "Test3: Expecting $expected{$k}  for key $k got $r{$k} ");
}

# Test dense ranking

# Test 1
@data = ({"id" =>"A", "val" => 1.0},
    {"id" => "B", "val" => 42.1},
    {"id" => "C", "val" => 42.1}
);

%expected = ("A", 1, "B",2,"C",2);
%r = Utils::dense_ranking(@data);

for my $k(keys(%expected)) {
    ok (  $r{$k} eq $expected{$k} , "Test1: Expecting $expected{$k}  for key $k");
}

# Test 2
@data = (
    {"id" => "A", "val" =>1.0},
    {"id" => "B", "val" =>1.0},
    {"id" => "C", "val" => 1.0},
    {"id" => "D", "val" =>1.0}
);
%expected = ("A",1,"B",1,"C",1,"D",1);
%r = Utils::dense_ranking(@data);

for my $k(keys(%expected)) {
    ok (  $r{$k} eq $expected{$k} , "Test1: Expecting $expected{$k}  for key $k");
}

# Test 3
@data = ({"id" => "A", "val" => 1.0},
    {"id" => "B", "val" => 1.0},
    {"id" => "C" ,"val" => 1.0},
    {"id" => "D", "val" => 2.5},
    {"id" => "E","val" =>3.5});

%expected = ("A",1,"B",1,"C",1,"D",2,"E", 3);
%r = Utils::dense_ranking(@data);

for my $k(keys(%expected)) {
    ok ($r{$k} eq $expected{$k} , "Test3: Expecting $expected{$k}  for key $k got $r{$k} ");
}

# Test reverse hash

# Test 1
my %toReverse = ("A",1,"B",1,"C",2,"D",2,"E", 3);
%expected = (1, ["A", "B"], 2, ["C", "D"], 3, ["E"]);
%r = Utils::revhash(%toReverse);

for my $k(keys(%expected)) {
    ok (  @{$r{$k}} eq @{$expected{$k}} , "Test1: Expecting @{$expected{$k}}  for key $k, was @{$r{$k}}");
}


# Test infer_sex
my %anatSexInfo = ();
# there will be no sex inference for this one
my @sexes = ($Utils::MALE_SEX, $Utils::FEMALE_SEX);
$anatSexInfo{'organId1'} = \@sexes;
# there will be sex inference because only one possibility
my @sexes2 = ($Utils::FEMALE_SEX);
$anatSexInfo{'organId2'} = \@sexes2;
# there will be sex inference because ne hermaphrodite in the species
my @sexes3 = ($Utils::MALE_SEX, $Utils::HERMAPHRODITE_SEX);
$anatSexInfo{'organId3'} = \@sexes3;
my %speciesSexInfo = ();
my @speciesSexes = ($Utils::MALE_SEX, $Utils::FEMALE_SEX);
$speciesSexInfo{'species1'} = \@speciesSexes;
ok(!defined Utils::infer_sex(\%anatSexInfo, \%speciesSexInfo, 'organId1', 'species1'),
    'Inference while there should be none');
ok(Utils::infer_sex(\%anatSexInfo, \%speciesSexInfo, 'organId2', 'species1') eq $Utils::FEMALE_SEX,
    'No inference for a sex-specific organ');
ok(Utils::infer_sex(\%anatSexInfo, \%speciesSexInfo, 'organId3', 'species1') eq $Utils::MALE_SEX,
    'No inference for a species with no hermaphrodite');



# Test summarizeExperimentCallAndQuality
my %calls = ();
# First, let's try with only one experiment
$calls{'expId1'}->{'libId1'}->{'call'}    = $Utils::PRESENT_CALL;
$calls{'expId1'}->{'libId1'}->{'quality'} = $Utils::HIGH_QUAL;
$calls{'expId1'}->{'libId2'}->{'call'}    = $Utils::ABSENT_CALL;
$calls{'expId1'}->{'libId2'}->{'quality'} = $Utils::LOW_QUAL;
$calls{'expId1'}->{'libId3'}->{'call'}    = $Utils::PRESENT_CALL;
$calls{'expId1'}->{'libId3'}->{'quality'} = $Utils::HIGH_QUAL;
$calls{'expId1'}->{'libId4'}->{'call'}    = $Utils::PRESENT_CALL;
$calls{'expId1'}->{'libId4'}->{'quality'} = $Utils::LOW_QUAL;

my %expectedSummary = ();
$expectedSummary{'expId1'}->{'pstHighEvidenceCount'} = 2;
$expectedSummary{'expId1'}->{'pstLowEvidenceCount'}  = 1;
$expectedSummary{'expId1'}->{'absHighEvidenceCount'} = 0;
$expectedSummary{'expId1'}->{'absLowEvidenceCount'}  = 1;
$expectedSummary{'expId1'}->{'expCall'}              = $Utils::PRESENT_CALL;
$expectedSummary{'expId1'}->{'expCallQuality'}       = $Utils::HIGH_QUAL;

my ($exclusion, $actualSummary) = Utils::summarizeExperimentCallAndQuality(\%calls);
ok($exclusion eq $Utils::CALL_NOT_EXCLUDED, 'incorrect reason for exclusion');
is_deeply($actualSummary, \%expectedSummary, 'Summaries should be equal');

# now, with several experiments
%calls = ();
# present experiments
$calls{'expId1'}->{'libId1'}->{'call'}    = $Utils::PRESENT_CALL;
$calls{'expId1'}->{'libId1'}->{'quality'} = $Utils::HIGH_QUAL;
$calls{'expId1'}->{'libId2'}->{'call'}    = $Utils::ABSENT_CALL;
$calls{'expId1'}->{'libId2'}->{'quality'} = $Utils::HIGH_QUAL;

$calls{'expId2'}->{'libId3'}->{'call'}    = $Utils::PRESENT_CALL;
$calls{'expId2'}->{'libId3'}->{'quality'} = $Utils::HIGH_QUAL;
$calls{'expId2'}->{'libId4'}->{'call'}    = $Utils::PRESENT_CALL;
$calls{'expId2'}->{'libId4'}->{'quality'} = $Utils::LOW_QUAL;
$calls{'expId2'}->{'libId5'}->{'call'}    = $Utils::ABSENT_CALL;
$calls{'expId2'}->{'libId5'}->{'quality'} = $Utils::HIGH_QUAL;
$calls{'expId2'}->{'libId6'}->{'call'}    = $Utils::ABSENT_CALL;
$calls{'expId2'}->{'libId6'}->{'quality'} = $Utils::LOW_QUAL;

$calls{'expId3'}->{'libId3'}->{'call'}    = $Utils::PRESENT_CALL;
$calls{'expId3'}->{'libId3'}->{'quality'} = $Utils::LOW_QUAL;
$calls{'expId3'}->{'libId4'}->{'call'}    = $Utils::PRESENT_CALL;
$calls{'expId3'}->{'libId4'}->{'quality'} = $Utils::LOW_QUAL;
$calls{'expId3'}->{'libId5'}->{'call'}    = $Utils::UNDEFINED_CALL;
$calls{'expId3'}->{'libId5'}->{'quality'} = undef;

#absent experiments
$calls{'expId4'}->{'libId1'}->{'call'}    = $Utils::ABSENT_CALL;
$calls{'expId4'}->{'libId1'}->{'quality'} = $Utils::HIGH_QUAL;
$calls{'expId4'}->{'libId2'}->{'call'}    = $Utils::ABSENT_CALL;
$calls{'expId4'}->{'libId2'}->{'quality'} = $Utils::LOW_QUAL;

$calls{'expId5'}->{'libId1'}->{'call'}    = $Utils::ABSENT_CALL;
$calls{'expId5'}->{'libId1'}->{'quality'} = $Utils::LOW_QUAL;
$calls{'expId5'}->{'libId2'}->{'call'}    = $Utils::ABSENT_CALL;
$calls{'expId5'}->{'libId2'}->{'quality'} = $Utils::LOW_QUAL;


%expectedSummary = ();
$expectedSummary{'expId1'}->{'pstHighEvidenceCount'} = 1;
$expectedSummary{'expId1'}->{'pstLowEvidenceCount'}  = 0;
$expectedSummary{'expId1'}->{'absHighEvidenceCount'} = 1;
$expectedSummary{'expId1'}->{'absLowEvidenceCount'}  = 0;
$expectedSummary{'expId1'}->{'expCall'}              = $Utils::PRESENT_CALL;
$expectedSummary{'expId1'}->{'expCallQuality'}       = $Utils::HIGH_QUAL;

$expectedSummary{'expId2'}->{'pstHighEvidenceCount'} = 1;
$expectedSummary{'expId2'}->{'pstLowEvidenceCount'}  = 1;
$expectedSummary{'expId2'}->{'absHighEvidenceCount'} = 1;
$expectedSummary{'expId2'}->{'absLowEvidenceCount'}  = 1;
$expectedSummary{'expId2'}->{'expCall'}              = $Utils::PRESENT_CALL;
$expectedSummary{'expId2'}->{'expCallQuality'}       = $Utils::HIGH_QUAL;

$expectedSummary{'expId3'}->{'pstHighEvidenceCount'} = 0;
$expectedSummary{'expId3'}->{'pstLowEvidenceCount'}  = 2;
$expectedSummary{'expId3'}->{'absHighEvidenceCount'} = 0;
$expectedSummary{'expId3'}->{'absLowEvidenceCount'}  = 0;
$expectedSummary{'expId3'}->{'expCall'}              = $Utils::PRESENT_CALL;
$expectedSummary{'expId3'}->{'expCallQuality'}       = $Utils::LOW_QUAL;

$expectedSummary{'expId4'}->{'pstHighEvidenceCount'} = 0;
$expectedSummary{'expId4'}->{'pstLowEvidenceCount'}  = 0;
$expectedSummary{'expId4'}->{'absHighEvidenceCount'} = 1;
$expectedSummary{'expId4'}->{'absLowEvidenceCount'}  = 1;
$expectedSummary{'expId4'}->{'expCall'}              = $Utils::ABSENT_CALL;
$expectedSummary{'expId4'}->{'expCallQuality'}       = $Utils::HIGH_QUAL;

$expectedSummary{'expId5'}->{'pstHighEvidenceCount'} = 0;
$expectedSummary{'expId5'}->{'pstLowEvidenceCount'}  = 0;
$expectedSummary{'expId5'}->{'absHighEvidenceCount'} = 0;
$expectedSummary{'expId5'}->{'absLowEvidenceCount'}  = 2;
$expectedSummary{'expId5'}->{'expCall'}              = $Utils::ABSENT_CALL;
$expectedSummary{'expId5'}->{'expCallQuality'}       = $Utils::LOW_QUAL;

($exclusion, $actualSummary) = Utils::summarizeExperimentCallAndQuality(\%calls);
ok($exclusion eq $Utils::CALL_NOT_EXCLUDED, 'incorrect reason for exclusion');
is_deeply($actualSummary, \%expectedSummary, 'Summaries should be equal');

# Now test when there exists only undefined results
%calls = ();
$calls{'expId1'}->{'libId1'}->{'call'}    = $Utils::UNDEFINED_CALL;
$calls{'expId1'}->{'libId1'}->{'quality'} = undef;
$calls{'expId1'}->{'libId2'}->{'call'}    = $Utils::UNDEFINED_CALL;
$calls{'expId1'}->{'libId2'}->{'quality'} = undef;

$calls{'expId2'}->{'libId3'}->{'call'}    = $Utils::UNDEFINED_CALL;
$calls{'expId2'}->{'libId3'}->{'quality'} = undef;
$calls{'expId2'}->{'libId4'}->{'call'}    = $Utils::UNDEFINED_CALL;
$calls{'expId2'}->{'libId4'}->{'quality'} = undef;

%expectedSummary = ();

($exclusion, $actualSummary) = Utils::summarizeExperimentCallAndQuality(\%calls);
ok($exclusion eq $Utils::EXCLUDED_FOR_UNDEFINED, 'incorrect reason for exclusion');
is_deeply($actualSummary, \%expectedSummary, 'Summaries should be equal');

done_testing();
