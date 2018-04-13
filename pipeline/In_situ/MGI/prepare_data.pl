#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Prepare data from MGI-MouseMine before insertion
#
#############################################################

$| = 1; # stdout not put in memory buffer

use Getopt::Long;
use File::Slurp;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;


# Define arguments & their default value
my ($bgee_connector) = ('');
my ($data)           = ('');
my ($Aport, $Sport)  = (0, 0);
my %opts = ('bgee=s'          => \$bgee_connector,     # Bgee connector string
            'data=s'          => \$data,
            'Aport=i'         => \$Aport,              # Anatomy mapper socket port
            'Sport=i'         => \$Sport,              # Stage mapper socket port
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $Aport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)  -data=map_mgi -Aport=\$(IDMAPPINGPORT) -Sport=\$(INBETWEENSTAGESPORT)
\t-bgee            Bgee connector string
\t-data            MGI/MouseMine tsv data file
\t-Aport           Anatomy mapper socket port
\t-Sport           Stage   mapper socket port
\n";
    exit 1;
}


# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

# Get MGI source id
my $selSrc = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName = ?');
$selSrc->execute('MGI')  or die $selSrc->errstr;
my $data_source_id = $selSrc->fetchrow_array;
$selSrc->finish;
$dbh->disconnect;
die "Data source MGI not found\n"  if ( !defined $data_source_id );


# Read MGI data
my %tsv = %{ Utils::read_spreadsheet("$data", "\t", 'csv', '', 1) };
#NOTE if Ensembl ids are provided in the data file, it means we skip MGI genes without Ensembl id !!!

# Get in between start-end stages + anatomy mapping
my @Stages     = map { ($_, $_ ) }
                 @{$tsv{'stage'}};
#warn @Stages, "\n";
my $doneStages = Utils::get_in_between_stages(\@Stages, $Sport);

my @Anat       = @{ $tsv{'structure.id'} };
#warn @Anat, "\n";
my $doneAnat   = Utils::get_anatomy_mapping(\@Anat, $Aport);


#Assay type   feature.symbol feature.id  feature.xref.id     stage   age                 structure.name                                    structure.id  strength   pattern        assayId      probe        image         specimen  sex             detected   taxonId
#RNA in situ  0610005C13Rik  MGI:1918911 ENSMUSG00000109644  TS23    embryonic day 15.5  anlage of loop of Henle of cortical renal tubule  EMAPA:31281   Ambiguous  Not Specified  MGI:5541511  MGI:5000953  GUDMAP:13955  1         Not Specified              10090
#NOTE keep Assay type column to be sure to have a column without empty values
print join("\t", '#data_source', qw(inSituExperimentId  inSituEvidenceId  organId  stageId  geneId  detectionFlag  inSituData  linked  speciesId  strain  sex)), "\n";
SPOT:
for my $line ( 0..$#{$tsv{'Assay type'}} ){
    next  if ( $tsv{'Assay type'}[$line] ne 'RNA in situ' );                    # Only In situ
    next  if ( !$tsv{'taxonId'}[$line] || $tsv{'taxonId'}[$line]    != 10090 ); # Only Mus musculus

    # stage mapping
    my $inbetween = $tsv{'stage'}[$line].','.$tsv{'stage'}[$line];
    my $stage = $doneStages->{ $inbetween } || '';
    if ( $stage eq '' || $stage =~ /Could not find any OWLClass corresponding to/ ){
        warn "No stage reference for [$tsv{'stage'}[$line]] in [$tsv{'assayId'}[$line]]\n";
        next SPOT;
    }

    # anatId mapping
    my $anatId = $doneAnat->{ $tsv{'structure.id'}[$line] } || '';
    if ( $anatId eq '' || $anatId =~ /Could not find any OWLClass corresponding to/ ){
        warn "No anatomyId reference for [$tsv{'structure.id'}[$line]] in [$tsv{'assayId'}[$line]]\n";
        next SPOT;
    }

    # quality assignment
    my $quality  = get_quality( lc $tsv{'strength'}[$line] );
    if ( $quality eq '' ){
        warn "No defined quality strength for [$tsv{'strength'}[$line]] in [$tsv{'assayId'}[$line]]\n";
        next SPOT;
    }

    # presence
    #TODO Additional/Complementary info with $tsv{'detected'}[$line] ???
    my $presence = $quality eq 'high quality' ? 'present'
                 : $quality eq 'poor quality' ? 'present'
                 : $quality eq 'absent'       ? 'absent'
                 :                              '';
    if ( $presence eq '' ){
        warn "No defined presence for [$tsv{'strength'}[$line]] in [$tsv{'assayId'}[$line]]\n";
        next SPOT;
    }
    $quality = 'high quality'  if ( $presence eq 'absent' );

    # image
    #NOTE image regexp for URL transformation
    my $image = $tsv{'image'}[$line];
    $image   =~ s{[+-]\/[+-]}{}g;
    $image   =~ s{ }{_}g;
    $image   =~ s{/}{_}g;
    $image   =~ s{-}{_}g;
    $image   =~ s{\W+}{}g;
    $image   =~ s{__+}{_}g;
    $image   .= '_id'  if ( $image ne '' );

    # Sex
    my $sex = $tsv{'sex'}[$line] eq 'Female' ? 'female'
            : $tsv{'sex'}[$line] eq 'Male'   ? 'male'
            :                                  'not annotated'; # sex field is empty or "Not Specified" in that case

    #FIXME Loop over $stage as in other prepare_data.pl scripts???
    # output
    print join("\t", $data_source_id,
                     $tsv{'assayId'}[$line],
                     $tsv{'assayId'}[$line].$tsv{'specimen'}[$line],
                     $anatId,
                     $stage,
                     $tsv{'feature.xref.id'}[$line],
                     $presence,
                     $quality,
                     $image,
                     $tsv{'taxonId'}[$line],
                     $Utils::WILD_TYPE_STRAIN, #FIXME Don't know how to get this info in MouseMine yet
                     $sex,
              ), "\n";
}

exit 0;


sub get_quality {
    my ($strength) = @_;

    # define levels of quality
    my %quality = (
        'not applicable' => '',             # -2
        'not specified'  => 'high quality', # -1
        'absent'         => 'absent',       #  1
        'present'        => 'high quality', #  2  # this is debatable... but with zfin we attribute high quality when there is no star
        'ambiguous'      => 'poor quality', #  3
        'trace'          => 'poor quality', #  4
        'weak'           => 'poor quality', #  5
        'moderate'       => 'high quality', #  6
        'strong'         => 'high quality', #  7
        'very strong'    => 'high quality', #  8
    );

    if ( exists $quality{$strength} ){
        return $quality{$strength};
    }
    else {
        die "Invalid quality strength: [$strength]\n";
        return;
    }
}

