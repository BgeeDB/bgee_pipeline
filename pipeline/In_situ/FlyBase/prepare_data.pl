#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Prepare data from FlyBase before insertion
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
my ($Sport, $Aport)  = (0, 0);
my %opts = ('bgee=s'    => \$bgee_connector,   # Bgee connector string
            'data=s'    => \$data,
            'Sport=i'   => \$Sport,            # Stage mapper socket port
            'Aport=i'   => \$Aport,            # Anatomy mapper socket port
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $data eq '' || $Sport == 0 || $Aport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)  -data=map_flybase -Sport=\$(INBETWEENSTAGESPORT) -Aport=\$(IDMAPPINGPORT)
\t-bgee        Bgee connector string
\t-data        FlyBase tsv data file
\t-Sport       Stage   mapper socket port
\t-Aport       Anatomy mapper socket port
\n";
    exit 1;
}


# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

# Get MGI source id
my $selSrc = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName = ?');
$selSrc->execute('FlyBase')  or die $selSrc->errstr;
my $data_source_id = $selSrc->fetchrow_array;
$selSrc->finish;
$dbh->disconnect;
die "Data source FlyBase not found\n"  if ( !defined $data_source_id );


# Read FlyBase data
my %tsv = %{ Utils::read_spreadsheet("$data", "\t", 'csv', '', 1) };


# Read all start_stage_id/end_stage_id to run the in between stages socket once
my @Stages;
for my $id ( 0..$#{$tsv{'start_stage_id'}} ){
    next  if ( $tsv{'start_stage_id'}[$id] eq '' || $tsv{'end_stage_id'}[$id] eq '' );
    push @Stages, $tsv{'start_stage_id'}[$id], $tsv{'end_stage_id'}[$id];
}
my $doneStages = Utils::get_in_between_stages(\@Stages, $Sport);
my @Anat = @{ $tsv{'anatomy_id'} };
my $doneAnat   = Utils::get_anatomy_mapping(\@Anat, $Aport);


my %seen;
#feid    gene_name  gene_id      pub_ref      pub_desc                                            organism_id  assay    molecule_assayed  stage        sex            anatomy                                expression_uniquename  anatomy_desc                           anatomy_id     start_stage_name   start_stage_id  end_stage_name  end_stage_id
#100080  Or42b      FBgn0033043  FBrf0208960  Bai et al., 2009, J. Neurosci. 29(41): 12940--12947 1            in situ  mRNA              adult stage  not specified  adult olfactory receptor neuron Or42b  FBex0014101            adult olfactory receptor neuron Or42b  FBbt:00067046  adult stage        FBdv:00005369   adult stage     FBdv:00005369
#NOTE expression_uniquename is not used anymore, replaced by gene report == gene_id
print join("\t", '#data_source', qw(inSituExperimentId  inSituEvidenceId  organId  stageId  geneId  detectionFlag  inSituData  linked  speciesId  strain  sex)), "\n";
SPOT:
for my $line ( 0..$#{$tsv{'feid'}} ){
    next  if ( $tsv{'assay'}[$line]            ne 'in situ' );
    next  if ( $tsv{'molecule_assayed'}[$line] ne 'mRNA' ); # Remove molecule_assayed = protein for now because protein can be extracted, translocated in the whole organism

    $seen{ $tsv{'pub_ref'}[$line] }++;

    # stage mapping
    my $inbetweenstages = $tsv{'start_stage_id'}[$line].','.$tsv{'end_stage_id'}[$line];
    my $stages = $doneStages->{$inbetweenstages} || '';
    if ( $stages eq '' || $inbetweenstages eq ',' || $stages =~ /Could not find any OWLClass corresponding to/ ){
        warn "Problem with in between stages for [$inbetweenstages] in [$tsv{'pub_ref'}[$line]]\n";
        next SPOT;
    }

    # Organ mapping
    my $anatId = $doneAnat->{$tsv{'anatomy_id'}[$line]} || '';
    if ( $anatId eq '' || $anatId =~ /Could not find any OWLClass corresponding to/ ){
        warn "Problem with anatId for [$tsv{'anatomy_id'}[$line]] in [$tsv{'pub_ref'}[$line]]\n";
        next SPOT;
    }

    # Sex
    my $sex = $tsv{'sex'}[$line] =~ /female/i ? 'female'
            : $tsv{'sex'}[$line] =~ /male/i   ? 'male'
            :                                   'not annotated'; # Mostly "not specified"

    # Species tax id
    my $speciesId = $tsv{'organism_id'}[$line] == 1   ? 7227 # Mostly Drosophila melanogaster                experiments
                  : $tsv{'organism_id'}[$line] == 214 ? 7237 # Few    Drosophila pseudoobscura pseudoobscura experiments
                  : $tsv{'organism_id'}[$line] == 78  ? 7244 # No     Drosophila virilis                     experiments yet
                  :                                     '';
    if ( $speciesId eq '' ){
        warn "Problem with organism_id [$tsv{'organism_id'}[$line]]\n";
        next SPOT;
    }

#TODO A single evidence per assay, duplicate only spots

    # quality assignment
    #NOTE FlyBase does not provide in situ quality information (available only for RNA Seq).
    my $quality  = 'high quality';

    # presence
    #NOTE FlyBase provides very few absent annotation, in a non-structured way.
    # So we don't catch them for now and let them mixed with "regular" present annotation.
    my $presence = 'present';


    STAGE_RANGE:
    for my $stage ( split(/\t/, $stages) ){
        # Skip this experiment if in between stage problem
        if ( $stage =~ /Exception/ || $stage =~ /Incorrect/ ){
            chomp $stages;
            warn "Problem for [$inbetweenstages] [$stages] in [$tsv{'pub_ref'}[$line]]";
            next STAGE_RANGE;#  if ( $stage =~ /Could not find any OWLClass corresponding/ );
        }

        # output
        print join("\t", $data_source_id,
                         $tsv{'pub_ref'}[$line],
                         $tsv{'pub_ref'}[$line].$seen{ $tsv{'pub_ref'}[$line] },
                         $anatId,
                         $stage,
                         $tsv{'gene_id'}[$line],
                         $presence,
                         $quality,
                         '', # No linked info
                         $speciesId,
                         $Utils::WILD_TYPE_STRAIN, #FIXME Not really strain info in FlyBase
                         $sex,
              ), "\n";
    }
}

exit 0;

