#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Roux, created 09/05/08, last modified 25/11/10
# USAGE: perl insert_in_situ_zfin.pl ...
########################################################

use Getopt::Long;
use LWP::Simple;
use File::Slurp;

use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

$| = 1;

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($data)           = ('');
my ($Sport, $Aport)  = (0, 0);
my %opts = ('bgee=s'   => \$bgee_connector,   # Bgee connector string
            'data=s'   => \$data,
            'Sport=i'  => \$Sport,            # Stage mapper socket port
            'Aport=i'  => \$Aport,            # Anatomy mapper socket port
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $data eq '' || $Sport == 0 || $Aport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -data=\$(SOURCE_FILES_DIR)\$(DIR_NAME)ZebrafishMine.data -Sport=\$(INBETWEENSTAGESPORT) -Aport=\$(IDMAPPINGPORT)
\t-bgee      Bgee    connector string
\t-data      map_zfin tsv data file
\t-Sport     Stage   mapper socket port
\t-Aport     Anatomy mapper socket port
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


# Retrieve data source id for ZFIN
my $selSrc = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName = ?');
$selSrc->execute('ZFIN')  or die $selSrc->errstr;
my $data_source_id = $selSrc->fetchrow_array;
$selSrc->finish;
if ( !defined $data_source_id ){
    die "Data source ZFIN not found\n";
}


# Correspondence ZFIN -> Ensembl gene Ids (one to many)
#NOTE for Zebrafish forces speciesId=7955 to avoid insert of other fishes data (e.g. Tetraodon)
my $selGeneXref = $dbh->prepare('SELECT x.XRefId, g.geneId FROM geneXRef x, gene g WHERE x.dataSourceId = ? AND x.bgeeGeneId = g.bgeeGeneId AND g.speciesId=7955');
$selGeneXref->execute($data_source_id)  or die $selGeneXref->errstr;
my %ensembl_gene;
while ( my @data = $selGeneXref->fetchrow_array ){
    # zfin geneId -> ensembl geneId
    $ensembl_gene{$data[0]} = $data[1];
}
$selGeneXref->finish; # Let this finish() because causes warnings
$dbh->disconnect;


# Read ZebrafishMine data
my %tsv = %{ Utils::read_spreadsheet("$data", "\t", 'csv', '', 1) };


# map_zfin file header
# primaryIdentifier   expressions.expressionFound   symbol   expressions.anatomy.name   expressions.anatomy.identifier   expressions.startStage.name   expressions.startStage.identifier   expressions.endStage.name   expressions.endStage.identifier   expressions.figures.images.primaryIdentifier   expressions.figures.images.label   expressions.figures.primaryIdentifier   expressions.publication.primaryIdentifier   organism.taxonId   expressions.assay   probe.quality
# Probe Quality (optional 0 - 5 rating)

# Read all start_stage_id/end_stage_id to run the in between stages socket once
my @Stages;
for my $id ( 0..$#{$tsv{'expressions.startStage.identifier'}} ){
    next  if ( $tsv{'expressions.startStage.identifier'}[$id] eq '' || $tsv{'expressions.endStage.identifier'}[$id] eq '' || $tsv{'expressions.startStage.identifier'}[$id] eq 'None' || $tsv{'expressions.endStage.identifier'}[$id] eq 'None' );
    push @Stages, $tsv{'expressions.startStage.identifier'}[$id], $tsv{'expressions.endStage.identifier'}[$id];
}
my $doneStages = Utils::get_in_between_stages(\@Stages, $Sport);
my @Anat = @{ $tsv{'expressions.anatomy.identifier'} };
my $doneAnat   = Utils::get_anatomy_mapping(\@Anat, $Aport);


# Output TSV
print join("\t", '#data_source', qw(inSituExperimentId  inSituEvidenceId  organId  stageId  geneId  detectionFlag  inSituData  linked  speciesId  strain  sex)), "\n";
for my $id ( 0..$#{$tsv{'primaryIdentifier'}} ){
    # Some other constrains are set in the script querying ZebraFishMine!: StandardEnvironment/WildType/mRNA in situ hybridization
    next  if ( $tsv{'organism.taxonId'}[$id]  != 7955 );                         # Only Tetraodon in ZFIN currently
    next  if ( $tsv{'primaryIdentifier'}[$id] eq 'None' || $tsv{'primaryIdentifier'}[$id] eq '' );

    if ( !exists $ensembl_gene{ $tsv{'primaryIdentifier'}[$id] } ){
        warn "No gene mapping for [$tsv{'primaryIdentifier'}[$id]]\n";
        next;        # Ensembl correspondance MUST exist
    }

    my $inbetweenstages = $tsv{'expressions.startStage.identifier'}[$id].','.$tsv{'expressions.endStage.identifier'}[$id];

    my $stage_id = $doneStages->{ $inbetweenstages }                          || '';
    my $organ_id = $doneAnat->{ $tsv{'expressions.anatomy.identifier'}[$id] } || '';
    if ( $stage_id eq '' || $stage_id =~ /Could not find any OWLClass corresponding to/ ){
        warn "No stage reference for [$inbetweenstages] in [$id / $tsv{'expressions.publication.primaryIdentifier'}[$id] | $tsv{'expressions.figures.primaryIdentifier'}[$id]]\n";
        next;
    }
    if ( $organ_id eq '' || $organ_id =~ /Could not find any OWLClass corresponding to/ ){
        warn "Problem with anatId for [$tsv{'expressions.anatomy.identifier'}[$id]] in [$id / $tsv{'expressions.publication.primaryIdentifier'}[$id] | $tsv{'expressions.figures.primaryIdentifier'}[$id]]\n";
        next
    }

    #NOTE "None" is the value for an empty field (Python code only?)
    my $quality = defined $tsv{'probe.quality'}[$id] && $tsv{'probe.quality'}[$id] eq 'None' ? 'high quality'
                : defined $tsv{'probe.quality'}[$id] && $tsv{'probe.quality'}[$id] >= 2      ? 'high quality'
                : !defined $tsv{'probe.quality'}[$id]                                        ? 'high quality'
                :                                                                              'poor quality';

    # Few direct sex annotations available
    my $sex = $tsv{'expressions.anatomy.name'}[$id] eq 'female organism' ? 'female'
            : $tsv{'expressions.anatomy.name'}[$id] eq 'male organism'   ? 'male'
            : 'not annotated';


    STAGE_RANGE:
    for my $stage ( split(/\t/, $stage_id) ){
        # Skip this experiment if in between stage problem
        if ( $stage =~ /Exception/ || $stage =~ /Incorrect/ ){
            chomp $stage_id;
            warn "Problem for [$inbetweenstages] [$stage_id] in [$tsv{'expressions.anatomy.identifier'}[$id]]";
            next STAGE_RANGE;#  if ( $stage =~ /Could not find any OWLClass corresponding/ );
        }

        print join("\t", $data_source_id,
                         $tsv{'expressions.publication.primaryIdentifier'}[$id] eq 'None' ? '' : $tsv{'expressions.publication.primaryIdentifier'}[$id],
                         $tsv{'expressions.figures.primaryIdentifier'}[$id]     eq 'None' ? '' : $tsv{'expressions.figures.primaryIdentifier'}[$id],
                         $organ_id,
                         $stage,
                         $ensembl_gene{ $tsv{'primaryIdentifier'}[$id] },
                         $tsv{'expressions.expressionFound'}[$id] eq 'true' ? 'present' : 'absent',
                         $quality,
                         $tsv{'expressions.figures.images.primaryIdentifier'}[$id] eq 'None' ? '' : $tsv{'expressions.figures.images.primaryIdentifier'}[$id],
                         $tsv{'organism.taxonId'}[$id],
                         #NOTE Of wild type lines defined in ZFIN (condition forced in zebra_query.py/pl)
                         # because some wild type lines are not really "wild type" according to Leyla Ruzicka <leyla@zfin.org>
                         # They are "wild type" according to papers' authors!
                         # Should be better investigatd looking at https://zfin.org/action/feature/wildtype-list
                         $Utils::WILD_TYPE_STRAIN, #NOTE and "normal conditions"
                         $sex,
                  ), "\n";
    }
}

exit 0;

