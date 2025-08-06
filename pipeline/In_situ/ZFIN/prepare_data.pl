#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

$| = 1;

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($wt_exp, $exp_assay, $exp_env, $exp_stage_anat, $mapping) = ('', '', '', '', '');
my ($Sport, $Aport)  = (0, 0);
my %opts = ('bgee=s'   => \$bgee_connector,   # Bgee connector string
            'wt_exp=s' => \$wt_exp,
            'assay=s'  => \$exp_assay,
            'env=s'    => \$exp_env,
            'cond=s'   => \$exp_stage_anat,
            'map=s'    => \$mapping,
            'Sport=i'  => \$Sport,            # Stage mapper socket port
            'Aport=i'  => \$Aport,            # Anatomy mapper socket port
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $wt_exp eq '' || $exp_assay eq '' || $exp_env eq '' || $exp_stage_anat eq '' || $mapping eq '' || $Sport == 0 || $Aport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -wt_exp=\$(SOURCE_FILES_DIR)\$(DIR_NAME)wildtype-expression_fish.txt -assay=\$(SOURCE_FILES_DIR)\$(DIR_NAME)xpat_fish.txt -env=\$(SOURCE_FILES_DIR)\$(DIR_NAME)xpat_environment_fish.txt -cond=\$(SOURCE_FILES_DIR)\$(DIR_NAME)xpat_stage_anatomy.txt -map=\$(SOURCE_FILES_DIR)\$(DIR_NAME)stage_ontology.txt -Sport=\$(INBETWEENSTAGESPORT) -Aport=\$(IDMAPPINGPORT)
\t-bgee      Bgee    connector string
\t-wt_exp    wildtype-expression_fish.txt file
\t-assay     xpat_fish.txt file
\t-env       xpat_environment_fish.txt file
\t-cond      xpat_stage_anatomy.txt file
\t-map       stage_ontology.txt file
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


# Read stage ontology mapping file
my %stages;
for my $line ( `cat $mapping` ){
    chomp $line;
    my @col = split(/\t/, $line);
    $stages{ $col[0] } = $col[1];
}


# Read ZFIN raw data: wildtype-expression_fish.txt
my $genes;
# Only Danio rerio in ZFIN currently && ONLY in situ
#Gene ID    Gene Symbol    Fish Name    Super Structure ID    Super Structure Name    Sub Structure ID    Sub Structure Name    Start Stage    End Stage    Assay    Assay MMO ID    Publication ID    Probe ID    Antibody ID    Fish ID
for my $line ( `cat $wt_exp | grep -v '^#' | grep -P '\tmRNA in situ hybridization\t'` ){
    chomp $line;
    my @col = split(/\t/, $line);
    $genes->{$col[0]}->{$col[11]}->{$col[14]}->{'genotype'} = $col[2]; # Fish Name: genotype from a wild-type list: https://zfin.org/action/feature/wildtype-list
    $genes->{$col[0]}->{$col[11]}->{$col[14]}->{'sex'}      = $col[4];
}

# Read ZFIN raw data: xpat_environment_fish.txt
# ONLY standard conditions
#Environment ID    ZECO Term Name    ZECO Term ID (ZECO:ID)    ......
my $env;
for my $line ( `cat $exp_env | grep -v '^#' | grep -P '\tstandard conditions\t'` ){
    chomp $line;
    my @col = split(/\t/, $line);
    $env->{$col[0]} = 1;
}

# Read ZFIN raw data: xpat_stage_anatomy.txt
#Expression Result ID    Expression ID    Start Stage ID    End Stage ID    Anatomy Super Term ID    Anatomy Sub Term ID    Expression Found
my $exp;
for my $line ( `cat $exp_stage_anat | grep -v '^#'` ){
    chomp $line;
    my @col = split(/\t/, $line);
    my $start = $stages{ $col[2] };
    my $end   = $stages{ $col[3] };
    my $anat  = ($col[5] ne '' && $col[5] !~ /^GO:/ and $col[5] !~ /^BSPO:/) ? $col[5] : $col[4]; # BSPO and GO are not often all available in Uberon for the mapping!
    $exp->{$col[1]}->{ "$start-$end-$anat" } = $col[6];
}

# Read ZFIN raw data: xpat_fish.txt
# ONLY in situ
#Gene ID    Gene Symbol    EST ID (Optional)    EST Symbol (Optional)    Expression Type    Expression Type MMO ID    Expression ID    Publication ID    Fish ID    Environment ID    Probe Quality (optional 0 - 5 rating)
my @Stages;
my @Anat;
for my $line ( `cat $exp_assay | grep -v '^#' | grep -P '\tmRNA in situ hybridization\t'` ){
    chomp $line;
    my @col = split(/\t/, $line);
    if ( exists $genes->{$col[0]}->{$col[7]}->{$col[8]} && exists $env->{$col[9]} && exists $exp->{$col[6]} ){
        $genes->{$col[0]}->{$col[7]}->{$col[8]}->{'quality'}     = $col[10];
        $genes->{$col[0]}->{$col[7]}->{$col[8]}->{'environment'} = $col[9];
        $genes->{$col[0]}->{$col[7]}->{$col[8]}->{'expression'}  = $col[6];
        $genes->{$col[0]}->{$col[7]}->{$col[8]}->{'image'}       = $col[3];
        for my $cond ( keys %{ $exp->{$col[6]} } ){
            my ($start, $end, $anat) = split(/-/, $cond);
            push @Anat,   $anat;
            push @Stages, $start, $end;

            push @{ $genes->{$col[0]}->{$col[7]}->{$col[8]}->{'condition'} }, $cond;
            $genes->{$col[0]}->{$col[7]}->{$col[8]}->{'presence'} = $exp->{ $col[6] }->{ "$cond" };
        }
    }
}


# map_zfin file header
my $doneStages = Utils::get_in_between_stages(\@Stages, $Sport);
my $doneAnat   = Utils::get_anatomy_mapping(\@Anat, $Aport);


# Output TSV
my $count = 0;
print join("\t", '#data_source', qw(inSituExperimentId  inSituEvidenceId  organId  stageId  geneId  detectionFlag  inSituData  linked  speciesId  strain  sex)), "\n";
for my $gene_id ( keys %{$genes} ){
    if ( !exists $ensembl_gene{ $gene_id } ){
        warn "No gene mapping for [$gene_id]\n";
        next;        # Ensembl correspondance MUST exist
    }

    for my $exp_id ( sort keys %{ $genes->{$gene_id} } ){
        for my $fish_id ( keys %{ $genes->{$gene_id}->{$exp_id} } ){
            next  if ( !exists $genes->{$gene_id}->{$exp_id}->{$fish_id}->{'condition'} );

            for my $cond ( @{ $genes->{$gene_id}->{$exp_id}->{$fish_id}->{'condition'} } ){
                my ($start, $end, $anat) = split(/-/, $cond);

                my $inbetweenstages = $start.','.$end;
                my $stage_id = $doneStages->{ $inbetweenstages }  || '';
                my $organ_id = $doneAnat->{ $anat }               || '';
                if ( $stage_id eq '' || $stage_id =~ /Could not find any OWLClass corresponding to/ ){
                    warn "No stage reference for [$start/$end] in [$gene_id / $ensembl_gene{ $gene_id } | $exp_id]\n";
                    next;
                }
                if ( $organ_id eq '' || $organ_id =~ /Could not find any OWLClass corresponding to/ ){
                    warn "Problem with anatId for [$anat] in [$gene_id / $ensembl_gene{ $gene_id } | $exp_id]\n";
                    next;
                }


                #NOTE "None" is the value for an empty field (Python code only?)
                my $quality = ! $genes->{$gene_id}->{$exp_id}->{$fish_id}->{'quality'}         ? 'high quality'
                            : $genes->{$gene_id}->{$exp_id}->{$fish_id}->{'quality'} eq 'None' ? 'high quality'
                            : $genes->{$gene_id}->{$exp_id}->{$fish_id}->{'quality'} eq ''     ? 'high quality'
                            : $genes->{$gene_id}->{$exp_id}->{$fish_id}->{'quality'} >= 2      ? 'high quality'
                            :                                                                    'poor quality';


                # Few direct sex annotations available
                my $sex = $genes->{$gene_id}->{$exp_id}->{$fish_id}->{'sex'} eq 'female organism' ? 'female'
                        : $genes->{$gene_id}->{$exp_id}->{$fish_id}->{'sex'} eq 'male organism'   ? 'male'
                        : 'not annotated';

                #TODO what to use as inSituEvidenceId, was ZDB-FIG-...-... before

                STAGE_RANGE:
                for my $stage ( sort split(/\t/, $stage_id) ){
                    # Skip this experiment if in between stage problem
                    if ( $stage =~ /Exception/ || $stage =~ /Incorrect/ ){
                        chomp $stage_id;
                        warn "Problem for [$inbetweenstages] [$stage_id] in [$anat]\n";
                        next STAGE_RANGE;
                    }

                    print join("\t", $data_source_id,
                                     $exp_id,
                                     'ZFIN'.++$count,
                                     $organ_id,
                                     $stage,
                                     $ensembl_gene{ $gene_id },
                                     $genes->{$gene_id}->{$exp_id}->{$fish_id}->{'presence'} eq 't' ? 'present' : 'absent',
                                     $quality,
                                     'ZFIN'.$count,
                                     7955,
                                     #NOTE Of wild type lines defined in ZFIN (condition forced in zebra_query.py/pl)
                                     # because some wild type lines are not really "wild type" according to Leyla Ruzicka <leyla@zfin.org>
                                     # They are "wild type" according to papers' authors!
                                     # Should be better investigated looking at https://zfin.org/action/feature/wildtype-list
                                     $genes->{$gene_id}->{$exp_id}->{$fish_id}->{'genotype'} eq 'WT' ? $Utils::WILD_TYPE_STRAIN : $genes->{$gene_id}->{$exp_id}->{$fish_id}->{'genotype'},
                                     $sex,
                              ), "\n";
                }
            }
        }
    }
}
exit 0;

