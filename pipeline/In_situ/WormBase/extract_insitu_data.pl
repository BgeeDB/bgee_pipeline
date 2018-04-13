#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use File::Slurp;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

#NOTE if several anatomy terms, life_stage is uncertain: maybe linked to some anatomy terms not all of them
# so life_stage(s) is mapped only to ID root stage: UBERON:0000104
my $id_root_stage = 'UBERON:0000104'; # "life cycle"
# Same for stages
my $id_root_organ = 'UBERON:0000465'; # "material anatomical entity"

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($wormb_data)     = ('');
my ($Aport, $Sport)  = (0, 0);
my %opts = ('bgee=s'        => \$bgee_connector,     # Bgee connector string
            'wormb_data=s'  => \$wormb_data,         # from ftp://caltech.wormbase.org/pub/wormbase/expr_dump/
            'Aport=i'       => \$Aport,              # Anatomy mapper socket port
            'Sport=i'       => \$Sport,              # Stage mapper socket port
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $wormb_data eq '' || $Aport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)  -wormb_data=expr_pattern.ace.20160829 -Aport=\$(IDMAPPINGPORT) -Sport=\$(INBETWEENSTAGESPORT)
\t-bgee            Bgee connector string
\t-wormb_data      WormBase/WormMine tsv data file
\t-Aport           Anatomy mapper socket port
\t-Sport           Stage   mapper socket port
\n";
    exit 1;
}

#NOTE a way to filter and get only In Situ data would be added in WormMine. When???


# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);
# Get WormBase source id
my $selSrc = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName = ?');
$selSrc->execute('WormBase')  or die $selSrc->errstr;
my $data_source_id = $selSrc->fetchrow_array;
$selSrc->finish;
$dbh->disconnect;
die "Data source WormBase not found\n"  if ( !defined $data_source_id );


# Read data
my $expression;
my ($Expr_pattern, $gene, $ref, $strain, $In_Situ) = ('', '', '', '', 0);
my (@Anatomy_term, @Life_stage);
my @Anat;
my @Stages;
EXP:
for my $line ( read_file("$wormb_data", chomp => 1) ){
    # Reset values
	if ( $line =~ /^$/ ){
        if ( $In_Situ == 1 ){
            $expression->{$Expr_pattern}->{'id'}           = $Expr_pattern;
            $expression->{$Expr_pattern}->{'gene'}         = $gene;
            $expression->{$Expr_pattern}->{'ref'}          = $ref;
            $expression->{$Expr_pattern}->{'strain'}       = $strain;
            @{ $expression->{$Expr_pattern}->{'anatomy'} } = @Anatomy_term;
            @{ $expression->{$Expr_pattern}->{'stage'} }   = @Life_stage;
            push @Anat,   @Anatomy_term;
            push @Stages, @Life_stage;
        }

        $Expr_pattern = '';
        $gene         = '';
        $ref          = '';
        $strain       = '';
        @Anatomy_term = ();
        @Life_stage   = ();
        $In_Situ      = 0;
    }
    # Get fields
    elsif ( $line =~ /^Expr_pattern\s*:\s*"([^"]+)"/ ){
        $Expr_pattern = $1;
    }
    #e.g. Anatomy_term    "WBbt:0005772" Certain
    # OR  Anatomy_term    "WBbt:0008217" Life_stage "WBls:0000041"
    elsif ( $line =~ /^Anatomy_term\s*"([^"]+)"/ ){
        my $anat = $1;
        if ( $line =~ /^Anatomy_term\s*"[^"]+"\s*(\w+.*)$/ ){
            $anat .= "__$1";
            if ( $1 =~ /^Life_stage\s*"([^"]+)"/ ){
                push @Stages, $1;
            }
        }
        push @Anatomy_term, $anat;
    }
    #e.g. Not_in_Anatomy_term     "WBbt:0008598"
    # OR  Not_in_Anatomy_term     "WBbt:0008598" Life_stage "WBls:0000024"
    elsif ( $line =~ /^Not_in_Anatomy_term\s*"([^"]+)"/ ){
        my $anat = 'NOT'.$1;
        if ( $line =~ /^Not_in_Anatomy_term\s*"[^"]+"\s*(\w+.*)$/ ){
            $anat .= "__$1";
            if ( $1 =~ /^Life_stage\s*"([^"]+)"/ ){
                push @Stages, $1;
            }
        }
        push @Anatomy_term, $anat;
    }
    elsif ( $line =~ /^Gene\s*"([^"]+)"/ ){
        $gene = $1;
    }
    elsif ( $line =~ /^In_Situ/ ){
        $In_Situ = 1;
    }
    #e.g. Life_stage  "WBls:0000023"
    # OR  Life_stage  "WBls:0000038" Anatomy_term "WBbt:0008217"
    elsif ( $line =~ /^Life_stage\s*"([^"]+)"/ ){
        my $stage = $1;
        if ( $line =~ /^Life_stage\s*"[^"]+"\s*(\w+.*)$/ ){
            $stage .= "__$1";
            if ( $1 =~ /^Anatomy_term\s*"([^"]+)"/ ){
                push @Anat, $1;
            }
        }
        push @Life_stage, $stage;
    }
    #e.g. Not_in_Life_stage       "WBls:0000003"
    # OR  Not_in_Life_stage       "WBls:0000003" Anatomy_term "WBbt:0008217"   NOTE Not yet seen
    elsif ( $line =~ /^Not_in_Life_stage\s*"([^"]+)"/ ){
        my $stage = 'NOT'.$1;
        if ( $line =~ /^Not_in_Life_stage\s*"[^"]+"\s*(\w+.*)$/ ){
            $stage .= "__$1";
            if ( $1 =~ /^Anatomy_term\s*"([^"]+)"/ ){
                push @Anat, $1;
            }
        }
        push @Life_stage, $stage;
    }
    elsif ( $line =~ /^Reference\s*"([^"]+)"/ ){
        $ref = $1;
    }
    elsif ( $line =~ /^Strain\s*"([^"]+)"/ ){
        $strain = $1;
    }
    else {
        # ...
    }
}


# Get anatomy mapping
@Anat = map  { $_ = ' '.$_; $_ }        #NOTE WormBase anatomy id are maybe too short (?!?) and if this char extension is not done, the socket gets stuck
        grep { defined $_ && $_ ne '' }
        map  { s/^NOT//; $_ }           # Remove NOT / absence of expression qualifier
        map  { s/__.+$//; $_ }          # Remove quality info if any
        @Anat;
#warn @Anat, "\n";
my $doneAnat   = Utils::get_anatomy_mapping(\@Anat, $Aport);
#warn @Stages, "\n";
@Stages        = map { ($_, $_) }
                 map { s/^NOT//; $_ }   # Remove NOT / absence of expression qualifier
                 map { s/__.+$//; $_ }  # Remove quality info if any
                 @Stages;               # Emulate start-end stages by doubling it
my $doneStages = Utils::get_in_between_stages(\@Stages, $Sport);


my ($empty, $certain, $partial, $uncertain, $others) = (0, 0, 0, 0, 0);
# Output TSV
print join("\t", '#data_source', qw(inSituExperimentId  inSituEvidenceId  organId  stageId  geneId  detectionFlag  inSituData  linked  speciesId  strain  sex  [ExprDesc])), "\n";
for my $id ( sort keys %$expression ){
    # Only wild type (N2 is C. elegans wild type)
    next  if ( $expression->{$id}->{'strain'} ne 'N2' && $expression->{$id}->{'strain'} ne '' );
    # Get out if no info for gene!
    next  if ( $expression->{$id}->{'gene'} eq '' );


    # RESOLVED anatomy /stage pairs (one-to-one)
    # e.g. Anatomy_term    "WBbt:0008217" Life_stage "WBls:0000041"
    ANAT_RESOLVED:
    for my $annotation ( sort grep { /__Life_stage/ } @{ $expression->{$id}->{'anatomy'} } ){
        my ($anatomy_id, $stage_id) = $annotation =~ /^(.+?)__Life_stage\s*"(.+?)"/;
        print_TSV($expression, $id, $annotation, $anatomy_id, $stage_id);
    }

    # e.g. Life_stage  "WBls:0000038" Anatomy_term "WBbt:0008217"
    STAGE_RESOLVED:
    for my $annotation ( sort grep { /__Anatomy_term/ } @{ $expression->{$id}->{'stage'} } ){
        my ($stage_id, $anatomy_id) = $annotation =~ /^(.+?)__Anatomy_term\s*"(.+?)"/;
        print_TSV($expression, $id, $annotation, $anatomy_id, $stage_id);
    }


    # UNRESOLVED anatomy /stage pairs (one-to-many or many-to-many)
    my @unresolved_stage = grep { !/__Anatomy_term/ } @{ $expression->{$id}->{'stage'} };
    my @unresolved_organ = grep { !/__Life_stage/ }   @{ $expression->{$id}->{'anatomy'} };
    # At least one unresolved organ remain WITH at least one unresolved stage !!!
    if ( scalar @unresolved_organ >= 1 || scalar @unresolved_stage >= 1 ){
        ANAT1:
        for my $anat ( sort @unresolved_organ ){
            #NOTE if several anatomy terms, life_stage is uncertain: maybe linked to some anatomy terms but not all of them
            # so life_stage(s) is mapped only to ID root stage: UBERON:0000104
            # See https://gitlab.isb-sib.ch/Bgee/bgee_pipeline/issues/72
            my @stage = @unresolved_stage;
            if ( scalar @unresolved_organ > 1 ){
                @stage = ($id_root_stage);
            }
            STAGE1:
            for my $stage ( sort @stage ){
                my $anatomy_id = $anat =~ /^(.+?)__/ ? $1 : $anat;
                print_TSV($expression, $id, $anat, $anatomy_id, $stage);
            }
        }

        STAGE2:
        for my $stage ( sort @unresolved_stage ){
            #NOTE if several stage terms, anatomy is uncertain: maybe linked to some stage terms but not all of them
            # so anatomy(ies) is mapped only to ID root anatomy: UBERON:0000465 "material anatomical entity"
            # See https://gitlab.isb-sib.ch/Bgee/bgee_pipeline/issues/72
            my @anat = @unresolved_organ;
            if ( scalar @unresolved_stage > 1 ){
                @anat = ($id_root_organ);
            }
            ANAT2:
            for my $anat ( sort @anat ){
                my $stage_id   = $stage =~ /^(.+?)__/ ? $1 : $stage;
                my $anatomy_id = $anat =~ /^(.+?)__/  ? $1 : $anat;
                print_TSV($expression, $id, $anat, $anatomy_id, $stage_id);
            }
        }
    }
}

exit 0;


sub get_quality {
    my ($term) = @_;

    my %quality = ( ''          => 'high quality',
                    'certain'   => 'high quality',
                    'partial'   => 'high quality',
                    'uncertain' => 'poor quality',
                  );

    return $quality{$term} || die "Invalid quality term [$term]\n";
}

sub print_TSV {
    my ($expression, $id, $annotation, $anatomy_id, $stage_id) = @_;

    # Get anatomy and stage ids + quality info
    my $qual_term = '';
    if ( $annotation =~ /__(\w+)/ && $1 ne 'Life_stage' && $1 ne 'Anatomy_term' ){
        $qual_term = $1;
    }

    # Get quality
    my $quality = get_quality(lc $qual_term);
    # Deal with NOT qualifier, available now
    if ( $anatomy_id =~ /^NOT/ || $stage_id =~ /^NOT/ ){
        $quality    = 'absent';
        $anatomy_id =~ s{^NOT}{};
        $stage_id   =~ s{^NOT}{};
    }

    # Build stage from .. to
    if ( $stage_id ne $id_root_stage ){
        #FIXME $doneStages->{ $id_root_stage.','.$id_root_stage } returns always empty ???
        $stage_id = $doneStages->{ $stage_id.','.$stage_id } || '';
    }
    if ( $stage_id eq '' || $stage_id =~ /Could not find any OWLClass corresponding to/ ){
        warn "No stage reference for [$stage_id] in [$id / $expression->{$id}->{'ref'} | @{$expression->{$id}->{'anatomy'}} | @{$expression->{$id}->{'stage'}}] {$annotation} {$quality}\n";
        return;
    }

    $expression->{$id}->{'ref'} = $expression->{$id}->{'ref'} ne '' ? $expression->{$id}->{'ref'} : 'WormBase:NO_REF';

    # Anatomy to Uberon
    if ( $anatomy_id ne $id_root_organ ){
        #FIXME $doneAnat->{$id_root_organ} returns always empty ???
        $anatomy_id = $doneAnat->{$anatomy_id} || '';
    }
    if ( $anatomy_id eq '' || $anatomy_id =~ /Could not find any OWLClass corresponding to/ ){
        warn "Problem with anatId for [$anatomy_id] in [$id / $expression->{$id}->{'ref'} | @{$expression->{$id}->{'anatomy'}} | @{$expression->{$id}->{'stage'}}] {$annotation} {$quality}\n";
        return;
    }

    # Presence
    # Not_in == absent now, but independant of quality that should be high for most absent of expresssion reported!
    my $presence = $quality eq 'high quality' ? 'present'
                 : $quality eq 'poor quality' ? 'present'
                 : $quality eq 'absent'       ? 'absent'
                 :                              '';
    if ( $presence eq '' ){
        warn "No defined presence for [$quality] in [$id / $expression->{$id}->{'ref'}]\n";
        return;
    }
    $quality = 'high quality'  if ( $presence eq 'absent' );

    #data_source  inSituExperimentId  inSituEvidenceId  organId  stageId  geneId  detectionFlag  inSituData  linked  speciesId  strain  sex  [ExprDesc]
    print join("\t", $data_source_id,
                     $expression->{$id}->{'ref'},
                     $id,
                     $anatomy_id,
                     $stage_id,
                     $expression->{$id}->{'gene'},
                     $presence,
                     $quality,
                     '',
                     6239,                            #FIXME C. elegans only ????
                     $expression->{$id}->{'strain'} eq '' ? $Utils::NOT_ANNOTATED_STRAIN : $expression->{$id}->{'strain'},
                     'not annotated',
              ), "\n";
    return;
}


=pod
   Data example from expr_pattern.ace.20140526 file extracted by WormBase helper Daniela Raciti <draciti@caltech.edu>

Expr_pattern : "Expr10000"
Anatomy_term    "WBbt:0003679" Certain
Anatomy_term    "WBbt:0005772" Certain
Anatomy_term    "WBbt:0005785" Certain
Gene    "WBGene00018488"
In_Situ
In_situ "Primers for the generation of the acs-1 cDNA PCR fragment and subsequent sense and anti-sense hybridization probes were: F, atgtcacaagtggccgcaatggacc and R, aagctcggagaaatgagagac."
Life_stage      "WBls:0000041"
Pattern "In addition to the neuronal and intestinal expression previously detected using an acs-1 promoter-driven GFP expression construct (Kniazeva et al. 2004), prominent GFP fluorescence was observed in the somatic gonad of adults but not in that of larvae. Results from RNA in situ hybridization confirmed acs-1 expression in the somatic gonad, while no expression was detected in the germline."
Reference       "WBPaper00040893"
Reflects_endogenous_expression_of       "WBGene00018488"
Reporter_gene
Transgene       "WBTransgene00015406"
=cut

