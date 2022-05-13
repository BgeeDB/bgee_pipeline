#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Roelli Patrick, created 19.06.12
# Insertion of in-situ data from BDGP
# USAGE: perl insert_in_situ_bdgp.pl <BDGP_mapping_to_FBbt.tsv> <BDGP annotations CSV file>
#
#NOTE Write 'tsv_BDGP' by itself at the end
#####################################

use Getopt::Long;
use File::Slurp;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

$| = 1;

# Define arguments & their default value
my ($bgee_connector, $bdgp_connector) = ('', '');
my ($debug, $mapping, $annotation, $stagecorresp)    = (0, '', '', '');
my ($Aport, $Sport) = (0, 0);
my %opts = ('debug'          => \$debug,            # more verbose
            'bgee=s'         => \$bgee_connector,   # Bgee connector string
            'bdgp=s'         => \$bdgp_connector,   # BDGP connector string
            'mapping=s'      => \$mapping,          # BDGP_terms_to_FBbt_terms.xls
            'annotation=s'   => \$annotation,       # insitu_annot.csv from BDGP web site
            'stagecorresp=s' => \$stagecorresp,     # Stage correspondance
            'Aport=i'        => \$Aport,            # Anatomy mapper socket port
            'Sport=i'        => \$Sport,            # Stage mapper socket port
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $bdgp_connector eq '' || $mapping eq '' || $annotation eq '' || $stagecorresp eq '' || $Aport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -bdgp=\$(BDGPCMD) -mapping=\$(BDGP2FBBT_MAPPING_FILE) -annotation=insitu_annot.csv -stagecorresp=\$(STAGECORRESP_FILE) -Aport=\$(IDMAPPINGPORT) -Sport=\$(INBETWEENSTAGESPORT)
\t-bgee          Bgee   connector string
\t-bdgp          BDGP   connector string
\t-debug         More verbose
\t-mapping       BDGP_terms_to_FBbt_terms.xls
\t-annotation    insitu_annot.csv from BDGP web site
\t-stagecorresp  Stage correspondance
\t-Aport         Anatomy mapper socket port
\t-Sport         Stage   mapper socket port
\n";
    exit 1;
}

# Bgee db connection
my $dbh  = Utils::connect_bgee_db($bgee_connector);
# BDGP db connection
my $bdgp = Utils::connect_bgee_db($bdgp_connector);

## Omitted_terms are those that are not in term table on "purpose". The list was given by Erwin Frise.
## just for the record... this hash is never used!
#my %omitted_terms = ('130' => 1,  '131' => 1,  '207' => 1,  '232' => 1,  '233' => 1,  '265' => 1,
#                     '502' => 1,  '509' => 1,  '515' => 1,  '532' => 1,  '551' => 1,  '552' => 1,
#                     '553' => 1,  '555' => 1,  '556' => 1,  '576' => 1,  '595' => 1,  '596' => 1,
#                    );

# Retrieve all genes in Bgee (some genes may be out dated in BDGP)
my %bgee_genes;
# In the current release there are no cgname found in bdgp matching a geneNameSynonym in bgee
my $selGene = $dbh->prepare('SELECT geneId, geneName FROM gene WHERE speciesId = ?');
#FIXME Other Drosophila species but melanogaster?
$selGene->execute(7227)  or die $selGene->errstr;
while ( my @data = $selGene->fetchrow_array ){
    $bgee_genes{'FlyBase'}{$data[0]} = $data[1];
    $bgee_genes{'Name'}{$data[1]}    = $data[0];
}
$selGene->finish;

# Retrieve experiment infos in the database from main table in bdgp with id as key
my %experiments; # experiments is the 'main' in the new database and we will use est_id as the exp id on bgee

# CONTROL: Take only experiments which have a gene, an annot, an annot_term and a term linked to it.
# We need an organ to create a spot in inSituSpot.
# The cgname allows us to retrieve FBgn from bgee when the flybase_id is missing.
# We can have either cgname or flybase_id being empty, not both
# M.id:         id of the line in the 'main' table, will be used to check for duplicated est_id
# M.cgname:     name of the gene
# M.est_id:     will be used as experimentId in Bgee, can be used in URLs to retrieve data from BDGP
# M.flybase_id: ID of the gene
# annot table:      evidence. One 'annot' groups different terms for a same stage
# annot_term table: links the annot table (with stage information) to the term table
# term table:       contains anatomical structures with their name in the column 'go_term'
my $selExp = $bdgp->prepare('SELECT DISTINCT M.id, M.cgname, M.est_id, M.flybase_id, M.species, M.fly_strain, M.assay
FROM main AS M
INNER JOIN annot      AS A  ON A.main_id   = M.id
INNER JOIN annot_term AS AT ON AT.annot_id = A.id
INNER JOIN term       AS T  ON AT.term_id  = T.id
WHERE (M.cgname!="" OR M.flybase_id!="") AND
(M.flybase_id IS NOT NULL AND M.cgname IS NOT NULL)'); # We take experiments where both are NOT NULL, not the case in the current release
$selExp->execute()  or die $selExp->errstr;

# Store the complete list of cg names and flybase IDs used for the mapping to bgee genes
my %cgNames;
my %flybaseIds;
while ( my @data = $selExp->fetchrow_array ){
    @data = map { Utils::trim($_) } @data;
    next  if ( $data[0] eq '' || $data[0] eq 'NULL' );

    $experiments{$data[0]}{'CG'}      = $data[1];
    $experiments{$data[0]}{'species'} = $data[4];
    $experiments{$data[0]}{'strain'}  = $data[5];
    $experiments{$data[0]}{'assay'}   = $data[6];
    if ( $experiments{$data[0]}{'CG'} eq 'NULL' ){
        $experiments{$data[0]}{'CG'} = '';
    }

    $experiments{$data[0]}{'est_id'}  = $data[2];
    $experiments{$data[0]}{'FlyBase'} = $data[3]; # /!\ Attention, there are two upper caracter in the word /!\
    if ( $experiments{$data[0]}{'FlyBase'} eq 'NULL' ){
        $experiments{$data[0]}{'FlyBase'} = '';
    }

    # if cg name defined
    if ( $experiments{$data[0]}{'CG'} ne '' ){
        $cgNames{$experiments{$data[0]}{'CG'}} = $experiments{$data[0]}{'FlyBase'};
    }
    # if flybase ID defined
    if ( $experiments{$data[0]}{'FlyBase'} ne '' ){
        $flybaseIds{$experiments{$data[0]}{'FlyBase'}} = $experiments{$data[0]}{'CG'};
    }
}
$selExp->finish;


my $experiments_count = keys(%experiments);
# Since BDGP doesn't use the same stages ids as us, we had to use a mapping stage_id -> dev ontology id.
# This mapping is contained in the file ../../../curation/expression_data/in_situ/bdgp/stages_correspondence_new.txt
my %stage_correspondence;
my $getId = $dbh->prepare('SELECT stageId FROM stageXRef WHERE stageXRefId = ?');
my %tsv = %{ Utils::read_spreadsheet("$stagecorresp", "\t", 'csv', '', 1)}; # stages_correspondence_new.txt

# Get dev stages mapping
my @Stages     = map { ($_, $_) } @{ $tsv{'stageId'} }; # Pseudo start-end by doubling stage
my $doneStages = Utils::get_in_between_stages(\@Stages, $Sport);

# stageId  BDGPId  stageName
for my $line ( 0..$#{$tsv{'stageId'}} ){
    my $real_stage_id = $tsv{'stageId'}[$line];
    $getId->execute($real_stage_id)  or die $getId->errstr;
    if ( my @resultId = $getId->fetchrow_array ){
        $real_stage_id = $resultId[0];
    }

    $stage_correspondence{ $tsv{'BDGPId'}[$line] }{'FBdv'}       = $real_stage_id;
    $stage_correspondence{ $tsv{'BDGPId'}[$line] }{'stage_name'} = $tsv{'stageName'}[$line];
}
$getId->finish;


# We now need to retrieve annotations from the BDGP CSV file,
# as it is the only way to retrieve correct "filtered" anatomical terms
#csv file: cgname, flybase_name, flybase_id, stage_id, go_term (=name of the anatomical term)
my %tsv2 = %{ Utils::read_spreadsheet("$annotation", ',', 'csv', '"', 1)}; # insitu_annot.csv
#csvAnnot: cgname -> stage_id -> go_term
my %csvAnnot = ();
# e.g. "a10","CG6642","FBgn0011293",1,"no staining"
for my $line ( 0..$#{$tsv2{0}} ){
    $csvAnnot{ Utils::trim($tsv2{0}[$line]) }{ Utils::trim($tsv2{3}[$line]) }{ Utils::trim($tsv2{4}[$line]) } = 1;
}


# Retrieve the evidence which are in the annot table in bdgp database along with the stage.
# We use a hash of hash for performance purpose in the iterations.
my %evidences;
# CONTROL: Do not take evidence without links to experiments. Not the case in current release.
my $selStage = $bdgp->prepare('SELECT A.id, M.id, A.stage FROM annot AS A INNER JOIN main AS M ON M.id=A.main_id');
$selStage->execute()  or die $selStage->errstr;
while ( my @data = $selStage->fetchrow_array ){
    $evidences{$data[1]}{$data[0]} = $data[2];
    # always check that the stage actually exists in our mapping
    if ( !defined $stage_correspondence{$data[2]}{'FBdv'} ){
        warn 'Warning: no correspondence found for BDGP stage ['.$data[2].
              "], you should fix that in the stages_correspondence file\n";
    }
}
$selStage->finish;
my $evidences_count;
for my $main_id ( keys %evidences ){
    $evidences_count += keys(%{$evidences{$main_id}});
}


# Retrieve the annot_term entries as spots, with the column term_id linking to the term table.
my %spots;
# CONTROL: Do not take spots that have no evidence linked to them or no organs.

# field "annotator": when equals to "unified", only these terms should be considered for the current annotation
# ("win" over all other annotations)

# We also use only annotations present in the CSV annotations file (filtered terms)
# this is why we need A.stage, M.cgname, and T.go_term
my $selAnnot = $bdgp->prepare('SELECT AT.id, A.id, AT.term_id, AT.annotator, A.stage, M.cgname, T.go_term
FROM annot_term AS AT
INNER JOIN annot AS A ON A.id       = AT.annot_id
INNER JOIN main  AS M ON A.main_id  = M.id
INNER JOIN term  AS T ON AT.term_id = T.id
ORDER BY A.id');

# The ids of each annot_term will be used as the inSituExpressionPatternId in bgee.inSituSpot
$selAnnot->execute()  or die $selAnnot->errstr;
my $previousAnnotId = undef;
my $annotId         = undef;
my $isUnified       = 0;
my %tempSpots       = ();
my $toDo            = 1;
while ( $toDo ){
    my @data = undef;
    if ( @data = $selAnnot->fetchrow_array ){
        # We need first to retrieve all annotations for an annot_id,
        # in order to see if there are some "unified" annotator used, subsuming all other annotations
        $annotId = $data[1];
    }
    else {
        $annotId = undef;
        # one last additional loop to examine the last annot_id
        $toDo = 0;
    }

    if ( defined $previousAnnotId && ((defined $annotId && $previousAnnotId ne $annotId) || !$toDo) ){ # one additional iteration to examine the last annot_id
        # add the temp spots of the previous annot_id to the real spots
        for my $annotTermId ( keys %tempSpots ){
            # add a spot only if there were no "unified" annotator for this annot_id,
            # or if it is currently a "unified" spot
            if ( !$isUnified || $tempSpots{$annotTermId}{'annotator'} eq 'unified' ){
                # and add a spot only if the term is present in the CSV annotations file
                # (acts as a filter)
                # csvAnnot: cgname -> stage_id -> go_term
                my $cgname  = $tempSpots{$annotTermId}{'cgname'};
                my $stage   = $tempSpots{$annotTermId}{'stage'};
                my $go_term = $tempSpots{$annotTermId}{'go_term'};
                if ( defined $csvAnnot{$cgname}{$stage}{$go_term} ){
                    $spots{$previousAnnotId}{$annotTermId}{'term_id'} = $tempSpots{$annotTermId}{'term_id'};
                }
#                else {
#                    warn "Unmapped stage: [$stage]\n";
#                }
            }
        }

        # reinitialize
        $isUnified = 0;
        %tempSpots = ();
    }
    if ( $toDo ){
        $tempSpots{$data[0]}{'term_id'}   = $data[2];
        $tempSpots{$data[0]}{'annotator'} = $data[3];
        if ( $tempSpots{$data[0]}{'annotator'} eq 'unified' ){
            $isUnified = 1;
        }
        $tempSpots{$data[0]}{'go_term'} = $data[6];
        $tempSpots{$data[0]}{'cgname'}  = $data[5];
        $tempSpots{$data[0]}{'stage'}   = $data[4];
    }
    $previousAnnotId = $annotId;
}
$selAnnot->finish;
my $spots_count;
for my $annot_id( keys %spots ){
    $spots_count += keys(%{$spots{$annot_id}});
}



# Retrieve anatomical ontology ids from the manual mapping
my %tsv3 = %{ Utils::read_spreadsheet("$mapping", "\t", 'csv', '', 1)}; # BDGP_terms_to_FBbt_terms.xls

# Get anatomy mapping
my @Anat       = grep { $_ ne 'none' } @{ $tsv3{'FBbt_id'} };
my $doneAnat   = Utils::get_anatomy_mapping(\@Anat, $Aport);

my %organs;
# #stage_name   BGDP_term   BDGP_id   FBbt_id   FBbt_name   matches_count   match_source
for my $line ( 0..$#{$tsv3{'stage_name'}} ){
    # this term list will give us the well defined terms in BDGP hence, those who are healthy and usable in bgee
    if ( $tsv3{'FBbt_id'}[$line] ne 'none' ){
        if ( defined $organs{ $tsv3{'BDGP_id'}[$line] }{ $tsv3{'stage_name'}[$line] } ){
            warn "Warning, ambiguous mapping in the mapping file for BDGP term ID [$tsv3{'BDGP_id'}[$line]]\n";
        }
        $organs{ $tsv3{'BDGP_id'}[$line] }{ $tsv3{'stage_name'}[$line] } = $doneAnat->{ $tsv3{'FBbt_id'}[$line] } || '';
    }
}

my $organs_count = 0;
for my $organ ( keys %organs ){
    $organs_count += keys(%{$organs{$organ}});
}
$bdgp->disconnect;

# SUMMARY
printf("%9d  experiments where found\n%9d  evidences   where found\n%9d  spots       where found\n%9d  BDGP organs are mapped to FBbt organs\n", $experiments_count, $evidences_count, $spots_count, $organs_count);



################################
# Filtering data
################################
print "Filtering data\n"  if ( $debug );

# Filter double est_ids entries.
# We have found that some est_ids are in more than one experiment.
# As we use the est_id as key for differentiating the experiments we have to delete those with the same est_id
my %none_uniq;
my %est_ids;
# First, store each est_id in a hash, and count the number of experiment
for my $exp ( keys %experiments ){
    $est_ids{$experiments{$exp}{'est_id'}}++;
}
# We have to iterate %experiments twice,
# otherwise we would realize that an est_id is associated to two experiments only after iterating the first experiment...
# so we would delete experiments only starting at the second one, letting pass the first one...
for my $exp ( keys %experiments ){
    if ( $est_ids{$experiments{$exp}{'est_id'}} > 1 ){
        # If est_id associated to more than one experiment-> store the exp id
        $none_uniq{$exp} = 1;
    }
}
for my $double_entry_exp_id ( keys %none_uniq ){ # Delete each doubled entry experiments
    warn "Removing main_id [$double_entry_exp_id] because of est_id used several times: ".$experiments{$double_entry_exp_id}{'est_id'}."\n";
    delete $experiments{$double_entry_exp_id};
}

# Filter missing organs in spots and evidences linked to them.
# There are two modes to run the filtering of missing spots and evidences.
# default mode = mode1: If one spot from an evidence is missing, delete the whole evidence.
# mode2: Delete only the spots.
# You can choose how you want to filter in changing the $filter_mode variable to mode1 or mode2.
my $filter_mode = 'mode1';

# Do the actual filtering.
my %expToDelete;
my %evidenceToDelete;
my %spotToDelete;
my %termNotFound;
# First we have to filter the experiments based on cgname and flybase_id.
for my $exp ( keys %experiments ){
    # There are 4 cases:
    # 1) cgname is empty and flybase_id is not-> delete it
    # 2) cgname is not empty and flybase_id is-> attempt to get flybase from bgee, if not possible -> delete
    # 3) cgname and Flybase_id are empty-> delete it
    # 4) cgname and flybase_id are not empty-> OK but check for ambiguity with bgee
    # case 1:
    if ( $experiments{$exp}{'CG'} eq '' && $experiments{$exp}{'FlyBase'} ne '' ){
        # try to get the gene by the flybase ID
        if ( defined $bgee_genes{'FlyBase'}{$experiments{$exp}{'FlyBase'}} ){
            my $bgeeName = $bgee_genes{'FlyBase'}{$experiments{$exp}{'FlyBase'}};
            # check whether the gene name in bgee corresponds to a cg name in bdgp, but used for another flybase ID
            if ( defined $cgNames{$bgeeName} && $cgNames{$bgeeName} ne '' && $cgNames{$bgeeName} ne $experiments{$exp}{'FlyBase'} ){
                $expToDelete{$exp} = 1;
                warn "Removing main_id [$exp] because of inconsistent CG name for a gene retrieved by flybase ID\n";
                next;
            }
#            else {
#                #print "Ye, it works\n";
#            }
        }
    }
    #case 2:
    elsif( $experiments{$exp}{'CG'} ne '' && $experiments{$exp}{'FlyBase'} eq '' ){
        if ( defined $bgee_genes{'Name'}{$experiments{$exp}{'CG'}} ){
            my $bgeeId = $bgee_genes{'Name'}{$experiments{$exp}{'CG'}};
            # check whether the ID in bgee corresponds to a flybase ID in bdgp, but used for another CG name
            if ( defined $flybaseIds{$bgeeId} && $flybaseIds{$bgeeId} ne '' && $flybaseIds{$bgeeId} ne $experiments{$exp}{'CG'} ){
                $expToDelete{$exp}=1;
                warn "Removing main_id [$exp] because of inconsistent flybase ID for a gene retrieved by CG name\n";
                next;
            }
            $experiments{$exp}{'FlyBase'} = $bgeeId;
        }
        else{
            $expToDelete{$exp} = 1;
            warn "Removing main_id [$exp] because of flybase_id corresponding to cgname not found in bgee\n";
            next;
        }
    }
    #case 3:
    elsif( $experiments{$exp}{'CG'} eq '' && $experiments{$exp}{'FlyBase'} eq '' ){
        $expToDelete{$exp} = 1;
        warn "Removing main_id [$exp] because of no cgname and no flybase_id\n";
        next;
    }
    #case 4:
    elsif( $experiments{$exp}{'CG'} ne '' && $experiments{$exp}{'FlyBase'} ne '' ){
        # If same cgname but different flybase_id->ambigous->delete it
        if ( defined $bgee_genes{'Name'}{$experiments{$exp}{'CG'}} && $experiments{$exp}{'FlyBase'} ne $bgee_genes{'Name'}{$experiments{$exp}{'CG'}} ){
            warn "Removing main_id [$exp] because of ambiguous cgname and flybase_id\n";
            $expToDelete{$exp} = 1;
            next;
        }
        # If same Flybase but different cgname -> ambigous -> delete it
        if( defined $bgee_genes{'FlyBase'}{$experiments{$exp}{'FlyBase'}} && $experiments{$exp}{'CG'} ne $bgee_genes{'FlyBase'}{$experiments{$exp}{'FlyBase'}} ){
            warn "Removing main_id [$exp] because of ambiguous cgname and flybase_id\n";
            $expToDelete{$exp} = 1;
            next;
        }
        if ( !defined $bgee_genes{'Name'}{$experiments{$exp}{'CG'}} ){
            $expToDelete{$exp} = 1;
            warn "Removing main_id [$exp] because gene not present in Bgee\n";
            next;
        }
    }
    elsif ( defined $bgee_genes{'CG'}{$experiments{$exp}{'FlyBase'}} ){
        # Check if we have the geneId and let it pass
        warn "My developer is shity, assertion error (only 4 possible states with 2 variables with 2 states).\n";
    }
    else {
        warn "My developer is shity, assertion error (only 4 possible states with 2 variables with 2 states)...\n";
        # If no match found -> delete the experiment
        $expToDelete{$exp} = 1;
        warn "Removing main_id [$exp] because we could not find mapping for";
        if ( defined $experiments{$exp}{'FlyBase'} ){
            warn ' Flybase ID: '.$experiments{$exp}{'FlyBase'};
        }
        if ( defined defined $experiments{$exp}{'CG'} ){
            warn ' CG name: '.$experiments{$exp}{'CG'};
        }
        warn "\n";
        next;
    }

    # From now on the hash of hash will be useful.
    # Since we do not want to iterate over all spots of all evidences,
    # we use the primary keys first and hence gain time.
    # It is also a CONTROL point of link between spots, evidence and experiments.

    # count the number of evidence for each experiment.
    my $evidenceCount = 0;
    for my $evidence ( keys %{$evidences{$exp}} ){ # For each evidence with exp as key do...
        $evidenceCount++;
        # count the number of spots for each evidence.
        my $spotCount = 0;
        my $evidenceRemoved = 0;
        for my $spot ( keys %{$spots{$evidence}} ){ #For each spot from one SPECIFIC evidence
            $spotCount++;
            # If the organ is not defined in our mapping...
            if ( !defined $organs{$spots{$evidence}{$spot}{'term_id'}}{$stage_correspondence{$evidences{$exp}{$evidence}}{'stage_name'}} ){
                $termNotFound{$spots{$evidence}{$spot}{'term_id'}}{$stage_correspondence{$evidences{$exp}{$evidence}}{'stage_name'}} = 1;
                # delete spot or evidence based on $filter_mode
                if ( $filter_mode eq 'mode2' ){
                    $spotToDelete{$evidence}{$spot} = 1;
                    $spotCount--;
                    warn "Removing annot_term id [$spot] for annot id [$evidence] because no mapping found for ".$spots{$evidence}{$spot}{'term_id'}."\n";
                }
                elsif ( !$evidenceRemoved ){
                    warn "Removing annot_term id [$spot] for annot id [$evidence] because no mapping found for ".$spots{$evidence}{$spot}{'term_id'}."\n";
                    $evidenceToDelete{$exp}{$evidence} = 1;
                    $evidenceCount--;
                    $evidenceRemoved = 1;
                }
            }
        }
        # If there are no more spots in one evidence, we have to delete it.
        if ( $spotCount == 0 ){
            $evidenceToDelete{$exp}{$evidence} = 1;
            $evidenceCount--;
            warn "Removing annot id [$evidence] because all spots from this evidence were removed\n";
        }
        elsif ( $spotCount < 0 ){
            die "My developer is shity, check script near line 314 for debug\n";
        }
    }
    # If there are no more evidence in one experiment, we have to delete it.
    if ( $evidenceCount == 0 )  {
        $expToDelete{$exp} = 1;
        warn "Removing main_id [$exp] because all evidences from this experiment were removed\n";
    }
    elsif ( $evidenceCount < 0 ){
        die "My developer is shity, check script near line 322 for debug\n";
    }
}

for my $exp ( keys %expToDelete ){
    delete $experiments{$exp};
}
for my $exp ( keys %evidenceToDelete ){
    for my $evidence ( keys %{$evidenceToDelete{$exp}} ){
        delete $evidences{$exp}{$evidence};
    }
}
for my $evidence ( keys %spotToDelete ){
    for my $spot ( keys %{$spotToDelete{$evidence}} ){
        delete $spots{$evidence}{$spot};
    }
}

# SUMMARY post filtering
my $experiments_left = keys(%experiments);
my $evidences_left   = 0;
my $spots_left       = 0;
for my $exp ( keys %experiments ){
    $evidences_left += keys(%{$evidences{$exp}});
    for my $evidence ( keys %{$evidences{$exp}} ){
        $spots_left += keys(%{$spots{$evidence}});
    }

}
printf("After filtering:\n%9d  experiments are left\n%9d  evidences   are left\n%9d  spots       are left\n", $experiments_left, $evidences_left, $spots_left);


###########################
# Prepare data for insertion into bgee
###########################
# BDGP IDs of the 'no staining' terms, to insert 'no expression' data
my %no_staining = ('489' => 1,  '490' => 1,  '491' => 1,  '492' => 1,  '493' => 1,  '494' => 1,);

# Prepare insertion of data source id for BDGP
my $selSrc = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName = ?');
$selSrc->execute('BDGP')  or die $selSrc->errstr;
my $data_source_id = $selSrc->fetchrow_array;
$selSrc->finish;
$dbh->disconnect;
# Since Bgee_v10, dataSources are all inserted at once at the beginning of the pipeline
die "Data source not found\n"  if ( !defined $data_source_id );


my $spot_id = 1; # Initiate the spot id
my $output = join("\t", '#data_source', qw(inSituExperimentId  inSituEvidenceId  organId  stageId  geneId  detectionFlag  inSituData  linked  speciesId  strain  sex))."\n";
EXP:
for my $exp ( keys %experiments ){
    next EXP  if ( $experiments{$exp}{'assay'}  ne 'in situ hybridization' ); # In situ only,        NOT 'antibody staining'
    next EXP  if ( $experiments{$exp}{'strain'} ne 'wt' );                    # Keep only wild-type, NOT 'transgenic line'

    # Get species taxid
    my $speciesId = $experiments{$exp}{'species'} eq 'D. melanogaster'  ? 7227 # Only D. melanogaster seen in bdgp database for now, but others are declared
                  : $experiments{$exp}{'species'} eq 'D. pseudoobscura' ? 7237
                  : $experiments{$exp}{'species'} eq 'D. virilis'       ? 7244
                  : $experiments{$exp}{'species'} eq 'D. simulans'      ? 7240
                  :                                                       '';
    if ( $speciesId eq '' ){
        warn "Invalid/not registered species [$experiments{$exp}{'species'}]\n";
        next EXP;
    }

    # Get Strain info
    my $strain = $Utils::WILD_TYPE_STRAIN; # Forced to be wild type only above

    # No sex info in BDGP, because mostly on embryo
    my $sex = 'not annotated';

    EVID:
    for my $evidence ( keys %{$evidences{$exp}} ){
        SPOT:
        for my $spot ( keys %{$spots{$evidence}} ){
            if ( exists $no_staining{$spots{$evidence}{$spot}{'term_id'}} ){
                # we insert this 'no staining' term only if there are no other terms beside 'no staining' in this evidence,
                # (e.g. it is not an evidence with both a "brain" term and a "no staining" term)
                # otherwise we cannot just map 'no staining' to the root of the anatomy.
                # We cannot just count the size of the spot hash, sometimes there are several 'no staining' terms
                my $insert = 1;
                for my $spot2 ( keys %{$spots{$evidence}} ){
                    if ( !exists $no_staining{$spots{$evidence}{$spot2}{'term_id'}} ){
                        $insert = 0;
                        last;
                    }
                }
                if ( $insert ){
                    my $upperStage = $stage_correspondence{$evidences{$exp}{$evidence}}{'FBdv'}.','.$stage_correspondence{$evidences{$exp}{$evidence}}{'FBdv'};
                    $output .= join("\t", $data_source_id,
                                          'BDGP_'.$experiments{$exp}{'est_id'},
                                          'BDGP_'.$evidence,
                                          $organs{$spots{$evidence}{$spot}{'term_id'}}{$stage_correspondence{$evidences{$exp}{$evidence}}{'stage_name'}},
                                          $doneStages->{$upperStage},
                                          $experiments{$exp}{'FlyBase'},
                                          'absent',
                                          'high quality',
                                          '',
                                          $speciesId,
                                          $strain,
                                          $sex,
                                   )."\n";
                    $spot_id++;
                }
            }
            else {
                my $upperStage = $stage_correspondence{$evidences{$exp}{$evidence}}{'FBdv'}.','.$stage_correspondence{$evidences{$exp}{$evidence}}{'FBdv'};
                $output .= join("\t", $data_source_id,
                                      'BDGP_'.$experiments{$exp}{'est_id'},
                                      'BDGP_'.$evidence,
                                      $organs{$spots{$evidence}{$spot}{'term_id'}}{$stage_correspondence{$evidences{$exp}{$evidence}}{'stage_name'}},
                                      $doneStages->{$upperStage},
                                      $experiments{$exp}{'FlyBase'},
                                      'present',
                                      'high quality',
                                      '',
                                      $speciesId,
                                      $strain,
                                      $sex,
                               )."\n";
                $spot_id++;
            }
        }
    }
}
write_file('tsv_BDGP', $output);

exit 0;

#TODO
# - check URL dataSource
#   Better to link to a gene summary page like:
#   http://www.fruitfly.org/cgi-bin/ex/bquery.pl?qtype=report&find=CG3187&
# - Quality -> not defined so always high quality
# - some anatomical terms are obsolete: a new ID can be found sometimes in the OBO file (field "consider:")

