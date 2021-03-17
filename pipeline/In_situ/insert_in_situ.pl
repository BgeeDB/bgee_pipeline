#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use List::MoreUtils qw{uniq};

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

$| = 1;

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($debug)          = (0);
my ($tsv)            = ('');
my ($sex_info)       = ('');
my %opts = ('debug'      => \$debug,            # more verbose
            'bgee=s'     => \$bgee_connector,   # Bgee connector string
            'tsv=s'      => \$tsv,              # tsv mapping file
            'sex_info=s' => \$sex_info,         # generated_files/uberon/uberon_sex_info.tsv
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $tsv eq '' || $sex_info eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -tsv=tsv_flybase -sex_info=\$(UBERON_SEX_INFO_FILE_PATH)
\t-bgee      Bgee    connector string
\t-debug     More verbose
\t-tsv       TSV mapping file
\t-sex_info  file containing sex-related info about anatomical terms
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


# Read TSV data
my %tsv = %{ Utils::read_spreadsheet("$tsv", "\t", 'csv', '', 1) };
#data_source  inSituExperimentId  inSituEvidenceId  organId  stageId  geneId  detectionFlag  inSituData  linked  speciesId  strain  sex  [ExprDesc]


# Fill inSituExperiment table
my $data_source_id   = $tsv{'data_source'}[0]  || die "Problem with data in [$tsv] file\n";
my $species          = $tsv{'speciesId'}[0]    || die "Problem with speciesId in [$tsv] file\n";


# Get gene_mapping to bgeeGeneId
my $gene_mapping = Utils::query_bgeeGene($dbh, $species);

# Get already known conditions
my $conditions         = Utils::query_conditions($dbh);

# Get simpler (upper level) stage equivalences
my $stage_equivalences = Utils::get_stage_equivalences($dbh);

# Get sex-related information needed for sub 'insert_get_condition'
my $anatSexInfo    = Utils::get_anat_sex_info($sex_info);
my $speciesSexInfo = Utils::get_species_sex_info($dbh);

my %expr_seen;
my $insExpr = $dbh->prepare('INSERT INTO inSituExperiment (inSituExperimentId, inSituExperimentDescription, dataSourceId) VALUES (?, ?, ?)');
for my $line ( 0..$#{$tsv{'inSituExperimentId'}} ){
    if ( ! exists $expr_seen{ $tsv{'inSituExperimentId'}[$line] } && exists $gene_mapping->{ $tsv{'geneId'}[$line] } ){
        my $desc = $tsv{'ExprDesc'}[$line] || '';
        $insExpr->execute($tsv{'inSituExperimentId'}[$line], $desc, $data_source_id)  or die "[$tsv{'inSituExperimentId'}[$line]] ".$insExpr->errstr;
        print "[$tsv{'inSituExperimentId'}[$line]]\t$desc\t[$data_source_id]\n"  if ( $debug );
        $expr_seen{ $tsv{'inSituExperimentId'}[$line] }++;
    }
}
$insExpr->finish;


# Fill inSituEvidence/inSituSpot tables
#NOTE for FlyBase evidenceDistinguishable is currently 0

my $inSituSpotId = 0;
my %evid_seen;
my $insEvid = $dbh->prepare('INSERT INTO inSituEvidence (inSituEvidenceId, inSituExperimentId, evidenceDistinguishable, inSituEvidenceUrlPart) VALUES (?, ?, ?, ?)');
my $insSpot = $dbh->prepare('INSERT INTO inSituSpot (inSituSpotId, inSituEvidenceId, conditionId, bgeeGeneId, detectionFlag, inSituData) VALUES (?, ?, ?, ?, ?, ?)');
for my $line ( 0..$#{$tsv{'data_source'}} ){
    if ( !exists $gene_mapping->{ $tsv{'geneId'}[$line] } ){
        warn "Unknown gene id (too old or too up-to-date) [$tsv{'geneId'}[$line]]\n";
        next;
    }
    # inSituEvidence
    if ( ! exists $evid_seen{ $tsv{'inSituExperimentId'}[$line].'-'.$tsv{'inSituEvidenceId'}[$line] } ){
        my $evidenceDistinguishable = $tsv =~ /flybase$/ ? 0 : 1;
        my $linked                  = $tsv{'linked'}[$line] || '';

        $insEvid->execute($tsv{'inSituEvidenceId'}[$line], $tsv{'inSituExperimentId'}[$line], $evidenceDistinguishable, $linked)  or die "[$tsv{'inSituEvidenceId'}[$line]] ".$insEvid->errstr;
        print "[$tsv{'inSituEvidenceId'}[$line]]\t[$tsv{'inSituExperimentId'}[$line]]\t[$evidenceDistinguishable]\t[$linked]\n"  if ( $debug );
        $evid_seen{ $tsv{'inSituExperimentId'}[$line].'-'.$tsv{'inSituEvidenceId'}[$line] }++;
    }

    # inSituSpot
    # Insert and get conditionId instead of anatEntityId/stageId to insert in inSituSpot table
    my $condKeyMap;
    ($condKeyMap, $conditions) = Utils::insert_get_condition($dbh, $conditions, $stage_equivalences,
                                                             $tsv{'organId'}[$line], $tsv{'stageId'}[$line],
                                                             $species, $tsv{'sex'}[$line], $tsv{'strain'}[$line],
                                                             $anatSexInfo, $speciesSexInfo,
                                                             $tsv{'inSituExperimentId'}[$line], $gene_mapping->{ $tsv{'geneId'}[$line] });
    #FIXME What to do with one to many XRefId / bgeeGeneId relationship

    $inSituSpotId++;
    my (undef, $spot) = split(/_/, $tsv);
    $spot .= "-$inSituSpotId";
    $insSpot->execute($spot, $tsv{'inSituEvidenceId'}[$line], $condKeyMap->{'conditionId'}, $gene_mapping->{ $tsv{'geneId'}[$line] }, $tsv{'detectionFlag'}[$line], $tsv{'inSituData'}[$line])
        or die "[$spot][$tsv{'inSituEvidenceId'}[$line]] ".$insSpot->errstr;
    print "[$spot]\t[$tsv{'inSituEvidenceId'}[$line]]\t[$condKeyMap->{'conditionId'}]\t[$gene_mapping->{ $tsv{'geneId'}[$line] }]\t[$tsv{'detectionFlag'}[$line]]\t[$tsv{'inSituData'}[$line]]\n"  if ( $debug );
}
$insEvid->finish;
$insSpot->finish;

exit 0;

