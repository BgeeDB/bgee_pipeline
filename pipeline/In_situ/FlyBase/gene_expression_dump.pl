#!/usr/bin/env perl
use strict;
use warnings;
use diagnostics;

use DBI;
use ChadoConnect;
require TAPget;

=pod

Makes dump of all gene expression data as pub, feature, anatomy, stage-range

usage ./gene_expression_dump.pl <DB name> > output_file.tsv

where DB name is a string identifying a database instance to connect to.  See doc for ChadoConnect.pm for details of available connections.

=cut

my $debug = 0;
unless ( $ARGV[0] ){
    die "Please specify a DB. See ChadoConnect doc for details of available DBs\n";
}


# make chado connection
my $dbh = ChadoConnect::make_con($ARGV[0]);

# First loop, modified from annotation_ref_gen_finder.pl, makes two hashes: $hash{fe_id}=fbex; $hash{fe_id}=out and array @fbex.

warn "*** Retrieving FBgn to FBex associations.***\n"  if ( $debug );
# "I split the query up into two SQL queries, depending on whether the expression statement specifies sex or not. I took that approach because when I request "sex", my query excluded expression statements where no sex is specified. So, I wrote a second query to report expression statements where no sex is specified. There might be a way to get both types of statements in one query, but that's beyond me. In any case, the expression statements returned by the two queries are mutually exclusive."
#
#NOTE gene.is_obsolete = 'f'  forces to use only current genes
#NOTE strain info not stored there (2016/10/17)!
#NOTE really in situ only now!
#
# "This SQL query returns ONLY expression statements where the sex is NOT specified. Requesting "sex" in the first query excluded expression statements where no sex was specified."
my $sth_sex_unspecified = $dbh->prepare("SELECT DISTINCT gene.name AS gene_name,
gene.uniquename AS fbgn,
gene.organism_id AS organism_id,
fe.feature_expression_id AS feid,
e.uniquename AS fbex,
assay.name AS assay,
molecule.name AS molecule_assayed,
stage.name AS stage,
'not specified' AS sex,
anatomy.name AS anatomy,
pub.uniquename AS fbrf,
pub.miniref AS miniref
FROM feature_expression fe
        JOIN expression e ON (e.expression_id = fe.expression_id)
        JOIN pub ON (fe.pub_id = pub.pub_id)
        JOIN feature_relationship fr ON (fr.subject_id = fe.feature_id)
        JOIN cvterm rel ON (fr.type_id = rel.cvterm_id)
        JOIN feature gene ON (fr.object_id = gene.feature_id)
        JOIN expression_cvterm ec1 ON (ec1.expression_id = e.expression_id)
        JOIN cvterm cvta ON (cvta.cvterm_id = ec1.cvterm_type_id)
        JOIN cv cves ON (cves.cv_id = cvta.cv_id)
        JOIN cvterm assay ON (assay.cvterm_id = ec1.cvterm_id)
        JOIN cv cva ON (cva.cv_id = assay.cv_id)
        JOIN expression_cvterm ec2 ON (ec2.expression_id = e.expression_id)
        JOIN cvterm cvtdv ON (cvtdv.cvterm_id = ec2.cvterm_type_id)
        JOIN cvterm stage ON (stage.cvterm_id = ec2.cvterm_id)
        JOIN cv cvdv ON (cvdv.cv_id = stage.cv_id)
        JOIN expression_cvterm ec4 ON (ec4.expression_id = e.expression_id)
        JOIN cvterm cvtbt ON (cvtbt.cvterm_id = ec4.cvterm_type_id)
        JOIN cvterm anatomy ON (anatomy.cvterm_id = ec4.cvterm_id)
        JOIN cv cvbt ON (cvbt.cv_id = anatomy.cv_id)
        JOIN feature gene_product ON (gene_product.feature_id = fe.feature_id)
        JOIN cvterm molecule ON (molecule.cvterm_id = gene_product.type_id)
        WHERE rel.name = 'associated_with'
        AND not e.expression_id in
                (SELECT DISTINCT e.expression_id
                FROM expression e, expression_cvterm ec, cvterm cvt
                WHERE e.expression_id = ec.expression_id and ec.cvterm_id = cvt.cvterm_id and cvt.name in ('male','virgin male','mated male','female','virgin female','mated female'))
        AND cves.name = 'expression slots'
        AND cves.cv_id = cvtdv.cv_id
        AND cves.cv_id = cvtbt.cv_id
        AND cvta.name = 'assay'
        AND cva.name = 'experimental assays'
        AND cvtdv.name = 'stage'
        AND cvdv.name = 'FlyBase development CV'
        AND cvtbt.name = 'anatomy'
        AND cvbt.name = 'FlyBase anatomy CV'
        AND gene.uniquename LIKE 'FBgn%'
        AND gene.is_obsolete = 'f'
        AND assay.name = 'in situ' -- LIMIT 15"); # uncomment limit for debugging

# "This SQL query returns ONLY expression statements where the sex IS specified. Note that in some atypical cases, we specify if the male/female is mated/virgin."
my $sth_sex_specified = $dbh->prepare("SELECT DISTINCT gene.name AS gene_name,
gene.uniquename AS fbgn,
gene.organism_id AS organism_id,
fe.feature_expression_id AS feid,
e.uniquename AS fbex,
assay.name AS assay,
molecule.name AS molecule_assayed,
stage.name AS stage,
sex.name AS sex,
anatomy.name AS anatomy,
pub.uniquename AS fbrf,
pub.miniref AS miniref
FROM feature_expression fe
        JOIN expression e ON (e.expression_id = fe.expression_id)
        JOIN pub ON (fe.pub_id = pub.pub_id)
        JOIN feature_relationship fr ON (fr.subject_id = fe.feature_id)
        JOIN cvterm rel ON (fr.type_id = rel.cvterm_id)
        JOIN feature gene ON (fr.object_id = gene.feature_id)
        JOIN expression_cvterm ec1 ON (ec1.expression_id = e.expression_id)
        JOIN cvterm cvta ON (cvta.cvterm_id = ec1.cvterm_type_id)
        JOIN cv cves ON (cves.cv_id = cvta.cv_id)
        JOIN cvterm assay ON (assay.cvterm_id = ec1.cvterm_id)
        JOIN cv cva ON (cva.cv_id = assay.cv_id)
        JOIN expression_cvterm ec2 ON (ec2.expression_id = e.expression_id)
        JOIN cvterm cvtdv ON (cvtdv.cvterm_id = ec2.cvterm_type_id)
        JOIN cvterm stage ON (stage.cvterm_id = ec2.cvterm_id)
        JOIN cv cvdv ON (cvdv.cv_id = stage.cv_id)
        JOIN expression_cvterm ec3 ON (ec3.expression_id = e.expression_id)
        JOIN cvterm cvtsx ON (cvtsx.cvterm_id = ec3.cvterm_type_id)
        JOIN cvterm sex ON (sex.cvterm_id = ec3.cvterm_id)
        JOIN cv cvsx ON (cvsx.cv_id = sex.cv_id)
        JOIN expression_cvterm ec4 ON (ec4.expression_id = e.expression_id)
        JOIN cvterm cvtbt ON (cvtbt.cvterm_id = ec4.cvterm_type_id)
        JOIN cvterm anatomy ON (anatomy.cvterm_id = ec4.cvterm_id)
        JOIN cv cvbt ON (cvbt.cv_id = anatomy.cv_id)
        JOIN feature gene_product ON (gene_product.feature_id = fe.feature_id)
        JOIN cvterm molecule ON (molecule.cvterm_id = gene_product.type_id)
        WHERE rel.name = 'associated_with'
        AND cves.name = 'expression slots'
        AND cves.cv_id = cvtdv.cv_id
        AND cves.cv_id = cvtsx.cv_id
        AND cves.cv_id = cvtbt.cv_id
        AND cvta.name = 'assay'
        AND cva.name = 'experimental assays'
        AND cvtdv.name = 'stage'
        AND cvdv.name = 'FlyBase development CV'
        AND cvtsx.name = 'stage'
        AND cvsx.name = 'FlyBase miscellaneous CV'
        AND sex.name in ('male','virgin male','mated male','female','virgin female','mated female')
        AND cvtbt.name = 'anatomy'
        AND cvbt.name = 'FlyBase anatomy CV'
        AND gene.uniquename LIKE 'FBgn%'
        AND gene.is_obsolete = 'f'
        AND assay.name = 'in situ' -- LIMIT 15"); # uncomment limit for debugging

# gene_name |     fbgn    | organism_id |  feid  |    fbex     |  assay  | molecule_assayed |       stage        |  sex   | anatomy |     fbrf    |                     miniref
# Nop60B    | FBgn0259937 |           1 |  75075 | FBex0008727 | in situ | mRNA             | oogenesis stage S8 | female | oocyte  | FBrf0105301 | Phillips et al., 1998, Mol. Gen. Genet. 260(1): 20--29

my $feid_out;
my $feid_fbex;
my $fbexr;
$sth_sex_unspecified->execute()  or die "\tUnable to execute query sex_unspecified\n";
($feid_out, $feid_fbex, $fbexr) = parse_query($sth_sex_unspecified, $feid_out, $feid_fbex, $fbexr);

$sth_sex_specified->execute()  or die "\tUnable to execute query sex_specified\n";
($feid_out, $feid_fbex, $fbexr) = parse_query($sth_sex_specified, $feid_out, $feid_fbex, $fbexr);
# $sth_sex_specified->finish;   ### Probably not necessary as should set connection to inactive once all data is collected by while loop.  Commented out for now, as may be more useful to get warning of incomplete queries on disconnection.
# $sth_sex_unspecified->finish; ### Probably not necessary as should set connection to inactive once all data is collected by while loop.  Commented out for now, as may be more useful to get warning of incomplete queries on disconnection.


# Make hash to send to TAPget.pm
my @fbex;
while ( my ($key, $value) = each %$fbexr ){
    push @fbex, $key;
}
# Take hash of fbex return datamodel
my $number_exp_statements = scalar(@fbex);
warn "*** Retrieving details of [$number_exp_statements] expression assertions (TAP statements). This may take a long time. Unparsable statements will be mentioned in STDERR ***\n"  if ( $debug );

my $TAP_dm = TAPget::get_TAP_data_model($dbh, \@fbex);  # But note - this works on names!

my $FBbtdv_name_id = roll_anatomy_stage_lookup($dbh);


# Output loop
warn "*** Printing gene expression details to STDOUT ***\n"  if ( $debug );
print join("\t", '#feid', qw(gene_name gene_id pub_ref pub_desc organism_id assay molecule_assayed stage sex anatomy expression_uniquename anatomy_desc anatomy_id start_stage_name start_stage_id end_stage_name end_stage_id)), "\n";
FEID:
for my $feid ( sort keys %$feid_fbex ){ # Force always the same order for diff
    my $fbex = $feid_fbex->{$feid};
    FBEX:
    for ( @{$fbex} ){
        print $feid_out->{$feid}, "\t$_";
        if ( $TAP_dm->{$_}->{'anatomy'} && exists $FBbtdv_name_id->{$TAP_dm->{$_}->{'anatomy'}} ){
            print "\t".$TAP_dm->{$_}->{'anatomy'}."\t".$FBbtdv_name_id->{$TAP_dm->{$_}->{'anatomy'}};
        }
        else {
            print "\torganism\t".$FBbtdv_name_id->{'organism'};  # If no anatomy specified - use 'organism'
        }
        if ( exists $TAP_dm->{$_}->{'amp_stages'} ){  # If only one stage is recorded - put this in both start and end slots.
            if ( @{$TAP_dm->{$_}->{'amp_stages'}} == 1 ){
                my $s = pop(@{$TAP_dm->{$_}->{'amp_stages'}});
                $TAP_dm->{$_}->{'start_stage'} = $s;
                $TAP_dm->{$_}->{'end_stage'}   = $s;
            }
        }
        if ( $TAP_dm->{$_}->{'start_stage'} && exists $FBbtdv_name_id->{$TAP_dm->{$_}->{'start_stage'}} ){
            print "\t".$TAP_dm->{$_}->{'start_stage'}."\t".$FBbtdv_name_id->{$TAP_dm->{$_}->{'start_stage'}};
        }
        else {
            print "\t\t";
        }
        if ( $TAP_dm->{$_}->{'end_stage'} && exists $FBbtdv_name_id->{$TAP_dm->{$_}->{'end_stage'}} ){
            print "\t".$TAP_dm->{$_}->{'end_stage'}."\t".$FBbtdv_name_id->{$TAP_dm->{$_}->{'end_stage'}};
        }
        else {
            print "\t\t";
        }
        print "\n";
    }
}

my $rc = $dbh->disconnect  or warn $dbh->errstr; # May not be necessary, but could provide useful warning of incomplete SQL queries.

exit 0;


sub roll_anatomy_stage_lookup {
    #NOTE - RELIES ON UNIQUE NAME ASSUMPTION FOR anatomy+stage!
    my %FBbtdv_name_id;
    my $dbh =$_[0];
    my $sth = $dbh->prepare('SELECT c.name, dbx.accession, db.name AS idp
 FROM cvterm c
 JOIN dbxref dbx ON (c.dbxref_id = dbx.dbxref_id)
 JOIN db ON (dbx.db_id = db.db_id)
 WHERE c.is_obsolete = \'0\' -- NOT IMPLEMENTED AS BOOLEAN IN CHADO!
 AND db.name IN (\'FBbt\', \'FBdv\')');
    $sth->execute()  or die "WARNING: ERR2: Unable to execute query\n";
    # Make perl data structure
    while ( my $hash_ref = $sth->fetchrow_hashref ){
        $FBbtdv_name_id{$hash_ref->{'name'}} = "$hash_ref->{'idp'}:$hash_ref->{'accession'}";
    }
    return \%FBbtdv_name_id;
}

sub parse_query {
    my ($sth, $feid_out, $feid_fbex, $fbexr) = @_;

    while ( my $hash_ref = $sth->fetchrow_hashref ){
        my $feid = $hash_ref->{'feid'};
        my $fbex = $hash_ref->{'fbex'}; # expression uniquename
        $feid_out->{$feid} = join("\t", $hash_ref->{'feid'},
                                        $hash_ref->{'gene_name'},
                                        $hash_ref->{'fbgn'},
                                        $hash_ref->{'fbrf'},
                                        $hash_ref->{'miniref'},
                                        $hash_ref->{'organism_id'},
                                        $hash_ref->{'assay'},
                                        $hash_ref->{'molecule_assayed'},
                                        $hash_ref->{'stage'},
                                        $hash_ref->{'sex'},
                                        $hash_ref->{'anatomy'},
                             );
        push (@{$feid_fbex->{$feid}}, $fbex); # Array of FBex for each
        $fbexr->{$fbex} = 1; # Unique list of expression statements, to avoid looking up anything multiple times in chado.
    }

    return ($feid_out, $feid_fbex, $fbexr);
}

