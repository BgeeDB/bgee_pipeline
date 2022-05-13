#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use List::MoreUtils qw{uniq any};
use LWP::Simple;
#use JSON::XS; # See http://blogs.perl.org/users/e_choroba/2018/03/numbers-and-strings-in-json.html
use Cpanel::JSON::XS;
use Data::Dumper;

use FindBin;
use lib "$FindBin::Bin/."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my ($gene_id, $bgee_connector) = ('', '');
my ($verbose)                  = (0);
my %opts = ('gene=s'    => \$gene_id,            # gene id to search in Ensembl to insert/update in the db
            'bgee=s'    => \$bgee_connector,     # Bgee connector string
            'verbose'   => \$verbose,            # verbose mode
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t     OR
\t     $0  -bgee=\$(BGEECMD)  -gene=wbgene00001127/ENSGALG00000039923
\t-gene     specific gene id to search in Ensembl to insert/update in the db
\t-bgee     Bgee    connector string
\t-verbose  verbose mode
\n";
    exit 1;
}


# EnsEMBL REST URLs
#doc: http://rest.ensembl.org/
#e.g. http://rest.ensembl.org/xrefs/id/ENSTGEP00000031995?content-type=application/json
my $ensembl_lookup = 'https://rest.ensembl.org/lookup/id/';
my $ensembl_xrefs  = 'https://rest.ensembl.org/xrefs/id/';
my $ensembl_type   = '?content-type=application/json';
# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

my $up_xref  = $dbh->prepare('INSERT INTO geneXRef               (bgeeGeneId, XRefId, XRefName, dataSourceId) VALUES (?, ?, ?, ?)');
my $up_syn   = $dbh->prepare('INSERT INTO geneNameSynonym        (bgeeGeneId, geneNameSynonym)                VALUES (?, ?)');
my $up_term  = $dbh->prepare('INSERT INTO geneToTerm             (bgeeGeneId, term)                           VALUES (?, ?)');
my $up_go    = $dbh->prepare('INSERT INTO geneToGeneOntologyTerm (bgeeGeneId, goId, goEvidenceCode)           VALUES (?, ?, ?)');
my $up_gene1 = $dbh->prepare('UPDATE gene SET geneName=?         WHERE bgeeGeneId=?');
my $up_gene2 = $dbh->prepare('UPDATE gene SET geneDescription=?  WHERE bgeeGeneId=?');

## DataSources
# Get used dataSources
my $sourceDB = $dbh->prepare('SELECT dataSourceId, dataSourceName FROM dataSource');
$sourceDB->execute()  or die $sourceDB->errstr;
my %InsertedDataSources = map { $_->[1] = lc $_->[1]; $_->[1] => $_->[0] } @{$sourceDB->fetchall_arrayref}; # List already inserted dataSources and return it in a hash dataSourceName (lowercase) => dataSourceId
$sourceDB->finish;
# Add extra dataSource aliases
# MUST be in lowercase to ease comparison
#TODO Add other species specific dataSource variant names:
#TODO Add CGNC (Chichen Gene Nomenclature Consortium)?
$InsertedDataSources{'flybasename_gene'} = $InsertedDataSources{'flybase'};
$InsertedDataSources{'flybasecgid_gene'} = $InsertedDataSources{'flybase'};
$InsertedDataSources{'flybase_symbol'}   = $InsertedDataSources{'flybase'};
$InsertedDataSources{'wormbase_gene'}    = $InsertedDataSources{'wormbase'};
$InsertedDataSources{'xenopus_jamboree'} = $InsertedDataSources{'xenbase'};
$InsertedDataSources{'zfin_id'}          = $InsertedDataSources{'zfin'};


#Test gene table for empty fields
#NOTE BLxxxx ids from Amphioxus are not properly inserted in Ensembl BUT similar ids may exist there!
my $sth_genes = $dbh->prepare('SELECT bgeeGeneId, geneId, geneName, geneDescription FROM gene WHERE (geneDescription="" OR geneName="") AND geneId NOT LIKE "BL%"');
if ( $gene_id ne '' ){
    $sth_genes = $dbh->prepare("SELECT bgeeGeneId, geneId, geneName, geneDescription FROM gene  WHERE geneId='$gene_id' AND geneId NOT LIKE 'BL%'");
}
$sth_genes->execute()  or die $sth_genes->errstr;
GENE:
while ( my ($st_own_id, $st_gene_id, $st_gene_name, $st_description) = ($sth_genes->fetchrow_array) ){
    next  if ( $st_gene_id !~ /^[A-Za-z]/ ); #Should be NCBI RefSeq gene id is starts by digits
    sleep 1;
    ## Get info from Ensembl/Ensembl Metazoa

	# Name and description
	my ($name, $desc) = ('', '');
    my $contentg = get("$ensembl_lookup$st_gene_id$ensembl_type");
    if ( defined $contentg ){
        my $json = decode_json( $contentg );
        $name = $json->{'display_name'} // '';
        $desc = $json->{'description'}  // '';
    }

    if ( $name =~ / \[provisional:(.+?)\]$/ ){ #e.g. XB5961369 [provisional:plpp3]
        $name = $1;
    }

#	# Xrefs
#    my $xrefs;
#    my $contentx = get("$ensembl_xrefs$st_gene_id$ensembl_type");
#    if ( defined $contentx ){
#        my $json = decode_json( $contentx );
#        XREF:
#        for my $xref ( @{ $json }){
#            next XREF  if ( $xref->{'dbname'} eq 'ArrayExpress' ); #May refer to another species ensembl id
#            next XREF  if ( $xref->{'dbname'} =~ /^Ens_Ga_/ );     #May refer to another species ensembl id
#            next XREF  if ( $xref->{'dbname'} eq 'KEGG_Enzyme' );  #Weird EC number syntax
#            next XREF  if ( $xref->{'dbname'} eq 'MEROPS' );       #Useful?
#
#            my $id = $xref->{'dbname'} eq 'GO'         ? $xref->{'display_id'}
#                   : $xref->{'dbname'} =~ /^ZFIN_ID_/  ? 'ZFIN_ID:'.$xref->{'display_id'}
#                   : $xref->{'dbname'} =~ /^Uniprot\// ? 'Uniprot:'.$xref->{'display_id'}
#                   : $xref->{'dbname'} =~ /^RefSeq_/   ? 'RefSeq:'.$xref->{'display_id'}
#                   : $xref->{'dbname'} =~ /^HGNC_/     ? 'HGNC:'.$xref->{'display_id'}
#                   : $xref->{'dbname'} =~ /^Vega_/     ? 'Vega:'.$xref->{'display_id'}
#                   :                                     $xref->{'dbname'}.':'.$xref->{'display_id'};
#            $xrefs->{ $id } = {};
#            if ( scalar @{ $xref->{'synonyms'} } > 0 ){
#                $xrefs->{ $id }->{'synonyms'} = join(';', @{ $xref->{'synonyms'} });
#            }
#        }
#    }


    # Update gene table
    if ( $st_gene_name     eq '' && $name     ne '' ){
        $up_gene1->execute($name,     $st_own_id)  or die $up_gene1->errstr;
        print "\tInsert name for $st_gene_id\n";
    }
    if ( $st_description   eq '' && $desc     ne '' && $desc ne 'Uncharacterized protein' ){
        $up_gene2->execute($desc,     $st_own_id)  or die $up_gene2->errstr;
        print "\tInsert desc for $st_gene_id\n";
    }

#    # Update synonym table
#    if ( scalar @{ $annotations->{'synonyms'} } > 0 ){
#        for my $syn ( @{ $annotations->{'synonyms'} } ){
#            $up_syn->execute($st_own_id, $syn);#  or die $up_syn->errstr;
#        }
#    }
#    # Update GO table
#    if ( scalar @{ $annotations->{'go'} } > 0 ){
#        for my $go_full ( @{ $annotations->{'go'} } ){
#            my ($go, $evidence) = split(/___/, $go_full);
#            $up_go->execute($st_own_id, $go, $evidence || '');#  or die $up_go->errstr;
#        }
#    }
#    # Update xref table
#    if ( scalar @{ $annotations->{'xrefs'} } > 0 ){
#        for my $xf ( @{ $annotations->{'xrefs'} } ){
#            next  if ( $xf =~ /;/ );
#            my ($db_source, $pid) = split(/:/, $xf, 2);
#            $up_xref->execute($st_own_id, $pid, '', $InsertedDataSources{lc $db_source});#  or die $up_xref->errstr;
#        }
#    }
#
#    # Update term tables
#    my @all;
#    push @all, $annotations->{'prot_name'}                         if ( $annotations->{'prot_name'} );
#    push @all, @{ $annotations->{'synonyms'} }                     if ( scalar @{ $annotations->{'synonyms'} } > 0 );
#    push @all, map { s/___.*$//; $_ } @{ $annotations->{'go'} }    if ( scalar @{ $annotations->{'go'} } > 0 );
#    push @all, map { s/^.+?://; $_ }  @{ $annotations->{'xrefs'} } if ( scalar @{ $annotations->{'xrefs'} } > 0 );
#    #without version digit
#    push @all, map { s/^.+?://; $_ } map  { s/\.\d+$//; $_ } grep { /\.\d+$/ } @{ $annotations->{'xrefs'} } if ( scalar @{ $annotations->{'xrefs'} } > 0 );
#    for my $term ( uniq sort @all ){
#        $up_term->execute($st_own_id, $term);#  or die $up_term->errstr;
#    }
}

$sth_genes->finish;
$up_gene1->finish;
$up_gene2->finish;
$up_xref->finish;
$up_syn->finish;
$up_term->finish;
$up_go->finish;
# Close db connections
$dbh->disconnect;

exit 0;

