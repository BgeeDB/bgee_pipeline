#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Data::Dumper;
$| = 1;

# Julien Wollbrett, created December 2021

# USAGE: perl insert_uniprot_xrefs.pl -uniprot_xrefs=<...> -bgee=<...> <OPTIONAL: -debug=1>
# This script inserts the uniprot xrefs retrieved from uniprot using the R script
# update_genes_from_UniProt.R.
# This script was created for Bgee 15.0. Some uniprot XRefs was already inserted but some were
# missing. 
# The script checks xrefs already inserted because inserting with a query "...on duplicate update" is not
# possible. We first have to check if same uniprot ID with different letter case was already 
# inserted. If this verification is not done first, it is then not possible to know if the table 
# geneToTerm has to be updated or not.
#
# -debug=1: if provided, run in verbose mode (print the update/insert SQL queries, not executing them)
#############################################################


# Define arguments & their default value
my ($bgee_connector) = ('');
my ($uniprot_xrefs)    = ('');
my ($debug)          = (0);
my %opts = ('bgee=s'          => \$bgee_connector,   # Bgee connector string
            'uniprot_xrefs=s' => \$uniprot_xrefs,    # generated_files/genes/xrefs_from_uniprot.txt file
            'debug'           => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $uniprot_xrefs eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD) -uniprot_xrefs=\$(XREFS_FROM_UNIPROT_FILE)
\t-bgee             Bgee connector string
\t-uniprot_xrefs    xrefs_from_uniprot.txt file
\t-debug            printing the update/insert SQL queries, not executing them
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

# it is safer to map the datasource on the name rather than on the ID
my $SWISSPROT_DATASOURCE_NAME = 'UniProtKB/Swiss-Prot';
my $TREMBL_DATASOURCE_NAME    = 'UniProtKB/TrEMBL';

print "Retrieves potentially new uniprot xrefs...\n";
my %new_xrefs;
open(my $IN, '<', $uniprot_xrefs)  or die "Could not read file [$uniprot_xrefs]\n";
while ( defined (my $line = <$IN>) ){
    next  if ( $line =~ /^GeneID\t/ );

    chomp $line;
    my @tmp = split(/\t/, $line);
    my $speciesId   = $tmp[5];
    my $geneId      = $tmp[0];
    my $uniprotId   = $tmp[1];
    my $uniprotName = $tmp[2];
    my $reviewed    = $tmp[3];

    # record the xrefs
    $new_xrefs{$speciesId}->{$geneId}->{$uniprotId} -> {'reviewed'} = $reviewed;
    $new_xrefs{$speciesId}->{$geneId}->{$uniprotName} -> {'reviewed'} = $reviewed;
}

print "retrieves uniprot xrefs already inserted in the database...\n";

my $retrieveQuery = $bgee->prepare("SELECT t1.geneId, t1.speciesId, t2.XRefId, t2.XRefName, t3.dataSourceID FROM gene AS t1 INNER JOIN geneXRef AS t2 ON t1.bgeegeneId = t2.bgeegeneId INNER JOIN dataSource AS t3 on t2.dataSourceId = t3.dataSourceId where t3.dataSourceName IN (?, ?)");
my %existing_xrefs;
$retrieveQuery->execute($SWISSPROT_DATASOURCE_NAME, $TREMBL_DATASOURCE_NAME)  or die $retrieveQuery->errstr;
while(my @values = $retrieveQuery -> fetchrow_array()){
    my $geneId       = $values[0];
    my $speciesId    = $values[1];
    my $uniprotId    = $values[2];
    my $uniprotName  = $values[3];
    my $dataSourceId = $values[4];
    #in this hash uniprotIds are always lower case in order to easily map to new uniprot XRefs
    $existing_xrefs{$speciesId}->{$geneId}->{lc $uniprotId} -> {'uniprotName'} = $uniprotName;
    $existing_xrefs{$speciesId}->{$geneId}->{lc $uniprotId} -> {'dataSourceId'} = $dataSourceId;
    $existing_xrefs{$speciesId}->{$geneId}->{lc $uniprotId} -> {'uniprotId'} = $uniprotId;
}

print "retrieves geneId to bgeegeneId mappings from the database...\n";

my $retrieveGenes = $bgee->prepare("SELECT speciesId, geneId, bgeegeneId FROM gene");
my %bgee_genes;
$retrieveGenes->execute()  or die $retrieveGenes->errstr;
while(my @values = $retrieveGenes -> fetchrow_array()){
    my $geneId       = $values[1];
    my $bgeegeneId   = $values[2];
    my $speciesId    = $values[0];
    $bgee_genes{$speciesId}->{$geneId} = $bgeegeneId;
}

print "retrieves internal IDs of swissprot and trembl datasources...\n";

my $retrieveDataSourceIds = $bgee->prepare("SELECT dataSourceId, dataSourceName from dataSource");
my %type_of_datasource;
$retrieveDataSourceIds->execute()  or die $retrieveDataSourceIds->errstr;
while(my @values = $retrieveDataSourceIds -> fetchrow_array()){
    my $dataSourceId   = $values[0];
    my $dataSourceName = $values[1];
    if($dataSourceName eq $SWISSPROT_DATASOURCE_NAME) {
        $type_of_datasource{"reviewed"} = $dataSourceId;
    }
    if($dataSourceName eq $TREMBL_DATASOURCE_NAME) {
        $type_of_datasource{"unreviewed"} = $dataSourceId;
    }
}

print "Inserting data into the database...\n";

my $insert_geneXref     = $bgee->prepare("INSERT INTO geneXRef (bgeegeneId, XRefId, XRefName, dataSourceId) VALUES (?, ?, ?, ?)");
my $insert_geneToTerm   = $bgee->prepare("INSERT INTO geneToTerm (bgeegeneId, term) VALUES (?, ?)");
my $update_datasourceId = $bgee->prepare("UPDATE FROM geneXRef SET dataSourceId = ? where bgeegeneId = ? and XRefId = ?");

for my $speciesId ( keys %new_xrefs ){
    for my $geneId ( keys %{$new_xrefs{$speciesId}} ){
        my $bgeegeneId = $bgee_genes{$speciesId}{$geneId};
        for my $uniprotId ( keys %{$new_xrefs{$speciesId}{$geneId}} ){
            my $existing_uniprot_id = $existing_xrefs{$speciesId}{$geneId}{lc $uniprotId}{'uniprotId'};
            my $to_insert_datasource = $type_of_datasource{$new_xrefs{$speciesId}{$geneId}{$uniprotId}{'reviewed'}};
            if (defined( $existing_xrefs{$speciesId}{$geneId}{lc $uniprotId}{'dataSourceId'} ) ) {
                # xref already exists. Check if should update geneXRef.dataSourceId. No insertion in geneToTerm
                my $inserted_datasource = $existing_xrefs{$speciesId}{$geneId}{lc $uniprotId}{'dataSourceId'};
                if (defined($to_insert_datasource) && $inserted_datasource != $to_insert_datasource) {
                    if ( $debug ){
                        print "UPDATE FROM geneXRef SET dataSourceId = $to_insert_datasource WHERE bgeegeneId = ", 
                            "$bgeegeneId AND XRefId = $existing_uniprot_id\n";
                    }
                    else {
                        $update_datasourceId->execute($to_insert_datasource, $bgeegeneId, $uniprotId)  or die $insLength->errstr;
                    }
                }
                else {
                    if ( $debug ){
                        print "already exist $bgeegeneId, $existing_uniprot_id, $inserted_datasource. Nothing to insert\n";
                    }
                }
            }
            # xref does not exist. Insert in both geneXRef and geneToTerm tables
            else {
                my $to_insert_datasource = $type_of_datasource{$new_xrefs{$speciesId}{$geneId}{$uniprotId}{'reviewed'}};
                if (!defined($to_insert_datasource)) {
                    $to_insert_datasource = $type_of_datasource{"unreviewed"};
                }
                if ($debug) {
                    print "INSERT INTO geneXRef (bgeegeneId, XRefId, XRefName, dataSourceId) VALUES ($bgeegeneId, $uniprotId, \"\", $to_insert_datasource)\n";
                    print "INSERT INTO geneToTerm (bgeegeneId, term) VALUES ($bgeegeneId, $uniprotId)\n";
                } else {
                    $insert_geneXref <- execute($bgeegeneId, $uniprotId, "", $to_insert_datasource);
                    $insert_geneToTerm <- execute($bgeegeneId, $uniprotId);
                }
            }
        }
    }
}

$insert_geneXref -> finish();
$insert_geneToTerm -> finish();
$update_datasourceId -> finish();

for my $speciesId ( keys %existing_xrefs ){
    for my $geneId ( keys %{$existing_xrefs{$speciesId}} ){
        for my $uniprotId ( keys %{$existing_xrefs{$speciesId}{$geneId}} ){
            if(!defined($new_xrefs{$speciesId}{$geneId}{lc $uniprotId}) && !defined($new_xrefs{$speciesId}{$geneId}{uc $uniprotId})) {
                print "exists in the database but not found in new XRefs $speciesId, $geneId, $uniprotId.\n";
            }
        }
    }
}
print "Done\n";
exit 0;

