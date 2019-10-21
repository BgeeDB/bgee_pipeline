#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Data::Dumper;
use Text::CSV;

my ($bgee_connector, $oncoMX_ids_file, $oncoMX_datasource_name, $uniProt_datasource_name) = ('', '', '','');
my $debug = 0;
my %opts = ('bgee=s'                     => \$bgee_connector,     		 # Bgee connector string
			'oncoMX_ids_file=s'   		 => \$oncoMX_ids_file,           # path to file containing all oncoMX ids
			'oncoMX_datasource_name=s'   => \$oncoMX_datasource_name,    # name of the OncoMX datasource in the Bgee database
			'uniProt_datasource_name=s'   => \$uniProt_datasource_name,  # name of the UniProt datasource in the Bgee database
            'debug'                      => \$debug,                     # debug mode, do not insert/update in database
           );
           
# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( $bgee_connector eq '' || $oncoMX_ids_file eq '' || $oncoMX_datasource_name eq '' || $uniProt_datasource_name eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. perl $0 -bgee=\$(BGEECMD) -oncoMX_ids_file=path/to/file.tsv oncoMX_datasource_name=OncoMX uniProt_datasource_name=Uniprot/SWISSPROT
\t-bgee      				Bgee connector string
\t-oncoMX_ids_file      	path to file containing all oncoMX ids
\t-oncoMX_datasource_name   name of the OncoMX datasource in the Bgee database
\t-uniProt_datasource_name  name of the UniProt/SWISSPROT datasource in the Bgee database
\t-debug                	Debug mode, do not insert/update in database
\n";
    exit 1;
}

my $unwanted_xref_pattern = "%\\_%";

##### LOAD OncoMX IDs ######

my %oncoMX_ids;
my $csv = Text::CSV->new({ sep_char => ',' });
open(my $lines, $oncoMX_ids_file) || die "failed to read input file: $!";
while (my $line = <$lines>) {
	chomp $line;
	my ($key, $value) = ('','');
	if ($csv->parse($line)) {
		my @fields = $csv->fields();
		$key = $fields[0];
		$value = $fields[1];
	}
	$oncoMX_ids{$key} = $value;
}

#print Dumper(\%oncoMX_ids);

##### MYSQL QUERIES ######

# retrieve bgeeGeneId and XRefId of all genes from Bgee having a uniprot id XRef
my $retrieve_uniprot_xrefs = "SELECT t1.bgeeGeneId, t2.XRefId FROM gene AS t1 INNER JOIN geneXRef AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId INNER JOIN dataSource AS t3 ON t2.dataSourceId = t3.dataSourceId where t3.datasourceName = ? and t2.XRefId NOT LIKE ? ";
# retrieve oncoMX datasource ID
my $retrieve_oncoMX_id = "SELECT dataSourceId from dataSource where datasourceName = ? ";
# retrieve oncoMX datasource ID
my $insert_oncoMX_xrefs = "INSERT INTO geneXRef(bgeeGeneId, XRefId, XRefName, dataSourceId) VALUES(?, ?, ?, ?)";

##### MAP OncoMX IDs and bgeeGeneIds using uniprot XRefs ######


# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

my %gene_to_xrefs;
my $sth = $dbh->prepare($retrieve_uniprot_xrefs);
$sth->execute($uniProt_datasource_name, $unwanted_xref_pattern)  or die $sth->errstr;
while(my @values = $sth -> fetchrow_array()){
	my $bgeeGeneId = $values[0];
	my $uniprotId = $values[1];
	if (exists($oncoMX_ids{$uniprotId})) {
		$gene_to_xrefs{$bgeeGeneId}{'xref_id'} = $uniprotId;
		$gene_to_xrefs{$bgeeGeneId}{'xref_name'} = $oncoMX_ids{$uniprotId};
	}
}
$sth->finish();

##### Insert oncoMX XRefs in Bgee ######
my $oncoMX_id = $dbh->selectrow_array($retrieve_oncoMX_id, undef, $oncoMX_datasource_name);

$sth = $dbh->prepare($insert_oncoMX_xrefs);
foreach my $bgeeGeneId (keys %gene_to_xrefs){
	if(!$debug) {
		$sth->execute($bgeeGeneId, $gene_to_xrefs{$bgeeGeneId}{'xref_id'}, $gene_to_xrefs{$bgeeGeneId}{'xref_name'}, $oncoMX_id)  or die $sth->errstr;
	} else {
		print "$bgeeGeneId, $gene_to_xrefs{$bgeeGeneId}{'xref_id'}, $gene_to_xrefs{$bgeeGeneId}{'xref_name'}, $oncoMX_id\n";
	}
}
$sth->finish();
$dbh->disconnect;
