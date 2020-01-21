#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Text::CSV;

my ($bgee_connector, $oncoMX_ids_file, $oncoMX_datasource_name) = ('', '', '');
my $debug = 0;
my %opts = ('bgee=s'                     => \$bgee_connector,     		 # Bgee connector string
			'oncoMX_ids_file=s'   		 => \$oncoMX_ids_file,           # path to file containing all oncoMX ids
			'oncoMX_datasource_name=s'   => \$oncoMX_datasource_name,    # name of the OncoMX datasource in the Bgee database
            'debug'                      => \$debug,                     # debug mode, do not insert/update in database
           );
           
# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( $bgee_connector eq '' || $oncoMX_ids_file eq '' || $oncoMX_datasource_name eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. perl $0 -bgee=\$(BGEECMD) -oncoMX_ids_file=path/to/file.tsv oncoMX_datasource_name=OncoMX
\t-bgee      				Bgee connector string
\t-oncoMX_ids_file      	path to file containing all oncoMX ids
\t-oncoMX_datasource_name   name of the OncoMX datasource in the Bgee database
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

# retrieve internal bgeeGeneId, ensembl gene Id and gene name for all human genes in Bgee
my $retrieve_human_genes = "SELECT bgeeGeneId, geneId, geneName FROM gene where speciesId = 9606";
# retrieve oncoMX datasource ID
my $retrieve_oncoMX_id = "SELECT dataSourceId from dataSource where datasourceName = ? ";
# retrieve oncoMX datasource ID
my $insert_oncoMX_xrefs = "INSERT INTO geneXRef(bgeeGeneId, XRefId, XRefName, dataSourceId) VALUES(?, ?, ?, ?)";

##### MAP OncoMX IDs and bgeeGeneIds using uniprot XRefs ######


# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);
my %bgeeId_to_xrefId;
my $sth = $dbh->prepare($retrieve_human_genes);
$sth->execute()  or die $sth->errstr;
while(my @values = $sth -> fetchrow_array()){
	my $bgeeGeneId = $values[0];
	my $ensemblId = $values[1];
	my $geneName = $values[2];
	if (exists($oncoMX_ids{$geneName})) {
		$bgeeId_to_xrefId{$bgeeGeneId}{'xref_id'} = $geneName;
		$bgeeId_to_xrefId{$bgeeGeneId}{'xref_name'} = $geneName;
	} elsif (exists($oncoMX_ids{$ensemblId})) {
		$bgeeId_to_xrefId{$bgeeGeneId}{'xref_id'} = $ensemblId;
		$bgeeId_to_xrefId{$bgeeGeneId}{'xref_name'} = $geneName;
	}
}
$sth->finish();

##### Insert oncoMX XRefs in Bgee ######
my $oncoMX_id = $dbh->selectrow_array($retrieve_oncoMX_id, undef, $oncoMX_datasource_name);

$sth = $dbh->prepare($insert_oncoMX_xrefs);
foreach my $bgeeGeneId (keys %bgeeId_to_xrefId){
	if(!$debug) {
		$sth->execute($bgeeGeneId, $bgeeId_to_xrefId{$bgeeGeneId}{'xref_id'}, $bgeeId_to_xrefId{$bgeeGeneId}{'xref_name'}, $oncoMX_id)  or die $sth->errstr;
	} else {
		print "$bgeeGeneId, $bgeeId_to_xrefId{$bgeeGeneId}{'xref_id'}, $bgeeId_to_xrefId{$bgeeGeneId}{'xref_name'}, $oncoMX_id\n";
	}
}
if(!$debug) {
	my $size = keys %bgeeId_to_xrefId;
	print "$size OncoMX XRefs have been inserted successfully in database $dbh->{Name}\n";
}
$sth->finish();
$dbh->disconnect;
