package Utils_insert_genes;

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Array::Utils qw(:all);
use Bio::EnsEMBL::Registry; # Require Ensembl API
use DBI;
use File::Basename;
use File::Slurp;
use IO::Socket;
use List::MoreUtils qw(uniq);
use Spreadsheet::Read qw{ReadData row};

###############################################################################
# Patrick Tran Van, created April 2019
#
# Common functions between `insert_genes.pl` and `insert_genes_nonEnsembl.pl`.
###############################################################################

# This sub checks if the species id can be mapped to another species id.

sub map_species {
	
	my ($dbh, $species)= @_;

	# Need to map to another genomeSpeciesId?
	my ($speciesBgee, $newSpecies, $scientific_name, $ensSource) = split('__', $species, -1);
	
	if ( $speciesBgee == $newSpecies || $newSpecies == 0 ){
	    # No mapping to another species
	    $species = $speciesBgee;
	}
	else {
	    $species = $newSpecies;
	}

	return($species, $speciesBgee, $newSpecies, $scientific_name, $ensSource);

}

# This sub get the alternatives GO id

sub get_alt_go {
	
	my ($dbh)= @_;

	my $altgoDB = $dbh->prepare('SELECT goAltId, goId FROM geneOntologyTermAltId');
	$altgoDB->execute()  or die $altgoDB->errstr;
	my %altid_go = map { lc $_->[0] => lc $_->[1] }
	               @{$altgoDB->fetchall_arrayref};
	
	$altgoDB->finish;

	return(%altid_go);

}

# This sub get obsolete GO to not insert them later on

sub get_obs_go {
	
	my ($dbh, $obsGO)= @_;
	
	print "obsGO : [$obsGO]\n";
	die "Missing [$obsGO] file\n"  if ( !-e "$obsGO" || -z "$obsGO" );
	my %obs_go = map { lc $_ => 1 }
             read_file("$obsGO", chomp => 1);

	return(%obs_go);

}

# This sub return id of data sources

sub get_source_id {
	
	my ($dbh)= @_;

	my $sourceDB = $dbh->prepare('SELECT dataSourceId, dataSourceName FROM dataSource');
	$sourceDB->execute()  or die $sourceDB->errstr;
	my %InsertedDataSources = map { $_->[1] = lc $_->[1]; $_->[1] => $_->[0] } @{$sourceDB->fetchall_arrayref}; # List already inserted dataSources and return it in a hash dataSourceName (lowercase) => dataSourceId
	$sourceDB->finish;

	return(%InsertedDataSources);

}

1;

