#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;


my $species = $ARGV[0]  or die "\n\t$0 <genus_species>    -> e.g. $0 homo_sapiens\n\n";


# Ensembl connection via Ensembl API/Registry
#NOTE local Ensembl:
#my $ensembl_connector = 'user=bgee__pass=XXX__host=annotbioinfo__port=3306';
# Remote up-to-date Ensembl
my $ensembl_connector = 'user=anonymous__pass=__host=ensembldb.ensembl.org__port=5306';
my $reg = Utils::connect_ensembl_registry($ensembl_connector, 0);


my $array_adaptor = $reg->get_adaptor($species, 'Funcgen', 'Array');
my @ac_ids = @{ $array_adaptor->fetch_all() };


print join("\t", '#Name', 'Vendor', 'Type', 'Format', 'Class', 'Description'), "\n";
map { print join("\t", $_->name, $_->vendor, $_->type, $_->format, $_->class, ($_->description || '' )), "\n" } @ac_ids;

#print Dumper @ac_ids;

exit 0;

