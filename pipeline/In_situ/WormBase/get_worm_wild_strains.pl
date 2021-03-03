#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

$, = "\t";
binmode(STDOUT, 'utf8');
no warnings ('uninitialized');
use Webservice::InterMine 'http://intermine.wormbase.org/tools/wormmine';

my $query = new_query(class => 'Strain');
$query->add_view(qw/
    primaryIdentifier
    name
    species
    remark
    genotype
    otherName
/);

$query->add_constraint(
    path   => 'Strain.species',
    op     => 'ONE OF',
    values => [
        'Caenorhabditis elegans',
        ,
    ],
    code  => 'A',
);
$query->add_constraint(
    path  => 'Strain.genotype',
    op    => 'CONTAINS',
    value => 'wild',
    code  => 'B',
);

my $it = $query->iterator();
while ( my $row = <$it> ){
    print join("\t", $row->{'primaryIdentifier'}, # strain code
                     $row->{'name'},              # human readable strain name
                     $row->{'species'},
                     $row->{'remark'},
                     $row->{'genotype'},
                     $row->{'otherName'},
              ), "\n";
}

exit 0;

