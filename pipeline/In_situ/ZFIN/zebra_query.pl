#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

binmode(STDOUT, 'utf8');
use Webservice::InterMine 1.0405 'http://zebrafishmine.org';

my $query = new_query(class => 'Gene');
$query->add_view(qw/
    primaryIdentifier
    expressions.expressionFound
    symbol
    expressions.anatomy.name
    expressions.anatomy.identifier
    expressions.startStage.name
    expressions.startStage.identifier
    expressions.endStage.name
    expressions.endStage.identifier
    expressions.figures.images.primaryIdentifier
    expressions.figures.images.label
    organism.taxonId
    expressions.assay
    expressions.probe.ThisseCloneRating
    expressions.figures.primaryIdentifier
    expressions.publication.primaryIdentifier
/);
#$query->add_outer_join('expressions.startStage');
#$query->add_outer_join('expressions.endStage');
$query->add_outer_join('expressions.figures.images');
$query->add_outer_join('expressions.probe');
$query->add_outer_join('expressions.publication');

$query->add_sort_order('primaryIdentifier', 'ASC');

$query->add_constraint(
    path  => 'Gene.expressions.fish.wildtype',
    op    => '=',
    value => 'true',
    code  => 'A',
);
$query->add_constraint(
    path  => 'Gene.expressions.environment.StandardEnvironment',
    op    => '=',
    value => 'true',
    code  => 'B',
);
$query->add_constraint(
    path  => 'Gene.expressions.assay',
    op    => '=',
    value => 'mRNA in situ hybridization',
    code  => 'C',
);


my $it = $query->iterator();
while (my $row = <$it>) {
    print join("\t", $row->{'primaryIdentifier'},
                     $row->{'expressions.expressionFound'},
                     $row->{'symbol'},
                     $row->{'expressions.anatomy.name'},
                     $row->{'expressions.anatomy.identifier'},
                     $row->{'expressions.startStage.name'},
                     $row->{'expressions.startStage.identifier'},
                     $row->{'expressions.endStage.name'},
                     $row->{'expressions.endStage.identifier'},
                     $row->{'expressions.figures.images.primaryIdentifier'},
                     $row->{'expressions.figures.images.label'},
                     $row->{'organism.taxonId'},
                     $row->{'expressions.assay'},
              ),"\n";
}

exit 0;

