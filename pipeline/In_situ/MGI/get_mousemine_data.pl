#!/usr/bin/env perl

######################################################################
# This is an automatically generated script to run your query.
# To use it you will require the InterMine Perl client libraries.
# These can be installed from CPAN, using your preferred client, eg:
#
#    sudo cpan Webservice::InterMine
#
# For help using these modules, please see these resources:
#
#  * https://metacpan.org/pod/Webservice::InterMine
#       - API reference
#  * https://metacpan.org/pod/Webservice::InterMine::Cookbook
#       - A How-To manual
#  * http://www.intermine.org/wiki/PerlWebServiceAPI
#       - General Usage
#  * http://www.intermine.org/wiki/WebService
#       - Reference documentation for the underlying REST API
#
######################################################################

use strict;
use warnings;
use diagnostics;

use List::Compare;

# Print unicode to standard out
binmode(STDOUT, 'utf8');

# This code makes use of the Webservice::InterMine library.
# The following import statement sets MouseMine as your default
use Webservice::InterMine 1.0405 'http://www.mousemine.org/mousemine';


# Query to get all "RNA in situ" assays with isWildType field == false
my $query_false = new_query(class => 'GXDExpression');

#NOTE Do not display/view the GXDExpression.genotype.alleles.isWildType field!
#     It makes disappear most rows if displayed!
# The view specifies the output columns
$query_false->add_view(qw/
    assayType
    feature.symbol
    feature.primaryIdentifier
    feature.crossReferences.identifier
    stage
    age
    structure.name
    structure.identifier
    strength
    pattern
    assayId
    probe
    image
    specimenNum
    sex
    detected
    feature.organism.taxonId
/);

# Your custom sort order is specified with the following code:
$query_false->add_sort_order('assayId', 'ASC');

#NOTE Exists also "In situ reporter (knock in)"
$query_false->add_constraint(
    path  => 'GXDExpression.assayType',
    op    => '=',
    value => 'RNA in situ',
    code  => 'A',
);
# GXDExpression.genotype.alleles.isWildType=true does not what it should!
# Get list without GXDExpression.genotype.alleles.isWildType constraint,
# then substract list with GXDExpression.genotype.alleles.isWildType=false from first list
# to retrieve what we expect!!!
$query_false->add_constraint(
    path  => 'GXDExpression.genotype.alleles.isWildType',
    op    => '=',
    value => 'false',
    code  => 'B',
);
$query_false->add_constraint(
    path  => 'GXDExpression.feature.crossReferences.source.name',
    op    => '=',
    value => 'Ensembl Gene Model',
    code  => 'C',
);

# Edit the code below to specify your own custom logic:
#$query->set_logic('B and A and D and C');


my $query = new_query(class => 'GXDExpression');
# The view specifies the output columns
$query->add_view(qw/
    assayType
    feature.symbol
    feature.primaryIdentifier
    feature.crossReferences.identifier
    stage
    age
    structure.name
    structure.identifier
    strength
    pattern
    assayId
    probe
    image
    specimenNum
    sex
    detected
    feature.organism.taxonId
/);

# Your custom sort order is specified with the following code:
$query->add_sort_order('assayId', 'ASC');

$query->add_constraint(
    path  => 'GXDExpression.assayType',
    op    => '=',
    value => 'RNA in situ',
    code  => 'A',
);
$query->add_constraint(
    path  => 'GXDExpression.feature.crossReferences.source.name',
    op    => '=',
    value => 'Ensembl Gene Model',
    code  => 'B',
);

my @with_false = $query_false->results(['as' => 'strings',]);
my @without    = $query->results(['as' => 'strings',]);
#warn scalar @without, ' ... ', scalar @with_false, "\n";

my $lc = List::Compare->new(\@without, \@with_false);
my @Lonly = $lc->get_unique;
#warn scalar @Lonly, "\n";


# Header
print "#Assay type\tfeature.symbol\tfeature.id\tfeature.xref.id\tstage\tage\tstructure.name\tstructure.id\tstrength\tpattern\tassayId\tprobe\timage\tspecimen\tsex\tdetected\ttaxonId\n";
for my $row ( @Lonly ){
    my %fields = map  { /^([^:]+): (.*)$/; $1 => trim($2) }
                 grep { !/^GXDExpression$/ } # Do not use class name - useless here - and in a different data format
                 split(/\t/, $row, -1);
    print join("\t", $fields{'assayType'},
                     $fields{'feature.symbol'},
                     $fields{'feature.primaryIdentifier'},
                     $fields{'feature.crossReferences.identifier'},
                     $fields{'stage'},
                     $fields{'age'},
                     $fields{'structure.name'},
                     $fields{'structure.identifier'},
                     $fields{'strength'},
                     $fields{'pattern'},
                     $fields{'assayId'},
                     $fields{'probe'},
                     $fields{'image'},
                     $fields{'specimenNum'},
                     $fields{'sex'},
                     $fields{'detected'},
                     $fields{'feature.organism.taxonId'},
              ), "\n";
}

exit 0;


sub trim {
    my ($string) = @_;

    $string =~ s{^\s+}{};
    $string =~ s{\s+$}{};

    return $string eq 'undef' ? '' : $string;
}

__END__
# Use an iterator to avoid having all rows in memory at once.
my $it = $query->iterator();
while (my $row = <$it>) {
#    print $row->{'assayType'}, $row->{'feature.symbol'}, $row->{'feature.primaryIdentifier'},
#        $row->{'feature.crossReferences.identifier'}, $row->{'feature.crossReferences.source.name'},
#        $row->{'stage'}, $row->{'age'}, $row->{'structure.name'}, $row->{'structure.identifier'},
#        $row->{'strength'}, $row->{'pattern'}, $row->{'assayId'}, $row->{'probe'}, $row->{'image'},
#        $row->{'genotype.alleles.isWildType'}, "\n";
    #SEB print only columns not set/frozen by a constraint
    print $row->{'assayType'}, $row->{'feature.symbol'}, $row->{'feature.primaryIdentifier'},
          $row->{'feature.crossReferences.identifier'},
          $row->{'stage'}, $row->{'age'}, $row->{'structure.name'}, $row->{'structure.identifier'},
          $row->{'strength'}, $row->{'pattern'}, $row->{'assayId'}, $row->{'probe'}, $row->{'image'}, $row->{'specimenNum'},
          "\n";
}

exit 0;

