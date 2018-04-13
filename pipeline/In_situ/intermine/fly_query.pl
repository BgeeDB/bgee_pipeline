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
#  * http://intermine.readthedocs.org/en/latest/web-services
#       - General Usage
#  * http://iodoc.labs.intermine.org
#       - Reference documentation for the underlying REST API
#
######################################################################

use strict;
use warnings;

# List of supported in situ data sources
my @in_situ_src = ('BDGP', 'fly-FISH');

# Print unicode to standard out
binmode(STDOUT, 'utf8');

# This code makes use of the Webservice::InterMine library.
# The following import statement sets FlyMine as your default
use Webservice::InterMine 1.0405 'http://www.flymine.org/flymine';

my $query = new_query(class => 'Gene');

# The view specifies the output columns
$query->add_view(qw/
    primaryIdentifier
    secondaryIdentifier
    symbol
    mRNAExpressionResults.expressed
    mRNAExpressionResults.stageRange
    mRNAExpressionResults.mRNAExpressionTerms.name
    mRNAExpressionResults.stages.name
    mRNAExpressionResults.stages.identifier
    mRNAExpressionResults.mRNAExpressionTerms.identifier
    organism.taxonId
    mRNAExpressionResults.dataSet.dataSource.name
/);

# edit the line below to change the sort order:
$query->add_sort_order('primaryIdentifier', 'ASC');

#$query->set_logic('A and B');

#TODO publication page to link to ? instead of mRNAExpressionResults.images.url

# Use an iterator to avoid having all rows in memory at once.
my $it = $query->iterator();
while ( my $row = <$it> ){
    # Should only be in situ but force in situ anyway for future releases
    next  if ( ! grep { /^$row->{'mRNAExpressionResults.dataSet.dataSource.name'}$/ } @in_situ_src );

    print join("\t",  $row->{'primaryIdentifier'},
                      $row->{'secondaryIdentifier'}                                  || '',
                      $row->{'symbol'}                                               || '',
                      $row->{'mRNAExpressionResults.expressed'},                     #NOTE some false exist!
                      $row->{'mRNAExpressionResults.stageRange'},
                      $row->{'mRNAExpressionResults.mRNAExpressionTerms.name'},
                      $row->{'mRNAExpressionResults.stages.name'},
                      $row->{'mRNAExpressionResults.stages.identifier'},
                      $row->{'mRNAExpressionResults.mRNAExpressionTerms.identifier'} || '', #FIXME NEED organ id !!! always empty in release v43.0 Jul 2016
#                     $row->{'mRNAExpressionResults.images.url'}, #FIXME useful ???
                      $row->{'organism.taxonId'},
                      $row->{'mRNAExpressionResults.dataSet.dataSource.name'},
                      'wild type',     # Rachel Lyne <rachel@intermine.org>: "As far as I can see, the BDGP in situ project used wild type Canton S flies while the flyFish project used Wild-type Oregon R flies."
                      'not annotated', # No Sex info
              ), "\n";
}

exit 0;

