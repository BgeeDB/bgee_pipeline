#!/usr/bin/env perl
use strict;
use warnings;

# Print unicode to standard out
binmode(STDOUT, 'utf8');

# The following import statement sets WormMine as your default
use Webservice::InterMine 1.0405 'http://intermine.wormbase.org/tools/wormmine';

my $query = new_query(class => 'Gene');

#FIXME How to filter for in situ only ?
#      "in WormMine you cannot select by ‘type of experiment’ -i.e. ISH, IHC, reporter etc..- yet."  Daniela Raciti <draciti@caltech.edu>
#TODO strain and sex ???
#      "Wild type strain in C. elegans is N2"  Daniela Raciti <draciti@caltech.edu>
#      "We generally do not capture sex information. C. elegans has 2 natural sexes, males and hermaphrodites. The vast majority of the experiments are done in hermaphrodites. We do capture ‘male’ in the rare instances authors analyze expression specifically in males."  Daniela Raciti <draciti@caltech.edu>
#TODO How to get Wild Type strains only in WormMine?

# The view specifies the output columns
$query->add_view(qw/
    primaryIdentifier
    symbol
    expressionPatterns.primaryIdentifier
    expressionPatterns.lifeStages.anatomyTerms.name
    expressionPatterns.lifeStages.anatomyTerms.primaryIdentifier
    expressionPatterns.lifeStages.publicName
    expressionPatterns.lifeStages.primaryIdentifier
    strains.primaryIdentifier
    organism.taxonId
/);

# edit the line below to change the sort order:
$query->add_sort_order('primaryIdentifier', 'ASC');

# Edit the code below to specify your own custom logic:
# $query->set_logic('A and B');

# Use an iterator to avoid having all rows in memory at once.
my $it = $query->iterator();
while (my $row = <$it>){
    print join("\t", $row->{'primaryIdentifier'},
                     $row->{'symbol'},
                     $row->{'expressionPatterns.primaryIdentifier'},
                     $row->{'expressionPatterns.lifeStages.anatomyTerms.name'}                || '', #NOTE May contain useful sex info such as "adult hermaphrodite Ce"
                     $row->{'expressionPatterns.lifeStages.anatomyTerms.primaryIdentifier'},
                     $row->{'expressionPatterns.lifeStages.publicName'},
                     $row->{'expressionPatterns.lifeStages.primaryIdentifier'},
                     $row->{'strains.primaryIdentifier'},
                     $row->{'organism.taxonId'}
              ),"\n";
}

exit 0;

