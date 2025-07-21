#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use List::Compare;
use File::Slurp;

my $with_false = $ARGV[0]  || die "\n\t$0 <with_false> <without>\n\n";
my $without    = $ARGV[1]  || die "\n\t$0 <with_false> <without>\n\n";

my @with_false = read_file($with_false, chomp => 1);
my @without    = read_file($without,    chomp => 1);


#NOTE Do not display/view the GXDExpression.genotype.alleles.isWildType field!
#     It makes disappear most rows if displayed!
#NOTE Exists also "In situ reporter (knock in)"

#NOTE GXDExpression.genotype.alleles.isWildType=true does not what it should!
# Get list without GXDExpression.genotype.alleles.isWildType constraint,
# then substract list with GXDExpression.genotype.alleles.isWildType=false from first list
# to retrieve what we expect!!!
my $lc = List::Compare->new(\@without, \@with_false);
my @Lonly = $lc->get_unique; # Get those items which appear (at least once) only in the first list


# Header
# Print unicode to standard out
binmode(STDOUT, 'utf8');
print "#Assay type\tfeature.symbol\tfeature.id\tfeature.xref.id\tstage\tage\tstructure.name\tstructure.id\tstrength\tpattern\tassayId\tprobe\timage\tspecimen\tsex\tdetected\ttaxonId\tgenotype.symbol\n";
ROW:
for my $row ( @Lonly ){
    $row =~ s{"}{}g;
    my @items = map { trim($_) } split(/\t/, $row, -1);

    print join("\t", @items), "\n"
}

exit 0;


sub trim {
    my ($string) = @_;

    $string =~ s{^\s+}{};
    $string =~ s{\s+$}{};

    return $string eq 'undef' ? '' : $string;
}

