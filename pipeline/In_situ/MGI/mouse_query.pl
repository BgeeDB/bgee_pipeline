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


# Query to obtain all "structures" for "RNA in situ" assays
######################################################
my $structure_query = new_query(class => 'GXDExpression');

# The view specifies the output columns
$structure_query->add_view(qw/
    structure.identifier
/);

$structure_query->add_constraint(
    path  => 'GXDExpression.assayType',
    op    => '=',
    value => 'RNA in situ',
    code  => 'A',
);

my %structures;
my $result = $structure_query->iterator();
while ( my $row = <$result> ){
        $structures{ $row->{'structure.identifier'} }++;
}

#warn 'Total uniq structures: ', scalar keys %structures, "\n";


# Query to get all "RNA in situ" assays with isWildType field == false
######################################################
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
    detected
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

my @result_false = $query_false->results_iterator(as => 'string');

my $it = $query_false->results_iterator(as => 'string');
    while ( my $row = <$it> ){
        # handle row as string
        push @result_false, $row;
    }

#warn 'Result with isWildType == False: ', scalar @result_false, "\n";


# Query to get all "RNA in situ" assays
# Have to iterate with "structure" because complete query it's too big
############################################################################
my @array;
my @result_all;
for ( sort keys %structures ){
    push @array, $_;
    next  if ( scalar @array < 10 );

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
        detected
    /);

    # Your custom sort order is specified with the following code:
    $query->add_sort_order('assayId', 'ASC');

    $query->add_constraint(
        path  => 'GXDExpression.structure.identifier',
        op    => '=',
        value => "$array[0]",
        code  => 'A',
    );
    $query->add_constraint(
        path  => 'GXDExpression.structure.identifier',
        op    => '=',
        value => "$array[1]",
        code  => 'B',
    );
    $query->add_constraint(
        path  => 'GXDExpression.structure.identifier',
        op    => '=',
        value => "$array[2]",
        code  => 'C',
    );
    $query->add_constraint(
        path  => 'GXDExpression.structure.identifier',
        op    => '=',
        value => "$array[3]",
        code  => 'D',
    );
    $query->add_constraint(
        path  => 'GXDExpression.structure.identifier',
        op    => '=',
        value => "$array[4]",
        code  => 'E',
    );
    $query->add_constraint(
        path  => 'GXDExpression.structure.identifier',
        op    => '=',
        value => "$array[5]",
        code  => 'F',
    );
    $query->add_constraint(
        path  => 'GXDExpression.structure.identifier',
        op    => '=',
        value => "$array[6]",
        code  => 'G',
    );
    $query->add_constraint(
        path  => 'GXDExpression.structure.identifier',
        op    => '=',
        value => "$array[7]",
        code  => 'H',
    );
    $query->add_constraint(
        path  => 'GXDExpression.structure.identifier',
        op    => '=',
        value => "$array[8]",
        code  => 'I',
    );
    $query->add_constraint(
        path  => 'GXDExpression.structure.identifier',
        op    => '=',
        value => "$array[9]",
        code  => 'J',
    );
    $query->add_constraint(
        path  => 'GXDExpression.assayType',
        op    => '=',
        value => 'RNA in situ',
        code  => 'K',
    );
    $query->add_constraint(
        path  => 'GXDExpression.feature.crossReferences.source.name',
        op    => '=',
        value => 'Ensembl Gene Model',
        code  => 'L',
    );

    # Your custom logic is specified with the code below:
    $query->set_logic('(A or B or C or D or E or F or G or H or I or J) and K and L');

#    # Print out number of results
#    my $query_count = $query->count;
#    warn "Results in structures: ", $query_count, "\n";

    my $it = $query->results_iterator(as => 'string');
    while (my $row = <$it>) {
        # handle row as string
        push @result_all, $row;
    }

    @array = ();
}


# Query to print the rest of strutures
# If the number of "strutures" is not multiple of 10, you need this query
##########################################################################
if ( scalar @array > 0 ){
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
        detected
    /);

    # Your custom sort order is specified with the following code:
    $query->add_sort_order('assayId', 'ASC');

    my $code = 'A';
    my @logic;
    for my $constraint ( @array ){
        push @logic, $code;
        $query->add_constraint(
            path  => 'GXDExpression.structure.identifier',
            op    => '=',
            value => "$constraint",
            code  => $code++,
        );
    }

    $query->add_constraint(
        path  => 'GXDExpression.assayType',
        op    => '=',
        value => 'RNA in situ',
        code  => 'K',
    );
    $query->add_constraint(
        path  => 'GXDExpression.feature.crossReferences.source.name',
        op    => '=',
        value => 'Ensembl Gene Model',
        code  => 'L',
    );

    # Your custom logic is specified with the code below:
    $query->set_logic('('.join(' or ', @logic).') and K and L');

    ## Print out number of results
    #my $query_count = $query->count;
    #warn "Results in rest structures: ", $query_count, "\n";

    my $result_rest = $query->results_iterator(as => 'string');
    while ( my $row = <$result_rest> ){
        # handle row as string
        push @result_all, $row;
    }
}

#warn scalar @result_all, ' ... ', scalar @result_false, "\n";

my $lc = List::Compare->new(\@result_all, \@result_false);
my @lonly = $lc->get_unique;
#warn scalar @lonly, "\n";

#print join("\n", @lonly);

# Header
print "#Assay type\tfeature.symbol\tfeature.id\tfeature.xref.id\tstage\tage\tstructure.name\tstructure.id\tstrength\tpattern\tassayId\tprobe\timage\tspecimen\tsex\tdetected\ttaxonId\n";
for my $row ( @lonly ){
    print join("\t", map { trim($_) }
                     split(/\t/, $row, -1)
              ), "\n";
}

exit 0;


sub trim {
    my ($string) = @_;

    $string =~ s{^\s+}{};
    $string =~ s{\s+$}{};

    return   $string eq 'undef' ? ''
           : $string eq '""'    ? ''
           :                      $string;
}

