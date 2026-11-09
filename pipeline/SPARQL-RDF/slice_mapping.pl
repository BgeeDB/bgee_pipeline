#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;

# Created Nov. 2026
# Splits the Ontop mapping of Bgee into parts that independent Ontop processes can materialize
# in parallel (see bgee_rdf_parallel.sh).
# One Ontop process serializes the RDF triples with a single thread: the materialization of the
# globalExpression table of easyBgee, with more than a billion rows, takes days while MySQL stays
# almost idle. The mappings reading this table (aliased "ge" in their SQL query) are thus copied
# once per slice of globalExpressionId, a condition on ge.globalExpressionId being appended to their
# query. The slices partition the range of identifiers: each row of globalExpression is read by
# exactly one slice. The other mappings, reading small tables, go together in a part named "dim".
# The triples produced by a mapping depend only on the rows its query returns, so the union of the
# triples of all parts is the one of the whole mapping, asserted or inferred from the ontology
# (the test mode of bgee_rdf_parallel.sh checks it).
# The path of each mapping file written is printed on its own line.
#############################################################

# Define arguments & their default value
my ($mapping, $prefix, $slice_count, $min_id, $max_id, $id_limit) = ('', '', 0, 0, 0, 0);
my %opts = ('mapping=s'  => \$mapping,     # Path to the Ontop mapping file
            'prefix=s'   => \$prefix,      # Prefix of the paths of the mapping files written
            'slices=i'   => \$slice_count, # Number of slices of the globalExpression table
            'min_id=i'   => \$min_id,      # Lowest globalExpressionId
            'max_id=i'   => \$max_id,      # Highest globalExpressionId
            'id_limit=i' => \$id_limit     # Test mode: only the rows with a lower globalExpressionId
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $mapping eq '' || $prefix eq '' || $slice_count < 1 || $max_id < 1 ||
        $min_id > $max_id || ($id_limit && $id_limit <= $min_id) ) {
    print "\n\tInvalid or missing argument:
\te.g. $0 -mapping=/path/to/mapping.obda -prefix=/path/to/parts/mapping -slices=12 -min_id=1 -max_id=1385000000
\t-mapping    Path to the Ontop mapping file
\t-prefix     Prefix of the paths of the mapping files written
\t-slices     Number of slices of the globalExpression table
\t-min_id     Lowest globalExpressionId
\t-max_id     Highest globalExpressionId
\t-id_limit   Test mode (optional): keep only the rows with a globalExpressionId lower than this
\t            value, the mappings of the other tables are not written
\n";
    exit 1;
}

##########################################
# READ THE MAPPING
##########################################
open(my $fh, '<', $mapping) || die "failed to read mapping file: $!";
my @lines = <$fh>;
close($fh);
chomp(@lines);

my ($start, $end);
for my $i (0 .. $#lines) {
    $start = $i if !defined $start && $lines[$i] =~ /^\[MappingDeclaration\]/;
    $end = $i if $lines[$i] =~ /^\]\]\s*$/;
}
die "no [MappingDeclaration] section found in $mapping\n"
    if !defined $start || !defined $end || $end <= $start;
my @header = @lines[0 .. $start];
my @footer = @lines[$end .. $#lines];

# One entry per mapping: its mappingId, target and source lines, until the next blank line
my @entries = ();
my @current = ();
for my $line (@lines[$start + 1 .. $end - 1]) {
    if ($line =~ /^\s*$/) {
        push @entries, [@current] if @current;
        @current = ();
    } else {
        push @current, $line;
    }
}
push @entries, [@current] if @current;

my (@dim_entries, @expression_entries);
for my $entry (@entries) {
    my ($mapping_id) = $entry->[0] =~ /^mappingId\s+(\S+)/;
    die "unexpected line, not a mappingId: $entry->[0]\n" if !defined $mapping_id;
    my ($source_idx) = grep { $entry->[$_] =~ /^source\s/ } 0 .. $#$entry;
    die "no source for the mapping $mapping_id\n" if !defined $source_idx;
    my $sql = join(' ', @{$entry}[$source_idx .. $#$entry]);
    if ($sql =~ /\bglobalExpression\s+as\s+ge\b/i) {
        die "no WHERE clause in the source of the mapping $mapping_id\n" if $sql !~ /\bWHERE\b/i;
        push @expression_entries, $entry;
    } elsif ($sql =~ /\bglobalExpression\b/i) {
        die "the mapping $mapping_id reads globalExpression without the alias ge, ",
            "the condition on the identifiers cannot be added\n";
    } else {
        push @dim_entries, $entry;
    }
}
die "no mapping reads globalExpression, nothing to split\n" if !@expression_entries;
print STDERR scalar(@expression_entries), " mappings read globalExpression, ",
             scalar(@dim_entries), " do not\n";

##########################################
# WRITE THE PARTS
##########################################
# Returns a copy of $entry whose source query is restricted by $condition
sub restrict_entry {
    my ($entry, $condition) = @_;
    my @restricted = @$entry;
    $restricted[-1] =~ s/\s+$//;
    $restricted[-1] .= " AND $condition";
    return \@restricted;
}
sub write_mapping {
    my ($file, @part_entries) = @_;
    open(my $out, '>', $file) || die "failed to write $file: $!";
    print $out join("\n", @header), "\n";
    print $out join("\n\n", map { join("\n", @$_) } @part_entries), "\n";
    print $out join("\n", @footer), "\n";
    close($out);
    print "$file\n";
}

# In test mode the identifiers are limited, and the reference part reads them all in one go
my $last_id = $id_limit? $id_limit - 1: $max_id;
my $total = $last_id - $min_id + 1;
my $width = int(($total + $slice_count - 1) / $slice_count);

if ($id_limit) {
    write_mapping("$prefix.ref.obda", map { restrict_entry($_,
        "ge.globalExpressionId >= $min_id AND ge.globalExpressionId < $id_limit") } @expression_entries);
} else {
    write_mapping("$prefix.dim.obda", @dim_entries);
}
for my $slice (0 .. $slice_count - 1) {
    my $low = $min_id + $slice * $width;
    last if $low > $last_id;
    my $high = $low + $width > $last_id? $last_id + 1: $low + $width;
    printf STDERR "Slice %02d: globalExpressionId from %d to %d\n", $slice, $low, $high - 1;
    write_mapping(sprintf('%s.s%02d.obda', $prefix, $slice), map { restrict_entry($_,
        "ge.globalExpressionId >= $low AND ge.globalExpressionId < $high") } @expression_entries);
}

exit 0;
