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
# once per slice of genes, a condition on ge.bgeeGeneId being appended to their query.
# The slices are ranges of bgeeGeneId, the first column of the primary key of globalExpression, so
# that MySQL reads only the rows of the slice (with a condition on globalExpressionId, it read the
# whole table for each slice). The first slice has no lower limit and the last one no upper limit:
# each row of globalExpression is read by exactly one slice, whatever the bounds given.
# The other mappings, reading small tables, go together in a part named "dim".
# The triples produced by a mapping depend only on the rows its query returns, so the union of the
# triples of all parts is the one of the whole mapping, asserted or inferred from the ontology
# (the test mode of bgee_rdf_parallel.sh checks it).
# The path of each mapping file written is printed on its own line.
#############################################################

# Define arguments & their default value
my ($mapping, $prefix, $bounds, $test_window) = ('', '', '', '');
my %opts = ('mapping=s'     => \$mapping,     # Path to the Ontop mapping file
            'prefix=s'      => \$prefix,      # Prefix of the paths of the mapping files written
            'bounds=s'      => \$bounds,      # First bgeeGeneId of each slice, comma-separated and ascending
            'test_window=s' => \$test_window  # Test mode: "low high", only the genes with low <= bgeeGeneId < high
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
my @first_genes = split(/,/, $bounds);
my ($test_low, $test_high) = $test_window =~ /^\s*(\d+)\s+(\d+)\s*$/;
my $valid_bounds = @first_genes && !grep { !/^\d+$/ } @first_genes;
$valid_bounds &&= !grep { $first_genes[$_] <= $first_genes[$_ - 1] } 1 .. $#first_genes;
# In test mode the slices must stay within the window, to be compared with the whole window
my $valid_window = $test_window eq '' || ($valid_bounds && defined $test_high && $test_low < $test_high &&
                                          $first_genes[0] >= $test_low && $first_genes[-1] < $test_high);
if ( !$test_options || $mapping eq '' || $prefix eq '' || !$valid_bounds || !$valid_window ) {
    print "\n\tInvalid or missing argument:
\te.g. $0 -mapping=/path/to/mapping.obda -prefix=/path/to/parts/mapping -bounds=257972,410322,523117
\t-mapping      Path to the Ontop mapping file
\t-prefix       Prefix of the paths of the mapping files written
\t-bounds       First bgeeGeneId of each slice, comma-separated and ascending
\t-test_window  Test mode (optional): \"low high\", keep only the genes with low <= bgeeGeneId < high,
\t              the bounds being within this window. The mappings of the other tables are not written
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

# In test mode only the genes of the window are read, and the reference part reads them all in one go
if ($test_window ne '') {
    write_mapping("$prefix.ref.obda", map { restrict_entry($_,
        "ge.bgeeGeneId >= $test_low AND ge.bgeeGeneId < $test_high") } @expression_entries);
} else {
    write_mapping("$prefix.dim.obda", @dim_entries);
}
for my $slice (0 .. $#first_genes) {
    # Outside the test mode, the first slice has no lower limit and the last one no upper limit
    my @conditions = ();
    if ($slice > 0) {
        push @conditions, "ge.bgeeGeneId >= $first_genes[$slice]";
    } elsif ($test_window ne '') {
        push @conditions, "ge.bgeeGeneId >= $test_low";
    }
    if ($slice < $#first_genes) {
        push @conditions, "ge.bgeeGeneId < $first_genes[$slice + 1]";
    } elsif ($test_window ne '') {
        push @conditions, "ge.bgeeGeneId < $test_high";
    }
    my $condition = join(' AND ', @conditions);
    printf STDERR "Slice %02d: %s\n", $slice, $condition || 'all the rows';
    write_mapping(sprintf('%s.s%02d.obda', $prefix, $slice),
        map { $condition? restrict_entry($_, $condition): $_ } @expression_entries);
}

exit 0;
