#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use JSON::XS;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

# ONE-OFF FIX FOR Bgee 16.0 -- DELETE THIS SCRIPT ONCE IT HAS BEEN RUN.
#
# The genes of Bgee 16.0 were inserted before the gene.geneLength column existed, so the
# column has to be filled afterwards. This script does it for one species, from the same
# Ensembl JSON dump insert_ensembl_genes.pl uses, and with the very same calculation
# (Utils::median_isoform_length).
#
# From Bgee 16.1 onwards insert_ensembl_genes.pl fills geneLength on its own, so this
# script must NOT be kept around: it would only be a second, divergent way of writing a
# column the genes pipeline already handles. See the fix_gene_length rule of the Makefile.
#
# Only the genes of the species given by -species are updated, matched on geneId, which
# mirrors how insert_ensembl_genes.pl inserted them. Nothing else in the gene table is
# touched, so the script can be re-run: it is idempotent.
#############################################################


# Define arguments & their default value
my ($species, $bgee_connector) = ('', '');
my ($specific_gene) = ('');
my ($debug) = (0);
my %opts = ('species=s'    => \$species,            # speciesId from Bgee db, as for insert_ensembl_genes.pl
            'bgee=s'       => \$bgee_connector,     # Bgee connector string
            'debug'        => \$debug,              # debug mode, do not insert/update in database
            'gene=s'       => \$specific_gene,      # Query a single gene
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $species eq '' || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -species=9606__0__Homo_sapiens__Ensembl  -bgee=\$(BGEECMD)  <Ensembl species JSON>
\t-species   speciesId from Bgee db with the genomeSpeciesId concatenated
\t-bgee      Bgee    connector string
\t-debug     Debug mode, do not insert/update in database
\t
\t-gene      Query a single specific gene (optional)
\n";
    exit 1;
}


## Load the JSON
my $input_json_file = $ARGV[0]  or die "\n\tMissing input JSON file\n\n";
my $json_text = do {
    open(my $json_fh, '<:encoding(UTF-8)', $input_json_file)  or die("Can't open \"$input_json_file\": $!\n");
    local $/;
    <$json_fh>;
};
my $ensembl_json = decode_json($json_text);
#NOTE keep only the genes section of the JSON!
map { delete( $ensembl_json->{$_} ) } grep { $_ ne 'genes' } keys %$ensembl_json;
my @genes = @{ $ensembl_json->{'genes'} };


# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


# Same -species format as insert_ensembl_genes.pl, so that the same string can be reused.
# Only the Bgee species ID matters here: for a species whose genome is borrowed from
# another one, the genes were inserted under the Bgee species ID, not the genome one.
my ($speciesBgee, $newSpecies, $scientific_name, $ensSource) = split('__', $species, -1);
die "\n\tInvalid speciesId [$speciesBgee] in -species\n\n"  if ( $speciesBgee !~ /^\d+$/ );


## Genes already in Bgee for that species, to tell apart what the JSON does not cover
## from what Bgee does not know about
my $geneDB = $dbh->prepare('SELECT geneId FROM gene WHERE speciesId = ?');
$geneDB->execute($speciesBgee)  or die $geneDB->errstr;
my %bgeeGenes = map { $_->[0] => 1 } @{$geneDB->fetchall_arrayref};
$geneDB->finish;
die "\n\tNo gene in Bgee for species [$speciesBgee]\n\n"  if ( !%bgeeGenes );
print "Species $speciesBgee: ", scalar(keys %bgeeGenes), " genes in Bgee, ", scalar(@genes), " in the JSON\n";


## Update the gene lengths
my $lengthDB = $dbh->prepare('UPDATE gene SET geneLength = ? WHERE geneId = ? AND speciesId = ?');
my ($updated, $no_length, @not_in_bgee) = (0, 0);
my %seen;
for my $gene ( @genes ){
    my $stable_id = $gene->{'id'};
    if ( $specific_gene ){
        next  if ( $stable_id ne $specific_gene );
    }

    if ( !exists $bgeeGenes{$stable_id} ){
        push @not_in_bgee, $stable_id;
        next;
    }
    $seen{$stable_id} = 1;

    my $gene_length = Utils::median_isoform_length($gene);
    if ( !defined $gene_length ){
        $no_length++;
        next;
    }

    if ( ! $debug ){
        $lengthDB->execute($gene_length, $stable_id, $speciesBgee)  or die $lengthDB->errstr;
    }
    else {
        print "UPDATE gene SET geneLength = $gene_length WHERE geneId = '$stable_id' AND speciesId = $speciesBgee\n";
    }
    $updated++;
}
$lengthDB->finish;
$dbh->disconnect;


## Report, so that a run that silently did almost nothing cannot go unnoticed
my @not_in_json = grep { !exists $seen{$_} } sort keys %bgeeGenes;
print "Updated geneLength for $updated gene(s)", ($debug ? ' (debug mode, nothing written)' : ''), "\n";
if ( $no_length ){
    warn "Warning: $no_length gene(s) of the JSON carry no exon, their geneLength is left untouched\n";
}
if ( @not_in_bgee ){
    warn "Warning: ", scalar(@not_in_bgee), " gene(s) of the JSON are not in Bgee for species $speciesBgee",
         " (e.g. ", join(', ', @not_in_bgee[0..($#not_in_bgee < 4 ? $#not_in_bgee : 4)]), ")\n";
}
if ( @not_in_json && !$specific_gene ){
    warn "Warning: ", scalar(@not_in_json), " gene(s) of Bgee are absent from the JSON, their geneLength stays NULL",
         " (e.g. ", join(', ', @not_in_json[0..($#not_in_json < 4 ? $#not_in_json : 4)]), ")\n";
}
# A run matching almost nothing usually means the wrong JSON was given for that species
die "\n\tOnly $updated gene(s) matched out of ", scalar(keys %bgeeGenes), ": wrong JSON for species $speciesBgee?\n\n"
    if ( !$specific_gene && $updated < scalar(keys %bgeeGenes) / 2 );

exit 0;

