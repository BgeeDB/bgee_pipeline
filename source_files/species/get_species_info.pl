#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use LWP::Simple;
use List::MoreUtils qw( uniq );

# UniProt REST taxonomy URL
my $uniprot_taxonomy = 'https://www.uniprot.org/taxonomy/';
my $uniprot_type     = '.rdf';

my $taxid   = $ARGV[0]  or die "\n\t$0 <taxid> [comment]\n\n";
my $comment = $ARGV[1]  // '';


# Request
my $content = get("$uniprot_taxonomy$taxid$uniprot_type");
# Parsing
if ( defined $content ){
    my @aliases;
    my ($species, $genus, $common_name) = ('', '','');

    if ( $content =~ /<scientificName>(.+?) (.+?)<\/scientificName>/ ){
        $genus   = $1;
        $species = $2;
    }
    if ( $content =~ /<commonName>(.+?)<\/commonName>/ ){
        $common_name = $1;
    }

    if ( $content =~ /<mnemonic>(.+?)<\/mnemonic>/ ){
        push @aliases, $1;
    }
    while ( $content =~ /<otherName>(.+?)<\/otherName>/g ){
        push @aliases, $1;
    }

# Print
#speciesId    genus    species    speciesCommonName    displayOrder    genomeFilePath    genomeVersion    dataSourceId    genomeSpeciesId    fakeGeneIdPrefix    keywords    comment
    print join("\t", $taxid,
                     $genus,
                     $species,
                     $common_name,
                     '',
                     lc($genus).'_'.lc($species).'/TOCOMPLETE',
                     '',
                     '',
                     $taxid,
                     '',
                     join('|', uniq sort @aliases),
                     $comment,
              ), "\n";
}
else {
    warn "Cannot get [$taxid] info\n";
}

exit 0;

