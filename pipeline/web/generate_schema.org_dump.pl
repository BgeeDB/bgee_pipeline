#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use File::Slurp;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;


# Define arguments & their default value
my ($bgee_connector) = ('');
my %opts = ('bgee=s'     => \$bgee_connector,   # Bgee connector string
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)  >bgee.schema_org.jsonld
\t-bgee     Bgee connector string
\n";
    exit 1;
}


# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);
my ($bgee_version) = $bgee_connector =~ /^.+__name=bgee_v(\d+_\d+)$/;
$bgee_version =~ s{_}{.}g;

$| = 1;


# Homepage
my $schema_default = '';

# species pages
my @schema_species;

# gene pages
my @schema_gene;

$bgee->disconnect;


print '[', join(',', map { s|</?script.*||g; $_ } $schema_default, @schema_species, @schema_gene), ']';

exit 0;



sub get_schema_default {
    return '{
    "@context": "https://schema.org/",
    "@id": "https://bgee.org/",
    "@graph": [
        {
            "@type": "Organization",
            "@id": "https://bgee.org/",
            "name": "Bgee - Bring Gene Expression Expertise",
            "url": "https://bgee.org/",
            "description": "The aim of Bgee is to help biologists to use and understand gene expression",
            "logo": "https://bgee.org/img/logo/bgee13_hp_logo.png",
            "sameAs": [
                "https://twitter.com/Bgeedb",
                "https://genomic.social/@bgeedb",
                "https://bgeedb.wordpress.com/"
            ],
            "parentOrganization": [
                {
                    "@type": "Organization",
                    "@id": "https://www.sib.swiss",
                    "name": "SIB Swiss Institute of Bioinformatics",
                    "url": "https://www.sib.swiss",
                    "sameAs": "https://en.wikipedia.org/wiki/Swiss_Institute_of_Bioinformatics"
                },
                {
                    "@type": "CollegeOrUniversity",
                    "@id": "https://unil.ch",
                    "name": "UNIL University of Lausanne",
                    "url": "https://unil.ch",
                    "sameAs": "https://en.wikipedia.org/wiki/University_of_Lausanne"
                },
                {
                    "@type": "EducationalOrganization",
                    "@id": "https://www.unil.ch/dee/robinson-rechavi-group",
                    "name": "Evolutionary Bioinformatics group",
                    "url": "https://www.unil.ch/dee/robinson-rechavi-group"
                }
            ]
        },
        {
            "@type": "Dataset",
            "@id": "https://bgee.org/",
            "http://purl.org/dc/terms/conformsTo": {
                "@id": "https://bioschemas.org/profiles/Dataset/1.0-RELEASE",
                "@type": "CreativeWork"
            },
            "url": "https://bgee.org/",
            "name": "Bgee gene expression data",
            "description": "Bgee is a database for retrieval and comparison of gene expression patterns across multiple animal species. It provides an intuitive answer to the question -where is a gene expressed?- and supports research in cancer and agriculture as well as evolutionary biology.",
            "keywords": [
                "bgee",
                "gene expression",
                "evolution",
                "ontology",
                "anatomy",
                "development",
                "evo-devo database",
                "anatomical ontology",
                "developmental ontology",
                "gene expression evolution"
            ],
            "creator": {
                "@id": "https://bgee.org/"
            },
            "license": "https://creativecommons.org/publicdomain/zero/1.0/",
            "version": "'.$bgee_version.'"
        }
    ]
}';
}


