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


#### constants ####
my $dateModified = '2021-07-01';



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
my $bgee_db_version = $bgee_version;
$bgee_version =~ s{_}{.}g;

$| = 1;


# Homepage
my $schema_default = get_schema_default();

# species pages
my @schema_species = '{['.join(',', @{ get_schema_species($bgee) }).']}';

# gene pages
my @schema_gene;
my @schema_gene_homolog;
my @schema_gene_expression;

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

sub get_schema_species {
    my ($db) = @_;

    #TODO Loop over data types with *Processed expression value files*, separated by *,*
    my $template = '{
    "@context": "https://schema.org/",
    "@id": "https://bgee.org/species/__TAXID__",
    "@type": "Taxon",
    "http://purl.org/dc/terms/conformsTo": {
        "@id": "https://bioschemas.org/profiles/Taxon/0.6-RELEASE",
        "@type": "CreativeWork"
    },
    "name": "__SPECIES NAME__",
    "identifier": __TAXID__,
    "sameAs": [
        "http://purl.obolibrary.org/obo/NCBITaxon___TAXID__",
        "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lvl=0&id=__TAXID__"__DBSRC_SPECIES_LINK__
    ],
    "taxonRank": [
        "http://rs.tdwg.org/ontology/voc/TaxonRank#Species",
        "http://purl.uniprot.org/core/Species",
        "http://purl.obolibrary.org/obo/NCBITaxon_species",
        "http://www.wikidata.org/entity/Q7432",
        "species"
    ],
    "subjectOf": [
        {
            "@type": "Dataset",
            "dateModified": "'.$dateModified.'",
            "citation": "https://doi.org/10.1093/nar/gkaa793",
            "description": "__SPECIES NAME__ calls of presence/absence of expression. Each call corresponds to a unique combination of a gene, an anatomical entity, a life stage, a sex, and a strain, with reported presence or absence of expression.",
            "includedInDataCatalog": {
                "@id": "https://bgee.org",
                "@type": "DataCatalog",
                "name": "Bgee"
            },
            "keywords": [
                "gene expression",
                "call",
                "__SPECIES NAME__"__COMMON NAME1__
            ],
            "creator": {
                "@type": "Organization",
                "url": "https://bgee.org/",
                "name": "The Bgee Team"
            },
            "license": "https://creativecommons.org/publicdomain/zero/1.0/",
            "name": "__SPECIES NAME__ gene expression calls",
            "url": "https://bgee.org/species/__TAXID__#expr-calls",
            "version": "'.$bgee_version.'",
            "hasPart": [
                {
                    "@type": "Dataset",
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "https://creativecommons.org/publicdomain/zero/1.0/",
                    "name": "__SPECIES NAME__ gene expression simple",
                    "description": "Anatomical entities only, file without advanced columns.",
                    "url": "https://bgee.org/species/__TAXID__#expr-calls-anat-simple",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/calls/expr_calls/__SPECIES_NAME___expr_simple.tsv.gz"
                        }
                    ]
                },
                {
                    "@type": "Dataset",
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "https://creativecommons.org/publicdomain/zero/1.0/",
                    "name": "__SPECIES NAME__ gene expression advanced",
                    "description": "Anatomical entities only, file with advanced columns.",
                    "url": "https://bgee.org/species/__TAXID__#expr-calls-anat-advanced",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/calls/expr_calls/__SPECIES_NAME___expr_advanced.tsv.gz"
                        }
                    ]
                },
                {
                    "@type": "Dataset",
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "https://creativecommons.org/publicdomain/zero/1.0/",
                    "name": "__SPECIES NAME__ gene expression simple with all conditions",
                    "description": "Anatomical entities, developmental stages, sexes and strains. File without advanced columns.",
                    "url": "https://bgee.org/species/__TAXID__#expr-calls-cond-simple",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/calls/expr_calls/__SPECIES_NAME___expr_simple_all_conditions.tsv.gz"
                        }
                    ]
                },
                {
                    "@type": "Dataset",
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "https://creativecommons.org/publicdomain/zero/1.0/",
                    "name": "__SPECIES NAME__ gene expression advanced with all conditions",
                    "description": "Anatomical entities, developmental stages, sexes and strains. File with advanced columns.",
                    "url": "https://bgee.org/species/__TAXID__#expr-calls-cond-advanced",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/calls/expr_calls/__SPECIES_NAME___expr_advanced_all_conditions.tsv.gz"
                        }
                    ]
                }
            ]
        },
        {
            "@type": "Dataset",
            "dateModified": "'.$dateModified.'",
            "creator": {
                "@type": "Organization",
                "url": "https://bgee.org/",
                "name": "The Bgee Team"
            },
            "citation": "https://doi.org/10.1093/nar/gkaa793",
            "description": "Annotations and experiment information (e.g., annotations to anatomy and development, quality scores used in QCs, library information), and processed expression values (e.g., read counts, TPM and FPKM values) for __SPECIES NAME__",
            "includedInDataCatalog": {
                "@id": "https://bgee.org/bgee'.$bgee_db_version.'",
                "@type": "DataCatalog",
                "name": "Bgee"
            },
            "keywords": [
                "annotations",
                "experiment information",
                "processed expression values",
                "__SPECIES NAME__"__COMMON NAME1__
            ],
            "license": "https://creativecommons.org/publicdomain/zero/1.0/",
            "name": "__SPECIES NAME__ processed expression values",
            "url": "https://bgee.org/species/__TAXID__#proc-values",
            "version": "'.$bgee_version.'",
            "hasPart": [
                {
                    "@type": "Dataset",
                    "http://purl.org/dc/terms/conformsTo": {
                        "@id": "https://bioschemas.org/profiles/Dataset/1.0-RELEASE",
                        "@type": "CreativeWork"
                    },
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "https://creativecommons.org/publicdomain/zero/1.0/",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/processed_expr_values/affymetrix/__SPECIES_NAME__/__SPECIES_NAME___Affymetrix_experiments_chips.tar.gz"
                        }
                    ],
                    "name": "__SPECIES NAME__ Affymetrix experiments chips",
                    "keywords": [
                        "Affymetrix"
                    ],
                    "description": "Affymetrix experiments/chips annotations and metadata.",
                    "url": "https://bgee.org/species/__TAXID__#proc-values-affymetrix"
                },
                {
                    "@type": "Dataset",
                    "http://purl.org/dc/terms/conformsTo": {
                        "@id": "https://bioschemas.org/profiles/Dataset/1.0-RELEASE",
                        "@type": "CreativeWork"
                    },
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "https://creativecommons.org/publicdomain/zero/1.0/",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/processed_expr_values/affymetrix/__SPECIES_NAME__/__SPECIES_NAME___Affymetrix_probesets.tar.gz"
                        }
                    ],
                    "name": "__SPECIES NAME__ Affymetrix probesets",
                    "description": "__SPECIES NAME__ Affymetrix probesets, data (signal intensities).",
                    "url": "https://bgee.org/species/__TAXID__#proc-values-affymetrix"
                },
                {
                    "@type": "Dataset",
                    "http://purl.org/dc/terms/conformsTo": {
                        "@id": "https://bioschemas.org/profiles/Dataset/1.0-RELEASE",
                        "@type": "CreativeWork"
                    },
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "https://creativecommons.org/publicdomain/zero/1.0/",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/processed_expr_values/rna_seq/__SPECIES_NAME__/__SPECIES_NAME___RNA-Seq_experiments_libraries.tar.gz"
                        }
                    ],
                    "name": "__SPECIES NAME__ RNA-Seq experiment libraries",
                    "keywords": [
                        "RNA-Seq"
                    ],
                    "description": "__SPECIES NAME__ RNA-Seq experiments/libraries annotations and metadata.",
                    "url": "https://bgee.org/species/__TAXID__#proc-values-rna-seq"
                },
                {
                    "@type": "Dataset",
                    "http://purl.org/dc/terms/conformsTo": {
                        "@id": "https://bioschemas.org/profiles/Dataset/1.0-RELEASE",
                        "@type": "CreativeWork"
                    },
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "https://creativecommons.org/publicdomain/zero/1.0/",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/processed_expr_values/rna_seq/__SPECIES_NAME__/__SPECIES_NAME___RNA-Seq_read_counts_TPM_FPKM.tar.gz"
                        }
                    ],
                    "name": "__SPECIES NAME__ RNA-Seq read counts, TPM and FPKM",
                    "description": "__SPECIES NAME__ RNA-Seq read counts, TPM (Transcript Per Million) and FPKM (Fragments Per Kilobase of transcript per Million mapped reads).",
                    "keywords": [
                        "RNA-Seq"
                    ],
                    "url": "https://bgee.org/species/__TAXID__#proc-values-rna-seq"
                }
            ]
        },
        {
            "@type": "WebPage",
            "url": "https://bgee.org/species/__TAXID__",
            "name": "Species: __SPECIES NAME____COMMON NAME2__"
        }
    ]
}';

    my $speciesdbh = $bgee->prepare('SELECT DISTINCT s.speciesId, s.genus, s.species, s.speciesCommonName, s.dataSourceId, d.baseUrl FROM species s, dataSource d WHERE s.dataSourceId=d.dataSourceId ORDER BY s.speciesId');
    $speciesdbh->execute()  or die $speciesdbh->errstr;
    my @species_json;
    while ( my ($speciesId, $genus, $species, $speciesCommonName, $dataSourceId, $baseUrl) = $speciesdbh->fetchrow_array() ){
        my $temp = $template;
        $temp =~ s{__TAXID__}{$speciesId}g;
        $temp =~ s{__SPECIES NAME__}{$genus $species}g;
        $temp =~ s{__SPECIES_NAME__}{${genus}_$species}g;
        if ( $speciesCommonName ne '' ){
            $temp =~ s{__COMMON NAME1__}{,\n                "$speciesCommonName"}g;
            $temp =~ s{__COMMON NAME2__}{ ($speciesCommonName)}g;
        }
        else {
            $temp =~ s{__COMMON NAME1__}{}g;
            $temp =~ s{__COMMON NAME2__}{}g;
        }
        if ( $baseUrl =~ /ensembl\.org/ ){
            $temp =~ s{__DBSRC_SPECIES_LINK__}{,\n        "$baseUrl${genus}_$species/"}g;
        }
        else {
            $temp =~ s{__DBSRC_SPECIES_LINK__}{}g;
        }
        push @species_json, $temp;
    }
    $speciesdbh->finish;

    return \@species_json;
}

sub get_schema_genes {
    my ($db) = @_;

    #TODO ALL data types, db fields: taxid, species name, ensembl/refseq links, ...
    #TODO Loop over alternateName, separated by *,*
    #TODO Loop over sameAs links, separated by *,*
    my $template = '{
    "@context": "https://schema.org/",
    "@type": "Gene",
    "@id": "https://bgee.org/gene/__GENEID__",
    "http://purl.org/dc/terms/conformsTo": {
        "@id": "https://bioschemas.org/profiles/Gene/1.0-RELEASE",
        "@type": "CreativeWork"
    },
    "description": "__GENEDESC__",
    "alternateName": [
        "wu:fb58g10",
        "apoc1l",
        "fb58g10"
    ],
    "identifier": "__GENEID__",
    "name": "__GENENAME__",
    "subjectOf": {
        "@type": "WebPage",
        "url": "https://bgee.org/bgee'.$bgee_db_version.'/gene/__GENEID__",
        "name": "Gene: __GENENAME__ - __GENEID__ - __SPECIES NAME__ (__COMMON NAME__)"
    },
    "taxonomicRange": {
        "@type": "Taxon",
        "@id": "https://bgee.org/bgee'.$bgee_db_version.'/species/__TAXID__",
        "name": "__SPECIES NAME__ (__COMMON NAME__)",
        "identifier": __TAXID__,
        "sameAs": "http://purl.obolibrary.org/obo/NCBITaxon___TAXID__"
    },
    "sameAs": [
        "__DBSRC_SPECIES_LINK__/Gene/Summary?g=__GENEID__",
        "https://zfin.org/ZDB-GENE-030131-1074",
        "https://www.ebi.ac.uk/ena/data/view/BX004983"
    ]
}';

    #TODO 7955                                                      <-> __TAXID__
    #TODO Danio rerio                                               <-> __SPECIES NAME__
    #TODO Danio_rerio                                               <-> __SPECIES_NAME__
    #TODO zebrafish                                                 <-> __COMMON NAME__   if no common name, do not display *"",* and * ()*
    #TODO ENSDARG00000092170                                        <-> __GENEID__
    #TODO apolipoprotein C-I [Source:ZFIN;Acc:ZDB-GENE-030131-1074] <-> __GENEDESC__
    #TODO apoc1                                                     <-> __GENENAME__
    #TODO      <-> __DBSRC_SPECIES_LINK__

    return;
}

sub get_schema_gene_homologs {
    my ($db) = @_;

    #TODO ALL data types, db fields: taxid, species name, ensembl/refseq links, ...
    #TODO Loop over taxon names with orthologs, separated by *,*
    my $template = '[{
        "@context": "https://schema.org/",
        "@type": "https://schema.org/Taxon",
        "@id": "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=__TAXIDANCESTOR__",
        "http://purl.org/dc/terms/conformsTo": {
            "@id": "https://bioschemas.org/profiles/Taxon/0.6-RELEASE",
            "@type": "CreativeWork"
        },
        "identifier": __TAXIDANCESTOR__,
        "name": "__SPECIESNAMEANCESTOR__"
    }
]';

    #TODO 186625        <-> __TAXIDANCESTOR__
    #TODO Clupeocephala <-> __SPECIESNAMEANCESTOR__

    return;
}

sub get_schema_gene_expression {
    my ($db) = @_;

    #TODO ALL data types, db fields: taxid, species name, ensembl/refseq links, ...
    #TODO Loop over expressedIn, separated by *,*
    my $template = '[
    {
        "@context": "https://schema.org/",
        "@type": "Gene",
        "@id": "https://bgee.org/gene/__GENEID__",
        "expressedIn": {
            "@type": "AnatomicalStructure",
            "@id": "__EXTIDURL__",
            "identifier": "__EXTID__",
            "name": "__EXTNAME__"
        }
    }
]';

    #TODO ENSDARG00000092170                            <-> __GENEID__
    #TODO http://purl.obolibrary.org/obo/UBERON_0002107 <-> __EXTIDURL__
    #TODO UBERON:0002107                                <-> __EXTID__
    #TODO liver                                         <-> __EXTNAME__

    return;
}

