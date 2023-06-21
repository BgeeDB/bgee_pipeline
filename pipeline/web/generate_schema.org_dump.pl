#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use File::Slurp;

use Cpanel::JSON::XS;
use LWP::Simple;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;


#### constants ####
my $dateModified  = '2021-07-01';
my $bgeeCitation  = 'https://doi.org/10.1093/nar/gkaa793';
my $bgeeLicense   = 'https://creativecommons.org/publicdomain/zero/1.0/';
my $bioschGene    = 'https://bioschemas.org/profiles/Gene/1.0-RELEASE';
my $bioschDataset = 'https://bioschemas.org/profiles/Dataset/1.0-RELEASE';
my $bioschTaxon   = 'https://bioschemas.org/profiles/Taxon/0.6-RELEASE';


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
my @schema_species = '['.join(',', @{ get_schema_species($bgee) }).']';

# gene pages
my @schema_gene    = '['.join(',', @{ get_schema_genes($bgee) }).']';

$bgee->disconnect;


print '[', join(',', map { s|</?script.*||g; $_ } $schema_default, @schema_species, @schema_gene), ']';

exit 0;



sub get_schema_default {
    return '{
    "@context": "https://schema.org/",
    "@id": "https://www.bgee.org/",
    "@graph": [
        {
            "@type": "Organization",
            "@id": "https://www.bgee.org/",
            "name": "Bgee - Bring Gene Expression Expertise",
            "url": "https://www.bgee.org/",
            "description": "The aim of Bgee is to help biologists to use and understand gene expression",
            "logo": "https://www.bgee.org/img/logo/bgee13_hp_logo.png",
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
            "@id": "https://www.bgee.org/",
            "http://purl.org/dc/terms/conformsTo": {
                "@id": "'.$bioschDataset.'",
                "@type": "CreativeWork"
            },
            "url": "https://www.bgee.org/",
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
                "@id": "https://www.bgee.org/"
            },
            "license": "'.$bgeeLicense.'",
            "version": "'.$bgee_version.'"
        }
    ]
}';
}

sub get_schema_species {
    my ($bgee) = @_;

    my $affy_template   = '                {
                    "@type": "Dataset",
                    "http://purl.org/dc/terms/conformsTo": {
                        "@id": "'.$bioschDataset.'",
                        "@type": "CreativeWork"
                    },
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://www.bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "'.$bgeeLicense.'",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://www.bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/processed_expr_values/affymetrix/__SPECIES_NAME__/__SPECIES_NAME___Affymetrix_experiments_chips.tar.gz"
                        }
                    ],
                    "name": "__SPECIES NAME__ Affymetrix experiments chips",
                    "keywords": [
                        "Affymetrix"
                    ],
                    "description": "Affymetrix experiments/chips annotations and metadata.",
                    "url": "https://www.bgee.org/species/__TAXID__#proc-values-affymetrix"
                },
                {
                    "@type": "Dataset",
                    "http://purl.org/dc/terms/conformsTo": {
                        "@id": "'.$bioschDataset.'",
                        "@type": "CreativeWork"
                    },
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://www.bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "'.$bgeeLicense.'",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://www.bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/processed_expr_values/affymetrix/__SPECIES_NAME__/__SPECIES_NAME___Affymetrix_probesets.tar.gz"
                        }
                    ],
                    "name": "__SPECIES NAME__ Affymetrix probesets",
                    "description": "__SPECIES NAME__ Affymetrix probesets, data (signal intensities).",
                    "url": "https://www.bgee.org/species/__TAXID__#proc-values-affymetrix"
                }';
    my $rnaseq_template = '                {
                    "@type": "Dataset",
                    "http://purl.org/dc/terms/conformsTo": {
                        "@id": "'.$bioschDataset.'",
                        "@type": "CreativeWork"
                    },
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://www.bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "'.$bgeeLicense.'",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://www.bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/processed_expr_values/rna_seq/__SPECIES_NAME__/__SPECIES_NAME___RNA-Seq_experiments_libraries.tar.gz"
                        }
                    ],
                    "name": "__SPECIES NAME__ RNA-Seq experiment libraries",
                    "keywords": [
                        "RNA-Seq"
                    ],
                    "description": "__SPECIES NAME__ RNA-Seq experiments/libraries annotations and metadata.",
                    "url": "https://www.bgee.org/species/__TAXID__#proc-values-rna-seq"
                },
                {
                    "@type": "Dataset",
                    "http://purl.org/dc/terms/conformsTo": {
                        "@id": "'.$bioschDataset.'",
                        "@type": "CreativeWork"
                    },
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://www.bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "'.$bgeeLicense.'",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://www.bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/processed_expr_values/rna_seq/__SPECIES_NAME__/__SPECIES_NAME___RNA-Seq_read_counts_TPM_FPKM.tar.gz"
                        }
                    ],
                    "name": "__SPECIES NAME__ RNA-Seq read counts, TPM and FPKM",
                    "description": "__SPECIES NAME__ RNA-Seq read counts, TPM (Transcript Per Million) and FPKM (Fragments Per Kilobase of transcript per Million mapped reads).",
                    "keywords": [
                        "RNA-Seq"
                    ],
                    "url": "https://www.bgee.org/species/__TAXID__#proc-values-rna-seq"
                }';
    my $fulll_template  = '                {
                    "@type": "Dataset",
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://www.bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "'.$bgeeLicense.'",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://www.bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/processed_expr_values/sc_full_length/__SPECIES_NAME__/__SPECIES_NAME___Full-Length_SC_RNA-Seq_experiments_libraries.tar.gz"
                        }
                    ],
                    "name": "__SPECIES NAME__ full-length Single cell RNA-Seq experiment libraries",
                    "description": "__SPECIES NAME__ full-length Single cell RNA-Seq experiments/ libraries annotations and metadata.",
                    "keywords": [
                        "Single cell full length RNA-Seq",
                        "Single cell RNA-Seq"
                    ],
                    "url": "https://www.bgee.org/species/__TAXID__#proc-values-fl-scrna-seq"
                },
                {
                    "@type": "Dataset",
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://www.bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "'.$bgeeLicense.'",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://www.bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/processed_expr_values/sc_full_length/__SPECIES_NAME__/__SPECIES_NAME___Full-Length_SC_RNA-Seq_read_counts_TPM_FPKM.tar.gz"
                        }
                    ],
                    "name": "__SPECIES NAME__ Full-Length Single Cell RNA-Seq read counts, TPM and FPKM",
                    "description": "__SPECIES NAME__ Full-Length Single Cell RNA-Seq read counts, TPM (Transcript Per Million) and FPKM (Fragments Per Kilobase of transcript per Million mapped reads).",
                    "keywords": [
                        "Single cell full length RNA-Seq",
                        "Single cell RNA-Seq"
                    ],
                    "url": "https://www.bgee.org/species/__TAXID__#proc-values-fl-scrna-seq"
                }';
#TODO    my $target_template = '';

    my $template = '{
    "@context": "https://schema.org/",
    "@id": "https://www.bgee.org/species/__TAXID__",
    "@type": "Taxon",
    "http://purl.org/dc/terms/conformsTo": {
        "@id": "'.$bioschTaxon.'",
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
            "citation": "'.$bgeeCitation.'",
            "description": "__SPECIES NAME__ calls of presence/absence of expression. Each call corresponds to a unique combination of a gene, an anatomical entity, a life stage, a sex, and a strain, with reported presence or absence of expression.",
            "includedInDataCatalog": {
                "@id": "https://www.bgee.org",
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
                "url": "https://www.bgee.org/",
                "name": "The Bgee Team"
            },
            "license": "'.$bgeeLicense.'",
            "name": "__SPECIES NAME__ gene expression calls",
            "url": "https://www.bgee.org/species/__TAXID__#expr-calls",
            "version": "'.$bgee_version.'",
            "hasPart": [
                {
                    "@type": "Dataset",
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://www.bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "'.$bgeeLicense.'",
                    "name": "__SPECIES NAME__ gene expression simple",
                    "description": "Anatomical entities only, file without advanced columns.",
                    "url": "https://www.bgee.org/species/__TAXID__#expr-calls-anat-simple",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://www.bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/calls/expr_calls/__SPECIES_NAME___expr_simple.tsv.gz"
                        }
                    ]
                },
                {
                    "@type": "Dataset",
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://www.bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "'.$bgeeLicense.'",
                    "name": "__SPECIES NAME__ gene expression advanced",
                    "description": "Anatomical entities only, file with advanced columns.",
                    "url": "https://www.bgee.org/species/__TAXID__#expr-calls-anat-advanced",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://www.bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/calls/expr_calls/__SPECIES_NAME___expr_advanced.tsv.gz"
                        }
                    ]
                },
                {
                    "@type": "Dataset",
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://www.bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "'.$bgeeLicense.'",
                    "name": "__SPECIES NAME__ gene expression simple with all conditions",
                    "description": "Anatomical entities, developmental stages, sexes and strains. File without advanced columns.",
                    "url": "https://www.bgee.org/species/__TAXID__#expr-calls-cond-simple",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://www.bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/calls/expr_calls/__SPECIES_NAME___expr_simple_all_conditions.tsv.gz"
                        }
                    ]
                },
                {
                    "@type": "Dataset",
                    "dateModified": "'.$dateModified.'",
                    "creator": {
                        "@type": "Organization",
                        "url": "https://www.bgee.org/",
                        "name": "The Bgee Team"
                    },
                    "license": "'.$bgeeLicense.'",
                    "name": "__SPECIES NAME__ gene expression advanced with all conditions",
                    "description": "Anatomical entities, developmental stages, sexes and strains. File with advanced columns.",
                    "url": "https://www.bgee.org/species/__TAXID__#expr-calls-cond-advanced",
                    "distribution": [
                        {
                            "@type": "DataDownload",
                            "encodingFormat": "TSV",
                            "contentUrl": "https://www.bgee.org/ftp/bgee_v'.$bgee_db_version.'/download/calls/expr_calls/__SPECIES_NAME___expr_advanced_all_conditions.tsv.gz"
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
                "url": "https://www.bgee.org/",
                "name": "The Bgee Team"
            },
            "citation": "'.$bgeeCitation.'",
            "description": "Annotations and experiment information (e.g., annotations to anatomy and development, quality scores used in QCs, library information), and processed expression values (e.g., read counts, TPM and FPKM values) for __SPECIES NAME__",
            "includedInDataCatalog": {
                "@id": "https://www.bgee.org/bgee'.$bgee_db_version.'",
                "@type": "DataCatalog",
                "name": "Bgee"
            },
            "keywords": [
                "annotations",
                "experiment information",
                "processed expression values",
                "__SPECIES NAME__"__COMMON NAME1__
            ],
            "license": "'.$bgeeLicense.'",
            "name": "__SPECIES NAME__ processed expression values",
            "url": "https://www.bgee.org/species/__TAXID__#proc-values",
            "version": "'.$bgee_version.'",
            "hasPart": [
__DATATYPES__
            ]
        },
        {
            "@type": "WebPage",
            "url": "https://www.bgee.org/species/__TAXID__",
            "name": "Species: __SPECIES NAME____COMMON NAME2__"
        }
    ]
}';

    my $speciesdbh = $bgee->prepare('SELECT DISTINCT s.speciesId, s.genus, s.species, s.speciesCommonName, s.dataSourceId, d.baseUrl, (SELECT GROUP_CONCAT(DISTINCT dataType) FROM dataSourceToSpecies WHERE speciesId=s.speciesId AND infoType="data" AND dataType !="est" AND dataType !="in situ") AS datatypes FROM species s, dataSource d, dataSourceToSpecies t WHERE s.dataSourceId=d.dataSourceId AND s.speciesId=t.speciesId ORDER BY s.speciesId');
    $speciesdbh->execute()  or die $speciesdbh->errstr;
    my @species_json;
    while ( my ($speciesId, $genus, $species, $speciesCommonName, $dataSourceId, $baseUrl, $datatypes) = $speciesdbh->fetchrow_array() ){
        my $temp = $template;

        # Datatypes
        #Loop over data types with *Processed expression value files*, separated by *,*
        my @datatypes = split(',', $datatypes);
        my @dt_temp;
        # can be: affymetrix,rna-seq,full-length single-cell RNA-Seq
        for my $dt ( sort @datatypes ){
            if ( $dt eq 'rna-seq' ){
                push @dt_temp, $rnaseq_template;
            }
            elsif ( $dt eq 'affymetrix' ){
                push @dt_temp, $affy_template;
            }
            elsif ( $dt eq 'full-length single-cell RNA-Seq' ){
                push @dt_temp, $fulll_template;
            }
            #TODO target based single-cell RNA-Seq!
        }
        warn "No datatypes for $genus $species\n"  if ( scalar @dt_temp == 0 );
        my $dt_temp = join(",\n", @dt_temp);
        $temp =~ s{__DATATYPES__}{$dt_temp};

        # Taxon info
        $temp =~ s{__TAXID__}{$speciesId}g;
        $temp =~ s{__SPECIES NAME__}{$genus $species}g;
        $species =~ s{ }{_}g; # for Canis lupus familiaris
        $temp =~ s{__SPECIES_NAME__}{${genus}_$species}g;
        if ( $speciesCommonName && $speciesCommonName ne '' ){
            $temp =~ s{__COMMON NAME1__}{,\n                "$speciesCommonName"}g;
            $temp =~ s{__COMMON NAME2__}{ ($speciesCommonName)}g;
        }
        else {
            $temp =~ s{__COMMON NAME1__}{}g;
            $temp =~ s{__COMMON NAME2__}{}g;
        }
        #FIXME this is wrong on the Bgee web site for non www.ensembl.org: do not use ensembl metazoa link, and provide a link for refseq based genomes
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
    my ($bgee) = @_;

    my $template = '{
    "@context": "https://schema.org/",
    "@type": "Gene",
    "@id": "https://www.bgee.org/gene/__GENEID__",
    "http://purl.org/dc/terms/conformsTo": {
        "@id": "'.$bioschGene.'",
        "@type": "CreativeWork"
    },
    "description": "__GENEDESC__",__ALTNAME__
    "identifier": "__GENEID__",
    "name": "__GENENAME__",
    "subjectOf": {
        "@type": "WebPage",
        "url": "https://www.bgee.org/bgee'.$bgee_db_version.'/gene/__GENEID__",
        "name": "Gene: __GENENAME__ - __GENEID__ - __SPECIES NAME____COMMON NAME2__"
    },
    "taxonomicRange": {
        "@type": "Taxon",
        "@id": "https://www.bgee.org/bgee'.$bgee_db_version.'/species/__TAXID__",
        "name": "__SPECIES NAME____COMMON NAME2__",
        "identifier": __TAXID__,
        "sameAs": "http://purl.obolibrary.org/obo/NCBITaxon___TAXID__"
    },
    "sameAs": [
__SAMEAS__
    ]
}';

    my $genesdbh = $bgee->prepare('SELECT DISTINCT g.bgeeGeneId, g.geneId, g.geneName, g.geneDescription, t.speciesId, t.genus, t.species, t.speciesCommonName,d.baseUrl FROM gene g, species t, dataSource d WHERE t.dataSourceId=d.dataSourceId AND g.speciesId=t.speciesId AND g.bgeeGeneId IN (1611228, 129767, 257884, 62875, 873906) ORDER BY g.geneId');
    $genesdbh->execute()  or die $genesdbh->errstr;
    my $genesSyndbh  = $bgee->prepare('SELECT GROUP_CONCAT(DISTINCT geneNameSynonym) AS syn FROM geneNameSynonym WHERE bgeeGeneId=?');
    my $genesXrefdbh = $bgee->prepare('SELECT GROUP_CONCAT(DISTINCT REPLACE(REPLACE(REPLACE(d.XRefUrl, "[xref_id]" ,x.XRefId), "[gene_id]", x.XRefId), "[species_ensembl_link]", "__SPECIES_NAME__")) AS geneXrefLink FROM geneXRef x, dataSource d WHERE x.bgeeGeneId = ? AND d.dataSourceId=x.dataSourceId AND x.dataSourceId!=101');
    my @genes_json;
    while ( my ($bgeeGeneId, $geneId, $geneName, $geneDesc, $taxId, $genus, $species, $speciesCommonName, $baseUrl) = $genesdbh->fetchrow_array() ){
        my $temp = $template;

        # Gene info
        $temp =~ s{__GENEID__}{$geneId}g;
        $temp =~ s{__GENEDESC__}{$geneDesc}g;
        if ( $geneName =~ /\\/ ){
            $geneName =~ s{\\}{\\\\}g;
        }
        $temp =~ s{__GENENAME__}{$geneName}g;

        # Alt gene syn
        $genesSyndbh->execute($bgeeGeneId)  or die $genesSyndbh->errstr;
        my ($geneSyn) = $genesSyndbh->fetchrow_array();
        if ( $geneSyn && $geneSyn ne '' ){
            my $syns = "\n    \"alternateName\": [\n";
            $syns .= join(",\n", map { "        \"$_\"" } split(/,/, $geneSyn));
            $syns .= "\n    ],";
            $temp =~ s{__ALTNAME__}{$syns}g;
        }
        else {
            $temp =~ s{__ALTNAME__}{}g;
        }

        # Taxon info
        $temp =~ s{__TAXID__}{$taxId}g;
        $temp =~ s{__SPECIES NAME__}{$genus $species}g;
        if ( $speciesCommonName && $speciesCommonName ne '' ){
            $temp =~ s{__COMMON NAME2__}{ ($speciesCommonName)}g;
        }
        else {
            $temp =~ s{__COMMON NAME2__}{}g;
        }

        # Same as
        $genesXrefdbh->execute($bgeeGeneId)  or die $genesXrefdbh->errstr;
        my $sameAs = '';
        if ( $baseUrl =~ /ensembl\.org/ ){
            #NOTE direct links to Ensembl protein/transcript
            if ( $geneId =~ /^FBtr/ ){
                $sameAs .= "        \"$baseUrl${genus}_$species/Transcript/Summary?g=$geneId\",\n";
            }
            elsif ( $geneId =~ /^FBpp/ ){
                $sameAs .= "        \"$baseUrl${genus}_$species/Transcript/ProteinSummary?g=$geneId\",\n";
            }
            else {
                $sameAs .= "        \"$baseUrl${genus}_$species/Gene/Summary?g=$geneId\",\n";
            }
        }
        $sameAs .= join(",\n", map { "        \"$_\"" } split(/,/, $genesXrefdbh->fetchrow_array()));
        $species =~ s{ }{_}g; # for Canis lupus familiaris
        $sameAs =~ s{__SPECIES_NAME__}{${genus}_$species}g;
        #NOTE direct links to Ensembl protein/transcript
        if ( $sameAs =~ /ensembl\.org.+g=FBtr\d+/ ){
            $sameAs =~ s{Gene/Summary\?g=FBtr}{Transcript/Summary?g=FBtr}g;
        }
        if ( $sameAs =~ /ensembl\.org.+g=FBpp\d+/ ){
            $sameAs =~ s{Gene/Summary\?g=FBpp}{Transcript/ProteinSummary?g=FBpp}g;
        }
        $temp =~ s{__SAMEAS__}{$sameAs}g;

        push @genes_json, join(',', '['.$temp,
                                    '['.join(',', @{get_schema_gene_homologs($bgee, $bgeeGeneId)})."\n".']',
                                    '['.join(',', @{get_schema_gene_expression($bgee, $geneId, $bgeeGeneId, $taxId)})."\n".']]',
                              );
    }
    $genesSyndbh->finish;
    $genesXrefdbh->finish;
    $genesdbh->finish;

    return \@genes_json;
}

sub get_schema_gene_homologs {
    my ($bgee, $bgeeGeneId) = @_;

    my $template = '{
        "@context": "https://schema.org/",
        "@type": "https://schema.org/Taxon",
        "@id": "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=__TAXIDANCESTOR__",
        "http://purl.org/dc/terms/conformsTo": {
            "@id": "'.$bioschTaxon.'",
            "@type": "CreativeWork"
        },
        "identifier": __TAXIDANCESTOR__,
        "name": "__SPECIESNAMEANCESTOR__"__COMMON NAME1__
    }';

    my $genesHomologdbh = $bgee->prepare('(SELECT DISTINCT geneOrthologs.taxonId, taxon.taxonScientificName, taxon.taxonCommonName, taxon.taxonLevel FROM geneOrthologs  INNER JOIN taxon ON geneOrthologs.taxonId = taxon.taxonId WHERE geneOrthologs.bgeeGeneId IN (?))  UNION  (SELECT DISTINCT geneOrthologs.taxonId, taxon.taxonScientificName, taxon.taxonCommonName, taxon.taxonLevel FROM geneOrthologs INNER JOIN taxon ON geneOrthologs.taxonId = taxon.taxonId WHERE geneOrthologs.targetGeneId IN (?)) UNION (SELECT DISTINCT geneParalogs.taxonId, taxon.taxonScientificName, taxon.taxonCommonName, taxon.taxonLevel FROM geneParalogs  INNER JOIN taxon ON geneParalogs.taxonId = taxon.taxonId WHERE geneParalogs.bgeeGeneId IN (?)) UNION (SELECT DISTINCT geneParalogs.taxonId, taxon.taxonScientificName, taxon.taxonCommonName, taxon.taxonLevel FROM geneParalogs INNER JOIN taxon ON geneParalogs.taxonId = taxon.taxonId WHERE geneParalogs.targetGeneId IN (?)) ORDER BY taxonLevel DESC');
    $genesHomologdbh->execute($bgeeGeneId, $bgeeGeneId, $bgeeGeneId, $bgeeGeneId)  or die $genesHomologdbh->errstr;
    my @homologs_json;
    while ( my ($taxonId, $taxonScientificName, $taxonCommonName, $taxonLevel) = $genesHomologdbh->fetchrow_array() ){
        my $temp = $template;

        $temp =~ s{__TAXIDANCESTOR__}{$taxonId}g;
        $temp =~ s{__SPECIESNAMEANCESTOR__}{$taxonScientificName}g;
        if ( $taxonCommonName && $taxonCommonName ne '' ){
            $temp =~ s{__COMMON NAME1__}{,\n        "alternateName": "$taxonCommonName"}g;
        }
        else {
            $temp =~ s{__COMMON NAME1__}{}g;
        }

        push @homologs_json, $temp;
    }
    $genesHomologdbh->finish;

    return \@homologs_json;
}

sub get_schema_gene_expression {
    my ($bgee, $geneId, $bgeeGeneId, $taxId) = @_;

    my $template = '{
        "@context": "https://schema.org/",
        "@type": "Gene",
        "@id": "https://www.bgee.org/gene/__GENEID__",
        "expressedIn": {
            "@type": "AnatomicalStructure",
            "@id": "http://purl.obolibrary.org/obo/__EXTIDURL__",
            "identifier": "__EXTID__",
            "name": "__EXTNAME__"
        }
    }';

    #NOTE Use the HTTP API to catch only expressed in, without the parent terms with the same or lower score
    my $URL = 'https://www.bgee.org/api/?display_type=json&page=gene&action=expression&gene_id='.$geneId.'&species_id='.$taxId.'&cond_param=anat_entity&cond_param=cell_type&data_type=all';
    my $content = get($URL);
    my @expressed_json;
    if ( defined $content ){
        my $json_ref = decode_json($content);
        for my $call ( @{ $json_ref->{'data'}->{'calls'} } ){
            my $anatEntityId   = $call->{'condition'}->{'anatEntity'}->{'id'};
            my $anatEntityName = $call->{'condition'}->{'anatEntity'}->{'name'};

            my $temp = $template;
            $temp =~ s{__GENEID__}{$geneId}g;
            $temp =~ s{__EXTNAME__}{$anatEntityName}g;
            $temp =~ s{__EXTID__}{$anatEntityId}g;
            $anatEntityId =~ s{:}{_};
            $temp =~ s{__EXTIDURL__}{$anatEntityId}g;

            push @expressed_json, $temp;
        }
    }
    else {
        warn "Cannot get the API answer [$URL]\n";
    }

    return \@expressed_json;
}

