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
my ($debug)          = (0);
my %opts = ('bgee=s'     => \$bgee_connector,   # Bgee connector string
            'debug'      => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\t-debug            printing the update/insert SQL queries, not executing them
\n";
    exit 1;
}


# global variables
#NOTE see https://www.sitemaps.org/ and https://en.wikipedia.org/wiki/Site_map for spec and information
#NOTE The produced **sitemap.xml** can be validated on different web sites, but has to be submitted and validated in Google Dashboard
my $url_limit = 50_000; # If more than that, has to split in several sitemap files, indexed in a sitemapindex file!
my $homepage  = 'https://bgee.org';

my $sitemap_idx  = 'sitemap.xml';
my $sitemap_main = 'sitemap_main.xml';
my @index;

my $sitemap_header = "<?xml version='1.0' encoding='UTF-8'?>\n<urlset xmlns='http://www.sitemaps.org/schemas/sitemap/0.9'>\n";
my $sitemap_footer = "\n</urlset>";

#TODO Add <lastmod> related to <loc> to help indexing?


# "static" pages
print "Write main/static pages\n"  if ( $debug );
my @static_pages;
push @static_pages, "<loc>$homepage/</loc><priority>0.8</priority>";
push @static_pages, "<loc>$homepage/sparql/</loc><priority>0.7</priority>";
my @basic_about     = ('', 'news', 'collaborations', 'publications', 'source', 'privacy-policy');
for my $baseUrlName ( sort @basic_about ){
    push @static_pages, "<loc>$homepage/about/$baseUrlName</loc><priority>0.7</priority>";
}
my @basic_support   = ('data-sets', 'top-anat', 'gene-expression-calls', 'processed-expression-values', 'rnaseq-processed-expression-values', 'singlecell-rnaseq-fulllength-processed-expression-values', 'affymetrix-processed-expression-values', 'faq');
for my $baseUrlName ( sort @basic_support ){
    push @static_pages, "<loc>$homepage/support/$baseUrlName</loc><priority>0.7</priority>";
}
my @basic_download  = ('gene-expression-calls', 'processed-expression-values', 'data-dumps');
for my $baseUrlName ( sort @basic_download ){
    push @static_pages, "<loc>$homepage/download/$baseUrlName</loc><priority>0.7</priority>";
}
my @basic_resources = ('r-packages', 'sparql', 'annotations', 'ontologies', 'source-code');
for my $baseUrlName ( sort @basic_resources ){
    push @static_pages, "<loc>$homepage/resources/$baseUrlName</loc><priority>0.7</priority>";
}
my @basic_analysis  = ('top-anat', 'expr-comparison');
for my $baseUrlName ( sort @basic_analysis ){
    push @static_pages, "<loc>$homepage/analysis/$baseUrlName</loc><priority>0.7</priority>";
}
my @basic_search    = ('genes', 'anatomical-homology', 'species');
for my $baseUrlName ( sort @basic_search ){
    push @static_pages, "<loc>$homepage/search/$baseUrlName</loc><priority>0.7</priority>";
}

#NOTE FTP links are invalid because different namespace (not bgee.org)

write_file("$sitemap_main", $sitemap_header, join("\n", map { "<url>$_</url>" } @static_pages), $sitemap_footer);
push @index, $sitemap_main;


# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

$| = 1;


# species pages
print "Write species pages\n"  if ( $debug );
my $speciesId = $bgee->prepare('SELECT DISTINCT speciesId FROM species ORDER BY speciesId');
$speciesId->execute()  or die $speciesId->errstr;
my $sp_sitemap = 'sitemap_species.xml';
my @species;
while ( my @data = $speciesId->fetchrow_array ){
    push @species, "<loc>$homepage/species/$data[0]</loc>";
}
write_file("$sp_sitemap", $sitemap_header, join("\n", map { "<url>$_</url>" } @species), $sitemap_footer);
push @index, $sp_sitemap;


# gene pages
print "Write gene pages\n"  if ( $debug );
my $geneId = $bgee->prepare('SELECT DISTINCT geneId, speciesId FROM gene ORDER BY geneId');
$geneId->execute()  or die $geneId->errstr;
my $count = 0;
my $split = 0;
my @pages;
while ( my @data = $geneId->fetchrow_array ){
    $count++;
    push @pages, "<loc>$homepage/gene/$data[0]/$data[1]</loc>";
    if ( $count > ($url_limit - 1000) ){
        $split++;
        my $sitemap = "sitemap_gene$split.xml";
        write_file("$sitemap", $sitemap_header, join("\n", map { "<url>$_</url>" } @pages), $sitemap_footer);
        push @index, $sitemap;
        $count = 0;
        @pages = ();
    }
}
$split++;
my $sitemap = "sitemap_gene$split.xml";
write_file("$sitemap", $sitemap_header, join("\n", map { "<url>$_</url>" } @pages), $sitemap_footer);
push @index, $sitemap;



# Write sitemap index
my $index = "<?xml version='1.0' encoding='UTF-8'?>\n<sitemapindex xmlns='http://www.sitemaps.org/schemas/sitemap/0.9'>\n";
for my $idx ( @index ){
    $index .= "<sitemap><loc>$homepage/$idx</loc></sitemap>\n";
}
$index .= '</sitemapindex>';

write_file("$sitemap_idx", $index);


$bgee->disconnect;
print "Done\n"  if ( $debug );

exit 0;

