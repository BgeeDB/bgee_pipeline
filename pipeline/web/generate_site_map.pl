#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use File::Slurp;
use List::MoreUtils qw(uniq);
use Sort::Naturally;

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
my @basic_cgi      = ('about', 'anat_similarities', 'collaborations', 'doc', 'download', 'expression_comparison', 'gene', 'privacy_policy', 'source', 'sparql', 'top_anat');
for my $baseUrlName ( sort @basic_cgi ){
    push @static_pages, "<loc>$homepage/?page=$baseUrlName</loc><priority>0.7</priority>";
}
my @basic_doc      = ('access', 'call_files', 'data_sets', 'faq', 'top_anat');
for my $baseUrlName ( sort @basic_doc ){
    push @static_pages, "<loc>$homepage/?page=doc&amp;action=$baseUrlName</loc><priority>0.7</priority>";
}
my @basic_download = ('expr_calls', 'mysql_dumps', 'proc_values');
for my $baseUrlName ( sort @basic_download ){
    push @static_pages, "<loc>$homepage/?page=download&amp;action=$baseUrlName</loc><priority>0.7</priority>";
}
my @basic_resources = ('annotations', 'ontologies', 'r_packages', 'source_code');
for my $baseUrlName ( sort @basic_resources ){
    push @static_pages, "<loc>$homepage/?page=resources&amp;action=$baseUrlName</loc><priority>0.7</priority>";
}

#NOTE FTP links are invalid because different namespace (not bgee.org)

write_file("$sitemap_main", $sitemap_header, join("\n", map { "<url>$_</url>" } @static_pages), $sitemap_footer);
push @index, $sitemap_main;


# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

$| = 1;


# gene pages
print "Write gene pages\n"  if ( $debug );
my $geneId = $bgee->prepare('SELECT DISTINCT geneId FROM gene ORDER BY geneId');
$geneId->execute()  or die $geneId->errstr;
my $count = 0;
my $split = 0;
my @pages;
while ( my @data = $geneId->fetchrow_array ){
    $count++;
    push @pages, "<loc>$homepage/?page=gene&amp;gene_id=$data[0]</loc>";
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

