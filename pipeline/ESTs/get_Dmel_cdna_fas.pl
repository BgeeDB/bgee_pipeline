#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use HTTP::Request;
use LWP::UserAgent;

my $xml = '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.7" >
	<Dataset name = "dmelanogaster_gene_ensembl" interface = "default" >
		<Filter name = "biotype" value = "protein_coding"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "cdna" />
	</Dataset>
</Query>';

my $path  = 'http://www.ensembl.org/biomart/martservice?';
#my $path = 'http://www.biomart.org/biomart/martservice?';


my $request  = HTTP::Request->new('POST', $path, HTTP::Headers->new(), 'query='.$xml."\n");
my $ua       = LWP::UserAgent->new;
my $all_data = '';
my $response = $ua->request($request,
               sub{
                   my($data, $res) = @_;
                   if ( $res->is_success ){
                       $all_data .= $data;
                   }
                   else {
                       #print ('Problems with the web server: '.$res->status_line);
                   }
                }, 1000);

if ( !$response->is_success || $all_data =~ /^query error/i ){
    my $message = $response->status_line;
    if ( $response->is_success ){
        $message = $all_data;
    }
    die('ERROR with the web server: ', $message, "\n", 'query=', $xml, "\n");
}

print $all_data;

exit 0;

