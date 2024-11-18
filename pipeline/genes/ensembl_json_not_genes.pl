#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use JSON::XS;

# To see the Ensembl JSON without the large `genes` part!

## Load the JSON
my $input_json_file = $ARGV[0]  or die "\n\tMissing input JSON file\n\n";
my $json_text = do {
    open(my $json_fh, '<:encoding(UTF-8)', $input_json_file)  or die("Can't open \"$input_json_file\": $!\n");
    local $/;
    <$json_fh>;
};
my $ensembl_json = decode_json($json_text);
#NOTE keep only the NOT genes section of the JSON!
map { delete( $ensembl_json->{$_} ) } grep { $_ eq 'genes' } keys %$ensembl_json;

print encode_json($ensembl_json);

exit 0;

