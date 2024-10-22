#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use List::MoreUtils qw{uniq any};
use LWP::UserAgent;
use JSON;
use Try::Tiny;

# Define arguments & their default value
my ($ensembl_species_info_file) = ('');
my ($verbose)                  = (0);
my %opts = (
            'ensembl_species_info_file=s'     => \$ensembl_species_info_file
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $ensembl_species_info_file eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -_species_file=inputFile.tsv  -ensembl_species_info_file=outputFile.tsv
\t-ensembl_species_info_file       final file containing info for all species to insert in Bgee
\n";
    exit 1;
}

############# Functions #############
sub get_json {
    my ($user_agent, $server, $ext) = @_;
    my $response = $user_agent->get($server.$ext, 'Content-Type' => 'application/json');
    die "Error: ", $response->status_line unless $response->is_success;
    return decode_json($response->decoded_content);
}

sub write_ensembl_species_info {
    my ($response, $user_agent, $ensembl_server, $fh, $datasource_id) = @_;
    SPECIES:
    foreach my $species (@{$response->{species}}) {
        # have to manage different specific cases
        # 1. species corresponding to a strain (e.g. Mus musculus)
        # we do not want to retrieve them. To do so we do 2 checks. The first one check that the name of 
        # the species ends with the name of the strain and the second is that no taxonomy info can be
        # retrieved using the species name.
        # 2. species name with two _ (e.g Canis lupus familiaris)
        # In that case we check if it is the genus or species name that contain a space.
        my $name = $species->{name};
        my $speciesId = $species->{taxon_id};

        # if there are two '_' in the species name. We have to understand why (strain or names with 
        # multiple spaces)
        my ($genus_string, $species_string) = ('', '');
        my @name_parts = split("_", $name);
        if (scalar @name_parts < 2 || scalar @name_parts > 3) {
            # warn "$name is composed by less than 2 words or more than 3 words. The current script can not manage such cases. This species will not be inserted in Bgee\n";
        } elsif (scalar @name_parts == 3) {
            my $strain = $species->{strain};
            if (defined $strain) {
                $strain =~ s/\///g;
            }
            # it could be a strain specific species
            if (defined $strain && lc($name_parts[2]) eq lc($strain)) {
                next SPECIES;
            } else {
                # first we check that it is not heterocephalus_glaber_female. It is the only exemple
                # currently in Ensembl where male and female assembly are available AND the sex is part
                # of the name of the species....
                if ($name eq 'heterocephalus_glaber_female') {
                    $genus_string = $name_parts[0];
                    $species_string = $name_parts[1];
                } else {
                    #try to retrieve taxonomy info if an error accors it means it was a strain related species
                    # and we do not keep it
                    my $response = $user_agent->get($ensembl_server.'/taxonomy/classification/'.$name, 'Content-Type' => 'application/json');
                    if (! $response->is_success) {
                        warn "the species $name with speciesId $speciesId was detected as a strain/sex version of a species and is not going to be inserted in Bgee\n";
                        next SPECIES;
                    }
                    my $taxon_info = decode_json($response->decoded_content);

                    # if a taxonomy was retrieved then it means the presence of the two "_" in the taxon
                    # name was
                    # expected. We then need to know if it is the genus or species name that contains a space
                    # to do so we check the taxonomy info and check if the genus name contains a space
                    if (index($taxon_info->[0]->{name}, " ") eq -1 ) {
                        # if the genus name contains a space then it means the species name contains a space
                        $genus_string = $name_parts[0];
                        $species_string = $name_parts[1].' '.$name_parts[2];
                    } else {
                        $genus_string = $name_parts[0].' '.$name_parts[1];
                        $species_string = $name_parts[2];
                    }

                }
            }

        } else {
            $genus_string = $name_parts[0];
            $species_string = $name_parts[1];
        }
        my $speciesCommonName = $species->{common_name};
        my $genomeSpeciesId = $species->{taxon_id};
        my $fakeGeneIdPrefix = '';
        my $genus_with_underscore = ($genus_string =~ s/ /_/r);
        my $species_with_underscore = ($species_string =~ s/ /_/r);
        my $genomeVersion = $species->{assembly};
        my $keywords = join('|', @{$species->{aliases}});

        my $genomeFilePath = $genus_with_underscore."_".$species_with_underscore."/".
            ucfirst($genus_with_underscore)."_".$species_with_underscore.".".$genomeVersion;

        print $fh "$speciesId\t$genus_string\t$species_string\t$speciesCommonName\t\t$genomeFilePath\t$genomeVersion\t$datasource_id\t$speciesId\t\t$keywords\t\n";
    }
}

############# Main #############

my $ua = LWP::UserAgent->new;
$ua->agent("BgeeSpeciesGenerator/0.1");

my $ensembl_dataSource_id = 2;
my $ensembl_metazoa_datasource_id = 24;
 
# first retrieve info for all ensembl species
my $ensembl_server = 'https://rest.ensembl.org';
my $ensembl_species = '/info/species';




open(my $fh, '>', $ensembl_species_info_file) or die "Could not open file $ensembl_species_info_file $!";
print $fh "speciesId\tgenus\tspecies\tspeciesCommonName\tdisplayOrder\tgenomeFilePath\tgenomeVersion\tdataSourceId\tgenomeSpeciesId\tfakeGeneIdPrefix\tkeywords\tcomment\n";

# get the species info for ensembl
my $response = get_json($ua, $ensembl_server, $ensembl_species);
write_ensembl_species_info($response, $ua, $ensembl_server, $fh, $ensembl_dataSource_id);

# get the species info for ensembl metazoa
$response = get_json($ua, $ensembl_server, $ensembl_species.'?division=ensemblMetazoa');
write_ensembl_species_info($response, $ua, $ensembl_server, $fh, $ensembl_metazoa_datasource_id);

close $fh;