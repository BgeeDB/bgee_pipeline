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
my ($ensembl_species_info_file, $filter_species_ids) = ('', '');
my ($verbose)                  = (0);
my %opts = (
            'ensembl_species_info_file=s'     => \$ensembl_species_info_file,
            'filter_species_ids=s'            => \$filter_species_ids
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $ensembl_species_info_file eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0  -_species_file=inputFile.tsv  -ensembl_species_info_file=outputFile.tsv
\t-ensembl_species_info_file       final file containing info for all species to insert in Bgee
\tfilter_species_ids               comma separated list of species IDs to filter out;
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

sub create_filtered_species_hash {
    my ($filter_species_ids) = @_;
    my %filtered_species_hash;
    foreach my $species_id (split(',', $filter_species_ids)) {
        $filtered_species_hash{$species_id} = 1;
    }
    return %filtered_species_hash;
}

sub write_ensembl_species_info {
    my ($response, $user_agent, $ensembl_server, $fh, $datasource_id, $filtered_species_hash) = @_;

    SPECIES:
    foreach my $species (@{$response->{species}}) {
        my $speciesId = $species->{taxon_id};
        # check if the species is in the filtered species hash
        if (keys %$filtered_species_hash != 0 && !exists $filtered_species_hash->{$speciesId}) {
            print "species not part of filtered species: $species->{name} (taxon_id: $speciesId)\n";
            next SPECIES;
        }
        # have to manage different specific cases
        # 1. species corresponding to a strain (e.g. Mus musculus)
        # we do not want to retrieve them. To do so we do 2 checks. The first one check that the name of 
        # the species ends with the name of the strain and the second is that no taxonomy info can be
        # retrieved using the species name.
        # 2. species name with two _ (e.g Canis lupus familiaris)
        # In that case we check if it is the genus or species name that contain a space.
        my $name = $species->{name};

        # if there are two '_' in the species name. We have to understand why (strain or names with 
        # multiple spaces)
        my ($genus_string, $species_string) = ('', '');
        my @name_parts = split("_", $name);
        if (scalar @name_parts < 2 || scalar @name_parts > 3) {
            warn "$name is composed by less than 2 words or more than 3 words. The current script can not manage such cases. This species will not be inserted in Bgee\n";
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

        # For D. pseudoobscura, we have do not want to use the subspecies so we have to hardcode
        # a remapping to D. pseudoobscura... 
        if ($speciesId == 46245) {
            $speciesId = 7237;
            $genomeVersion = "UCI_Dpse_MV25";
            my $accession = "gca009870125v2rs";
            $genomeFilePath = $genus_with_underscore."_".$species_with_underscore."_".$accession."/".
            ucfirst($genus_with_underscore)."_".$species_with_underscore."_".$accession.".".$genomeVersion;

        }
        # For G. gorilla,  we do not want to use the subspecies ID for which come the reference genome
        # so we have to hardcode a remapping to G. gorilla. It will allow to keep the same taxon ID as
        # in previous releases of Bgee (15.2 and before)
        if ($speciesId == 9595) {
            $speciesId = 9593;
        }
        # For Gasterosteus aculeatus, we do not want to use the subspecies ID for which come the reference genome
        # so we have to hardcode a remapping to Gasterosteus aculeatus. It will allow to keep the same taxon ID as
        # in previous releases of Bgee (15.2 and before)
        if ($speciesId == 481459) {
            $speciesId = 69293;
        }
        
        # add an empty column corresponding to comments. This column will be manually filled.
        my$comments = '';
        print $fh "$speciesId\t$genus_string\t$species_string\t$speciesCommonName\t\t$genomeFilePath\t$genomeVersion\t$datasource_id\t$speciesId\t\t$keywords\t$comments\n";
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

my %filtered_species_hash = create_filtered_species_hash($filter_species_ids);

open(my $fh, '>', $ensembl_species_info_file) or die "Could not open file $ensembl_species_info_file $!";
print $fh "speciesId\tgenus\tspecies\tspeciesCommonName\tdisplayOrder\tgenomeFilePath\tgenomeVersion\tdataSourceId\tgenomeSpeciesId\tfakeGeneIdPrefix\tkeywords\tcomment\n";

# get the species info for ensembl
my $response = get_json($ua, $ensembl_server, $ensembl_species);
write_ensembl_species_info($response, $ua, $ensembl_server, $fh, $ensembl_dataSource_id, \%filtered_species_hash);

# get the species info for ensembl metazoa
$response = get_json($ua, $ensembl_server, $ensembl_species.'?division=ensemblMetazoa');
write_ensembl_species_info($response, $ua, $ensembl_server, $fh, $ensembl_metazoa_datasource_id, \%filtered_species_hash);

close $fh;