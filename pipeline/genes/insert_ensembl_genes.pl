#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use JSON::XS;
use List::MoreUtils qw{uniq};
use List::Compare;
use File::Slurp;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;


# Define arguments & their default value
my ($species, $bgee_connector) = ('', '');
my ($uniprot_mapping) = ('');
my ($specific_gene) = ('');
my ($debug) = (0);
my %opts = ('species=s'    => \$species,            # speciesCommonName from TSV for or Bgee db
            'bgee=s'       => \$bgee_connector,     # Bgee connector string
            'uni=s'        => \$uniprot_mapping,    # UniProt mapping AC -> ID
            'debug'        => \$debug,              # debug mode, do not insert/update in database
            'gene=s'       => \$specific_gene,      # Query a single gene
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $species eq '' || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -species=9606__0__Homo_sapiens__Ensembl  -bgee=\$(BGEECMD)  -uni=\$(OUTPUT_DIR)/uniprot.map  <Ensembl species JSON>
\t-species   speciesId from Bgee db with the genomeSpeciesId concatenated
\t-bgee      Bgee    connector string
\t-uni       UniProt mapping file AC -> ID
\t-debug     Debug mode, do not insert/update in database
\t
\t-gene      Query a single specific gene (optional)
\n";
    exit 1;
}


## Load the JSON
my $input_json_file = $ARGV[0]  or die "\n\tMissing input JSON file\n\n";
my $json_text = do {
    open(my $json_fh, '<:encoding(UTF-8)', $input_json_file)  or die("Can't open \"$input_json_file\": $!\n");
    local $/;
    <$json_fh>;
};
my $ensembl_json = decode_json($json_text);
#NOTE keep only the genes section of the JSON!
map { delete( $ensembl_json->{$_} ) } grep { $_ ne 'genes' } keys %$ensembl_json;
my @genes = @{ $ensembl_json->{'genes'} };


# UniProt mapping
my $uniprot_map;
if ( -e "$uniprot_mapping" ){
    #Format: AC\tID
    map { my @tmp = split(/\t/, $_); $uniprot_map->{$tmp[0]} = $tmp[1]; } read_file("$uniprot_mapping", chomp => 1);
}


# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


# Need to map to another genomeSpeciesId? #NOTE currently not used
my ($speciesBgee, $newSpecies, $scientific_name, $ensSource) = split('__', $species, -1);
#TODO female or male naked-mole rat???
if ( $scientific_name eq 'Heterocephalus_glaber' ){
    #Specific genome versions for female and male
    $scientific_name = 'Heterocephalus_glaber_female';
}


## Biotypes
# Get BioTypes for this species
my @specificBioTypes = uniq
                       map { $_->{'biotype'} }
                       @genes;
# Get previously inserted BioTypes
my $biotypeDB = $dbh->prepare('SELECT geneBioTypeName FROM geneBioType');
$biotypeDB->execute()  or die $biotypeDB->errstr;
my @InsertedBioTypes = uniq map { $_->[0] } @{$biotypeDB->fetchall_arrayref};
$biotypeDB->finish;


## Insert only new BioTypes
my $lc = List::Compare->new(\@specificBioTypes, \@InsertedBioTypes);
my @newBioTypes = $lc->get_unique; # Get entries in the 1st list not in the 2nd
if ( exists $newBioTypes[0] ){
    $biotypeDB = $dbh->prepare('INSERT INTO geneBioType (geneBioTypeName) VALUES (?)');
    for my $biotype ( @newBioTypes ){
        #NOTE Fix to avoid warning because different cases between Ensembl db
        # and we would like to keep the right case for names, not everything in lc or uc!
        next  if ( $biotype =~ /^3prime_overlapping_ncRNA$/i      && grep { /^3prime_overlapping_ncRNA$/i }      @InsertedBioTypes );
        next  if ( $biotype =~ /^bidirectional_promoter_lncRNA$/i && grep { /^bidirectional_promoter_lncRNA$/i } @InsertedBioTypes );
        next  if ( $biotype =~ /^miRNA$/i                         && grep { /^miRNA$/i }                         @InsertedBioTypes );
        if ( ! $debug ){
            $biotypeDB->execute($biotype)  or die $biotypeDB->errstr;
        }
    }
    $biotypeDB->finish;

    print "Inserting new BioTypes\n", join("\t", @newBioTypes), "\n";
}


## DataSources
# Get used dataSources
my $sourceDB = $dbh->prepare('SELECT dataSourceId, dataSourceName FROM dataSource');
$sourceDB->execute()  or die $sourceDB->errstr;
my %InsertedDataSources = map { $_->[1] = lc $_->[1]; $_->[1] => $_->[0] } @{$sourceDB->fetchall_arrayref}; # List already inserted dataSources and return it in a hash dataSourceName (lowercase) => dataSourceId
$sourceDB->finish;
my $datasourceId = $ensSource eq 'Ensembl'        ? $InsertedDataSources{'ensembl'}
                 : $ensSource eq 'EnsemblMetazoa' ? $InsertedDataSources{'ensemblmetazoa'}
                 : '';

# Add extra dataSource aliases
# MUST be in lowercase to ease comparison
#TODO Add other species specific dataSource variant names
#TODO Add CGNC (Chichen Gene Nomenclature Consortium)?
$InsertedDataSources{'flybasename_gene'}    = $InsertedDataSources{'flybase'};
$InsertedDataSources{'flybasecgid_gene'}    = $InsertedDataSources{'flybase'};
$InsertedDataSources{'flybase_symbol'}      = $InsertedDataSources{'flybase'};
$InsertedDataSources{'wormbase_gene'}       = $InsertedDataSources{'wormbase'};
$InsertedDataSources{'wormbase_cds'}        = $InsertedDataSources{'wormbase'};
$InsertedDataSources{'wormbase_transcript'} = $InsertedDataSources{'wormbase'};
$InsertedDataSources{'wormpep_id'}          = $InsertedDataSources{'wormbase'};
$InsertedDataSources{'xenopus_jamboree'}    = $InsertedDataSources{'xenbase'};
$InsertedDataSources{'zfin_id'}             = $InsertedDataSources{'zfin'};
# Not species specific
$InsertedDataSources{'ensembl_metazoa'}     = $InsertedDataSources{'ensemblmetazoa'};
$InsertedDataSources{'refseq_mrna'}         = $InsertedDataSources{'refseq nucleotide'};
$InsertedDataSources{'refseq_ncrna'}        = $InsertedDataSources{'refseq nucleotide'};
$InsertedDataSources{'refseq_peptide'}      = $InsertedDataSources{'refseq protein'};
$InsertedDataSources{'refseq_dna'}          = $InsertedDataSources{'refseq genomic'};
$InsertedDataSources{'protein_id'}          = $InsertedDataSources{'embl'};
$InsertedDataSources{'mirbase_trans_name'}  = $InsertedDataSources{'mirbase'};
$InsertedDataSources{'uniprot/swissprot'}   = $InsertedDataSources{'uniprotkb/swiss-prot'};
$InsertedDataSources{'uniprot/sptrembl'}    = $InsertedDataSources{'uniprotkb/trembl'};
my %UnknownDataSources;


## Gene info (id, description)
# Get individual gene info
my $geneDB    = $dbh->prepare('INSERT INTO gene (geneId, geneName, geneDescription, geneBioTypeId, speciesId, ensemblGene, seqRegionName)
                                  VALUES (?, ?, ?, (SELECT geneBioTypeId FROM geneBioType WHERE geneBioTypeName=?), ?, ?, ?)');
my $synonymDB = $dbh->prepare('INSERT INTO geneNameSynonym (bgeeGeneId, geneNameSynonym)
                                  VALUES (?, ?)');
my $xrefDB    = $dbh->prepare('INSERT INTO geneXRef (bgeeGeneId, XRefId, XRefName, dataSourceId)
                                  VALUES (?, ?, ?, ?)');
print "Inserting gene info...\n";
GENE:
for my $gene (sort {$a->{'id'} cmp $b->{'id'}} (@genes)) { #Sort to always get the same order
    if ( $specific_gene ){
        next  if ( $gene->{'id'} ne $specific_gene );
    }

    my $stable_id       = $gene->{'id'};
    my $external_name   = $gene->{'name'}             || '';
    my $description     = $gene->{'description'}      || '';
    my $biotype         = $gene->{'biotype'}          || die "Invalid BioType for $stable_id\n";
    my $seq_region_name = $gene->{'seq_region_name'}  || '';

    ## Cleaning
    # Remove useless whitespace(s)
    $description     =~ s{  +}{ }g;
    $description     =~ s{[\.,]+ *\[Source:}{ \[Source:};
    # Remove HTML tags in gene names
    $external_name   =~ s{<[^>]+?>}{}g;
    $external_name   =~ s{ \[provisional:(.+?)\]$}{$1}; #e.g. XB5961369 [provisional:plpp3

    ## Insert gene info
    my $bgeeGeneId;
    if ( ! $debug ){
        #NOTE ensemblGene = 1 for Ensembl & EnsemblMetazoa
        $geneDB->execute($stable_id, $external_name, $description, $biotype, $speciesBgee, 1, $seq_region_name)  or die $geneDB->errstr;
        $bgeeGeneId = $dbh->{'mysql_insertid'};
        die "Cannot get bgeeGeneId [$bgeeGeneId]\n"  if ( $bgeeGeneId !~ /^\d+$/ );
    }
    else {
        print "\n[$stable_id] [$external_name] [$description]   [$biotype] [$speciesBgee] [$seq_region_name]\n";
    }

    ## Get gene synonyms, if any
    #NOTE Synonyms shown on the web site appear to come from THE main gene xref synonyms
    #     But we need all of them, so no filtering !!!
#TODO Remove part of synonym within "{...}" and split the rest on "|" !!!
    #from xrefs
    my @synonyms = map  { s{^\s+}{}; s{\s+$}{}; lc $_ }              # Trim & lowercase
                   grep { $_ ne $stable_id && $_ ne $external_name } # Avoid putting main name/id as synonym
                   map  {
                          if ( exists $_->{'synonyms'}->[0] ){
                            (
                              $_->{'display_id'},
                              @{ $_->{'synonyms'} },
                            )
                          }
                          else {
                            (
                              $_->{'display_id'},
                            )
                          }
                        } # Official xref synonyms
                   @{ $gene->{'xrefs'} };
    #from transcript id and xrefs
    push @synonyms, map  { s{^\s+}{}; s{\s+$}{}; lc $_ }              # Trim & lowercase
                    grep { $_ ne $stable_id && $_ ne $external_name } # Avoid putting main name/id as synonym
                    map  { $_->{'display_id'} }
                    map  { @{ $_->{'xrefs'} } }
                    @{ $gene->{'transcripts'} };
    push @synonyms, map { s{^\s+}{}; s{\s+$}{}; lc $_ }              # Trim & lowercase
                    map { $_->{'id'} }
                    @{ $gene->{'transcripts'} };
    #from exon id
    push @synonyms, map { s{^\s+}{}; s{\s+$}{}; lc $_ }              # Trim & lowercase
                    map { $_->{'id'} }
                    map { @{ $_->{'exons'} } }
                    @{ $gene->{'transcripts'} };
    #from translation id and xrefs
    for my $transcript ( @{ $gene->{'transcripts'} } ){ # Do like this because some genes may have some coding AND non-coding transcripts at the same time, e.g. ENSG00000000003
        if ( exists $transcript->{'translations'} ){
            push @synonyms, map  { s{^\s+}{}; s{\s+$}{}; lc $_ }              # Trim & lowercase
                            grep { $_ ne $stable_id && $_ ne $external_name } # Avoid putting main name/id as synonym
                            map  { $_->{'display_id'} }
                            map  { @{ $_->{'xrefs'} } }
                            @{ $transcript->{'translations'} };
            push @synonyms, map { s{^\s+}{}; s{\s+$}{}; lc $_ }              # Trim & lowercase
                            map { $_->{'id'} }
                            @{ $transcript->{'translations'} };
        }
    }
    SYNONYM:
    for my $syn ( uniq sort @synonyms ){
        if ( ! $debug ){
            $synonymDB->execute($bgeeGeneId, $syn)  or die $synonymDB->errstr;
        }
        else {
            print "synonym: [$syn]\n";
        }
    }


    ## Get Xref (linked in dataSource table BUT not GO)
    my %xrefs = map  { my $dbname = $_->{'dbname'};
                       my $pid    = $_->{'primary_id'};
                       "$dbname##$pid" => $_->{'display_id'} }           # Remove duplicates
#                grep { exists $InsertedDataSources{lc $_->{'dbname'}} }  # Only external db in dataSource table
                @{ $gene->{'xrefs'} };
    #from transcript id and xrefs
    map { $xrefs{"Ensembl##$_"} = $_ }
        map { $_->{'id'} }
        @{ $gene->{'transcripts'} };
    map { my $dbname = $_->{'dbname'};
          my $pid    = $_->{'primary_id'};
          $xrefs{"$dbname##$pid"} = $_->{'display_id'}
        }
        map  { @{ $_->{'xrefs'} } }
        @{ $gene->{'transcripts'} };
    #from exon id
    map { $xrefs{"Ensembl##$_"} = $_ }
        map { $_->{'id'} }
        map { @{ $_->{'exons'} } }
        @{ $gene->{'transcripts'} };
    #from translation id and xrefs
    for my $transcript ( @{ $gene->{'transcripts'} } ){ # Do like this because some genes may have some coding AND non-coding transcripts at the same time, e.g. ENSG00000000003
        if ( exists $transcript->{'translations'} ){
            map { $xrefs{"Ensembl##$_"} = $_ }
                map { $_->{'id'} }
                @{ $transcript->{'translations'} };
            map { my $dbname = $_->{'dbname'};
                  my $pid    = $_->{'primary_id'};
                $xrefs{"$dbname##$pid"} = $_->{'display_id'}
                }
                map { @{ $_->{'xrefs'} } }
                @{ $transcript->{'translations'} };
        }
    }

    # Get UniProt ID, missing in Ensembl that contains only UniProt AC
    UNIPROT_ID:
    for my $uniprot ( sort grep { /^Uniprot\/SPTREMBL##/ || /^Uniprot\/SWISSPROT##/} keys %xrefs ){
        my ($dbname, $pid) = split('##', $uniprot);
        my $uid = $uniprot_map->{ $pid } || do { warn "\tCannot get UniProt ID for [$pid]\n"; next};
        $xrefs{"$dbname##$uid"} = $uid;
    }

    XREF:
    for my $xref ( uniq sort grep { !/^GO##/ } keys %xrefs ){
        next  if ( $xref =~ /;/ );
        my ($dbname, $pid) = split('##', $xref);
        $xrefs{$xref} = ''  if ( $xrefs{$xref} eq $pid );
        if ( ! $InsertedDataSources{lc $dbname} ){
            $UnknownDataSources{lc $dbname} = 1;
            next XREF;
        }

        if ( ! $debug ){
            $xrefDB->execute($bgeeGeneId, $pid, $xrefs{$xref}, $InsertedDataSources{lc $dbname})  or warn $xrefDB->errstr;
        }
        else {
            print "xref: [$stable_id] [$pid] [$xrefs{$xref}] [$dbname]\n";
        }
    }
}
$geneDB->finish;
$synonymDB->finish;
$xrefDB->finish;
print "Gene nbr for $scientific_name: ", scalar @genes, "\n\n";
warn join(', ', sort keys %UnknownDataSources), "\n";


# Close db connections
$dbh->disconnect;

exit 0;

