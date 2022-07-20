#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use List::MoreUtils qw{uniq};
use List::Compare;
use LWP::Simple;
use File::Slurp;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;


# Define arguments & their default value
my ($species, $bgee_connector, $ensembl_connector) = ('', '', '');
my ($obsGO) = ('');
my ($debug) = (0);
my %opts = ('species=s'    => \$species,            # speciesCommonName from TSV for or Bgee db
            'bgee=s'       => \$bgee_connector,     # Bgee connector string
            'ensembl=s'    => \$ensembl_connector,  # Ensembl connector string
            'obsGO=s'      => \$obsGO,              # go.obsolete file
            'debug'        => \$debug,              # debug mode, do not insert/update in database
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $species eq '' || $bgee_connector eq '' || $ensembl_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -species=9606__0__Ensembl  -bgee=\$(BGEECMD)  -ensembl=\$(ENSCMD)  -obsGO=go.obsolete
\t-species   speciesId from Bgee db with the genomeSpeciesId concatenated
\t-bgee      Bgee    connector string
\t-ensembl   Ensembl connector string
\t-obsGO     go.obsolete file
\t-debug     Debug mode, do not insert/update in database
\n";
    exit 1;
}

# Ensembl connection via Ensembl API/Registry
my $reg = Utils::connect_ensembl_registry($ensembl_connector, $debug);
# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


# Need to map to another genomeSpeciesId?
my ($speciesBgee, $newSpecies, $scientific_name, $ensSource) = split('__', $species, -1);
my @prefix;
if ( $speciesBgee == $newSpecies || $newSpecies == 0 ){
    # No mapping to another species
    $species = $speciesBgee;
}
else {
    $species = $newSpecies;
    my $selSpecies = $dbh->prepare('SELECT fakeGeneIdPrefix FROM species WHERE speciesId=? AND genomeSpeciesId=?');
    $selSpecies->execute($speciesBgee, $newSpecies)  or die $selSpecies->errstr;
    @prefix = map { $_->[0] } @{$selSpecies->fetchall_arrayref};
    $selSpecies->finish();
    die "Too many prefixes returned [@prefix]\n"  if ( exists $prefix[1] );
}


# Get a slice adaptor for the  $species  core database
my $gene_adaptor;
if ( $scientific_name eq 'Heterocephalus_glaber' ){
    #Specific genome versions for female and male
    $scientific_name = 'Heterocephalus_glaber_female';
}
if ( $ensSource eq 'Ensembl' || $ensSource eq 'EnsemblMetazoa' ){
    $gene_adaptor = $reg->get_adaptor( lc $scientific_name, 'Core', 'Gene' );
}
else {
    die "Unknown Ensembl db source [$ensSource]\n";
}
# Fetch all clones from a slice adaptor (returns a list reference)
my @genes = @{$gene_adaptor->fetch_all()};
#print scalar @genes, " genes for $scientific_name\n";
#exit;



## Biotypes
# Get previously inserted BioTypes
my $biotypeDB = $dbh->prepare('SELECT geneBioTypeName FROM geneBioType');
$biotypeDB->execute()  or die $biotypeDB->errstr;
my @InsertedBioTypes = uniq map { $_->[0] } @{$biotypeDB->fetchall_arrayref};
$biotypeDB->finish;

# Get BioTypes for this species
my @specificBioTypes = uniq
                       map { $_->biotype() }
                       @genes; # Get through genes list and return uniq (non-redundant) biotype from it

# Insert only new BioTypes
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


## Alt GO
# Get alt_id GO
my $altgoDB = $dbh->prepare('SELECT goAltId, goId FROM geneOntologyTermAltId');
$altgoDB->execute()  or die $altgoDB->errstr;
my %altid_go = map { lc $_->[0] => lc $_->[1] }
               @{$altgoDB->fetchall_arrayref};
$altgoDB->finish;


## Obsolete GO
# Get obsolete GO to not insert them later on
die "Missing [go.obsolete] file\n"  if ( !-e 'go.obsolete' || -z 'go.obsolete' );
my %obs_go = map { lc $_ => 1 }
             read_file("$obsGO", chomp => 1);


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
#TODO Add other species specific dataSource variant names:
#TODO Add CGNC (Chichen Gene Nomenclature Consortium)?
$InsertedDataSources{'flybasename_gene'} = $InsertedDataSources{'flybase'};
$InsertedDataSources{'flybasecgid_gene'} = $InsertedDataSources{'flybase'};
$InsertedDataSources{'flybase_symbol'}   = $InsertedDataSources{'flybase'};
$InsertedDataSources{'wormbase_gene'}    = $InsertedDataSources{'wormbase'};
$InsertedDataSources{'xenopus_jamboree'} = $InsertedDataSources{'xenbase'};
$InsertedDataSources{'zfin_id'}          = $InsertedDataSources{'zfin'};


## Gene info (id, description)
# Get individual gene info
my $geneDB       = $dbh->prepare('INSERT INTO gene (geneId, geneName, geneDescription, geneBioTypeId, speciesId, dataSourceId)
                                  VALUES (?, ?, ?, (SELECT geneBioTypeId FROM geneBioType WHERE geneBioTypeName=?), ?, ?)');
my $synonymDB    = $dbh->prepare('INSERT INTO geneNameSynonym (bgeeGeneId, geneNameSynonym)
                                  VALUES (?, ?)');
my $xrefDB       = $dbh->prepare('INSERT INTO geneXRef (bgeeGeneId, XRefId, XRefName, dataSourceId)
                                  VALUES (?, ?, ?, ?)');
my $goDB         = $dbh->prepare('INSERT INTO geneToGeneOntologyTerm (bgeeGeneId, goId, goEvidenceCode)
                                  VALUES (?, ?, ?)');
my $geneToTermDB = $dbh->prepare('INSERT INTO geneToTerm (bgeeGeneId, term)
                                  VALUES (?, ?)');
print "Inserting gene info...\n";
GENE:
for my $gene (sort {$a->stable_id cmp $b->stable_id} (@genes)) { #Sort to always get the same order
    my $display_id    = $gene->display_id();
    my $stable_id     = $gene->stable_id();
    my $external_name = $gene->external_name() || '';
    my $external_db   = $gene->external_db()   || '';
    my $description   = $gene->description()   || '';
    my $biotype       = $gene->biotype()       || die "Invalid BioType for $stable_id\n";

    ## Cleaning
    # Remove useless whitespace(s)
    $description     =~ s{  +}{ }g;
    $description     =~ s{[\.,]+ *\[Source:}{ \[Source:};
    # Remove HTML tags in gene names
    $external_name   =~ s{<[^>]+?>}{}g;

    warn "Different ids [$display_id] [$stable_id]\n"  if ( $display_id ne $stable_id);

    ## external genome mapping
    my $id2insert = $stable_id;
    if ( exists $prefix[0] && $prefix[0] ne '' ){
        #FIXME only for Ensembl ids now, AND non-human Ensembl ids
        $id2insert =~ s{^ENS[A-Z][A-Z][A-Z]G}{$prefix[0]};
    }

    ## Insert gene info
    my $bgeeGeneId;
    if ( ! $debug ){
        $geneDB->execute($id2insert, $external_name, $description, $biotype, $speciesBgee, $datasourceId)  or die $geneDB->errstr;
        $bgeeGeneId = $dbh->{'mysql_insertid'};
        die "Cannot get bgeeGeneId [$bgeeGeneId]\n"  if ( $bgeeGeneId !~ /^\d+$/ );
    }
    else {
        print "\n[$id2insert] [$external_name] [$description]   [$biotype] [$speciesBgee]\n";
    }


    ## Get gene synonyms, if any
    #NOTE Synonyms shown on the web site appear to come from THE main gene xref synonyms
    #     But we need all of them, so no filtering based on main external_db !!!
#TODO Remove part of synonym within "{...}" and split the rest on "|" !!!
    my @synonyms = uniq sort                                                              # non-redundant & sorted
                   map  { s{^\s+}{}; s{\s+$}{}; lc $_ }                                   # Trim & lowercase
                   grep { $_ ne $stable_id && $_ ne $display_id && $_ ne $external_name } # Avoid putting $display_id as synonym
                   map  { @{$_->get_all_synonyms} }                                       # Official xref synonyms
#                   grep { $_->dbname() eq $external_db }                                  # Only external_db that are main gene db source
                   @{$gene->get_all_xrefs()};
    SYNONYM:
    for my $syn ( @synonyms ){
        if ( ! $debug ){
            $synonymDB->execute($bgeeGeneId, $syn)  or die $synonymDB->errstr;
        }
        else {
            print "synonym: [$syn]\n";
        }
    }


    ## Get Xref (linked in dataSource table BUT not GO)
    # Show all datasources (but GO)
    if ( $debug ){
        print join(' | ', uniq sort
                          map  { $_->dbname() }
                          grep { $_->dbname() ne 'GO' }
                          @{$gene->get_all_xrefs()}
                  ), "\n";
    }
    my %xrefs = map  { my $dbname = $_->dbname();
                       my $pid    = $_->primary_id();
                       "$dbname##$pid" => $_->display_id() }           # Remove duplicates
                grep { exists $InsertedDataSources{lc $_->dbname()} }  # Only external db in dataSource table
                grep { $_->dbname() ne 'GO' }                          # GO xrefs have there own table, so not in XRefs
                @{$gene->get_all_xrefs()};

    # Get UniProt ID, missing in Ensembl that contains only UniProt AC
    UNIPROT_ID:
    for my $uniprot ( sort grep { /^Uniprot\/SPTREMBL##/ || /^Uniprot\/SWISSPROT##/} keys %xrefs ){
        my ($dbname, $pid) = split('##', $uniprot);
        my $content = get("http://www.uniprot.org/uniprot/?query=id:$pid&format=tab&columns=entry%20name");
        if ( defined $content ){
            my (undef, $uniprot_id) = split("\n", $content, -1);
            #TODO CHECK !!!
            for my $uid ( split(/\s*;\s*/, $uniprot_id) ){
                $xrefs{"$dbname##$uid"} = $uid;
            }
        }
        else {
            warn "\tCannot get UniProt ID for [$pid]\n";
        }
    }

    XREF:
    for my $xref ( uniq sort keys %xrefs ){
        #TODO CHECK !!!
        next  if ( $xref =~ /;/ );
        my ($dbname, $pid) = split('##', $xref);
        $xrefs{$xref} = ''  if ( $xrefs{$xref} eq $pid );
        if ( ! $debug ){
            $xrefDB->execute($bgeeGeneId, $pid, $xrefs{$xref}, $InsertedDataSources{lc $dbname})  or die $xrefDB->errstr;
        }
        else {
            print "xref: [$id2insert] [$pid] [$xrefs{$xref}] [$dbname]\n";
        }
    }


    ## Get GO xref
    #TODO alt_id GO test  GO:0007243
    my %GO = map  { lc $_->display_id => ${$_->{'linkage_types'}->[0]}[0] }
             grep { $_->dbname() eq 'goslim_goa' } # Extra goslim_goa terms
             @{$gene->get_all_xrefs()};
    %GO    = map  { lc $_->display_id => ${$_->{'linkage_types'}->[0]}[0] }
             grep { $_->dbname() eq 'GO' } # Replace by GO itself and its Evidence Code if any, so keep intersection between GO & goslim_goa in Ensembl
             @{$gene->get_all_xrefs()};
    %GO    = map  { exists $altid_go{$_} ? ($altid_go{$_} => $GO{$_}) : ($_ => $GO{$_}) } # Replace alt_id GO by main GO if any
             grep { !exists $obs_go{$_} }                                                 # Skip obsolete GO, not inserted in Bgee previously
             keys %GO;
    GO:
    for my $go ( uniq sort keys %GO ){
        if ( ! $debug ){
            $goDB->execute($bgeeGeneId, uc $go, $GO{$go})  or do {warn "[$id2insert] [$go] [$GO{$go}]\n"};# die $goDB->errstr};
        }
        else {
            print "go:   [$id2insert] [$go] [$GO{$go}]\n";
        }
    }


    ## Everything in geneToTerm
    # gene_id + version
    my @all;
    push @all, $display_id                      if ( $display_id );
    push @all, $id2insert                       if ( $id2insert );
    push @all, $external_name                   if ( $external_name );
    #NOTE version() does not seem to return something for species with non Ensembl gene ids such as C. elegans WBGene00000001
    push @all, $id2insert.'.'.$gene->version()  if ( $id2insert && $gene->version() );
    # transcript ids
    push @all, map  { if ( $_->version() ){ ($_->stable_id(), $_->stable_id().'.'.$_->version()) } else { $_->stable_id() } }
               grep { defined $_->stable_id() }
               @{$gene->get_all_Transcripts};
    # exon ids
    push @all, map  { if ( $_->version() ){ ($_->stable_id(), $_->stable_id().'.'.$_->version()) } else { $_->stable_id() } }
               grep { defined $_->stable_id() }
               @{$gene->get_all_Exons};
    # translation ids
    push @all, map  { if ( $_->version() ){ ($_->stable_id(), $_->stable_id().'.'.$_->version()) } else { $_->stable_id() } }
               grep { defined $_->stable_id() }
               map  { $_->translation() }
               grep { defined $_->translation() }
               @{$gene->get_all_Transcripts};
    # Xref display ids
    push @all, map  { $_->display_id() }
               @{$gene->get_all_xrefs()};
    # Xref primary ids
    push @all, map  { $_->primary_id() }
               @{$gene->get_all_xrefs()};
    # Xref synonyms
    push @all, map  { @{ $_->get_all_synonyms } }
               @{$gene->get_all_xrefs()};
    # Extra Xref in gene description
    if ( $description ne '' ){
        #EC number
        while ( $description =~ /E\.?C\.?\s*([1-6]\.[\d\-]+\.[\d\-]+\.[\d\-]+)/g ){
            push @all, $1;
        }
        #Entry source
        while ( $description =~ /\[Source:\s*.+?\s*;Acc:\s*(.+?)\s*\]/g ){
            push @all, $1;
        }
    }
    # Remove duplicates AND empty strings AND trim!!!
    @all = uniq sort                                          # non-redundant & sorted
           grep { $_ ne '' && defined $_ }                    # Non-empty & non-undef
           grep { !exists $obs_go{$_} }                       # Skip obsolete GO, not inserted in Bgee previously
           map  { exists $altid_go{$_} ? $altid_go{$_} : $_ } # Replace alt_id GO by main GO if any
           map  { s{^\s+}{}; s{\s+$}{}; lc $_ }               # Trim + lowercase because same entry in different cases
           @all;
    ALL:
    for my $term ( @all ){
        if ( ! $debug ){
            $geneToTermDB->execute($bgeeGeneId, $term)  or die $geneToTermDB->errstr;
        }
        else {
            print "term: [$id2insert] [$term]\n";
        }
    }
}
$geneDB->finish;
$synonymDB->finish;
$xrefDB->finish;
$goDB->finish;
$geneToTermDB->finish;
print "Gene nbr for $scientific_name: ", scalar @genes, "\n\n";


# Close db connections
$reg->clear();
$dbh->disconnect;

exit 0;

