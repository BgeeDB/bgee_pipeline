#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use File::Basename;
use File::Slurp;
use List::MoreUtils qw{uniq any};
use List::Compare;
use LWP::Simple;
use XML::Fast;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;


# Define arguments & their default value
my ($species, $bgee_connector, $bgee_species) = ('', '', '');
my ($obsGO) = ('');
my ($debug) = (0);
my %opts = ('species=s'     => \$species,            # speciesCommonName from TSV for or Bgee db
            'bgee=s'        => \$bgee_connector,     # Bgee connector string
            'bgeeSpecies=s' => \$bgee_species,       # Bgee species main file
            'obsGO=s'       => \$obsGO,              # go.obsolete file
            'debug'         => \$debug,              # debug mode, do not insert/update in database
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $species eq '' || $bgee_connector eq '' || $bgee_species eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -species=9606__0__NonEnsembl  -bgee=\$(BGEECMD)  -obsGO=go.obsolete  -bgeeSpecies=\$(SPECIESFILEPATH)
\t-species     speciesId from Bgee db with the genomeSpeciesId concatenated
\t-bgee        Bgee    connector string
\t-obsGO       go.obsolete file
\t-bgeeSpecies bgeeSpecies.tsv file
\t-debug       Debug mode, do not insert/update in database
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


# Need to map to another genomeSpeciesId?
my ($speciesBgee, $newSpecies, $scientific_name, $NonEnsSource) = split('__', $species, -1);
if ( $NonEnsSource eq 'Ensembl' || $NonEnsSource eq 'EnsemblMetazoa' ){
    die "This script is not for Ensembl sources\n";
}
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


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#NOTE currently this script is tested only with RefSeq GTF as source!
#NOTE and works only for  "protein_coding"  genes!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# Get species genome genomeFilePath to fetch GTF annotations
my ($species_info) = grep { /^$speciesBgee\t/ } read_file("$bgee_species", chomp => 1);
my @species_info = split(/\t/, $species_info);
my $genomeFilePath = $species_info[5];

# Fetch GTF annotation file for gene annotations
my $gtf_file = basename("$genomeFilePath.gtf.gz");
if ( $genomeFilePath =~ /^\w+\/\w+?_((GC[FA])_(\d\d\d)(\d\d\d)(\d\d\d).*)$/ ){
    # From NCBI, RefSeq or GenBank assembly annotations
    #See https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/ for help
    # e.g. macaca_fuscata/macaca_fuscata_GCA_003118495.1_macFus_1.0
    #      manis_javanica/manis_javanica_GCF_014570535.1_YNU_ManJav_2.0
    #if ( is_success( getstore("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/$2/$3/$4/$5/$1/${1}_genomic.gtf.gz", basename("$genomeFilePath.gtf.gz")) ) ){
    #FIXME don't knwo why LWP::Simple fails on that URL! Passive FTP ???
    if ( system("wget -O '$gtf_file' 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/$2/$3/$4/$5/$1/${1}_genomic.gtf.gz'")==0 ){
        # done
    }
    else {
        die "Download of GTF file [$genomeFilePath] failed\n";
    }
}
else {
    die "genomeFilePath [$genomeFilePath] not valid\n";
}


# Read gene annotations
open(my $GTF, "zcat $gtf_file|")  or die "Cannot read the GTF file [$gtf_file]\n";
my $annotations;
#NOTE only "protein_coding" gene_biotype is used here because information are got from UniProt
while(<$GTF>){
    next  if ( !/gene_biotype "protein_coding"/ );

    #NC_030727.1	Gnomon	gene	141537768	141585606	.	+	.	gene_id "LOC108709019"; db_xref "GeneID:108709019"; gbkey "Gene"; gene "LOC108709019"; gene_biotype "protein_coding"; 
    #NC_030735.1	BestRefSeq%2CGnomon	gene	66815561	66907873	.	+	.	gene_id "bcl2.S"; db_xref "GeneID:100271914"; db_xref "Xenbase:XB-GENE-6251434"; description "B-cell CLL/lymphoma 2 S homeolog"; gbkey "Gene"; gene "bcl2.S"; gene_biotype "protein_coding"; gene_synonym "bcl-2"; gene_synonym "bcl2"; gene_synonym "xBcl-2"; gene_synonym "xbcl2"; 
    my @gene_info = split(/\t/, $_);
    my $GeneID;
    my %info;
    # Current list of possible annotations:
    # db_xref, description, gbkey, gene, gene_biotype, gene_id, gene_synonym, partial, transcript_id
    # You may have several  db_xref  and  gene_synonym 
    for my $annot ( grep { /\w/ } split(/" *; +/, $gene_info[8]) ){
        if ( $annot =~ /^(\w+)\s+"([^"]+)$/ ){
            my $key   = $1;
            my $value = $2;
            if ( $key eq 'db_xref' || $key eq 'gene_synonym' ){
                push @{ $info{$key} }, $value;
                if ( $value =~ /^GeneID:(\d+)$/ ){
                    $GeneID = $1;
                }
            }
            else {
                $info{$key} = $value;
            }
        }
    }
    #NOTE some gene ids may not be uniq in the Bgee database
    $info{'gene_id'} = 'LOC'.$GeneID;

    $annotations->{ $info{'gene_id'} } = \%info;
}
close $GTF;
unlink "$gtf_file";


## Biotypes
# Get previously inserted BioTypes
my $biotypeDB = $dbh->prepare('SELECT geneBioTypeName FROM geneBioType');
$biotypeDB->execute()  or die $biotypeDB->errstr;
my @InsertedBioTypes = uniq map { $_->[0] } @{$biotypeDB->fetchall_arrayref};
$biotypeDB->finish;

# Get BioTypes for this species
my @specificBioTypes = uniq
                       map { $annotations->{$_}->{'gene_biotype'} }
                       keys %$annotations; # Get through genes list and return uniq (non-redundant) biotype from it

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
my $geneDB       = $dbh->prepare('INSERT INTO gene (geneId, geneName, geneDescription, geneBioTypeId, speciesId)
                                  VALUES (?, ?, ?, (SELECT geneBioTypeId FROM geneBioType WHERE geneBioTypeName=?), ?)');
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
    
}
__END__
#https://www.uniprot.org/uniprot/?query=ENSDARG00000024771&format=xml&force=true&columns=id,entry name,reviewed,protein names,genes,organism,database(CCDS),database(EMBL),database(RefSeq),database(PDB),database(Ensembl),database(EnsemblMetazoa),database(GeneID),go
#https://www.uniprot.org/uniprot/?query=ENSDARG00000024771&format=xml&force=true
#Cases with one ensembl -> several ensembl xrefs: ENSG00000139618
# Fetch all clones from a slice adaptor (returns a list reference)
my @genes = @{$gene_adaptor->fetch_all()};



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
        $geneDB->execute($id2insert, $external_name, $description, $biotype, $speciesBgee)  or die $geneDB->errstr;
        $bgeeGeneId = $dbh->{mysql_insertid};
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
    for my $xref ( sort keys %xrefs ){
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
    for my $go ( sort keys %GO ){
        if ( ! $debug ){
            $goDB->execute($bgeeGeneId, uc $go, $GO{$go})  or do {warn "[$id2insert] [$go] [$GO{$go}]\n"; die $goDB->errstr};
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

