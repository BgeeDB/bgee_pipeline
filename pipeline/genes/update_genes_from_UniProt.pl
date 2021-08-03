#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use List::MoreUtils qw{uniq any};
use LWP::Simple;
use XML::Fast;
use Data::Dumper;

use DBI;


# Define arguments & their default value
my ($gene_id, $dbname) = ('', '');
my ($verbose) = (0);
my %opts = ('gene=s'   => \$gene_id,   # gene id to search in UniProt to insert/update in the db
            'db=s'     => \$dbname,    # Bgee database name
            'verbose'  => \$verbose,   # verbose mode
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $dbname eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -db=bgee_v15_dev
\t     OR
\t     $0  -db=bgee_v15_dev  -gene=wbgene00001127
\t-gene     gene id to search in UniProt to insert/update in the db
\t-db       Bgee database name to search in UniProt all missing info
\t-verbose  verbose mode
\n";
    exit 1;
}


my $dbh = DBI->connect("dbi:mysql:database=$dbname;host=localhost;port=3306", 'root', 'BdD.pw')  or die $DBI::errstr;
my $up_xref  = $dbh->prepare('INSERT INTO geneXRef               (bgeeGeneId, XRefId, XRefName, dataSourceId) VALUES (?, ?, ?, ?)');
my $up_syn   = $dbh->prepare('INSERT INTO geneNameSynonym        (bgeeGeneId, geneNameSynonym)                VALUES (?, ?)');
my $up_term  = $dbh->prepare('INSERT INTO geneToTerm             (bgeeGeneId, term)                           VALUES (?, ?)');
my $up_go    = $dbh->prepare('INSERT INTO geneToGeneOntologyTerm (bgeeGeneId, goId, goEvidenceCode)           VALUES (?, ?, ?)');
my $up_gene1 = $dbh->prepare('UPDATE gene SET geneName=?         WHERE bgeeGeneId=?');
my $up_gene2 = $dbh->prepare('UPDATE gene SET geneDescription=?  WHERE bgeeGeneId=?');

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


#Test gene table for empty fields
#NOTE 454019 rows! => 454019 requests to UniProt
#NOTE BLxxxx ids from Amphioxus are not in UniProt BUT similar ids may exist there!
my $sth_genes = $dbh->prepare('SELECT bgeeGeneId, geneId, geneName, geneDescription FROM gene WHERE (geneDescription="" OR geneName="") AND geneId NOT LIKE "BL%"');
if ( $gene_id ne '' ){
    $sth_genes = $dbh->prepare("SELECT bgeeGeneId, geneId, geneName, geneDescription FROM gene  WHERE geneId='$gene_id' AND geneId NOT LIKE 'BL%'");
}
$sth_genes->execute()  or die $sth_genes->errstr;
GENE:
while ( my ($st_own_id, $st_gene_id, $st_gene_name, $st_description) = ($sth_genes->fetchrow_array) ){
    sleep 1;
    # Get info from UniProt
    my $annotations;
    my $content = get('https://www.uniprot.org/uniprot/?query=database%3A(type%3Aensembl+'.$st_gene_id.')+OR+database%3A(type%3Aensemblmetazoa+'.$st_gene_id.')&sort=score&format=xml&force=true'); #with Ensembl or EnsemblMetazoa xref
    #WARNING some cases with one ensembl -> several ensembl xrefs: ENSG00000139618
    if ( defined $content && $content ne '' && $content =~ /<\/uniprot>/ ){
        my $hash = xml2hash $content;
        #NOTE not easy to test if exists an array in hash ref. It works with eval!
        #NOTE may return several entries, keep the first one (the best one?)
        my $root = eval { exists $hash->{'uniprot'}->{'entry'}->[0] } ? $hash->{'uniprot'}->{'entry'}->[0] : $hash->{'uniprot'}->{'entry'};

        # Check the UniProt entry contains the xref used to query it, and is for the right species
        print "GeneID:$st_gene_id\n"  if ( $verbose );
        my @ensembl_xref = grep { $_->{'-type'} =~ /^Ensembl/ } @{ $root->{'dbReference'} }; #Ensembl, EnsemblMetazoa, ...
        if ( grep { $_->{'-value'} eq "$st_gene_id" && $_->{'-type'} eq 'gene ID' } map { @{ $_->{'property'} } } @ensembl_xref ){
            $annotations->{'gene_id'} = $st_gene_id;
#            print Dumper $root->{'dbReference'};
            #Get Ensembl protein ID
            ($annotations->{'prot_id'}) = map { $_->{'-value'} } grep { $_->{'-type'} eq 'protein sequence ID'} map { @{ $_->{'property'} } } @ensembl_xref;
            #Get Ensembl transcript ID
            ($annotations->{'transcript_id'}) = map { $_->{'-id'} } @ensembl_xref;

            #Get protein name
#            print Dumper $root->{'gene'};
            $annotations->{'prot_name'} = '';
            if ( ref $root->{'gene'} ne 'ARRAY' ){
                #NOTE to avoid some weird syntax such as GeneID:386601
                my @gene_names = eval { exists $root->{'gene'}->{'name'}->[0] } ? @{ $root->{'gene'}->{'name'} } : ($root->{'gene'}->{'name'});
                for my $gene_name ( sort @gene_names ){
                    next  if ( ref $gene_name eq 'ARRAY' );
                    if ( $gene_name->{'-type'} eq 'primary' ){
                        $annotations->{'prot_name'} = $gene_name->{'#text'};
                    }
                    else {
                        push @{ $annotations->{'synonyms'} }, $gene_name->{'#text'};
                    }
                }
            }
            push @{ $annotations->{'synonyms'} }, $annotations->{'prot_id'}, $annotations->{'transcript_id'};
            @{ $annotations->{'synonyms'} } = grep { $_ ne $annotations->{'prot_name'} } uniq @{ $annotations->{'synonyms'} };
            #Get protein description if any
            my @prot_desc_type = sort keys %{ $root->{'protein'} };
            my $prot_root;
            if ( scalar @prot_desc_type >= 1 ){
                $prot_root = $root->{'protein'}->{'recommendedName'} || $root->{'protein'}->{ $prot_desc_type[0] };
            }
            if ( $prot_root ){
                my $prot_name = eval { exists $prot_root->[0] } ? $prot_root->[0] : $prot_root;
                my $prot_desc = eval { exists $prot_name->{'fullName'}->{'#text'} } ? $prot_name->{'fullName'}->{'#text'} : $prot_name->{'fullName'};
                $annotations->{'prot_desc'} = $prot_desc  if ( $prot_desc !~ /LOC\d+/ );
            }

            #Get xrefs
#            print Dumper $root->{'dbReference'};
            push @{ $annotations->{'xrefs'} }, 'UniProtKB/TrEMBL:'.$root->{'name'};
            my @uniprotAC = eval { exists $root->{'accession'}->[0] } ? @{ $root->{'accession'} } : ($root->{'accession'});
            push @{ $annotations->{'xrefs'} }, map { 'UniProtKB/TrEMBL:'.$_ } @uniprotAC;
            # Xrefs
            my @used_xref_db = ('EMBL', 'CCDS', 'RefSeq', 'FlyBase', 'MGI', 'RGD', 'WormBase', 'Xenbase', 'ZFIN', 'HGNC', 'VGNC', 'Ensembl', 'EnsemblMetazoa');
            for my $dbref ( sort @{ $root->{'dbReference'} }){
                if ( any { $dbref->{'-type'} eq $_ } @used_xref_db ){
                    push @{ $annotations->{'xrefs'} }, $dbref->{'-type'}.':'.$dbref->{'-id'};
                    if ( exists $dbref->{'property'} ){
                        my @properties = eval { exists $dbref->{'property'}->[0] } ? @{ $dbref->{'property'} } : ($dbref->{'property'});
                        for my $property ( sort @properties ){
                            if ( $property->{'-type'} =~ /sequence ID$/ ){
                                push @{ $annotations->{'xrefs'} }, $dbref->{'-type'}.':'.$property->{'-value'};
                            }
                        }
                    }
                }
                elsif ( $dbref->{'-type'} eq 'GO' ){
                    for my $property ( sort @{ $dbref->{'property'} } ){
                        if ( $property->{'-type'} eq 'evidence' ){
                            push @{ $annotations->{'go'} }, $dbref->{'-id'}.'___'.$property->{'-value'};
                            last;
                        }
                    }
                }
            }
            @{ $annotations->{'xrefs'} } = map { s/GeneID:/GenBank:/; s/MGI:MGI:/MGI:/; $_ } uniq @{ $annotations->{'xrefs'} };
            @{ $annotations->{'go'} }    = uniq map { uc($_) }                                    @{ $annotations->{'go'} };
        }
        else {
            #FIXME Try to get info if the ensembl match is not in the main entry (check in all entries OR better query ???)
            warn "No info in the main entry [$st_gene_id]\n";
            next GENE;
        }
    }
    else {
        warn "No info found for [$st_gene_id]\n";
        next GENE;
    }

    #warn Dumper $annotations  if ( $verbose );
    # Update gene table
    if ( $st_gene_name     eq '' && exists $annotations->{'prot_name'}     && $annotations->{'prot_name'}     ne '' ){
        $up_gene1->execute($annotations->{'prot_name'},     $st_own_id)  or die $up_gene1->errstr;
    }
    if ( $st_description   eq '' && exists $annotations->{'prot_desc'}     && $annotations->{'prot_desc'}     ne '' && $annotations->{'prot_desc'} ne 'Uncharacterized protein' ){
        $up_gene2->execute($annotations->{'prot_desc'},     $st_own_id)  or die $up_gene2->errstr;
    }

    # Update synonym table
    if ( scalar @{ $annotations->{'synonyms'} } > 0 ){
        for my $syn ( @{ $annotations->{'synonyms'} } ){
            $up_syn->execute($st_own_id, $syn);#  or die $up_syn->errstr;
        }
    }
    # Update GO table
    if ( scalar @{ $annotations->{'go'} } > 0 ){
        for my $go_full ( @{ $annotations->{'go'} } ){
            my ($go, $evidence) = split(/___/, $go_full);
            $up_go->execute($st_own_id, $go, $evidence || '');#  or die $up_go->errstr;
        }
    }
    # Update xref table
    if ( scalar @{ $annotations->{'xrefs'} } > 0 ){
        for my $xf ( @{ $annotations->{'xrefs'} } ){
            next  if ( $xf =~ /;/ );
            my ($db_source, $pid) = split(/:/, $xf, 2);
            $up_xref->execute($st_own_id, $pid, '', $InsertedDataSources{lc $db_source});#  or die $up_xref->errstr;
        }
    }

    # Update term tables
    my @all;
    push @all, $annotations->{'prot_name'}                         if ( $annotations->{'prot_name'} );
    push @all, @{ $annotations->{'synonyms'} }                     if ( scalar @{ $annotations->{'synonyms'} } > 0 );
    push @all, map { s/___.*$//; $_ } @{ $annotations->{'go'} }    if ( scalar @{ $annotations->{'go'} } > 0 );
    push @all, map { s/^.+?://; $_ }  @{ $annotations->{'xrefs'} } if ( scalar @{ $annotations->{'xrefs'} } > 0 );
    #without version digit
    push @all, map { s/^.+?://; $_ } map  { s/\.\d+$//; $_ } grep { /\.\d+$/ } @{ $annotations->{'xrefs'} } if ( scalar @{ $annotations->{'xrefs'} } > 0 );
    for my $term ( uniq sort @all ){
        $up_term->execute($st_own_id, $term);#  or die $up_term->errstr;
    }
}

$sth_genes->finish;
$up_gene1->finish;
$up_gene2->finish;
$up_xref->finish;
$up_syn->finish;
$up_term->finish;
$up_go->finish;
# Close db connections
$dbh->disconnect;

exit 0;

