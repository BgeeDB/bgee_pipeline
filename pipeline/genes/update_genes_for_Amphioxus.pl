#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use List::MoreUtils qw{uniq};
use Data::Dumper;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my ($geneid, $bgee_connector) = ('', '');
my ($annotations)             = ('');
my ($verbose)                 = (0);
my %opts = ('annot=s'  => \$annotations,    # Amphioxus annotation file
            'gene=s'   => \$geneid,         # gene id to update in the db
            'verbose'  => \$verbose,        # verbose mode
            'bgee=s'   => \$bgee_connector, # Bgee connector string
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $annotations eq '' || !-r $annotations ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -annot=Bla_annot_final.gtf.gz  -bgee=\$(BGEECMD)
\t     OR
\t     $0  -annot=Bla_annot_final.gtf.gz  -bgee=\$(BGEECMD)  -gene=BL00000
\t-annot    Amphioxus annotation file
\t-gene     gene id to update in the db
\t-bgee     Bgee connector string
\t-verbose  verbose mode
\n";
    exit 1;
}


# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);
my $sel_id   = $dbh->prepare('SELECT bgeeGeneId FROM gene WHERE geneId=?');
my $up_xref  = $dbh->prepare('INSERT INTO geneXRef               (bgeeGeneId, XRefId, XRefName, dataSourceId) VALUES (?, ?, ?, ?)');
#my $up_syn   = $dbh->prepare('INSERT INTO geneNameSynonym        (bgeeGeneId, geneNameSynonym)                VALUES (?, ?)');
my $up_term  = $dbh->prepare('INSERT INTO geneToTerm             (bgeeGeneId, term)                           VALUES (?, ?)');
my $up_gene1 = $dbh->prepare('UPDATE gene SET geneName=?         WHERE bgeeGeneId=?');
#my $up_gene2 = $dbh->prepare('UPDATE gene SET geneDescription=?  WHERE bgeeGeneId=?');


# Read annotation file
#comes from https://amphiencode.github.io/Data/   -> Bla_annot_final.gtf
my $gene;
open(my $ANNOT, "zcat $annotations |")  or die ;
while(<$ANNOT>){
    next  if ( $_ =~ /^#/ );

    #Sc0000095	protein_coding	exon	778337	778580	.	-	.	gene_id "BL00000"; transcript_id "BL00000_evm0"; exon_number "1"; status "both"; oldID "Blg10491.0"; gene_name "FAM13A";
    my ($prot_id, undef, undef, undef, undef, undef, undef, undef, $rest) = split(/\t/, $_);
    my ($gene_id, $transcript_id, $oldID, $gene_name) = ('', '', '', '');
    for my $an ( split(/;\s*/, $rest) ){
        if ( $an =~ /^gene_id "(BL\d+)"$/ ){
            $gene_id = $1;
        }
        elsif ( $an =~ /^transcript_id "(.+?)"$/ ){
            $transcript_id = $1;
        }
        elsif ( $an =~ /^oldID "(.+?)"$/ ){
            $oldID = $1;
        }
        elsif ( $an =~ /^gene_name "(.+?)"$/ ){
            $gene_name = $1;
        }
    }

    next  if ( $gene_id eq '' );
    $gene->{$gene_id}->{'name'} = $gene_name;
    push @{ $gene->{$gene_id}->{'xrefs'} }, $prot_id, $transcript_id, $oldID;
}
close $ANNOT;


# Update db
for my $gene_id ( sort keys %$gene ){
    if ( $geneid ne '' ){
        next  if ( $gene_id ne $geneid );
    }

    # Get bgeeGeneId
    $sel_id->execute($gene_id) or die $sel_id->errstr;
    my @id = @{$sel_id->fetchall_arrayref};
    if ( scalar @id != 1 ){
        warn "Several or no bgeeGeneId for [$gene_id]\n";
        next;
    }

    # Update gene table
    if ( $id[0][0] && exists $gene->{$gene_id}->{'name'} && $gene->{$gene_id}->{'name'} ne '' ){
        $up_gene1->execute($gene->{$gene_id}->{'name'}, $id[0][0])  or die $up_gene1->errstr;
    }
    # Update xref table
    @{ $gene->{$gene_id}->{'xrefs'} } = uniq sort @{ $gene->{$gene_id}->{'xrefs'} };
    if ( scalar @{ $gene->{$gene_id}->{'xrefs'} } > 0 ){
        for my $xf ( @{ $gene->{$gene_id}->{'xrefs'} } ){
            # 39 is the dataSourceId for Amphiencode 
            $up_xref->execute($id[0][0], $xf, '', 39);#  or die $up_xref->errstr;
        }
    }

    # Update term tables
    my @all;
    push @all, $gene->{$gene_id}->{'name'}        if ( $gene->{$gene_id}->{'name'} ne '' );
    push @all, @{ $gene->{$gene_id}->{'xrefs'} }  if ( scalar @{ $gene->{$gene_id}->{'xrefs'} } > 0 );
    #without version digit
    push @all, map  { s/\.\d+$//; $_ } grep { /\.\d+$/ } @{ $gene->{$gene_id}->{'xrefs'} }  if ( scalar @{ $gene->{$gene_id}->{'xrefs'} } > 0 );
    for my $term ( uniq sort @all ){
        $up_term->execute($id[0][0], $term);#  or die $up_term->errstr;
    }
}

$sel_id->finish;
$up_gene1->finish;
#$up_gene2->finish;
$up_xref->finish;
#$up_syn->finish;
$up_term->finish;
# Close db connections
$dbh->disconnect;

exit 0;

