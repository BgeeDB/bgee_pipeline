#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

# mysql -u bgee -p -h bioinfo bgee_v15_2 -e "SELECT GROUP_CONCAT(geneId), bgeeGeneId, speciesId, geneName, geneDescription, COUNT(*) AS dupl_count FROM gene WHERE geneDescription != '' AND geneName != '' AND geneDescription NOT LIKE '%Source:RFAM;%' GROUP BY geneDescription, geneName, speciesId HAVING COUNT(*) > 1 ORDER BY dupl_count;"
# ./check_duplicated_genes_because_of_alt_assemblies.pl  -gene=ENSG00000204371  -bgee=user=bgee__pass=bgee__host=bioinfo.unil.ch__port=3306__name=bgee_v15_2
#TODO from "[Source:", get the source id to allow to merge the same Ensembl genes coming from different locations: "official" chromosome and its alternative assemblies
# https://bgee.atlassian.net/browse/BA-170

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;


# Define arguments & their default value
my ($gene, $bgee_connector) = ('', '');
my ($all) = (0);
my ($debug) = (0);
my %opts = ('gene=s' => \$gene,            # gene to check for duplicates
            'bgee=s' => \$bgee_connector,  # Bgee connector string
            'debug'  => \$debug,           # debug mode, do not insert/update in database
            'all'    => \$all,             # Return a global list of all potential duplicates
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || ($gene eq '' && $all == 0) || $bgee_connector eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -gene=ENSG00000204371  -bgee=\$(BGEECMD)
\t-gene     Gene to check for duplicates
\t-bgee     Bgee connector string
\t-debug    Debug mode, do not insert/update in database
\t-all      Return a global list of all potential duplicates
\n";
    exit 1;
}


# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


#TODO display also the field "seqRegionName" to help decide
# Return a global list of all potential duplicates
if ( $all ){
    my $allDB = $dbh->prepare("SELECT GROUP_CONCAT(geneId), bgeeGeneId, speciesId, geneName, geneDescription, COUNT(*) AS dupl_count FROM gene WHERE geneDescription != '' AND geneName != '' AND geneDescription NOT LIKE '%Source:RFAM;%' GROUP BY geneDescription, geneName, speciesId HAVING COUNT(*) > 1 ORDER BY dupl_count");
    $allDB->execute()  or die $allDB->errstr;
    for my $dupl ( @{ $allDB->fetchall_arrayref } ){
        print '[', join("]\t[", @{$dupl}), "]\n";
    }
    $allDB->finish;
    $dbh->disconnect;
    exit 0;
}


# Get gene info
my $geneDB = $dbh->prepare('SELECT geneName, geneDescription, geneBioTypeId, speciesId FROM gene WHERE geneId=?');
$geneDB->execute($gene)  or die $geneDB->errstr;
my ($geneName, $geneDescription, $geneBioTypeId, $speciesId) = $geneDB->fetchrow_array;
$geneDB->finish;


# Search for duplicates
my $duplicatesDB = $dbh->prepare('SELECT GROUP_CONCAT(geneId), COUNT(*) AS dupl_count, geneDescription, geneName FROM gene WHERE geneDescription=? AND geneName=? AND geneBioTypeId=? AND speciesId=? GROUP BY geneDescription, geneName, speciesId HAVING COUNT(*) > 1 ORDER BY dupl_count');
$duplicatesDB->execute($geneDescription, $geneName, $geneBioTypeId, $speciesId)  or die $duplicatesDB->errstr;
my ($dupl_geneIds, $dupl_count) = $duplicatesDB->fetchrow_array;
$duplicatesDB->finish;


# Are there duplicates?
my @duplicates = split(',', $dupl_geneIds || $gene);
if ( ! exists $duplicates[1] ){
    print "It does not look $gene has duplicates: ", '[', join("]\t[", ($dupl_geneIds || $gene), $geneDescription, $geneName), "]\n";
    $dbh->disconnect;
    exit 0;
}

print "It looks $gene has $dupl_count duplicates: ", '[', join("]\t[", $dupl_geneIds, $geneDescription, $geneName), "]\n";
#TODO => for species we are sure they match on alternative chromosome, simpler to remove genes on those alt chr!
#        So, mainly human, zebrahish and ???? With https://ftp.ensembl.org/pub/current_fasta/*/dna/*.dna.alt* file
#        human      EHMT2  gene on chr  *6* and Scaffold HSCHR*6*_MHC_QBL_CTG1
#        zebrafish  prune  gene on chr *16* and CHR_ALT_CTG*16*_1_20
#        It may happen for mouse, on genome patches!

#TODO Is there a cascade delete for genes from the gene table?
#     to delete all those gene xrefs, terms, ...
#TODO Before deleting, a check has to be done to see if there are some expression for the one(s) to delete!


#NOTE Take care that some duplicates may be due to bad assemblies, and may be solved in future assemblies.
#NOTE Take care that some duplicates may be linked to wrong gene models, and be splice variants. E.g. ENSMUSG00000007440 & ENSMUSG00000102206


# Close db connections
$dbh->disconnect;

exit 0;
