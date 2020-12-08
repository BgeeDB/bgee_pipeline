#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use FindBin;
use List::MoreUtils qw{uniq};
use File::Basename;
use File::Path qw(make_path);
use LWP::Simple;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

# Define arguments & their default value
my ($RNAseqSample, $outDir)          = ('', '');
my ($ensRelease, $ensMetazoaRelease) = (0, 0);
my %opts = ('RNAseqSample=s'        => \$RNAseqSample,
            'ensRelease=i'          => \$ensRelease,
            'ensMetazoaRelease=i'   => \$ensMetazoaRelease,
            'outDir=s'              => \$outDir,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $RNAseqSample eq '' || $ensRelease == 0 || $ensMetazoaRelease == 0 || $outDir eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -RNAseqSample=\$(RNASEQ_SAMPINFO_FILEPATH) -ensRelease=\$(ENSRELEASE) -ensMetazoaRelease=\$(ENSMETAZOARELEASE) -outDir=\$(RNASEQGTFDATAPATH)
\t-RNAseqSample         rnaseq_sample_info.txt file
\t-ensRelease           Ensembl Release number
\t-ensMetazoaRelease    Ensembl Metazoa Release number
\t-outDir               Output Directory to store GTF files
\n";
    exit 1;
}

if ( !-e $outDir ){
    make_path($outDir);
}
die "Invalid or missing [$outDir]: $?\n"  if ( !-d $outDir || !-w $outDir );


# Read RNAseqLib file & get GTF files for these species
my %tsv = %{ Utils::read_spreadsheet("$RNAseqSample", "\t", 'csv', '', 1) };

chdir $outDir; # Go into GTF folder
for my $genomeFilePath ( uniq sort @{$tsv{'genomeFilePath'}} ){
    next  if ( -e basename("$genomeFilePath.gtf.gz") && -s basename("$genomeFilePath.gtf.gz") );

    if ( is_success( getstore("ftp://ftp.ensembl.org/pub/release-$ensRelease/gtf/$genomeFilePath.$ensRelease.gtf.gz", basename("$genomeFilePath.gtf.gz")) ) ){
        # genomeFilePath exists in Ensembl FTP && GTF file downloaded
    }
    elsif ( is_success( getstore("ftp://ftp.ensemblgenomes.org/pub/release-$ensMetazoaRelease/metazoa/gtf/$genomeFilePath.$ensMetazoaRelease.gtf.gz", basename("$genomeFilePath.gtf.gz")) ) ){
        # genomeFilePath exists in Ensembl Metazoa FTP && GTF file downloaded
    }
    elsif ( $genomeFilePath =~ /^\w+\/\w+?_((GC[FA])_(\d\d\d)(\d\d\d)(\d\d\d).*)$/ ){
        # From NCBI, RefSeq or GenBank assembly annotations
        #See https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/ for help
        # e.g. macaca_fuscata/macaca_fuscata_GCA_003118495.1_macFus_1.0
        #      manis_javanica/manis_javanica_GCF_014570535.1_YNU_ManJav_2.0
        if ( is_success( getstore("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/$2/$3/$4/$5/${1}_genomic.gtf.gz", basename("$genomeFilePath.gtf.gz")) ) ){
            # genomeFilePath exists in NCBI FTP && GTF file downloaded
        }
    }
    else {
        #warn "\tftp://ftp.ensembl.org/pub/release-$ensRelease/gtf/$genomeFilePath.gtf.gz\n";
        warn "No GTF file found for [$genomeFilePath] in Ensembl $ensRelease or Ensembl Metazoa $ensMetazoaRelease\n";
    }
}

exit 0;

