#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use FindBin;
use List::MoreUtils qw{uniq};
use File::Basename;
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

die "Invalid or missing [$outDir]: $?\n"  if ( !-d $outDir || !-w $outDir );


# Read RNAseqLib file & get GTF files for these species
my %tsv = %{ Utils::read_spreadsheet("$RNAseqSample", "\t", 'csv', '', 1) };
chdir $outDir; # Go into GTF folder
for my $genomeFilePath ( uniq sort @{$tsv{'genomeFilePath'}} ){
    next  if ( -e basename("$genomeFilePath.$ensRelease.gtf.gz")        && -s basename("$genomeFilePath.$ensRelease.gtf.gz") );
    next  if ( -e basename("$genomeFilePath.$ensMetazoaRelease.gtf.gz") && -s basename("$genomeFilePath.$ensMetazoaRelease.gtf.gz") );

    if ( is_success( getstore("ftp://ftp.ensembl.org/pub/release-$ensRelease/gtf/$genomeFilePath.$ensRelease.gtf.gz", basename("$genomeFilePath.$ensRelease.gtf.gz")) ) ){
        # genomeFilePath exists in Ensembl FTP && GTF file downloaded
    }
    elsif ( is_success( getstore("ftp://ftp.ensemblgenomes.org/pub/release-$ensMetazoaRelease/metazoa/gtf/$genomeFilePath.$ensMetazoaRelease.gtf.gz", basename("$genomeFilePath.$ensMetazoaRelease.gtf.gz")) ) ){
        # genomeFilePath exists in Ensembl Metazoa FTP && GTF file downloaded
    }
    else {
        #warn "\tftp://ftp.ensembl.org/pub/release-$ensRelease/gtf/$genomeFilePath.gtf.gz\n";
        warn "No GTF file found for [$genomeFilePath] in Ensembl $ensRelease or Ensembl Metazoa $ensMetazoaRelease\n";
    }
}

exit 0;

