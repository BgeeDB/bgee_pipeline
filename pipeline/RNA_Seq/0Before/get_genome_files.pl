#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use FindBin;
use File::Basename;
use LWP::Simple;

#NOTE genome and indexed_genome files MUST be at the same location!

# Define arguments & their default value
my ($GTF_dir, $outDir)               = ('', '');
my ($ensRelease, $ensMetazoaRelease) = (0, 0);
my %opts = ('GTF_dir=s'             => \$GTF_dir,
            'ensRelease=i'          => \$ensRelease,
            'ensMetazoaRelease=i'   => \$ensMetazoaRelease,
            'outDir=s'              => \$outDir,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $GTF_dir eq '' || $ensRelease == 0 || $ensMetazoaRelease == 0 || $outDir eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -GTF_dir=\$(RNASEQ_DOWNLOAD_GTF) -ensRelease=\$(ENSRELEASE) -ensMetazoaRelease=\$(ENSMETAZOARELEASE) -outDir=\$(RNASEQ_DOWNLOAD_GTF)
\t-GTF_dir              GTF folder
\t-ensRelease           Ensembl Release number
\t-ensMetazoaRelease    Ensembl Metazoa Release number
\t-outDir               Output Directory to store genome/indexed_genome files
\n";
    exit 1;
}

die "Invalid or missing [$outDir]: $?\n"  if ( !-d $outDir || !-w $outDir );


chdir $outDir; # Go into output/genome folder
for my $gtf (glob($GTF_dir."/*.gtf.gz") ){
    die "Problem with GTF files [$gtf]\n"  if ( -z "$gtf" );

    my $species_name = basename($gtf);
    $species_name   =~ s{^(\w+).+$}{$1};
    $species_name    = lc $species_name;

    my $prefix   = basename($gtf);
    $prefix     =~ s{\.gtf.gz$}{};
    next  if ( -e "$prefix.genome.fa" && -s "$prefix.genome.fa" );
    if ( is_success( getstore("ftp://ftp.ensembl.org/pub/release-$ensRelease/fasta/$species_name/dna/$prefix.dna.toplevel.fa.gz", "$prefix.genome.fa.gz") ) ){
        # exists in Ensembl FTP && genome file downloaded
        # uncompress downloaded fasta file
        system("gunzip -f $prefix.genome.fa.gz")==0  or die "Failed to unzip genome files: $?\n";

    }
    #FIXME No more $ensMetazoaRelease in dna file name for Ensembl Metazoa 49, but is in gtf file name!
    elsif ( is_success( getstore("ftp://ftp.ensemblgenomes.org/pub/release-$ensMetazoaRelease/metazoa/fasta/$species_name/dna/$prefix.$ensMetazoaRelease.dna.toplevel.fa.gz", "$prefix.genome.fa.gz") ) ||
            is_success( getstore("ftp://ftp.ensemblgenomes.org/pub/release-$ensMetazoaRelease/metazoa/fasta/$species_name/dna/$prefix.dna.toplevel.fa.gz", "$prefix.genome.fa.gz") ) ){
        # exists in Ensembl Metazoa FTP && genome file downloaded
        # uncompress downloaded fasta file
        system("gunzip -f $prefix.genome.fa.gz")==0  or die "Failed to unzip genome files: $?\n";

    }
    elsif ( $prefix =~ /^\w+?_((GC[FA])_(\d\d\d)(\d\d\d)(\d\d\d).*)$/ ){
        # From NCBI, RefSeq or GenBank assembly annotations
        #See https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/ for help
        # e.g. macaca_fuscata_GCA_003118495.1_macFus_1.0
        #      manis_javanica_GCF_014570535.1_YNU_ManJav_2.0
        if ( is_success( getstore("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/$2/$3/$4/$5/$1/${1}_genomic.fna.gz", basename("$prefix.genome.fa.gz")) ) ){
            # exists in NCBI FTP && genome file downloaded
            # uncompress downloaded fasta file
            system("gunzip -f $prefix.genome.fa.gz")==0  or die "Failed to unzip genome files: $?\n";
        }
    }
    else {
        warn "No genome file found for [$species_name] in Ensembl $ensRelease, Ensembl Metazoa $ensMetazoaRelease or NCBI\n";
    }
}

exit 0;

