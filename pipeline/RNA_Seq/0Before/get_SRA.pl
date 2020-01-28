#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

use File::Slurp;

my $annotation_file = $ARGV[0]  // die "\n\t$0 <RNASeq library annotation file (rna_seq_sample_info.txt)> <RNASeq library already downloaded (rna_seq_sample_downloaded.txt)>\n\n";
my $downloaded_lib  = $ARGV[1]  // die "\n\t$0 <RNASeq library annotation file (rna_seq_sample_info.txt)> <RNASeq library already downloaded (rna_seq_sample_downloaded.txt)>\n\n";


my $SRATK_PATH         = $ENV{'SRATK_PATH'};
my $ASPERA_CONNECT_DIR = $ENV{'ASPERA_CONNECT_DIR'};

my $BASE               = "/scratch/wally/FAC/FBM/DEE/mrobinso/bgee";
my $SRA_PATH           = $BASE.'/sra';
my $FASTQ_PATH         = $BASE.'/FASTQ/RNAseq';


# Private experimentId to store encrypted
my @private_exp_id     = ('SRP012682'); # E.g. GTEx


if ( !$SRATK_PATH || !$ASPERA_CONNECT_DIR ){
    die "\n\tSRATK_PATH and/or ASPERA_CONNECT_DIR not defined! They are used to find NCBI SRA & ASPERA executables\n\n";
}


# Read already downloaded libraries
my %already_downloaded = read_file("$downloaded_lib", chomp=>1);
print scalar keys %already_downloaded, "\n";
exit;


open(my $ANNOTATION, '<', "$annotation_file")  or die "\n\tCannot read/open [$annotation_file]\n\n";
#libraryId   experimentId   speciesId   organism        genomeFilePath                        database   platform                       libraryType   libraryInfo   readLength   runIds
#SRX081869   GSE30352       9031        Gallus gallus   gallus_gallus/Gallus_gallus.Galgal4   Ensembl    Illumina Genome Analyzer IIx   SINGLE                      76           SRR306710
LIB:
while (<$ANNOTATION>){
    chomp $_;
    my @tmp = map { s{^"|"$}{}g; $_ } # remove quotes
              split(/\t/, $_);
    my $sra_list   = $tmp[10];
    my $library_id = $tmp[0];
    my $exp_id     = $tmp[1];
    next LIB  if ( $library_id =~ /^#/ || $sra_list =~ /^#/ ); # Header or commented line
    next LIB  if ( exists $already_downloaded{$library_id} );

    SRA:
    for my $sra_id ( sort split(/,/, $sra_list) ){
        if ( $sra_id =~ /^[SEDC]RR\d+/ ){ #S: SRA/NCBI; E: EBI; D: DDBJ; C: GSA_China
            # Check if FASTQ are there AND if SRA have been removed  -> already stored
            if ( !-e "$SRA_PATH/$sra_id.sra" && (-s "$FASTQ_PATH/$library_id/$sra_id.fastq.gz" || -s "$FASTQ_PATH/$library_id/$sra_id.fastq.gz.enc") ){
                warn "\t[$library_id/$sra_id] single-end already stored\n";
                next SRA;
            }
            if ( !-e "$SRA_PATH/$sra_id.sra" && (-s "$FASTQ_PATH/$library_id/${sra_id}_1.fastq.gz" || -s "$FASTQ_PATH/$library_id/${sra_id}_1.fastq.gz.enc")
                                             && (-s "$FASTQ_PATH/$library_id/${sra_id}_2.fastq.gz" || -s "$FASTQ_PATH/$library_id/${sra_id}_2.fastq.gz.enc")){
                warn "\t[$library_id/$sra_id] paired-end already stored\n";
                next SRA;
            }

            # Run prefetch to get SRA file
            #NOTE prefetch automatically checks what has already been downloaded and completes if needed
            #NOTE cd to the "project's workspace directory" the ONLY place where the SRA download works for private SRR
            system("cd $BASE; $SRATK_PATH/bin/prefetch --quiet --max-size 500G -t ascp -a \"$ASPERA_CONNECT_DIR/bin/ascp|$ASPERA_CONNECT_DIR/etc/asperaweb_id_dsa.openssh\" $sra_id")==0
                or do { warn "\tFailed to get [$library_id/$sra_id]\n"; next SRA };

            # Convert in FastQ
            #   Transforming ".sra" file into ".fastq" file using fastq-dump software, the option --split-3 cause splitting paired reads into separate files '_1.fastq' and '_2.fastq',
            #   output files go to sample subfolder in fastq folder
            #NOTE fastq-dump has problems with paths containing //
            #NOTE cd to the "project's workspace directory" the ONLY place where the SRA download works for private SRR
            mkdir "$FASTQ_PATH/$library_id";
            system("cd $BASE; $SRATK_PATH/bin/fastq-dump --split-3 --gzip --outdir $FASTQ_PATH/$library_id/  $SRA_PATH/$sra_id.sra")==0
                or do { warn "\tFailed to convert [$library_id/$sra_id]\n"; next SRA };

            # Check FastQ file size
            for my $fastq ( glob("$FASTQ_PATH/$library_id/*.gz") ){
                if ( -s $fastq < 1_000_000 ){
                    warn "$fastq file size looks very small!";
                }
            }

            # Run FastP (A quality control tool for high throughput sequence data) for ALL SRR (runs)
            # as well as basic read length statistics with R
            QC:
            for my $fastq ( glob("$FASTQ_PATH/$library_id/*.gz") ){
                my ($run_id) = $fastq =~ /([^\/]+)\.fastq\.gz/;
                system("fastp -i $fastq --json $FASTQ_PATH/$library_id/${run_id}.fastp.json --html $FASTQ_PATH/$library_id/${run_id}.fastp.html --thread 2  > $FASTQ_PATH/$library_id/${run_id}.fastp.log 2>&1")==0
                    or do { warn "\tfastp failed for [$FASTQ_PATH/$library_id/$run_id]\n"; next QC };
                system("xz -9 $FASTQ_PATH/$library_id/${run_id}.fastp.html $FASTQ_PATH/$library_id/${run_id}.fastp.json");
                #TODO Would be nice to have all basic stats from FastP
                # Basic read length stats
                system("echo -e \"#min\tmax\tmedian\tmean\" > $FASTQ_PATH/$library_id/${run_id}.R.stat");
                system("zcat $fastq | sed -n '2~4p' | awk '{print length($0)}' | Rscript -e 'd<-scan(\"stdin\", quiet=TRUE);cat(min(d), max(d), median(d), mean(d), sep=\"\t\");cat(\"\n\")' >> $FASTQ_PATH/$library_id/${run_id}.R.stat");
            }

            # If private (need encryption):
            if ( (scalar grep { /^$exp_id$/ } @private_exp_id) >= 1 ){
                for my $fastq ( glob("$FASTQ_PATH/$library_id/*.gz") ){
                    #NOTE Replace -salt by -d for decrypting and gz.enc as input and gz as output
                    system("openssl enc -aes-128-cbc -salt -in $fastq -out $fastq.enc -pass file:$BASE/.passw  &&  rm -f $fastq")==0
                        or do { warn "\tFailed to encrypt [$library_id/$sra_id]\n"; next SRA };
                }
            }

            # If fine, SRA cleaning
            unlink "$SRA_PATH/$sra_id.sra";
        }
        else {
            warn "\t[$sra_id] is not an SRA id\n";
        }
    }
}
close $ANNOTATION;

exit 0;

