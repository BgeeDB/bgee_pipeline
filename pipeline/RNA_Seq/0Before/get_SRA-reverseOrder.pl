#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

my $annotation_file = $ARGV[0]  // die "\n\t$0 <RNASeq library annotation file (rna_seq_sample_info.txt)> <odd|even>\n\n";
my $oddeven         = $ARGV[1]  // 'odd';


my $SRATK_PATH         = '/software/UHTS/Analysis/sratoolkit/2.5.2';
my $ASPERA_CONNECT_DIR = '/software/Utility/aspera_connect/3.6.1.110647';

my $BASE               = '/scratch/local/weekly/bbgee'; # == project's workspace directory
my $SRA_PATH           = $BASE.'/sra';
my $FASTQ_PATH         = $BASE.'/FASTQ';
my $REMOTE_BASE        = '/opt/gtexfile/FASTQ';
my $REMOTE_CONNECT     = 'adminbgee@bigbgee.unil.ch';


# Private experimentId to store encrypted
my @private_exp_id     = ('SRP012682'); # E.g. GTEx


open(my $REVERSE_ANNOTATION, "tac $annotation_file |")  or die "\n\tCannot read/open [$annotation_file]\n\n";
#libraryId   experimentId   speciesId   organism        genomeFilePath                        database   platform                       libraryType   libraryInfo   readLength   runIds
#SRX081869   GSE30352       9031        Gallus gallus   gallus_gallus/Gallus_gallus.Galgal4   Ensembl    Illumina Genome Analyzer IIx   SINGLE                      76           SRR306710
LIB:
while (<$REVERSE_ANNOTATION>){
    chomp $_;
    my @tmp = map { s{^"|"$}{}g; $_ } # remove quotes
              split(/\t/, $_);
    my $sra_list   = $tmp[10];
    my $library_id = $tmp[0];
    my $exp_id     = $tmp[1];
    next LIB  if ( $library_id =~ /^#/ || $sra_list =~ /^#/ ); # Header or commented line

    # Odd or even library_id only
    my ($digit) = $library_id =~ /(\d+)$/;
    if (    $oddeven eq 'odd'  && ($digit % 2)!=1 ){
        next;
    }
    elsif ( $oddeven eq 'even' && ($digit % 2)!=0 ){
        next;
    }

    # Already available remotely?
    if ( `ssh $REMOTE_CONNECT ls -1 $REMOTE_BASE/$library_id/` ){
        warn "\t[$library_id] already on bigbgee!\n";
        next;
    }

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

            # Compute FastQC (A quality control tool for high throughput sequence data) for ALL SRR (runs)
            mkdir "$FASTQ_PATH/$library_id/FASTQC";
            #TODO Test this step!
            QC:
            for my $fastq ( glob("$FASTQ_PATH/$library_id/*.gz") ){
                my ($run_id) = $fastq =~ /([^\/]+)\.fastq\.gz/;
                system("fastqc -o $FASTQ_PATH/$library_id/FASTQC $fastq > $FASTQ_PATH/$library_id/FASTQC/$run_id.fastqc.log 2>&1")==0
                    or do { warn "\tfastqc failed for [$FASTQ_PATH/$library_id/FASTQC/$run_id]\n"; next QC };
            }

            # If private (need encryption):
            if ( (scalar grep { /^$exp_id$/ } @private_exp_id) >= 1 ){
                for my $fastq ( glob("$FASTQ_PATH/$library_id/*.gz") ){
                    #NOTE Replace -salt by -d for decrypting and gz.enc as input and gz as output
                    system("openssl enc -aes-128-cbc -salt -in $fastq -out $fastq.enc -pass file:/home/bbgee/.passw  &&  rm -f $fastq")==0
                        or do { warn "\tFailed to encrypt [$library_id/$sra_id]\n"; next SRA };
                }
            }

            # Remote copy
            system("ssh $REMOTE_CONNECT mkdir -p $REMOTE_BASE/$library_id/")==0                      or die "\tCannot create remote FASTQ/$library_id/ folder\n";
            system("scp $FASTQ_PATH/$library_id/*  ${REMOTE_CONNECT}:$REMOTE_BASE/$library_id/")==0  or die "\tCannot scp [$library_id] to remote\n";


            # If fine, SRA cleaning
            unlink "$SRA_PATH/$sra_id.sra";
            system("rm -Rf $FASTQ_PATH/$library_id/");
        }
        else {
            warn "\t[$sra_id] is not an SRA id\n";
        }
    }
}
close $REVERSE_ANNOTATION;

exit 0;

