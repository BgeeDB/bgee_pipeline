#!/bin/bash

#SBATCH --partition=wally
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=10-00:00:00

#SBATCH --output=/users/smoretti/bgee_pipeline.git/pipeline/RNA_Seq/get_sra.out
#SBATCH --error=/users/smoretti/bgee_pipeline.git/pipeline/RNA_Seq/get_sra.err
#SBATCH --export=NONE
#SBATCH --job-name=Download_new_RNASeq_lib


module add Bioinformatics/Software/vital-it; module add UHTS/Analysis/sratoolkit/2.9.6.1;
module add Bioinformatics/Software/vital-it; module add Utility/aspera_connect/3.7.4.147727;
module add Bioinformatics/Software/vital-it; module add UHTS/Quality_control/fastp/0.19.5;
module add Bioinformatics/Software/vital-it; module add R/3.6.1;


export SCRIPT_PATH=/users/smoretti/bgee_pipeline.git/pipeline/RNA_Seq
export ANNOTATION_FILE=../../generated_files/RNA_Seq/rna_seq_sample_info.txt
export DONE_FILE=../../generated_files/RNA_Seq/rna_seq_sample_downloaded.txt

/software/bin/time -vo $SCRIPT_PATH/get_sra.time  perl $SCRIPT_PATH/0Before/get_SRA.pl  $SCRIPT_PATH/$ANNOTATION_FILE  $SCRIPT_PATH/$DONE_FILE

