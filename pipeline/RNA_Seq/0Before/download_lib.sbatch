#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=3-00:00:00

#SBATCH --output=/users/smoretti/bgee_pipeline.git/pipeline/RNA_Seq/get_sra.out
#SBATCH --error=/users/smoretti/bgee_pipeline.git/pipeline/RNA_Seq/get_sra.err
#SBATCH --export=NONE
#SBATCH --job-name=Download_new_RNASeq_lib
#SBATCH --mail-user smoretti@unil.ch

module load gcc/9.3.0;    module load sratoolkit/2.10.9;
module load gcc/9.3.0;    module load fastp/0.20.1;
module load gcc/9.3.0;    module load r/4.0.4;
module load gcc/9.3.0;    module load perl/5.32.1;

export PERL5LIB=/work/FAC/FBM/DEE/mrobinso/bgee/PERL_LIBS/lib/perl5:$PERL5LIB
export SCRIPT_PATH=/users/smoretti/bgee_pipeline.git/pipeline/RNA_Seq
export ANNOTATION_FILE=../../generated_files/RNA_Seq/rna_seq_sample_info.txt
export DONE_FILE=../../generated_files/RNA_Seq/rna_seq_sample_downloaded.txt

/software/bin/time -vo $SCRIPT_PATH/get_sra.time  perl $SCRIPT_PATH/0Before/get_SRA.pl  $SCRIPT_PATH/$ANNOTATION_FILE  $SCRIPT_PATH/$DONE_FILE

