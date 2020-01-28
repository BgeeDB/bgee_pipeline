#!/bin/bash

#SBATCH --partition=long
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/home/smoretti/QC_on_encrypted_FASTQ.out
#SBATCH --error=/home/smoretti/QC_on_encrypted_FASTQ.err
#SBATCH --export=NONE
#SBATCH --job-name=Download_new_RNASeq_lib


module load Bioinformatics/Software/vital-it

module add UHTS/Analysis/sratoolkit/2.9.6.1
module add Utility/aspera_connect/3.7.4.147727
module add UHTS/Quality_control/fastp/0.19.5
module add R/3.6.1


export SCRIPT_PATH=...
export ANNOTATION_FILE=...
export DONE_FILE=...

perl $SCRIPT_PATH/0Before/get_SRA.pl  $SCRIPT_PATH/$ANNOTATION_FILE  $SCRIPT_PATH/$DONE_FILE

