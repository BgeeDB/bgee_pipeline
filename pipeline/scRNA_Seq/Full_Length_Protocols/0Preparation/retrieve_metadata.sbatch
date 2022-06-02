#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=2-00:00:00

#SBATCH --output=/users/jwollbre/GIT/pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/retrieve_metadata.out
#SBATCH --error=/users/jwollbre/GIT/pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/retrieve_metadata.err
#SBATCH --export=NONE
#SBATCH --job-name=retrieve_metadata_information_fulllength

module load Bioinformatics/Software/vital-it

module add R/3.6.1


export SCRIPT_PATH=/users/jwollbre/GIT/pipeline/pipeline/scRNA_Seq/Full_Length_Protocols
export pass_annotationControl=../../../generated_files/scRNA_Seq/Full_length/passScRNASeqLibrary.tsv
export output_folder=../../../generated_files/scRNA_Seq/Full_length/
export ROUT=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/scRNASeq_libraries_FullLength/

mkdir -p $ROUT
R CMD BATCH --no-save --no-restore "--args pass_annotationControl=\"$pass_annotationControl\" output_folder=\"$output_folder\"" $SCRIPT_PATH/0Preparation/retrieve_metadata.R $ROUT/retrieve_metadata.Rout
