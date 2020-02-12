#!/bin/bash

#SBATCH --partition=long
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/retrieve_metadata.out
#SBATCH --error=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/retrieve_metadata.err
#SBATCH --export=NONE
#SBATCH --job-name=retrieve_metadata_information

module load Bioinformatics/Software/vital-it

module add R/3.6.1


export SCRIPT_PATH=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export NEW_scRNASeqLibrary=../../generated_files/scRNA_Seq/NEW_scRNASeqLibrary.tsv
export output_folder=../../generated_files/scRNA_Seq/
export ROUT=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA-Seq_all_results/

R CMD BATCH --no-save --no-restore "--args NEW_scRNASeqLibrary=\"$NEW_scRNASeqLibrary\" output_folder=\"$output_folder\"" $SCRIPT_PATH/0Preparation/retrieve_metadata.R $ROUT/retrieve_metadata.Rout
