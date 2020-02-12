#!/bin/bash

#SBATCH --partition=long
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/Sum_RawCounts_cellPopulation.out
#SBATCH --error=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/Sum_RawCounts_cellPopulation.err
#SBATCH --export=NONE
#SBATCH --job-name=Sum_RawCounts_cellPopulation

module load Bioinformatics/Software/vital-it

module add R/3.6.1


export SCRIPT_PATH=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export NEW_scRNASeq_sample_info=../../generated_files/scRNA_Seq/NEW_scRNASeq_sample_info.txt
export cells_folder=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNASeq_libraries/
export output_folder=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA-Seq_all_results/Sum_Raw_Counts_Cell_Population/
export ROUT=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA-Seq_all_results/


R CMD BATCH --no-save --no-restore "--args NEW_scRNASeq_sample_info=\"$NEW_scRNASeq_sample_info\" cells_folder=\"$cells_folder\" output_folder=\"$output_folder\"" $SCRIPT_PATH/1Run/Sum_RawCounts_cellPopulation.R $ROUT/Sum_RawCounts_cellPopulation.Rout
