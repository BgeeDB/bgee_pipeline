#!/bin/bash

#SBATCH --partition=long
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/QC_cellPopulation.out
#SBATCH --error=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/QC_cellPopulation.err
#SBATCH --export=NONE
#SBATCH --job-name=QC_cellPopulation

module load Bioinformatics/Software/vital-it

module add R/3.6.1


export SCRIPT_PATH=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export scrna_seq_sample_info=../../generated_files/scRNA_Seq/scrna_seq_sample_info.tsv
export cells_folder=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNASeq_libraries/
export output_folder=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA-Seq_all_results/Quality_Control_Cell_Population/
export plot=yes
export ROUT=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA-Seq_all_results/

R CMD BATCH --no-save --no-restore "--args scrna_seq_sample_info=\"$scrna_seq_sample_info\" cells_folder=\"$cells_folder\" output_folder=\"$output_folder\" plot=\"$plot\"" $SCRIPT_PATH/1Run/QC_cellPopulation.R $ROUT/QC_cellPopulation.Rout
