#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/Sum_RawCounts_cellPopulation.out
#SBATCH --error=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/Sum_RawCounts_cellPopulation.err
#SBATCH --export=NONE
#SBATCH --job-name=Sum_RawCounts_cellPopulation

module use /software/module/

module add R/3.6.1


export SCRIPT_PATH=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export NEW_scRNASeq_sample_info=/data/unil/dee/dee_robi/scRNA-Seq_all_results/Quality_Control_Cell_Population/NEW_scRNASeq_sample_info.txt
export cells_folder=/data/unil/dee/dee_robi/scRNA-Seq_all_results/Kallisto_All_Cells/
export output_folder=/data/unil/dee/dee_robi/scRNA-Seq_all_results/Sum_Raw_Counts_Cell_Population/
export ROUT=/data/unil/dee/dee_robi/scRNA-Seq_all_results/


R CMD BATCH --no-save --no-restore "--args NEW_scRNASeq_sample_info=\"$NEW_scRNASeq_sample_info\" cells_folder=\"$cells_folder\" output_folder=\"$output_folder\"" $SCRIPT_PATH/1Run/Sum_RawCounts_cellPopulation.R $ROUT/Sum_RawCounts_cellPopulation.Rout
