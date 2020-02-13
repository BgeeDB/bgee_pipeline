#!/bin/bash

#SBATCH --partition=long
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/Sum_Calls_cellPopulation.out
#SBATCH --error=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/Sum_Calls_cellPopulation.err
#SBATCH --export=NONE
#SBATCH --job-name=Sum_Calls_cellPopulation

module load Bioinformatics/Software/vital-it

module add R/3.6.1


export SCRIPT_PATH=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export NEW_scRNASeq_sample_info=/data/unil/dee/dee_robi/scRNA-Seq_all_results/Quality_Control_Cell_Population/NEW_scRNASeq_sample_info.txt
export cells_folder=/data/unil/dee/dee_robi/scRNA-Seq_all_results/Kallisto_All_Cells/
export output_folder=/data/unil/dee/dee_robi/scRNA-Seq_all_results/Calls_present_Cell_Population/
export ratioValue=0.05
export ROUT=/data/unil/dee/dee_robi/scRNA-Seq_all_results/


R CMD BATCH --no-save --no-restore "--args NEW_scRNASeq_sample_info=\"$NEW_scRNASeq_sample_info\" cells_folder=\"$cells_folder\" output_folder=\"$output_folder\" ratioValue=\"$ratioValue\"" $SCRIPT_PATH/1Run/Sum_Calls_cellPopulation.R $ROUT/Sum_Calls_cellPopulation.Rout

