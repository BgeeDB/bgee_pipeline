#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=1-00:00:00

#SBATCH --output=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/QC_cellPop_bimodality_TB.out
#SBATCH --error=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/QC_cellPop_bimodality_TB.err
#SBATCH --export=NONE
#SBATCH --job-name=QC_cellPop_bimodality_TB
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module use /software/module/
module add R/3.6.1;

export SCRIPT_PATH=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/
export scRNASeq_Info=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt
export folder_data=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/QC_CellType_identification/
export output_folder=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/QC_CellType_identification/
export ROUT=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/

R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" folder_data=\"$folder_data\" output_folder=\"$output_folder\"" $SCRIPT_PATH/1Run/QC_cellPop_bimodality.R $ROUT/QC_cellPop_bimodality.Rout
