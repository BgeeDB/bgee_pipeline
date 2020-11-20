#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/QC_cellPop_bimodality_TB.out
#SBATCH --error=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/QC_cellPop_bimodality_TB.err
#SBATCH --export=NONE
#SBATCH --job-name=QC_cellPop_bimodality_TB
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;

export SCRIPT_PATH=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols
export scRNASeq_Info=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt
export folder_data=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/DROPLET_10X/QC_CellType_identification/
export output_folder=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/DROPLET_10X/QC_CellType_identification/
export ROUT=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/DROPLET_10X/

R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" folder_data=\"$folder_data\" output_folder=\"$output_folder\"" $SCRIPT_PATH/1Run/QC_cellPop_bimodality.R $ROUT/QC_cellPop_bimodality.Rout
