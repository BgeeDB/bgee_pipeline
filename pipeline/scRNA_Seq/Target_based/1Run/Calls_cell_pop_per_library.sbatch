#!/bin/bash

#SBATCH --partition=urblauna
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=20:00:00

#SBATCH --output=/users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Calls_library_cellPop.out
#SBATCH --error=/users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Calls_library_cellPop.err
#SBATCH --job-name=calls

module use /software/module/

module load R/3.6.1;

export SCRIPT_PATH=/users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/
export scRNASeq_Info=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt
export bimodalityFile=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_1/scRNA-Seq_all_results_bgee_v15_1/DROPLET_10X/QC_CellType_identification//bimodality_targetBased.txt
export folder_data=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_1/scRNA-Seq_all_results_bgee_v15_1/DROPLET_10X/QC_CellType_identification/
export folder_gtf=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/rna_seq/GTF_15/
export desired_pValue_cutoff=0.05
export output_folder=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_1/scRNA-Seq_all_results_bgee_v15_1/DROPLET_10X/calls/
export ROUT=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_1/scRNA-Seq_all_results_bgee_v15_1/DROPLET_10X/

R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" bimodalityFile=\"$bimodalityFile\" folder_data=\"$folder_data\" folder_gtf=\"$folder_gtf\" desired_pValue_cutoff=\"$desired_pValue_cutoff\" output_folder=\"$output_folder\"" $SCRIPT_PATH/1Run/Calls_cell_pop_per_library.R $ROUT/Calls_cell_pop_per_library.Rout
