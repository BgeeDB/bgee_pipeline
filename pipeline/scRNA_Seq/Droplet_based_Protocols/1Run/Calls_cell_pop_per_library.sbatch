#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=4-00:00:00

#SBATCH --output=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Calls_library_cellPop.out
#SBATCH --error=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Calls_library_cellPop.err
#SBATCH --export=NONE
#SBATCH --job-name=Calls_cell_pop_per_library
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;

export SCRIPT_PATH=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/
export scRNASeq_Info=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt
export bimodalityFile=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/QC_CellType_identification//bimodality_targetBased.txt
export folder_data=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/QC_CellType_identification/
export folder_refIntergenic=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/ref_intergenic_bgee_v15_0/
export desired_pValue_cutoff=0.05
export output_folder=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/calls/
export ROUT=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/

R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" bimodalityFile=\"$bimodalityFile\" folder_data=\"$folder_data\" folder_refIntergenic=\"$folder_refIntergenic\" desired_pValue_cutoff=\"$desired_pValue_cutoff\" output_folder=\"$output_folder\"" $SCRIPT_PATH/1Run/Calls_cell_pop_per_library.R $ROUT/Calls_cell_pop_per_library.Rout
