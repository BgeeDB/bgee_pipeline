#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Calls_library_cellPop.out
#SBATCH --error=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Calls_library_cellPop.err
#SBATCH --export=NONE
#SBATCH --job-name=Calls_library_cellPop
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;


export SCRIPT_PATH=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols
export scRNASeq_Info=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt
export InformationAllLibraries=../../../source_files/scRNA_Seq/InformationAllLibraries.txt
export folder_data=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/FASTQ/scRNAseq/10X/
export folder_refIntergenic=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/Rerence_intergenic_folder/
export desired_pValue_cutoff=0.05
export output_folder=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/DROPLET_10X/Calls_IndCell_CellPopulation/
export ROUT=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/DROPLET_10X/


R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" InformationAllLibraries=\"$InformationAllLibraries\" folder_data=\"$folder_data\" folder_refIntergenic=\"$folder_refIntergenic\" desired_pValue_cutoff=\"$desired_pValue_cutoff\" output_folder=\"$output_folder\"" $SCRIPT_PATH/1Run/Calls_cell_pop_per_library.R $ROUT/Calls_cell_pop_per_library.Rout
