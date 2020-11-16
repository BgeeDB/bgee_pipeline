#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Call_PresentGenes_indivCell_cellPop.out
#SBATCH --error=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Call_PresentGenes_indivCell_cellPop.err
#SBATCH --export=NONE
#SBATCH --job-name=Call_PresentGenes_indivCell_cellPop
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;

export SCRIPT_PATH=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols
export scRNASeq_Info=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt
export InformationAllLibraries=/data/unil/dee/dee_robi/scRNA-Seq_all_results/DROPLET_10X/QC_CellType_identification/InformationAllLibraries.txt
export folder_data=/data/unil/dee/dee_robi/FASTQ/scRNAseq/10X/
export desired_r_cutoff=0.05
export output=/data/unil/dee/dee_robi/scRNA-Seq_all_results/DROPLET_10X/Calls_IndCell_CellPopulation/
export ROUT=/data/unil/dee/dee_robi/scRNA-Seq_all_results/DROPLET_10X/

R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" InformationAllLibraries=\"$InformationAllLibraries\" folder_data=\"$folder_data\" desired_r_cutoff=\"$desired_r_cutoff\" output=\"$output\"" $SCRIPT_PATH/1Run/Call_PresentGenes_indivCell_cellPop.R $ROUT/Call_PresentGenes_indivCell_cellPop.Rout
