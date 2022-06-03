#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=6-00:00:00

#SBATCH --output=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/QC_CelltypeIdentification.out
#SBATCH --error=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/QC_CelltypeIdentification.err
#SBATCH --export=NONE
#SBATCH --job-name=QC_CellTypeIdentification
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;

export SCRIPT_PATH=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/
export scRNASeq_Info=../../../generated_files/scRNA_Seq/Target_based/scRNA_Seq_info_TargetBased.txt
export infoFolder=../../../generated_files/scRNA_Seq/Target_based/cleaning_barcodes/
export output=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/QC_CellType_identification/
export kallisto_bus_results=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/kallisto_bus_results/
export folderSupport=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/rna_seq/GTF_15/
export ROUT=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/
export R_LIBS_USER=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/R_LIBS

R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" kallisto_bus_results=\"$kallisto_bus_results\" folderSupport=\"$folderSupport\" infoFolder=\"$infoFolder\" output=\"$output\"" $SCRIPT_PATH/1Run/QC_CelltypeIdentification.R $ROUT/QC_CelltypeIdentification.Rout
