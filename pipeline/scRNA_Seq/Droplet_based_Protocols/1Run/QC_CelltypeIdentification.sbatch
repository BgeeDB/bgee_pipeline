#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/QC_CelltypeIdentification.out
#SBATCH --error=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/QC_CelltypeIdentification.err
#SBATCH --export=NONE
#SBATCH --job-name=QC_CellTypeIdentification
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;

export SCRIPT_PATH=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols
export scRNASeq_Info=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt
export folder_data=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/FASTQ/scRNAseq/10X/
export infoFolder=../../../source_files/scRNA_Seq/
export output=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/DROPLET_10X/QC_CellType_identification/
export kallisto_bus_results=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/DROPLET_10X/QC_CellType_identification/
export ROUT=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/DROPLET_10X/
export R_LIBS_USER=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/R_LIBS

R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" folder_data=\"$folder_data\" kallisto_bus_results=\"$kallisto_bus_results\" infoFolder=\"$infoFolder\" output=\"$output\"" $SCRIPT_PATH/1Run/QC_CelltypeIdentification.R $ROUT/QC_CelltypeIdentification.Rout
