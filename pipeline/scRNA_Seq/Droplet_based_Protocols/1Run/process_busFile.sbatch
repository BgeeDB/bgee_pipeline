#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/process_busFile.out
#SBATCH --error=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/process_busFile.err
#SBATCH --export=NONE
#SBATCH --job-name=process_busFile
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;

export SCRIPT_PATH=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols
export scRNASeq_Info=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt
export folder_data=/data/unil/dee/dee_robi/FASTQ/scRNAseq/10X/
export folderSupport=/data/unil/dee/dee_robi/GIT/pipeline/rna_seq/GTF_15.0/
export whiteList_Path=../../../source_files/scRNA_Seq/
export ROUT=/data/unil/dee/dee_robi/scRNA-Seq_all_results/DROPLET_10X/

R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" folder_data=\"$folder_data\" folderSupport=\"$folderSupport\" whiteList_Path=\"$whiteList_Path\"" $SCRIPT_PATH/1Run/process_busFile.R $ROUT/process_busFile.Rout
