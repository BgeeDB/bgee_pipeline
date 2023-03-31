#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=4-00:00:00

#SBATCH --output=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/process_busFile_6.out
#SBATCH --error=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/process_busFile_6.err
#SBATCH --export=NONE
#SBATCH --job-name=process_busFile_6
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module use /software/module/;
module add R/3.6.1;
module add UHTS/Analysis/bustools/0.40.0;

export SCRIPT_PATH=/home/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/
export scRNASeq_Info=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased_6.txt
export kallisto_bus_results=/scratch/beegfs/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/kallisto_bus_results_bgee_v15_0/
export folderSupport=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/rna_seq/GTF_15/
export whiteList_Path=../../../source_files/scRNA_Seq/Target_based/
export ROUT=/scratch/beegfs/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/

R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" kallisto_bus_results=\"$kallisto_bus_results\" folderSupport=\"$folderSupport\" whiteList_Path=\"$whiteList_Path\"" $SCRIPT_PATH/1Run/process_busFile.R $ROUT/process_busFile_6.Rout
