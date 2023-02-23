#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=6-00:00:00

#SBATCH --output=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Kallisto_bus.out
#SBATCH --error=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Kallisto_bus.err
#SBATCH --export=NONE
#SBATCH --job-name=Kallisto_bus
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module use /software/module/
module add R/3.6.1;
module add UHTS/Analysis/kallisto/0.46.0;

export SCRIPT_PATH=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/
export metadata_file=../../../generated_files/scRNA_Seq/metadata_info_10X.txt
export annotation_file=../../../source_files/scRNA_Seq/scRNASeqLibrary.tsv
export folder_data=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/FASTQ/scRNAseq/10X/
export folderSupport=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/rna_seq/GTF_15/
export output=../../../generated_files/scRNA_Seq/
export kallisto_bus_results=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/kallisto_bus_results/
export ROUT=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/


R CMD BATCH --no-save --no-restore "--args metadata_file=\"$metadata_file\" annotation_file=\"$annotation_file\" folder_data=\"$folder_data\" folderSupport=\"$folderSupport\" output=\"$output\" kallisto_bus_results=\"$kallisto_bus_results\"" $SCRIPT_PATH/1Run/Kallisto_bus.R $ROUT/Kallisto_bus.Rout
