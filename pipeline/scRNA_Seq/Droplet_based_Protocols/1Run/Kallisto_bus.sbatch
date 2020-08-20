#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Kallisto_bus.out
#SBATCH --error=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Kallisto_bus.err
#SBATCH --export=NONE
#SBATCH --job-name=Kallisto_bus
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;


export SCRIPT_PATH=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols
export metadata_file=../../../generated_files/scRNA_Seq/metadata_info_10X.txt
export annotation_file=../../../source_files/scRNA_Seq/scRNASeqLibrary.tsv
export folder_data=/data/unil/dee/dee_robi/FASTQ/scRNAseq/10X/
export folderSupport=/data/unil/dee/dee_robi/rna_seq/GTF_bgee_v15_0/
export output=../../../generated_files/scRNA_Seq/
export ROUT=/data/unil/dee/dee_robi/scRNA-Seq_all_results/DROPLET_10X/


R CMD BATCH --no-save --no-restore "--args metadata_file=\"$metadata_file\" annotation_file=\"$annotation_file\" folder_data=\"$folder_data\" folderSupport=\"$folderSupport\" output=\"$output\"" $SCRIPT_PATH/1Run/Kallisto_bus.R $ROUT/Kallisto_bus.Rout
