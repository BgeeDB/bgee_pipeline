#!/bin/bash

#SBATCH --partition=wally
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=10-00:00:00

#SBATCH --output=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/download_SRA.out
#SBATCH --error=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/download_SRA.err
#SBATCH --export=NONE
#SBATCH --job-name=download_SRA
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;

export SCRIPT_PATH=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols
export metadata_info=../../../generated_files/scRNA_Seq/metadata_info_10X.txt
export librariesDownloadedJura=../../../generated_files/scRNA_Seq/juraLibrariesDownloaded.tsv
export output=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/
export ROUT=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/

R CMD BATCH --no-save --no-restore "--args metadata_info=\"$metadata_info\" librariesDownloadedJura=\"$librariesDownloadedJura\" output=\"$output\"" $SCRIPT_PATH/0Preparation/download_SRA.R $ROUT/download_SRA.Rout
