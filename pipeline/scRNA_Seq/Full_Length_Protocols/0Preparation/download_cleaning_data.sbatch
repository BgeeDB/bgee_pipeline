#!/bin/bash

#SBATCH --partition=long
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/download_cleaning_data.out
#SBATCH --error=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/download_cleaning_data.err
#SBATCH --export=NONE
#SBATCH --job-name=download_cleaning_data

module load Bioinformatics/Software/vital-it

module add R/3.6.1


export SCRIPT_PATH=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols
export metadata_info=../../../generated_files/scRNA_Seq/metadata_info.txt
export librariesDownloadedJura=../../../source_files/scRNA_Seq/juraLibrariesDownloaded.tsv
export output_folder=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA_Seq_All/scRNASeq_libraries_FullLength/
export ROUT=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA_Seq_All/


R CMD BATCH --no-save --no-restore "--args metadata_info=\"$metadata_info\" librariesDownloadedJura=\"$librariesDownloadedJura\" output_folder=\"$output_folder\"" $SCRIPT_PATH/0Preparation/download_cleaning_data.R $ROUT/download_cleaning_data.Rout
