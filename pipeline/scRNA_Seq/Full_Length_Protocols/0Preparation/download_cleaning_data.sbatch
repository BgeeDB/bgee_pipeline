#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00

#SBATCH --output=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Full_Length_Protocols/download_cleaning_data.out
#SBATCH --error=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Full_Length_Protocols/download_cleaning_data.err
#SBATCH --export=NONE
#SBATCH --job-name=download_cleaning_data

module use /software/module/

module add R/3.6.1


export SCRIPT_PATH=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Full_Length_Protocols
export metadata_info=../../../generated_files/scRNA_Seq/Full_length/metadata_info.txt
export librariesDownloadedJura=../../../generated_files/scRNA_Seq/juraLibrariesDownloaded.tsv
export output_folder=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/scRNASeq_libraries_FullLength/
export ROUT=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/


R CMD BATCH --no-save --no-restore "--args metadata_info=\"$metadata_info\" librariesDownloadedJura=\"$librariesDownloadedJura\" output_folder=\"$output_folder\"" $SCRIPT_PATH/0Preparation/download_cleaning_data.R $ROUT/download_cleaning_data.Rout
