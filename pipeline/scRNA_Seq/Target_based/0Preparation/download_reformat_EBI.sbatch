#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00

#SBATCH --output=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/download_reformat_EBI.out
#SBATCH --error=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/download_reformat_EBI.err
#SBATCH --export=NONE
#SBATCH --job-name=download_reformat_EBI
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module use /software/module/
module add R/3.6.1;

export SCRIPT_PATH=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols
export metadata_info=../../../generated_files/scRNA_Seq/Target_based/metadata_info_10X.txt
export librariesDownloadedJura=../../../generated_files/scRNA_Seq/juraLibrariesDownloaded.tsv
export output=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/
export bamtofastq=/users/sfonseca1/bin/bamtofastq-1.2.0
export ROUT=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/

R CMD BATCH --no-save --no-restore "--args metadata_info=\"$metadata_info\" librariesDownloadedJura=\"$librariesDownloadedJura\" output=\"$output\" bamtofastq=\"$bamtofastq\"" $SCRIPT_PATH/0Preparation/download_reformat_EBI.R $ROUT/download_reformat_EBI.Rout
