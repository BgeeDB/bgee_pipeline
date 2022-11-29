#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00

#SBATCH --output=/scratch/jwollbre/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/download_SRA.out
#SBATCH --error=/scratch/jwollbre/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/download_SRA.err
#SBATCH --export=NONE
#SBATCH --job-name=download_SRA
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module use /software/module/
module add R/3.6.1;
module add gcc/10.4.0;
module add sratoolkit/3.0.0;

export SCRIPT_PATH=/scratch/jwollbre/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols
export metadata_info=../../../generated_files/scRNA_Seq/Target_based/metadata_info_10X.txt
export librariesDownloadedJura=../../../generated_files/scRNA_Seq/juraLibrariesDownloaded.tsv
export output=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/
export ROUT=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/

R CMD BATCH --no-save --no-restore "--args metadata_info=\"$metadata_info\" librariesDownloadedJura=\"$librariesDownloadedJura\" output=\"$output\"" $SCRIPT_PATH/0Preparation/download_SRA.R $ROUT/download_SRA.Rout
