#!/bin/bash

#SBATCH --partition=urblauna
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=2-00:00:00

#SBATCH --output=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_2/scRNA-Seq_all_results_bgee_v15_2/DROPLET_10X//cleaning_barcodes/cleaning_barcode.out
#SBATCH --error=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_v15_2/scRNA-Seq_all_results_bgee_v15_2/DROPLET_10X//cleaning_barcodes/cleaning_barcode.err
#SBATCH --export=NONE
#SBATCH --job-name=cleaning_barcodes
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module use /software/module/
module add R/3.6.1;

export SCRIPT_PATH=/users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Target_based/
export barcodesFolder=/users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Target_based/../../../source_files/scRNA_Seq/Target_based/
export output=/users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Target_based/../../../generated_files/scRNA_Seq/Target_based/cleaned_barcodes/

Rscript $SCRIPT_PATH/0Preparation/cleaning_barcodes.R barcodesFolder=\"$barcodesFolder\" output=\"$output\"
