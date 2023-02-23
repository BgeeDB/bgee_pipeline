#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00

#SBATCH --output=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/cleaning_barcodes.out
#SBATCH --error=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/cleaning_barcodes.err
#SBATCH --export=NONE
#SBATCH --job-name=cleaning_barcodes
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module use /software/module/
module add R/3.6.1;

export SCRIPT_PATH=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols
export barcodesFolder=../../../source_files/scRNA_Seq/Target_based/
export output=../../../generated_files/scRNA_Seq/Target_based/cleaning_barcodes/
export ROUT=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/

R CMD BATCH --no-save --no-restore "--args barcodesFolder=\"$barcodesFolder\" output=\"$output\"" $SCRIPT_PATH/0Preparation/cleaning_barcodes.R $ROUT/cleaning_barcodes.Rout
