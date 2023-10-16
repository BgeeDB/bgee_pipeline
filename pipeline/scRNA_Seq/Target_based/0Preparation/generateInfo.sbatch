#!/bin/bash

#SBATCH --partition=urblauna
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00

#SBATCH --output=/scratch/jwollbre//bgee_v15_new_intergenic/scRNA-Seq_all_results_bgee_v15_new_intergenic/DROPLET_10X//generate_info/generateInfo.out
#SBATCH --error=/scratch/jwollbre//bgee_v15_new_intergenic/scRNA-Seq_all_results_bgee_v15_new_intergenic/DROPLET_10X//generate_info/generateInfo.err
#SBATCH --export=NONE
#SBATCH --job-name=generateInfo_TB

module use /software/module/
module add R/3.6.1;

export SCRIPT_PATH=/users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/
export folder_gtf=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/rna_seq/GTF_new_intergenic_deconvoluted/
export metadata_info_file=/users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/../../../generated_files/scRNA_Seq/Target_based/metadata_info_10X.txt

Rscript $SCRIPT_PATH/0Preparation/generateInfo.R folder_gtf=\"$folder_gtf\" metadata_info_file=\"$metadata_info_file\"
