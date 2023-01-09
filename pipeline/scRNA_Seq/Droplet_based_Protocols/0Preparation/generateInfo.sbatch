#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00

#SBATCH --output=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/generateInfo.out
#SBATCH --error=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/generateInfo.err
#SBATCH --export=NONE
#SBATCH --job-name=generateInfo_TB
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load /software/module/§
module add R/3.6.1;

export SCRIPT_PATH=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/
export folder_gtf=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/rna_seq/GTF_15/
export metadata_info_file=/home/jwollbre/Documents/git/bgee_pipeline_droplet_test/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/../../../generated_files/scRNA_Seq/Target_based/metadata_info_10X.txt
export ROUT=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/bgee_v15_0/scRNA-Seq_all_results_bgee_v15_0/DROPLET_10X/

R CMD BATCH --no-save --no-restore "--args folder_gtf=\"$folder_gtf\" metadata_info_file=\"$metadata_info_file\"" $SCRIPT_PATH/0Preparation/generateInfo.R $ROUT/generateInfo.Rout
