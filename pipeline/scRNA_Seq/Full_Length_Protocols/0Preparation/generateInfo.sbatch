#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/generateInfo_FL.out
#SBATCH --error=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/generateInfo_FL.err
#SBATCH --export=NONE
#SBATCH --job-name=generateInfo_FL
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;

export SCRIPT_PATH=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols
export folder_gtf=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/rna_seq/GTF_15/
export ROUT=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/FULL_LENGTH

R CMD BATCH --no-save --no-restore "--args folder_gtf=\"$folder_gtf\"" $SCRIPT_PATH/0Preparation/generateInfo.R $ROUT/generateInfo.Rout
