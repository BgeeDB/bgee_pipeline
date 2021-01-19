#!/bin/bash

#SBATCH --partition=wally
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=24:00:00

#SBATCH --output=/users/sfonseca1/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/ensembl_To_uniprot.out
#SBATCH --error=/users/sfonseca1/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/ensembl_To_uniprot.error
#SBATCH --export=NONE
#SBATCH --job-name=ensembl_To_uniprot
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;

export SCRIPT_PATH=/users/sfonseca1/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols
export scRNASeq_Info=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt
export output=../../../generated_files/scRNA_Seq/
export ROUT=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/
  
R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" output=\"$output\"" $SCRIPT_PATH/0Preparation/ensembl_To_uniprot.R $ROUT/ensembl_To_uniprot.Rout
