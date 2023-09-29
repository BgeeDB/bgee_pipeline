#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00

#SBATCH --output=/users/jwollbre/GIT/pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/ensembl_To_uniprot.out
#SBATCH --error=/users/jwollbre/GIT/pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/ensembl_To_uniprot.error
#SBATCH --export=NONE
#SBATCH --job-name=ensembl_To_uniprot
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module use /softwate/module/
module add R/3.6.1;


export SCRIPT_PATH=/users/jwollbre/GIT/pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/
export scRNASeq_Info=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt
export output=../../../generated_files/scRNA_Seq/
export ROUT=../../../generated_files/scRNA_Seq/
  
R CMD BATCH --no-save --no-restore "--args scRNASeq_Info=\"$scRNASeq_Info\" output=\"$output\"" $SCRIPT_PATH/0Preparation/ensembl_To_uniprot.R $ROUT/ensembl_To_uniprot.Rout
