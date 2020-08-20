#!/bin/bash

#SBATCH --partition=wally
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=10-00:00:00

#SBATCH --output=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/retrieve_metadata.out
#SBATCH --error=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/retrieve_metadata.err
#SBATCH --export=NONE
#SBATCH --job-name=retrieve_metadata
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;

export SCRIPT_PATH=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols
export scRNASeqExperiment=../../../source_files/scRNA_Seq/scRNASeqExperiment.tsv
export scRNASeqLibrary=../../../source_files/scRNA_Seq/scRNASeqLibrary.tsv
export output=../../../generated_files/scRNA_Seq/
export ROUT=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/

R CMD BATCH --no-save --no-restore "--args scRNASeqExperiment=\"$scRNASeqExperiment\" scRNASeqLibrary=\"$scRNASeqLibrary\" output=\"$output\"" $SCRIPT_PATH/0Preparation/retrieve_metadata.R $ROUT/retrieve_metadata.Rout
