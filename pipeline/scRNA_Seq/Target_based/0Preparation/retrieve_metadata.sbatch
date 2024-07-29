#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00

#SBATCH --output=/users/jwollbre/bgee_pipeline/pipeline/scRNA_Seq/Target_based/retrieve_metadata.out
#SBATCH --error=/users/jwollbre/bgee_pipeline/pipeline/scRNA_Seq/Target_based/retrieve_metadata.err
#SBATCH --export=NONE
#SBATCH --job-name=retrieve_metadata_targetbased
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module use /software/module/
module add R/3.6.1;

export SCRIPT_PATH=/users/jwollbre/bgee_pipeline/pipeline/scRNA_Seq/Target_based
export scRNASeqExperiment=../../../generated_files/scRNA_Seq/Target_based/scRNASeqExperiment.tsv
export scRNASeqTBLibrary=../../../generated_files/scRNA_Seq/Target_based/scRNASeqLibrary_merged.tsv
export metadata_file=../../../generated_files/scRNA_Seq/Target_based/metadata_info_10X.txt
export information_file=../../../generated_files/scRNA_Seq/scRNA_Seq_info_TargetBased.txt
export ROUT=/work/FAC/FBM/DEE/mrobinso/bgee//downloads/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/

mkdir -p $ROUT
export R_LIBS_USER=/work/FAC/FBM/DEE/mrobinso/bgee//R_LIBS
R CMD BATCH --no-save --no-restore "--args scRNASeqExperiment=\"$scRNASeqExperiment\" scRNASeqTBLibrary=\"$scRNASeqTBLibrary\" metadata_file=\"$metadata_file\" information_file=\"$information_file\"" $SCRIPT_PATH/0Preparation/retrieve_metadata.R $ROUT/retrieve_metadata.Rout
