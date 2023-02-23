#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00

#SBATCH --output=/scratch/jwollbre/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/retrieve_metadata.out
#SBATCH --error=/scratch/jwollbre/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/retrieve_metadata.err
#SBATCH --export=NONE
#SBATCH --job-name=retrieve_metadata_targetbased
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module use /software/module/
module add R/3.6.1;

export SCRIPT_PATH=/scratch/jwollbre/git/bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols
export scRNASeqExperiment=../../../source_files/scRNA_Seq/scRNASeqExperiment.tsv
export scRNASeqTBLibrary=../../../source_files/scRNA_Seq/Target_based/scRNASeqTBLibrary.tsv
export output_file=../../../generated_files/scRNA_Seq/Target_based/metadata_info_10X.txt
export ROUT=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/

mkdir -p $ROUT
export R_LIBS_USER=/work/FAC/FBM/DEE/mrobinso/bgee/R_LIBS
R CMD BATCH --no-save --no-restore "--args scRNASeqExperiment=\"$scRNASeqExperiment\" scRNASeqTBLibrary=\"$scRNASeqTBLibrary\" output_file=\"$output_file\"" $SCRIPT_PATH/0Preparation/retrieve_metadata.R $ROUT/retrieve_metadata.Rout
