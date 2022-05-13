#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00

#SBATCH --output=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/retrieve_metadata.out
#SBATCH --error=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/retrieve_metadata.err
#SBATCH --export=NONE
#SBATCH --job-name=retrieve_metadata_targetbased
#SBATCH --mail-user sara.fonsecacosta@unil.ch

module load Bioinformatics/Software/vital-it
module add R/3.6.1;

export SCRIPT_PATH=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols
export scRNASeqExperiment=../../../source_files/scRNA_Seq/scRNASeqExperiment.tsv
export scRNASeqLibrary=../../../source_files/scRNA_Seq/scRNASeqLibrary.tsv
export metadata_file_path=../../../generated_files/scRNA_Seq/metadata_info_10X.tsv
export metadata_notmatch_file_path=../../../generated_files/scRNA_Seq/metadata_notMatch_10X.tsv
export ROUT=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/

mkdir -p $ROUT
export R_LIBS_USER=/work/FAC/FBM/DEE/mrobinso/bgee/R_LIBS
R CMD BATCH --no-save --no-restore "--args scRNASeqExperiment=\"$scRNASeqExperiment\" scRNASeqLibrary=\"$scRNASeqLibrary\" output=\"$output\"" $SCRIPT_PATH/0Preparation/retrieve_metadata.R $ROUT/retrieve_metadata.Rout
