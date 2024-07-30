#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00

#SBATCH --output=/users/jwollbre/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/retrieve_metadata.out
#SBATCH --error=/users/jwollbre/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/retrieve_metadata.err
#SBATCH --export=NONE
#SBATCH --job-name=retrieve_metadata_information_fulllength

module use /software/module/;
module add R/3.6.1


export SCRIPT_PATH=/users/jwollbre/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols
export pass_annotationControl=../../../generated_files/scRNA_Seq/Full_length/passScRNASeqLibrary.tsv
export metadata_info_file=../../../generated_files/scRNA_Seq/Full_length/metadata_info.tsv
export metadata_info_not_match_file=../../../generated_files/scRNA_Seq/Full_length/metadata_info_not_match.tsv
export ROUT=/work/FAC/FBM/DEE/mrobinso/bgee_sensitive/downloads/scRNA_Seq_All/scRNASeq_libraries_FullLength/

mkdir -p $ROUT
R CMD BATCH --no-save --no-restore "--args pass_annotationControl=\"$pass_annotationControl\" metadata_info_file=\"$metadata_info_file\" metadata_info_not_match_file=\"$metadata_info_not_match_file\"" $SCRIPT_PATH/0Preparation/retrieve_metadata.R $ROUT/retrieve_metadata.Rout
