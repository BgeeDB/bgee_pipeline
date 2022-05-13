#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=2-00:00:00

#SBATCH --output=/users/jwollbre/GIT/pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/retrieve_metadata.out
#SBATCH --error=/users/jwollbre/GIT/pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/retrieve_metadata.err
#SBATCH --export=NONE
#SBATCH --job-name=retrieve_metadata_information_fulllength

module load Bioinformatics/Software/vital-it

module add R/3.6.1


export SCRIPT_PATH=/users/jwollbre/GIT/pipeline/pipeline/scRNA_Seq/Full_Length_Protocols
export NEW_scRNASeqLibrary=../../../generated_files/scRNA_Seq/NEW_scRNASeqLibrary.tsv
export metadata_info=../../../generated_files/scRNA_Seq/metadata_info.txt
export metadata_info_not_match=../../../generated_files/scRNA_Seq/
export ROUT=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/scRNASeq_libraries_FullLength/

mkdir -p $ROUT
R CMD BATCH --no-save --no-restore "--args NEW_scRNASeqLibrary=\"$NEW_scRNASeqLibrary\" metadata_info=\"$metadata_info\" metadata_info_not_match=\"metadata_info_not_match\"" $SCRIPT_PATH/0Preparation/retrieve_metadata.R $ROUT/retrieve_metadata.Rout
