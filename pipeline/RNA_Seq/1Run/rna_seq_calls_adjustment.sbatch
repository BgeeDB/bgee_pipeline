#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/unil/dee/dee_robi/GIT/pipeline/RNA_Seq/1Run/calls_adjustment.out
#SBATCH --error=/data/unil/dee/dee_robi/GIT/pipeline/RNA_Seq/1Run/calls_adjustment.err
#SBATCH --export=NONE
#SBATCH --job-name=calls_adjustment

module load Bioinformatics/Software/vital-it

module add R/3.6.1

export SCRIPT_PATH=/data/unil/dee/dee_robi/GIT/pipeline/RNA_Seq/1Run/
export rna_seq_sample_info=../generated_files/RNA_Seq/rna_seq_sample_info.tsv
export RNASeq_folder_calls=/scratch/beegfs/FAC/FBM/DEE/mrobinso/bgee_sensitive/presence_absence_$(DBNAME)/
export cutoff=0.01
export ROUT=/scratch/beegfs/FAC/FBM/DEE/mrobinso/bgee_sensitive/presence_absence_$(DBNAME)/

R CMD BATCH --no-save --no-restore '--args rna_seq_sample_info=\"$rna_seq_sample_info\" RNASeq_folder_calls=\"$RNASeq_folder_calls\" cutoff=\"$cutoff\"' $SCRIPT_PATH/rna_seq_calls_adjustment.R $ROUT/rna_seq_calls_adjustment.Rout
