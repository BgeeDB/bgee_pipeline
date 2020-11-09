#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/analysis.out
#SBATCH --error=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/analysis.err
#SBATCH --export=NONE
#SBATCH --job-name=analysis

module load Bioinformatics/Software/vital-it

module add R/3.6.1


export SCRIPT_PATH=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export scrna_seq_sample_info=../../generated_files/scRNA_Seq/scrna_seq_sample_info.tsv
export cells_folder=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/Kallisto_All_Cells/
export infoFolder=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/rna_seq/GTF_15/
export ROUT=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/FULL_LENGTH


R CMD BATCH --no-save --no-restore "--args scrna_seq_sample_info=\"$scrna_seq_sample_info\" cells_folder=\"$cells_folder\" infoFolder=\"$infoFolder\"" $SCRIPT_PATH/1Run/analysis.R $ROUT/analysis.Rout

