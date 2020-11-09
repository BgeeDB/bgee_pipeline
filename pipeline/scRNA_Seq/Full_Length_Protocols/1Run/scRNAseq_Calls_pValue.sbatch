#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/scRNAseq_Calls_pValue_FL.out
#SBATCH --error=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/scRNAseq_Calls_pValue_FL.err
#SBATCH --export=NONE
#SBATCH --job-name=scRNAseq_Calls_pValue_FL

module load Bioinformatics/Software/vital-it

module add R/3.6.1

export SCRIPT_PATH=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export NEW_scRNASeq_sample_info=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/FULL_LENGTH/Quality_Control_Cell_Population/NEW_scRNASeq_sample_info.txt
export cells_folder=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/FULL_LENGTH/Kallisto_All_Cells/
export refIntergenic_folder=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/Rerence_intergenic_folder/
export desired_pValue_cutoff=0.05
export output_folder=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/FULL_LENGTH/scRNAseq_Callpresent/
export ROUT=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/FULL_LENGTH

R CMD BATCH --no-save --no-restore "--args NEW_scRNASeq_sample_info=\"$NEW_scRNASeq_sample_info\" cells_folder=\"$cells_folder\" refIntergenic_folder=\"$refIntergenic_folder\" desired_pValue_cutoff=\"$desired_pValue_cutoff\" output_folder=\"$output_folder\"" $SCRIPT_PATH/1Run/scRNAseq_Calls_pValue.R $ROUT/scRNAseq_Calls_pValue.Rout
