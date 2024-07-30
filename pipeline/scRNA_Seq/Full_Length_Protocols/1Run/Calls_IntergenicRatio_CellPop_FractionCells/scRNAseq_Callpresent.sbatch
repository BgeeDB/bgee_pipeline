#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/scRNAseq_Callpresent_FL.out
#SBATCH --error=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/scRNAseq_Callpresent_FL.err
#SBATCH --export=NONE
#SBATCH --job-name=scRNAseq_Callpresent_FL

module use /software/module/

module add R/3.6.1


export SCRIPT_PATH=/data/unil/dee/dee_robi/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export NEW_scRNASeq_sample_info=/data/unil/dee/dee_robi/scRNA-Seq_all_results/Quality_Control_Cell_Population/NEW_scRNASeq_sample_info.txt
export cells_folder=/data/unil/dee/dee_robi/scRNA-Seq_all_results/Kallisto_All_Cells/
export sum_species=/PATH_TO_SUM_RNA/
export gaussian_choice=../../generated_files/RNA_Seq/gaussian_choice_by_species.txt
export ratioValue=0.05
export desired_pValue_cutoff=0.05
export output_folder=/data/unil/dee/dee_robi/scRNA-Seq_all_results/scRNAseq_Callpresent/
export ROUT=/data/unil/dee/dee_robi/scRNA-Seq_all_results/

R CMD BATCH --no-save --no-restore "--args NEW_scRNASeq_sample_info=\"$NEW_scRNASeq_sample_info\" cells_folder=\"$cells_folder\" sum_species=\"$sum_species\" gaussian_choice=\"$gaussian_choice\" ratioValue=\"$ratioValue\" desired_pValue_cutoff=\"$desired_pValue_cutoff\" output_folder=\"$output_folder\"" $SCRIPT_PATH/1Run/scRNAseq_Callpresent.R $ROUT/scRNAseq_Callpresent.Rout
