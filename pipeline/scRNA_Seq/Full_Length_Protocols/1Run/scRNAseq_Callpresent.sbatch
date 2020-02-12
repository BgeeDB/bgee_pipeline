#!/bin/bash

#SBATCH --partition=long
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/scRNAseq_Callpresent.out
#SBATCH --error=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/scRNAseq_Callpresent.err
#SBATCH --export=NONE
#SBATCH --job-name=scRNAseq_Callpresent

module load Bioinformatics/Software/vital-it

module add R/3.6.1


export SCRIPT_PATH=/users/sfonseca1/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export NEW_scRNASeq_sample_info=../../generated_files/scRNA_Seq/NEW_scRNASeq_sample_info.txt
export cells_folder=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNASeq_libraries/
export sum_species=/scratch/beegfs/FAC/FBM/DEE/bgee_sensitive/all_results_$(DBNAME)
export gaussian_choice=../../generated_files/RNA_Seq/gaussian_choice_by_species.txt
export ratioValue=0.05
export output_folder=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA-Seq_all_results/Sum_Raw_Counts_Cell_Population/
export ROUT=/scratch/wally/FAC/FBM/DEE/mrobinso/bgee/scRNA-Seq_all_results/


R CMD BATCH --no-save --no-restore "--args NEW_scRNASeq_sample_info=\"$NEW_scRNASeq_sample_info\" cells_folder=\"$cells_folder\" sum_species=\"$sum_species\" gaussian_choice=\"$gaussian_choice\" ratioValue=\"$ratioValue\" output_folder=\"$output_folder\"" $SCRIPT_PATH/1Run/scRNAseq_Callpresent.R $ROUT/scRNAseq_Callpresent.Rout
