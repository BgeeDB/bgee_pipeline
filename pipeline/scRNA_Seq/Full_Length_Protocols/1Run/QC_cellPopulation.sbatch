#!/bin/bash

#SBATCH --partition=urblauna
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=2-00:00:00

#SBATCH --output=/users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/QC_cellPopulation.out
#SBATCH --error=/users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/QC_cellPopulation.err
#SBATCH --export=NONE
#SBATCH --job-name=QC_cellPopulation

module use /software/module/

module add R/3.6.1


export SCRIPT_PATH=/users/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export scrna_seq_sample_info=../../../generated_files/scRNA_Seq/Full_length/scrna_seq_sample_info.tsv
export cells_folder=/scratch/jwollbre//bgee_v15_new_intergenic/scRNA-Seq_all_results_bgee_v15_new_intergenic/FULL_LENGTH/BgeeCall/
export sample_info_pass=../../../generated_files/scRNA_Seq/Full_length/NEW_scRNASeq_sample_info.tsv
export sample_info_discarded=../../../generated_files/scRNA_Seq/Full_length/Discard_scRNASeq_sample_info.tsv
export modality_info=/scratch/jwollbre//bgee_v15_new_intergenic/scRNA-Seq_all_results_bgee_v15_new_intergenic/FULL_LENGTH/Quality_Control_Cell_Population/Modality_Cell_type_per_experiment.tsv
export calls_file_name=gene_level_abundance+calls.tsv
export plot=yes
export ROUT=/scratch/jwollbre//bgee_v15_new_intergenic/scRNA-Seq_all_results_bgee_v15_new_intergenic/FULL_LENGTH/Quality_Control_Cell_Population/

mkdir -p $ROUT
R CMD BATCH --no-save --no-restore "--args scrna_seq_sample_info=\"$scrna_seq_sample_info\" cells_folder=\"$cells_folder\" sample_info_pass=\"$sample_info_pass\" sample_info_discarded=\"$sample_info_discarded\" modality_info=\"$modality_info\" calls_file_name=\"$calls_file_name\" plot=\"$plot\"" $SCRIPT_PATH/1Run/QC_cellPopulation.R $ROUT/QC_cellPopulation.Rout
