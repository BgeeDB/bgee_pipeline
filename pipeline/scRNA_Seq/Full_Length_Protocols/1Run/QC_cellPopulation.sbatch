#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=2-00:00:00

#SBATCH --output=/home/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/QC_cellPopulation.out
#SBATCH --error=/home/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/QC_cellPopulation.err
#SBATCH --export=NONE
#SBATCH --job-name=QC_cellPopulation

module load Bioinformatics/Software/vital-it

module add R/3.6.1


export SCRIPT_PATH=/home/jwollbre/Documents/git/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export scrna_seq_sample_info=../../../generated_files/scRNA_Seq/scrna_seq_sample_info.txt
export cells_folder=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/test_scRNA-Seq_all_results/
export sample_info_pass=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/test_scRNA-Seq_all_results/FULL_LENGTH/Quality_Control_Cell_Population/NEW_scRNASeq_sample_info.txt
export sample_info_discarded=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/test_scRNA-Seq_all_results/FULL_LENGTH/Quality_Control_Cell_Population/Discard_scRNASeq_sample_info.tsv
export modality_info=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/test_scRNA-Seq_all_results/FULL_LENGTH/Quality_Control_Cell_Population/Modality_Cell_type_per_experiment.tsv
export calls_file_name=gene_level_abundance+calls.tsv
export plot=yes
export ROUT=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/test_scRNA-Seq_all_results/

R CMD BATCH --no-save --no-restore "--args scrna_seq_sample_info=\"$scrna_seq_sample_info\" cells_folder=\"$cells_folder\" sample_info_pass=\"$sample_info_pass\" sample_info_discarded=\"$sample_info_discarded\" modality_info=\"$modality_info\" calls_file_name=\"$calls_file_name\" plot=\"$plot\"" $SCRIPT_PATH/1Run/QC_cellPopulation.R $ROUT/QC_cellPopulation.Rout
