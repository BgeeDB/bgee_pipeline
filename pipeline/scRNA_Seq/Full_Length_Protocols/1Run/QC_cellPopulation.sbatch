#!/bin/bash

#SBATCH --partition=normal
#SBATCH --account=mrobinso_bgee_sensitive

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00

#SBATCH --output=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/QC_cellPopulation.out
#SBATCH --error=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/QC_cellPopulation.err
#SBATCH --export=NONE
#SBATCH --job-name=QC_cellPopulation

module load Bioinformatics/Software/vital-it

module add R/3.6.1


export SCRIPT_PATH=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/
export scrna_seq_sample_info=../../generated_files/scRNA_Seq/scrna_seq_sample_info.tsv
export sample_info_pass=../../generated_files/scRNA_Seq/NEW_scRNASeq_sample_info.tsv
export sample_info_discarded=../../generated_files/scRNA_Seq/Discard_scRNASeq_sample_info.tsv
export modality_info=../../generated_files/scRNA_Seq/Modality_Cell_type_per_experiment.tsv
export calls_file_name=abundance+geneLevel+intergenic.tsv
export plot=yes
export ROUT=/data/FAC/FBM/DEE/mrobinso/bgee_sensitive/scRNA-Seq_all_results/FULL_LENGTH

R CMD BATCH --no-save --no-restore "--args scrna_seq_sample_info=\"$scrna_seq_sample_info\" cells_folder=\"$cells_folder\" output_folder=\"$output_folder\" plot=\"$plot\"" $SCRIPT_PATH/1Run/QC_cellPopulation.R $ROUT/QC_cellPopulation.Rout
