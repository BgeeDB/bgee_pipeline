PIPELINEROOT := ../
DIR_NAME := scRNA_Seq/

include $(PIPELINEROOT)Makefile.common

################## STEPS TO RUN ON AXIOM SERVER ##################

## Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations and store them in the source_files directory
## NOTE: this files should be already in the source_files (see whats happen with download)
get_annot:
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeqLibrary.tsv               && $(MV) scRNASeqLibrary.tsv                $(INPUT_DIR)/scRNASeqLibrary.tsv                2>/dev/null  || rm -f scRNASeqLibrary.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeqExperiment.tsv            && $(MV) scRNASeqExperiment.tsv             $(INPUT_DIR)/scRNASeqExperiment.tsv             2>/dev/null  || rm -f scRNASeqExperiment.tsv
	@touch $@

## Control the annotation (run rule in front)
control_annotation: get_annot 0Preparation/pre_process_control_annotation.R
	@echo --- Starting the control of the annotation  ---
	$(CLUSTER_R_CMD) R CMD BATCH --no-save --no-restore '--args scRNASeqLibrary="$(SC_RNASEQ_LIB_FILEPATH)" output_folder="${PWD}/$(OUTPUT_DIR)"' 0Preparation/pre_process_control_annotation.R 0Preparation/pre_process_control_annotation.Rout > $@.tmp 2>&1
	@echo --- DONE ---
	@$(MV) $@.tmp $@

## Retrieve metadata (run rule with sbatch)
retrieve_metadata: control_annotation 0Preparation/retrieve_metadata.R
	@echo --- Retrieve metadata information ---
	@sed -i 's@\(BASE *= \).*@\1"$(SC_RNASEQ_DOWNLOAD_PATH)";@'                                             0Preparation/retrieve_metadata.R
	@sed -i 's@--output=.*@--output=${PWD}/retrieve_metadata.out@'                                          0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/retrieve_metadata.err@'                                            0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                            0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export NEW_scRNASeqLibrary=.*@export NEW_scRNASeqLibrary=../$(SC_RNASEQ_NEW_LIB_FILEPATH)@'  0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export output_folder=.*@export output_folder=../$(OUTPUT_DIR)@'                              0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_FULL_LENGTH)@'                          0Preparation/retrieve_metadata.sbatch
	@sbatch 0Preparation/retrieve_metadata.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Download cleaning data that still is not downloaded/present in Jura server (run rule with sbatch)
download_cleaning_data: retrieve_metadata 0Preparation/download_cleaning_data.R
	@echo --- Starting downloading the data ---
	@sed -i 's@\(BASE *= \).*@\1"$(SC_RNASEQ_DOWNLOAD_PATH)";@'                                                   0Preparation/download_cleaning_data.R
	@sed -i 's@--output=.*@--output=${PWD}/download_cleaning_data.out@'                                           0Preparation/download_cleaning_data.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/download_cleaning_data.err@'                                             0Preparation/download_cleaning_data.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                  0Preparation/download_cleaning_data.sbatch
	@sed -i 's@export metadata_info=.*@export metadata_info=$(SC_RNASEQ_METADATA_FILEPATH)@'                      0Preparation/download_cleaning_data.sbatch
	@sed -i 's@export librariesDownloadedJura=.*@export librariesDownloadedJura=$(SC_RNASEQ_JURA_LIB_FILEPATH)@'  0Preparation/download_cleaning_data.sbatch
	@sed -i 's@export output_folder=.*@export output_folder=$(SC_RNASEQ_DOWNLOAD_LIB_DIR)@'                       0Preparation/download_cleaning_data.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_FULL_LENGTH)@'                                0Preparation/download_cleaning_data.sbatch
	@sbatch 0Preparation/download_cleaning_data.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Add new libraries download to the jura file list (run rule in front)
list_new_downloads: download_cleaning_data
	@find $(SC_RNASEQ_DOWNLOAD_LIB_DIR) -type f -name \*.fastq.gz\* | xargs -r dirname | sed -e 's@^.*/@@' | sort -u > /tmp/new_downloads
	@cat $(SC_RNASEQ_JURA_LIB_FILEPATH) >>/tmp/new_downloads
	@sort -u /tmp/new_downloads >$(SC_RNASEQ_JURA_LIB_FILEPATH)
	@rm -f /tmp/new_downloads
	@$(GIT) add $(SC_RNASEQ_JURA_LIB_FILEPATH)
	@$(GIT) commit -m 'Add new downloaded libraries' $(SC_RNASEQ_JURA_LIB_FILEPATH) || true
	@touch $@

commit_annotation_and_metadata: list_new_downloads
	# Commit the library information file that will be used for the rest of the pipeline
	@$(GIT) add $(INPUT_DIR)scRNASeqLibrary.tsv
	@$(GIT) add $(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv
	@$(GIT) add $(OUTPUT_DIR)metadata_info.txt
	@$(GIT) commit -m 'Update metadata for scRNASeq for $(DBNAME)' $(INPUT_DIR)scRNASeqLibrary.tsv $(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv $(OUTPUT_DIR)NOT_PASS_scRNASeqLibrary.tsv $(OUTPUT_DIR)metadata_info.txt $(OUTPUT_DIR)metadata_notMatch.txt || true
	@$(GIT) push
	@echo -e "All information is ready, you can go to JURA cluster to prepare the scRNA-Seq information file and continue the pipeline \n"
	@touch $@


## NOTE: Copy all git repository to JURA!
################## STEPS TO RUN ON JURA SERVER ##################

clusterJURA:
	@echo -e "\tBe sure everything is up-to-date before running single cell RNASeq pipeline\n"
	@touch $@

check_tools: clusterJURA
	@echo -e "\tGo to 'cd pipeline/scRNA_Seq/Full_Length_Protocols/' and be prepared to work\n"
	@echo -e "\n\tRun this command to give access to all modules installed on vital-it\n\tmodule add Bioinformatics/Software/vital-it\n"
	# Check if all required tools/libs are available
	@module add Bioinformatics/Software/vital-it || true
	@$(CLUSTER_R_CMD) which R                          > $@.tmp
	@$(CLUSTER_R_CMD) R -e 'library("dplyr")'         >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("stringr")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("data.table")'    >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("HelpersMG")'     >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("tools")'         >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("readr")'         >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("rjson")'         >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("plyr")'          >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("ggplot2")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("gridExtra")'     >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("gghighlight")'   >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("LaplacesDemon")' >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("reshape2")'      >> $@.tmp  2>/dev/null
	@which xz                                         >> $@.tmp
	@which sbatch                                     >> $@.tmp
	@$(CLUSTER_KALLISTO_CMD)   which kallisto         >> $@.tmp
	@$(CLUSTER_FASTP_CMD)      which fastp            >> $@.tmp
	@$(MV) $@.tmp $@


## Prepare single cell info file (run rule in front)
prepare_singlecell_info: check_tools $(SC_RNASEQ_NEW_LIB_FILEPATH) 0Preparation/prepare_scrna_seq_sample_info.R
	@$(RM) $(SC_RNASEQ_SAMPINFO_FILEPATH)
	@echo --- prepare single cell information file ---
	@$(CLUSTER_FASTP_CMD) $(CLUSTER_R_CMD) R CMD BATCH --no-save --no-restore '--args NEW_scRNASeqLibrary="$(SC_RNASEQ_NEW_LIB_FILEPATH)" raw_cells_folder="$(SC_RNASEQ_FASTQ_FULL_LENGTH)" output_folder="${PWD}/$(OUTPUT_DIR)"' 0Preparation/prepare_scrna_seq_sample_info.R 0Preparation/prepare_scrna_seq_sample_info.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@

## Run kallisto for all library/cells
kallisto: prepare_singlecell_info $(SC_RNASEQ_SAMPINFO_FILEPATH)
	mkdir -p $(SC_RNASEQ_CLUSTER_KALLISTO)
	@perl 1Run/slurm_scheduler_Kallisto_scRNASeq.pl -scrna_seq_sample_info=$(SC_RNASEQ_SAMPINFO_FILEPATH) -raw_cells_folder=$(SC_RNASEQ_FASTQ_FULL_LENGTH) -infoFolder=$(RNASEQ_CLUSTER_GTF) -output_folder=$(SC_RNASEQ_CLUSTER_KALLISTO) -cluster_kallisto_cmd='$(CLUSTER_KALLISTO_CMD)' -cluster_R_cmd='$(CLUSTER_R_CMD)'
	@echo "At the end it is a good idea to relaunch the kallisto rule step to be sure everything was run!"
	@touch $@

## Convert TranscriptID to GeneID and recalculations (TPM/FPKM)
analysis: kallisto 1Run/analysis.R
	@echo --- Starting the analysis ---
	@sed -i 's@\(BASE *= \).*@\1"$(SENSITIVE_GENERAL_PATH)";@'                                                1Run/analysis.R
	@sed -i 's@--output=.*@--output=${PWD}/analysis.out@'                                                     1Run/analysis.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/analysis.err@'                                                       1Run/analysis.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                              1Run/analysis.sbatch
	@sed -i 's@export scrna_seq_sample_info=.*@export scrna_seq_sample_info=$(SC_RNASEQ_SAMPINFO_FILEPATH)@'  1Run/analysis.sbatch
	@sed -i 's@export cells_folder=.*@export cells_folder=$(SC_RNASEQ_CLUSTER_KALLISTO)@'                     1Run/analysis.sbatch
	@sed -i 's@export infoFolder=.*@export infoFolder=$(RNASEQ_CLUSTER_GTF)@'                                 1Run/analysis.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES)@'                                      1Run/analysis.sbatch
	@sbatch 1Run/analysis.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Quality control by cell population
QC_cellPopulation: analysis 1Run/QC_cellPopulation.R
	@echo --- Starting the quality control per cell population ---
	@sed -i 's@\(BASE *= \).*@\1"$(SENSITIVE_GENERAL_PATH)";@'                                                1Run/QC_cellPopulation.R
	@sed -i 's@--output=.*@--output=${PWD}/QC_cellPopulation.out@'                                            1Run/QC_cellPopulation.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/QC_cellPopulation.err@'                                              1Run/QC_cellPopulation.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                              1Run/QC_cellPopulation.sbatch
	@sed -i 's@export scrna_seq_sample_info=.*@export scrna_seq_sample_info=$(SC_RNASEQ_SAMPINFO_FILEPATH)@'  1Run/QC_cellPopulation.sbatch
	@sed -i 's@export cells_folder=.*@export cells_folder=$(SC_RNASEQ_CLUSTER_KALLISTO)@'                     1Run/QC_cellPopulation.sbatch
	@sed -i 's@export output_folder=.*@export output_folder=$(SC_RNASEQ_CLUSTER_QC_CELLPOP_RES)@'             1Run/QC_cellPopulation.sbatch
	@sed -i 's@export plot=.*@export plot=yes@'                                                               1Run/QC_cellPopulation.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES)@'                                      1Run/QC_cellPopulation.sbatch
	@sbatch 1Run/QC_cellPopulation.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Summing the raw counts per cell-type, experiment and species
Sum_RawCounts_cellPopulation: QC_cellPopulation 1Run/Sum_RawCounts_cellPopulation.R
	@echo --- Start running the sum of raw counts of cell populations that pass QC ---
	@sed -i 's@\(BASE *= \).*@\1"$(SENSITIVE_GENERAL_PATH)";@'                                                          1Run/Sum_RawCounts_cellPopulation.R
	@sed -i 's@--output=.*@--output=${PWD}/Sum_RawCounts_cellPopulation.out@'                                           1Run/Sum_RawCounts_cellPopulation.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/Sum_RawCounts_cellPopulation.err@'                                             1Run/Sum_RawCounts_cellPopulation.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                        1Run/Sum_RawCounts_cellPopulation.sbatch
	@sed -i 's@export NEW_scRNASeq_sample_info=.*@export NEW_scRNASeq_sample_info=$(SC_RNASEQ_SAMPINFO_NEW_FILEPATH)@'  1Run/Sum_RawCounts_cellPopulation.sbatch
	@sed -i 's@export cells_folder=.*@export cells_folder=$(SC_RNASEQ_CLUSTER_KALLISTO)@'                               1Run/Sum_RawCounts_cellPopulation.sbatch
	@sed -i 's@export output_folder=.*@export output_folder=$(SC_RNASEQ_CLUSTER_SUM_RAW_COUNTS_CELLPOP_RES)@'           1Run/Sum_RawCounts_cellPopulation.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES)@'                                                1Run/Sum_RawCounts_cellPopulation.sbatch
	@sbatch 1Run/Sum_RawCounts_cellPopulation.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Call present genes per individual library/cell
scRNAseq_Callpresent: $(PIPELINEROOT)/RNASeq/sum_by_species Sum_RawCounts_cellPopulation 1Run/scRNAseq_Callpresent.R
	@echo --- Starting calling present genes per library ---
	@sed -i 's@\(BASE *= \).*@\1"$(SENSITIVE_GENERAL_PATH)";@'                                                          1Run/scRNAseq_Callpresent.R
	@sed -i 's@--output=.*@--output=${PWD}/scRNAseq_Callpresent.out@'                                                   1Run/scRNAseq_Callpresent.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/scRNAseq_Callpresent.err@'                                                     1Run/scRNAseq_Callpresent.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                        1Run/scRNAseq_Callpresent.sbatch
	@sed -i 's@export NEW_scRNASeq_sample_info=.*@export NEW_scRNASeq_sample_info=$(SC_RNASEQ_SAMPINFO_NEW_FILEPATH)@'  1Run/scRNAseq_Callpresent.sbatch
	@sed -i 's@export cells_folder=.*@export cells_folder=$(SC_RNASEQ_CLUSTER_KALLISTO)@'                               1Run/scRNAseq_Callpresent.sbatch
	@sed -i 's@export sum_species=.*@export sum_species=$(RNASEQ_CLUSTER_SUM_RES)@'                                     1Run/scRNAseq_Callpresent.sbatch
	@sed -i 's@export gaussian_choice=.*@export gaussian_choice=$(RNASEQ_CLUSTER_GAUSSIAN_CHOICE)@'                     1Run/scRNAseq_Callpresent.sbatch
	@sed -i 's@export ratioValue=.*@export ratioValue=$(INTERGENIC_N_PROPORTION)@'                                      1Run/scRNAseq_Callpresent.sbatch
	@sed -i 's@export output_folder=.*@export output_folder=$(SC_RNASEQ_CLUSTER_CALL_PRESENT)@'                         1Run/scRNAseq_Callpresent.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES)@'                                                1Run/scRNAseq_Callpresent.sbatch
	@sbatch 1Run/scRNAseq_Callpresent.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Call present genes per cell population
Sum_Calls_cellPopulation: scRNAseq_Callpresent 1Run/Sum_Calls_cellPopulation.R
	@echo --- Starting calling present genes per cell population ---
	@sed -i 's@\(BASE *= \).*@\1"$(SENSITIVE_GENERAL_PATH)";@'                                                          1Run/Sum_Calls_cellPopulation.R
	@sed -i 's@--output=.*@--output=${PWD}/Sum_Calls_cellPopulation.out@'                                               1Run/Sum_Calls_cellPopulation.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/Sum_Calls_cellPopulation.err@'                                                 1Run/Sum_Calls_cellPopulation.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                        1Run/Sum_Calls_cellPopulation.sbatch
	@sed -i 's@export NEW_scRNASeq_sample_info=.*@export NEW_scRNASeq_sample_info=$(SC_RNASEQ_SAMPINFO_NEW_FILEPATH)@'  1Run/Sum_Calls_cellPopulation.sbatch
	@sed -i 's@export cells_folder=.*@export cells_folder=$(SC_RNASEQ_CLUSTER_KALLISTO)@'                               1Run/Sum_Calls_cellPopulation.sbatch
	@sed -i 's@export output_folder=.*@export output_folder=$(SC_RNASEQ_CLUSTER_CALL_PRESENT_CELLPOP_RES)@'             1Run/Sum_Calls_cellPopulation.sbatch
	@sed -i 's@export ratioValue=.*@export ratioValue=$(SC_RNASEQ_RATIO_PROPORTION)@'                                   1Run/Sum_Calls_cellPopulation.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES)@'                                                1Run/Sum_Calls_cellPopulation.sbatch
	@sbatch 1Run/Sum_Calls_cellPopulation.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

final_status: Sum_Calls_cellPopulation
	@$(GIT) status


tar_all: final_status
	cd $(SC_RNASEQ_CLUSTER_ALL_RES)
	tar -cvfSp tarball_scRNASeq.tar .
	gzip -9 tarball_scRNASeq.tar
	# TODO: cp tarball_scRNASeq.tar.gz to archive

