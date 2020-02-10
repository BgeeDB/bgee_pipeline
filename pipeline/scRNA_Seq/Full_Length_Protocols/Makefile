PIPELINEROOT := ../../
DIR_NAME := scRNA_Seq/

include $(PIPELINEROOT)Makefile.common

################## STEPS TO RUN ON BGEE SERVER ##################

## Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations and store them in the source_files directory
## this files should be already in the source_files (we cannot download the file in the bgee server!)
get_annot:
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeqLibrary.tsv               && $(MV) scRNASeqLibrary.tsv                $(INPUT_DIR)/scRNASeqLibrary.tsv                2>/dev/null  || rm -f scRNASeqLibrary.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeqExperiment.tsv            && $(MV) scRNASeqExperiment.tsv             $(INPUT_DIR)/scRNASeqExperiment.tsv             2>/dev/null  || rm -f scRNASeqExperiment.tsv
	@touch $@

## Control the annotation 
control_annotation: get_annot 0Preparation/pre_process_control_annotation.R
	@echo --- Starting the control of the annotation  ---
	R CMD BATCH --no-save --no-restore '--args scRNASeqLibrary="$(SCRNASEQ_LIB_FILEPATH)" output_folder="${PWD}/$(OUTPUT_DIR)"' 0Preparation/pre_process_control_annotation.R 0Preparation/pre_process_control_annotation.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@

## Retrieve metadata
retrieve_metadata: control_annotation 0Preparation/retrieve_metadata.R
	@echo --- Retrieve metadata information ---
	R CMD BATCH --no-save --no-restore '--args NEW_scRNASeqLibrary="${PWD}/$(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv" output_folder="${PWD}/$(OUTPUT_DIR)"' 0Preparation/retrieve_metadata.R 0Preparation/retrieve_metadata.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@

commit_annotation_and_metadata: retrieve_metadata
	# Commit the library information file that will be used for the rest of the pipeline
	@$(GIT) add $(INPUT_DIR)scRNASeqLibrary.tsv
	@$(GIT) add $(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv
	@$(GIT) add $(OUTPUT_DIR)metadata_info.txt
	@$(GIT) commit -m 'Update metadata for scRNASeq for $(DBNAME)' $(INPUT_DIR)scRNASeqLibrary.tsv $(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv $(OUTPUT_DIR)metadata_info.txt
	@$(GIT) push
	@echo -e "\t$(RNASEQ_SAMPINFO_FILEPATH) is ready, you can go to regular cluster to download new SRA files with get_SRA.pl as *admin* user\n\tDo a 'git pull' before starting\n"
	@$(MV) $@.tmp $@
	
################## STEPS TO RUN ON AXIOM ##################
## First time before pull: just clone the github + Makefile.Config and then Pull axiom after just pull
pull_annotation_and_metadata: commit_annotation_and_metadata
	@$(GIT) pull
	@touch $@

## Download cleaning data that still is not downloaded and present in Jura server
download_cleaning_data: pull_annotation_and_metadata 0Preparation/download_cleaning_data.R
	@echo --- Starting downloading the data ---
	R CMD BATCH --no-save --no-restore '--args metadata_info="${PWD}/$(OUTPUT_DIR)metadata_info.txt" librariesDownloadedJura="$(SCRNASEQ_JURA_LIB)" output_folder="$(SC_RNASEQ_DOWNLOAD_LIB_DIR)"' 0Preparation/download_cleaning_data.R 0Preparation/download_cleaning_data.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@
	
## add new downloads to the list (see sebastien changes)
list_new_downloads: download_cleaning_data
	@find $(SC_RNASEQ_DOWNLOAD_LIB_DIR) -type f -name \*.fastq.gz\* | xargs dirname | sed -e 's@^.*/@@' | sort -u > /tmp/new_downloads
	@cat $(SC_RNASEQ_DOWNLOAD_LIB_DIR) >>/tmp/new_downloads
	@sort -u /tmp/new_downloads >$(SC_RNASEQ_DOWNLOAD_LIB_DIR)
	@rm -f /tmp/new_downloads
	@git status
	@touch $@

## Prepare single cell info file
prepare_singlecell_info: list_new_downloads $(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv 0Preparation/prepare_scrna_seq_sample_info.R
	@echo --- prepare single cell information file ---
	R CMD BATCH --no-save --no-restore '--args NEW_scRNASeqLibrary="$(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv" cells_folder="$(SC_RNASEQ_DOWNLOAD_LIB_DIR)" output_folder="${PWD}/$(OUTPUT_DIR)"' 0Preparation/prepare_scrna_seq_sample_info.R 0Preparation/prepare_scrna_seq_sample_info.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@
	
commit_singlecell_info: prepare_singlecell_info
	# Commit the single cell information file
	@$(GIT) add $(SC_RNASEQ_SAMPINFO_FILEPATH)
	@$(GIT) commit -m 'Create the scRNA-Seq info file for $(DBNAME)' $(SC_RNASEQ_SAMPINFO_FILEPATH)
	@$(GIT) push
	@echo -e "\t$(SC_RNASEQ_SAMPINFO_FILEPATH) is ready.\n"
	@$(MV) $@.tmp $@
	
	
## NOTE: Copy all git repository to JURA!
################## STEPS TO RUN ON JURA ##################

clusterJURA:

	@echo -e "\tBe sure everything is up-to-date before running single cell RNASeq pipeline\n"
	@touch $@

check_tools: clusterJURA
	@echo -e "\n\tFirst of all, go to '$(RNASEQ_CLUSTER_READONLY)GIT/' and do a 'git pull'"
	@echo -e "\tThen 'cd pipeline/scRNA_Seq/Full_Length_Protocols/' and be prepared to work\n"
	@echo -e "\n\tRun this command to give access to all modules installed on vital-it\n\tmodule add Bioinformatics/Software/vital-it\n"
	@echo -e "\n\tRun this command to prevent errors with Utils.pm\n\tmodule add Development/Ensembl_API/$(ENSRELEASE);\n"
	# Check if all required tools/libs are available 
	@module add Bioinformatics/Software/vital-it || true
	@$(CLUSTER_R_CMD) which R                         >> $@.tmp
	@$(CLUSTER_R_CMD) R -e 'library("dplyr")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("stringr")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("data.table")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("HelpersMG")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("tools")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("readr")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("rjson")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("plyr")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("ggplot2")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("gridExtra")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("gghighlight")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("LaplacesDemon")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("reshape2")'       >> $@.tmp  2>/dev/null
	@which xz                                         >> $@.tmp
	@which sbatch                                     >> $@.tmp
	@$(CLUSTER_KALLISTO_CMD)   which kallisto         >> $@.tmp
	@$(MV) $@.tmp $@

## !!!!!! Make script for kallisto in R (just function and loop)
## clean the slurm_scheduler.pl for single cell! slurm_scheduler.pl + my_R
KALLISTO_abundance_all_intergenic: prepare_indexed_transcriptome $(RNASEQ_SAMPINFO_FILEPATH)
	# Running the pipeline
	@perl 1Run/slurm_scheduler.pl -sample_info_file=$(RNASEQ_SAMPINFO_FILEPATH) -exclude_sample_file=$(RNASEQ_SAMPEXCLUDED_FILEPATH) -output_log_folder=$(RNASEQ_CLUSTER_LOG) -index_folder=$(RNASEQ_CLUSTER_GTF) -fastq_folder=$(RNASEQ_BIGBGEE_FASTQ) -kallisto_out_folder=$(RNASEQ_CLUSTER_ALL_RES) -ens_release=$(ENSRELEASE) -ens_metazoa_release=$(ENSMETAZOARELEASE) -data_host=$(DATAHOST) -data_login=$(DATALOGIN) -enc_passwd_file=$(ENCRYPT_PASSWD_FILE) -cluster_kallisto_cmd="$(CLUSTER_KALLISTO_CMD)" -cluster_R_cmd="$(CLUSTER_R_CMD)" >$@.tmp 2>&1
	@echo "TODO: At the end it is a good idea to relaunch the abundance_all_intergenic step to be sure everything was run!"
	@$(MV) $@.tmp $@
	
	
	
analysis:	

## Quality control by cell population
QC_cellPopulation: abundance_all_intergenic 1Run/QC_cellPopulation.R
	@echo --- Starting the quality control per cell population ---
	R CMD BATCH --no-save --no-restore '--args scrna_seq_sample_info="${PWD}/$(OUTPUT_DIR)scrna_seq_sample_info.txt" cells_folder="$(SC_RNASEQ_DOWNLOAD_LIB_DIR)" output_folder="$(SC_RNASEQ_CLUSTER_QC_CELLPOP_RES)" plot="yes"' 1Run/QC_cellPopulation.R 1Run/QC_cellPopulation.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@
	
## Sum raw counts of all libraries by cell population, experiment and species
Sum_RawCounts_cellPopulation: QC_cellPopulation 1Run/Sum_RawCounts_cellPopulation.R
	@echo --- Start running the sum of raw counts of cell populations that pass QC ---
	R CMD BATCH --no-save --no-restore '--args NEW_scRNASeq_sample_info="$(SC_RNASEQ_CLUSTER_QC_CELLPOP_RES)NEW_scRNASeq_sample_info.tsv" cells_folder="$(SC_RNASEQ_DOWNLOAD_LIB_DIR)" output_folder="$(SC_RNASEQ_CLUSTER_SUM_RAW_COUNTS_CELLPOP_RES)"' 1Run/Sum_RawCounts_cellPopulation.R 1Run/Sum_RawCounts_cellPopulation.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@
	
## NOTE: This step of single cell pipeline is dependent of the intergenic region from RNA-Seq pipeline
## Call present genes per individual library
scRNAseq_Callpresent: Sum_RawCounts_cellPopulation 1Run/scRNAseq_Callpresent.R
	@echo --- Starting calling present genes per library ---
	R CMD BATCH --no-save --no-restore '--args NEW_scRNASeq_sample_info="$(SC_RNASEQ_CLUSTER_QC_CELLPOP_RES)NEW_scRNASeq_sample_info.tsv" cells_folder="$(SC_RNASEQ_DOWNLOAD_LIB_DIR)" sum_species="$(RNASEQ_CLUSTER_SUM_RES)" gaussian_choice="$(RNASEQ_CLUSTER_GAUSSIAN_CHOICE)" ratioValue="$(INTERGENIC_N_PROPORTION)" output_folder="$(SC_RNASEQ_CLUSTER_CALL_PRESENT_CELLPOP_RES)"' 1Run/scRNAseq_Callpresent.R 1Run/scRNAseq_Callpresent.Rout
	@echo --- DONE ---
	@$(MV) $@.tmp $@

## Call present genes per cell population
Sum_Calls_cellPopulation: scRNAseq_Callpresent 1Run/Sum_Calls_cellPopulation.R
	@echo --- Starting calling present genes per cell population ---
	R CMD BATCH --no-save --no-restore '--args NEW_scRNASeq_sample_info="$(SC_RNASEQ_CLUSTER_QC_CELLPOP_RES)NEW_scRNASeq_sample_info.tsv" cells_folder="$(SC_RNASEQ_DOWNLOAD_LIB_DIR)" output_folder="$(SC_RNASEQ_CLUSTER_CALL_PRESENT_CELLPOP_RES)" ratioValue="$(SC_RNASEQ_RATIO_PROPORTION)"' 1Run/Sum_Calls_cellPopulation.R 1Run/Sum_Calls_cellPopulation.Rout
	@echo --- DONE ---
	@$(MV) $@.tmp $@

