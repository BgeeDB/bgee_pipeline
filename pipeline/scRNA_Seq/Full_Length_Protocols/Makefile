PIPELINEROOT := ../../
DIR_NAME := scRNA_Seq/

include $(PIPELINEROOT)Makefile.common

################## STEPS TO RUN ON AXIOM SERVER ##################

## Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations and store them in the source_files directory
## NOTE: this files should be already in the source_files (see whats happen with download)
get_annot:
	@$(WGET) $(ANNOTATION_GIT_URL)/Strains/StrainMapping.tsv                   && $(MV) StrainMapping.tsv                  $(STRAIN_MAPPING_FILE)
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeqLibrary.tsv               && $(MV) scRNASeqLibrary.tsv                $(INPUT_DIR)/scRNASeqLibrary.tsv                2>/dev/null  || rm -f scRNASeqLibrary.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeqExperiment.tsv            && $(MV) scRNASeqExperiment.tsv             $(INPUT_DIR)/scRNASeqExperiment.tsv             2>/dev/null  || rm -f scRNASeqExperiment.tsv
	@touch $@

## Control the annotation (run rule in front)
control_annotation: get_annot 0Preparation/pre_process_control_annotation.R
	@echo --- Starting the control of the annotation  ---
	$(CLUSTER_R_CMD) R CMD BATCH --no-save --no-restore '--args scRNASeqLibrary="$(SC_RNASEQ_LIB_FILEPATH)" passScRNASeqLibrary="$(SC_RNASEQ_LIB_PASS_FILEPATH)" notPassScRNASeqLibrary="$(SC_RNASEQ_LIB_NOT_PASS_FILEPATH)" cellsThreshold="$(MIN_CELLS_PER_LIBRARY_FULL_LENGTH)" ' 0Preparation/pre_process_control_annotation.R 0Preparation/pre_process_control_annotation.Rout > $@.tmp 2>&1
	@echo --- DONE ---
	# Map strain names
	#FIXME has to be done on other scRNAseq generated files? Or others come from that one?
	@perl -e 'use lib "../.."; use Utils; Utils::map_strain_names("../$(SC_RNASEQ_NEW_LIB_FILEPATH)", "$(STRAIN_MAPPING_FILE)")' >../$(SC_RNASEQ_NEW_LIB_FILEPATH).ori
	@$(MV) ../$(SC_RNASEQ_NEW_LIB_FILEPATH).ori ../$(SC_RNASEQ_NEW_LIB_FILEPATH)
	@$(MV) $@.tmp $@

## Retrieve metadata (run rule with sbatch)
retrieve_metadata: control_annotation 0Preparation/retrieve_metadata.R
	@echo --- Retrieve metadata information ---
	@sed -i 's@--output=.*@--output=${PWD}/retrieve_metadata.out@'                                               	0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/retrieve_metadata.err@'                                                 	0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                 	0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export NEW_scRNASeqLibrary=.*@export NEW_scRNASeqLibrary=$(SC_RNASEQ_LIB_PASS_FILEPATH)@'    		0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export metadata_info=.*@export metadata_info=$(SC_RNASEQ_METADATA_FILEPATH)@'                		0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export metadata_info_not_matchh=.*@export metadata_info=$(SC_RNASEQ_METADATA_NOT_MATCH_FILEPATH)@'   0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_FULL_LENGTH)@'                          		0Preparation/retrieve_metadata.sbatch
	@sbatch 0Preparation/retrieve_metadata.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Download cleaning data that still is not downloaded/present in Jura server (run rule with sbatch)
download_cleaning_data: retrieve_metadata 0Preparation/download_cleaning_data.R
	@echo --- Starting downloading the data ---
	@sed -i 's@--output=.*@--output=${PWD}/download_cleaning_data.out@'                                              0Preparation/download_cleaning_data.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/download_cleaning_data.err@'                                                0Preparation/download_cleaning_data.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                     0Preparation/download_cleaning_data.sbatch
	@sed -i 's@export metadata_info=.*@export metadata_info=$(SC_RNASEQ_METADATA_FILEPATH)@'                         0Preparation/download_cleaning_data.sbatch
	@sed -i 's@export librariesDownloadedJura=.*@export librariesDownloadedJura=$(SC_RNASEQ_JURA_LIB_FILEPATH)@'     0Preparation/download_cleaning_data.sbatch
	@sed -i 's@export output_folder=.*@export output_folder=$(SC_RNASEQ_DOWNLOAD_PATH_FULL_LENGTH)@'                 0Preparation/download_cleaning_data.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_LIB_DIR)@'                                       0Preparation/download_cleaning_data.sbatch
	@sbatch 0Preparation/download_cleaning_data.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Add new libraries download to the jura file list (run rule in front)
#XXX: are we sure all file names are *.fastq.gz ? Maybe add fq.gz/fastq.gzip/fq.gzip too ?
#XXX: why creating file in the /tmp directory?
list_new_downloads: download_cleaning_data
	@find $(SC_RNASEQ_DOWNLOAD_PATH_FULL_LENGTH) -type f -name \*.fastq.gz\* | xargs -r dirname | sed -e 's@^.*/@@' | sort -u > /tmp/new_downloads
	@cat $(SC_RNASEQ_JURA_LIB_FILEPATH) >>/tmp/new_downloads
	@sort -u /tmp/new_downloads >$(SC_RNASEQ_JURA_LIB_FILEPATH)
	@rm -f /tmp/new_downloads
	@$(GIT) add $(SC_RNASEQ_JURA_LIB_FILEPATH)
	@$(GIT) commit -m 'Add new downloaded libraries' $(SC_RNASEQ_JURA_LIB_FILEPATH) || true
	@touch $@

commit_annotation_and_metadata: list_new_downloads
	# Commit the library information file that will be used for the rest of the pipeline
	@$(GIT) add $(SC_RNASEQ_LIB_FILEPATH)
	@$(GIT) add $(SC_RNASEQ_LIB_PASS_FILEPATH)
	@$(GIT) add $(SC_RNASEQ_LIB_NOT_PASS_FILEPATH)
	@$(GIT) add $(SC_RNASEQ_METADATA_FILEPATH)
	@$(GIT) add $(SC_RNASEQ_METADATA_NOT_MATCH_FILEPATH)
	@$(GIT) commit -m 'Update metadata for scRNASeq for $(DBNAME)' $(SC_RNASEQ_LIB_FILEPATH) $(SC_RNASEQ_LIB_PASS_FILEPATH) $(SC_RNASEQ_LIB_NOT_PASS_FILEPATH) $(SC_RNASEQ_METADATA_FILEPATH) $(SC_RNASEQ_METADATA_NOT_MATCH_FILEPATH) || true
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
#TODO: should be written in perl to retrieve species from bgee DB (speciesName -> speciesID is hardcoded in this script)
prepare_singlecell_info: $(SC_RNASEQ_LIB_PASS_FILEPATH) 0Preparation/prepare_scrna_seq_sample_info.R
	@$(RM) $(SC_RNASEQ_SAMPINFO_FILEPATH)
	@echo --- prepare single cell information file ---
	@$(CLUSTER_FASTP_CMD) $(CLUSTER_R_CMD) R CMD BATCH --no-save --no-restore '--args NEW_scRNASeqLibrary="$(SC_RNASEQ_LIB_PASS_FILEPATH)" raw_cells_folder="$(SC_RNASEQ_FASTQ_FULL_LENGTH)" output_sample_info_file="$(SC_RNASEQ_SAMPINFO_FILEPATH)"' 0Preparation/prepare_scrna_seq_sample_info.R 0Preparation/prepare_scrna_seq_sample_info.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@

# generate kallisto indexes with BgeeCall. Indexes already generated during the bulk RNA-Seq pipeline can be reused.
indexes_bgeecall: prepare_singlecell_info
	# generate BgeeCall input file from the scrna_seq_sample_info.txt file
	# Run in front because it only parse a small file in order to create a new one.
	@perl 1Run/sc_sample_info_to_bgeecall.pl -sc_sample_info_file=$(SC_RNASEQ_SAMPINFO_FILEPATH) -transcriptome_dir=$(RNASEQ_CLUSTER_GTF) -annotation_dir=$(RNASEQ_CLUSTER_GTF) -fastq_dir=$(SC_RNASEQ_FASTQ_FULL_LENGTH) -output_dir=$(SC_RNASEQ_CLUSTER_ALL_RES) -bgeecall_file=$(SC_RNASEQ_BGEECALL_FILE_FULL_LENGTH) -ref_intergenic_dir=$(CLUSTER_REF_INTERGENIC_FOLDER) >>$@.tmp 2>&1
	# generate kallisto indexes with BgeeCall
	@$(CLUSTER_R_CMD) $(CLUSTER_KALLISTO_CMD) R CMD BATCH --vanilla --slave '--args bgeecall_input_file="$(SC_RNASEQ_BGEECALL_FILE_FULL_LENGTH)" account="mrobinso_bgee_sensitive" time="2:00:00" partition="normal" working_path="$(RNASEQ_CLUSTER_BGEECALL_OUTPUT)"' $(PIPELINEROOT)$(RNASEQPATH)1Run/bgeecall_index.R 1Run/bgeecall_index.Rout
	@$(MV) $@.tmp $@

# Use BgeeCall R package to generate present/absent calls.
presence_absence_bgeecall: indexes_bgeecall
	@$(CLUSTER_R_CMD) $(CLUSTER_KALLISTO_CMD) R CMD BATCH --vanilla --slave '--args bgeecall_input_file="$(SC_RNASEQ_BGEECALL_FILE_FULL_LENGTH)" account="mrobinso_bgee_sensitive" time="2:00:00" partition="normal" working_path="$(RNASEQ_CLUSTER_BGEECALL_OUTPUT)"' $(PIPELINEROOT)$(RNASEQPATH)1Run/bgeecall_calls.R 1Run/bgeecall_calls.Rout
	@touch $@

## Quality control by cell population to perform the calls (not pass if not bimodal)
QC_cellPopulation: presence_absence_bgeecall 1Run/QC_cellPopulation.R
	@echo --- Starting the quality control per cell population ---
	@sed -i 's@--output=.*@--output=${PWD}/QC_cellPopulation.out@'                                               		1Run/QC_cellPopulation.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/QC_cellPopulation.err@'                                                 		1Run/QC_cellPopulation.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                               		1Run/QC_cellPopulation.sbatch
	@sed -i 's@export scrna_seq_sample_info=.*@export scrna_seq_sample_info=$(SC_RNASEQ_SAMPINFO_FILEPATH)@'     		1Run/QC_cellPopulation.sbatch
	@sed -i 's@export cells_folder=.*@export cells_folder=$(SC_RNASEQ_CLUSTER_ALL_RES)@'                        		1Run/QC_cellPopulation.sbatch
	@sed -i 's@export sample_info_pass=.*@export sample_info_pass=$(SC_RNASEQ_SAMPINFO_PASS_FILEPATH)@'             	1Run/QC_cellPopulation.sbatch
	@sed -i 's@export sample_info_discarded=.*@export sample_info_discarded=$(SC_RNASEQ_SAMPINFO_NOT_PASS_FILEPATH)@'   1Run/QC_cellPopulation.sbatch
	@sed -i 's@export modality_info=.*@export modality_info=$(SC_RNASEQ_MODALITY_CELL_POP)@'                			1Run/QC_cellPopulation.sbatch
	@sed -i 's@export calls_file_name=.*@export calls_file_name=$(SC_RNASEQ_ABUNDANCEFILE)@'               				1Run/QC_cellPopulation.sbatch
	@sed -i 's@export plot=.*@export plot=yes@'                                                                  		1Run/QC_cellPopulation.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES)@'                                         		1Run/QC_cellPopulation.sbatch
	@sbatch 1Run/QC_cellPopulation.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

final_status: scRNAseq_Calls
	@$(GIT) status
	
final_commit: final_status
	# Commit the scRNASeq info files before and after QC, the bimodal file and information file about all libraries
	@$(GIT) add $(OUTPUT_DIR)scrna_seq_sample_info.txt
	@$(GIT) add $(OUTPUT_DIR)Modality_Cell_type_per_experiment.tsv
	@$(GIT) add $(OUTPUT_DIR)NEW_scRNASeq_sample_info.tsv
	@$(GIT) add $(OUTPUT_DIR)Discard_scRNASeq_sample_info.tsv
	@$(GIT) add $(OUTPUT_DIR)All_samples_stats_FL.tsv
	@$(GIT) commit -m 'Update info files for scRNASeq full-length for $(DBNAME)' $(OUTPUT_DIR)scrna_seq_sample_info.txt $(OUTPUT_DIR)Modality_Cell_type_per_experiment.tsv $(OUTPUT_DIR)NEW_scRNASeq_sample_info.tsv $(OUTPUT_DIR)Discard_scRNASeq_sample_info.tsv $(OUTPUT_DIR)All_samples_stats_FL.tsv || true
	@$(GIT) push
	@echo -e "All information is ready, you can make a tar of the results. \n"
	@touch $@

tar_all: final_commit
	cd $(SC_RNASEQ_CLUSTER_ALL_RES)
	tar -cvfSp tarball_scRNASeq.tar .
	gzip -9 tarball_scRNASeq.tar
	# TODO: cp tarball_scRNASeq.tar.gz to archive

