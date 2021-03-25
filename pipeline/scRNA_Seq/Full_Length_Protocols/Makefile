PIPELINEROOT := ../../
DIR_NAME := scRNA_Seq/

include $(PIPELINEROOT)Makefile.common

################## STEPS TO RUN ON AXIOM SERVER ##################

## Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations and store them in the source_files directory
## NOTE: this files should be already in the source_files (see whats happen with download)
## TODO : Should also upload worm annotations
get_annot:
	@$(WGET) $(ANNOTATION_GIT_URL)/Strains/StrainMapping.tsv                   && $(MV) StrainMapping.tsv                  $(STRAIN_MAPPING_FILE)
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeqLibrary.tsv               && $(MV) scRNASeqLibrary.tsv                $(INPUT_DIR)/scRNASeqLibrary.tsv                2>/dev/null  || rm -f scRNASeqLibrary.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeqExperiment.tsv            && $(MV) scRNASeqExperiment.tsv             $(INPUT_DIR)/scRNASeqExperiment.tsv             2>/dev/null  || rm -f scRNASeqExperiment.tsv
	@touch $@

## Control the annotation (run rule in front)
## XXX: Why no controls on the worm annotation file? Worm annotation are not used in the pipeline!!!!
control_annotation: get_annot 0Preparation/pre_process_control_annotation.R
	@echo --- Starting the control of the annotation  ---
	$(CLUSTER_R_CMD) R CMD BATCH --no-save --no-restore '--args scRNASeqLibrary="$(SC_RNASEQ_LIB_FILEPATH)" passScRNASeqLibrary="$(SC_RNASEQ_LIB_PASS_FILEPATH)" notPassScRNASeqLibrary="$(SC_RNASEQ_LIB_NOT_PASS_FILEPATH)" cellsThreshold="$(MIN_CELLS_PER_LIBRARY_FULL_LENGTH)" ' 0Preparation/pre_process_control_annotation.R 0Preparation/pre_process_control_annotation.Rout > $@.tmp 2>&1
	@echo --- DONE ---
	# Map strain names
	#FIXME has to be done on other scRNAseq generated files? Or others come from that one?
	@perl -e 'use lib "../.."; use Utils; Utils::map_strain_names("$(SC_RNASEQ_LIB_PASS_FILEPATH)", "$(STRAIN_MAPPING_FILE)")' > $(SC_RNASEQ_LIB_PASS_FILEPATH).ori
	@$(MV) $(SC_RNASEQ_LIB_PASS_FILEPATH).ori $(SC_RNASEQ_LIB_PASS_FILEPATH)
	@$(MV) $@.tmp $@

## Retrieve metadata (run rule with sbatch)
## NOTE: metadata here are information allowing to download fastq libraries. Should modify name of this file to be something like download_metadata.txt
retrieve_metadata: control_annotation 0Preparation/retrieve_metadata.R
	@echo --- Retrieve metadata information ---
	@sed -i 's@--output=.*@--output=${PWD}/retrieve_metadata.out@'                                                           0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/retrieve_metadata.err@'                                                             0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION}@'                                                             0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT}@'                                                                   0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                             0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export NEW_scRNASeqLibrary=.*@export NEW_scRNASeqLibrary=$(SC_RNASEQ_LIB_PASS_FILEPATH)@'                     0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export metadata_info=.*@export metadata_info=$(SC_RNASEQ_METADATA_FILEPATH)@'                                 0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export metadata_info_not_match=.*@export metadata_info_not_match=$(SC_RNASEQ_METADATA_NOT_MATCH_FILEPATH)@'   0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_FULL_LENGTH)@'                                           0Preparation/retrieve_metadata.sbatch
	@sbatch 0Preparation/retrieve_metadata.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Download cleaning data that still is not downloaded/present in Jura server (run rule with sbatch)
download_cleaning_data: retrieve_metadata 0Preparation/download_cleaning_data.R
	@echo --- Starting downloading the data ---
	@sed -i 's@--output=.*@--output=${PWD}/download_cleaning_data.out@'                                              0Preparation/download_cleaning_data.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/download_cleaning_data.err@'                                                0Preparation/download_cleaning_data.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION}@'                                                     0Preparation/download_cleaning_data.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT}@'                                                           0Preparation/download_cleaning_data.sbatch
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
#TODO: remove hardcoded mapping between speciesName and speciesID in the script
prepare_singlecell_info: $(SC_RNASEQ_LIB_PASS_FILEPATH) 0Preparation/prepare_scrna_seq_sample_info.R 
	@$(RM) $(SC_RNASEQ_SAMPINFO_FILEPATH)
	@echo --- prepare single cell information file ---
	@$(CLUSTER_FASTP_CMD) $(CLUSTER_R_CMD) R CMD BATCH --no-save --no-restore '--args NEW_scRNASeqLibrary="$(SC_RNASEQ_LIB_PASS_FILEPATH)" raw_cells_folder="$(SC_RNASEQ_FASTQ_FULL_LENGTH)" output_sample_info_file="$(SC_RNASEQ_SAMPINFO_FILEPATH)"' 0Preparation/prepare_scrna_seq_sample_info.R 0Preparation/prepare_scrna_seq_sample_info.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@

# generate kallisto indexes with BgeeCall. Indexes already generated during the bulk RNA-Seq pipeline can be reused.
# run script from bulk RNASeq pipeline
indexes_bgeecall: prepare_singlecell_info
	# generate BgeeCall input file from the scrna_seq_sample_info.txt file
	# Run in front because it only parse a small file in order to create a new one.
	@perl 1Run/sc_sample_info_to_bgeecall.pl -sc_sample_info_file=$(SC_RNASEQ_SAMPINFO_FILEPATH) -transcriptome_dir=$(RNASEQ_CLUSTER_GTF) -annotation_dir=$(RNASEQ_CLUSTER_GTF) -fastq_dir=$(SC_RNASEQ_FASTQ_FULL_LENGTH) -output_dir=$(SC_RNASEQ_CLUSTER_BGEECALL) -bgeecall_file=$(SC_RNASEQ_BGEECALL_FILE_FULL_LENGTH) -ref_intergenic_dir=$(CLUSTER_REF_INTERGENIC_FOLDER) >>$@.tmp 2>&1
	# generate kallisto indexes with BgeeCall. Use same working directory than bulk RNASeq as indexes are the same
	@$(CLUSTER_R_CMD) $(CLUSTER_KALLISTO_CMD) R CMD BATCH --vanilla --slave '--args bgeecall_input_file="$(SC_RNASEQ_BGEECALL_FILE_FULL_LENGTH)" account="$(CLUSTER_ACCOUNT_SENSITIVE)" time="10:00:00" partition="$(CLUSTER_PARTITION_SENSITIVE)" working_path="$(RNASEQ_CLUSTER_BGEECALL_OUTPUT)"' $(PIPELINEROOT)$(RNASEQPATH)1Run/bgeecall_index.R 1Run/bgeecall_index.Rout
	@$(MV) $@.tmp $@

# Use BgeeCall R package to generate present/absent calls.
# run script from bulk RNASeq pipeline
presence_absence_bgeecall: indexes_bgeecall
	#use same working_directory than bulk RNASeq as indexes are the same
	@$(CLUSTER_R_CMD) $(CLUSTER_KALLISTO_CMD) R CMD BATCH --vanilla --slave '--args bgeecall_input_file="$(SC_RNASEQ_BGEECALL_FILE_FULL_LENGTH)" account="$(CLUSTER_ACCOUNT_SENSITIVE)" time="2-00:00:00" partition="$(CLUSTER_PARTITION_SENSITIVE)" working_path="$(RNASEQ_CLUSTER_BGEECALL_OUTPUT)" decrypt_file_path=""' $(PIPELINEROOT)$(RNASEQPATH)1Run/bgeecall_calls.R 1Run/bgeecall_calls.Rout
	@touch $@

## Quality control by cell population to perform the calls (not pass if not bimodal)
QC_cellPopulation: presence_absence_bgeecall 1Run/QC_cellPopulation.R
	@echo --- Starting the quality control per cell population ---
	@sed -i 's@--output=.*@--output=${PWD}/QC_cellPopulation.out@'                                                      1Run/QC_cellPopulation.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/QC_cellPopulation.err@'                                                        1Run/QC_cellPopulation.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION_SENSITIVE}@'                                              1Run/QC_cellPopulation.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT_SENSITIVE}@'                                                    1Run/QC_cellPopulation.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                                       1Run/QC_cellPopulation.sbatch
	@sed -i 's@export scrna_seq_sample_info=.*@export scrna_seq_sample_info=$(SC_RNASEQ_SAMPINFO_FILEPATH)@'            1Run/QC_cellPopulation.sbatch
	@sed -i 's@export cells_folder=.*@export cells_folder=$(SC_RNASEQ_CLUSTER_BGEECALL)@'                               1Run/QC_cellPopulation.sbatch
	@sed -i 's@export sample_info_pass=.*@export sample_info_pass=$(SC_RNASEQ_SAMPINFO_PASS_FILEPATH)@'                 1Run/QC_cellPopulation.sbatch
	@sed -i 's@export sample_info_discarded=.*@export sample_info_discarded=$(SC_RNASEQ_SAMPINFO_NOT_PASS_FILEPATH)@'   1Run/QC_cellPopulation.sbatch
	@sed -i 's@export modality_info=.*@export modality_info=$(SC_RNASEQ_MODALITY_CELL_POP)@'                            1Run/QC_cellPopulation.sbatch
	@sed -i 's@export calls_file_name=.*@export calls_file_name=$(ABUNDANCEFILE)@'                                      1Run/QC_cellPopulation.sbatch
	@sed -i 's@export plot=.*@export plot=yes@'                                                                         1Run/QC_cellPopulation.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_QC_CELLPOP_RES)@'                                         1Run/QC_cellPopulation.sbatch
	@sbatch 1Run/QC_cellPopulation.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

# run script from bulk RNASeq pipeline
check_presence_absence: presence_absence_bgeecall 
	# check that presence absence calls have been generated for all libraries present in the rna_seq_sample_info file
	# generate file containing calls info for all libraries and use this file to generate different plots
	@$(CLUSTER_R_CMD) R CMD BATCH --vanilla --slave '--args bgeecall_sample_info="$(SC_RNASEQ_BGEECALL_FILE_FULL_LENGTH)" calls_dir="$(SC_RNASEQ_CLUSTER_BGEECALL)" presence_absence_report="$(SC_RNASEQ_CLUSTER_CALLS_STATS)" kallisto_report="$(SC_RNASEQ_CLUSTER_KALLISTO_STATS)"' $(PIPELINEROOT)$(RNASEQPATH)1Run/rna_seq_calls_plot.R $(SC_RNASEQ_FL_CLUSTER_R_LOG)rna_seq_calls_plot.Rout
	@touch $@

# Do not use slurm to run this light processing
# run script from bulk RNASeq pipeline
calculate_fpkm: presence_absence_bgeecall
	@$(CLUSTER_R_CMD) R CMD BATCH --vanilla --slave '--args all_results_dir="$(SC_RNASEQ_CLUSTER_BGEECALL)" calls_file_name="$(ABUNDANCEFILE)"' $(PIPELINEROOT)$(RNASEQPATH)1Run/calculate_fpkm.R $(SC_RNASEQ_FL_CLUSTER_R_LOG)calculate_fpkm.Rout
	@touch $@

final_status: scRNAseq_Calls
	@$(GIT) status
	
final_commit: final_status
	# Touch all files so that they are not removed from $(RNASEQ_CLUSTER_SCRATCH)
	@find $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH) -exec touch {} \; 2>warnings.$@
	# Back-up all data
	@$(RM) $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_all_results_$(DBNAME).tar* $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_presence_absence_$(DBNAME).tar* $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_qc_cell_population_$(DBNAME).tar*
	@tar -C $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH) -cf $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_presence_absence_$(DBNAME).tar        BgeeCall/                                   2>>warnings.$@
	@tar -C $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH) -cf $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_qc_cell_population_$(DBNAME).tar      Quality_Control_Cell_Population/            2>>warnings.$@
	@tar -C $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH) -cf $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_logs_R_$(DBNAME).tar                  Rout/                                       2>>warnings.$@
	@gzip -9 $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_presence_absence_$(DBNAME).tar             2>>warnings.$@
	@gzip -9 $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_qc_cell_population_$(DBNAME).tar           2>>warnings.$@
	@gzip -9 $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_logs_R_$(DBNAME).tar                       2>>warnings.$@
	# Move back-up data to $(RNASEQ_CLUSTER_ALL_RES_BACKUP)
	@$(MV) $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_presence_absence_$(DBNAME).tar.gz $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_qc_cell_population_$(DBNAME).tar.gz $(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_logs_R_$(DBNAME).tar.gz $(RNASEQ_CLUSTER_ALL_RES_BACKUP) 2>>warnings.$@
	# Whole archive is probably too big to be copied to our servers (>100Gb). Commit only the final gene-level expression + calls files:
	@$(CP) $(SC_RNASEQ_CLUSTER_CALLS_STATS)            $(SCRNASEQSAMPSTATS)   2>>warnings.$@
	@$(CP) $(SC_RNASEQ_CLUSTER_KALLISTO_STATS)         $(SCRNASEQREPORTINFO)  2>>warnings.$@	
	#To do to complete
	@echo -e "\tTODO: Commit/Push $(SCRNASEQREPORTINFO) $(SCRNASEQSAMPSTATS) $(SC_RNASEQ_SAMPINFO_FILEPATH) $(SC_RNASEQ_MODALITY_FILEPATH) $(SC_RNASEQ_SAMPINFO_PASS_FILEPATH) $(SC_RNASEQ_SAMPINFO_NOT_PASS_FILEPATH) files\n"
	@echo -e "\tTODO: Save tarballs in directory [$(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)] to $(CLUSTER_ARCHIVE_PATH)/full_length_rna_seq/all_results_$(DBNAME)/\n"
	@echo -e "\tTODO: Copy tarball [$(SC_RNASEQ_CLUSTER_ALL_RES_FULL_LENGTH)full_length_presence_absence_$(DBNAME).tar] to development server ($(PIPEHOST)) for db insertion"
	@echo -e "\t      For insertion only  '*$(ABUNDANCEFILE)' $(SCRNASEQSAMPSTATS) $(SCRNASEQREPORTINFO) look to be required\n"
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@echo -e "\tYou can go out of the sensitive cluster now\n"
	@touch $@

################################ INSERT CALLS IN DATABASE ################################

#NOTE For iterative updates die commands after insert/update/delete statement failure have to be changed to warn
#     because with iterative updates the database may contain duplicates in regard to what you are inserting/updating/deleting!
#NOTE In order to insert faster for Bgee 15 we first inserted everything except libraries and conditions. It allowed to disable the autocommit and commit results per library.
#     To do so we had to drop foreign key library ID in the result table
#     It could be nice to update the script to choose to insert everything, everything except libraries and conditions, or only libraries and conditions (e.g add script argument to manage that)
#     It is not possible to disable autocommit when inserting conditions. It creates an sql error of duplicated entries
insert_scRNA_seq: $(SCRNASEQALLRES) $(SCRNASEQSAMPSTATS) $(SCRNASEQREPORTINFO) $(SC_RNASEQ_EXPERIMENT_FILEPATH) $(UBERON_SEX_INFO_FILE_PATH) $(CUSTOM_UBERON_FILE_PATH) $(DEV_STAGE_ONT_FILE_PATH) $(SC_RNASEQ_SAMPINFO_PASS_FILEPATH) $(SC_RNASEQ_SAMPINFO_NOT_PASS_FILEPATH) $(EXTRAMAPPING_FILEPATH)
	# Launch the organ stage mapping tool (using $(CUSTOM_UBERON_FILE_PATH) and $(DEV_STAGE_ONT_FILE_PATH))
	@$(IDMAPPING)  $(IDMAPPINGPORT) &
	@$(STGMAPPING) $(STGMAPPINGPORT) &
	@sleep 50 # sleep because mappers need time to load Uberon
	# Insert RNA Seq data
	#@perl 3Insertion/insert_scrna_seq.pl -bgee=$(BGEECMD) -scRnaSeqExperiment=$(SC_RNASEQ_EXPERIMENT_FILEPATH) -library_info=$(SC_RNASEQ_SAMPINFO_PASS_FILEPATH) -excluded_libraries=$(SC_RNASEQ_SAMPINFO_NOT_PASS_FILEPATH) -library_stats=$(SCRNASEQSAMPSTATS) -report_info=$(SCRNASEQREPORTINFO) -all_results=$(SCRNASEQALLRES) -sex_info=$(UBERON_SEX_INFO_FILE_PATH) -extraMapping=$(EXTRAMAPPING_FILEPATH) -Aport=$(IDMAPPINGPORT) -Sport=$(STGMAPPINGPORT) > $@.tmp 2>warnings.$@
	@echo "Check inconsistencies between condition species and gene species (there should be none): " >> $@.tmp
	@$(MYSQL) -e "SELECT t1.* FROM scRnaSeqFullLengthResult AS t1 INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId INNER JOIN scRnaSeqFullLengthLibrary AS t3 ON t1.scRnaSeqFullLengthLibraryId = t3.scRnaSeqFullLengthLibraryId INNER JOIN cond AS t4 on t3.conditionId = t4.conditionId WHERE t2.speciesId != t4.speciesId" >> $@.tmp 2>>warnings.$@
	@echo >> $@.tmp
	@echo "Distinct strains in full length single cell RNA-Seq conditions, check that they are correct (e.g., no 'wild type' instead of 'wild-type')" >> $@.tmp
	@$(MYSQL) -e "SELECT DISTINCT t1.strain, t1.speciesId FROM cond AS t1 INNER JOIN scRnaSeqFullLengthLibrary AS t2 ON t1.conditionId = t2.conditionId ORDER BY t1.strain" >> $@.tmp 2>>warnings.$@
	@echo >> $@.tmp
	@echo "Statistics on libraries" >> $@.tmp
	@$(MYSQL) -e "SELECT t1.speciesId, COUNT(DISTINCT t2.scRnaSeqFullLengthLibraryId), MIN(t2.tpmThreshold), MAX(t2.tpmThreshold), AVG(t2.tpmThreshold), MIN(t2.allGenesPercentPresent), MAX(t2.allGenesPercentPresent), AVG(t2.allGenesPercentPresent), MIN(t2.proteinCodingGenesPercentPresent), MAX(t2.proteinCodingGenesPercentPresent), AVG(t2.proteinCodingGenesPercentPresent), MIN(t2.mappedReadsCount), MAX(t2.mappedReadsCount), AVG(t2.mappedReadsCount) FROM cond AS t1 INNER JOIN scRnaSeqFullLengthLibrary AS t2 ON t1.conditionId = t2.conditionId GROUP BY t1.speciesId ORDER BY COUNT(DISTINCT t2.scRnaSeqFullLengthLibraryId) DESC" >> $@.tmp 2>>warnings.$@
	@echo >> $@.tmp
	@echo "Statistics on conditions" >> $@.tmp
	@$(MYSQL) -e "SELECT t1.speciesId, COUNT(DISTINCT t2.scRnaSeqFullLengthLibraryId), COUNT(DISTINCT t2.conditionId), COUNT(DISTINCT t1.exprMappedConditionId), COUNT(DISTINCT t1.anatEntityId), COUNT(DISTINCT t1.cellTypeId), COUNT(DISTINCT t1.stageId), COUNT(DISTINCT t1.anatEntityId, t1.cellTypeId), COUNT(DISTINCT t1.anatEntityId, t1.cellTypeId, t1.stageId, t1.sex), GROUP_CONCAT(DISTINCT t1.sex ORDER BY t1.sex SEPARATOR ', '), GROUP_CONCAT(DISTINCT t1.strain ORDER BY t1.strain SEPARATOR ', ') FROM cond AS t1 INNER JOIN scRnaSeqFullLengthLibrary AS t2 ON t1.conditionId = t2.conditionId GROUP BY t1.speciesId ORDER BY COUNT(DISTINCT t2.scRnaSeqFullLengthLibraryId) DESC" >> $@.tmp 2>>warnings.$@
	@echo >> $@.tmp
	@echo "Same condition information, but for mapped conditions of expression tables" >> $@.tmp
	@$(MYSQL) -e "SELECT t3.speciesId, COUNT(DISTINCT t2.scRnaSeqFullLengthLibraryId), COUNT(DISTINCT t3.anatEntityId), COUNT(DISTINCT t3.cellTypeId), COUNT(DISTINCT t3.stageId), COUNT(DISTINCT t3.anatEntityId, t3.cellTypeId), COUNT(DISTINCT t3.anatEntityId, t3.cellTypeId, t3.stageId, t3.sex), GROUP_CONCAT(DISTINCT t3.sex ORDER BY t3.sex SEPARATOR ', '), GROUP_CONCAT(DISTINCT t3.strain ORDER BY t3.strain SEPARATOR ', ') FROM cond AS t1 INNER JOIN scRnaSeqFullLengthLibrary AS t2 ON t1.conditionId = t2.conditionId INNER JOIN cond AS t3 ON t1.exprMappedConditionId = t3.conditionId GROUP BY t3.speciesId ORDER BY COUNT(DISTINCT t2.scRnaSeqFullLengthLibraryId) DESC" >> $@.tmp 2>>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

check_conditions: insert_scRNA_seq
	@echo "Conditions with anat. entity not existing in related species:" > $@.tmp
	@$(MYSQL) -e "SELECT DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain FROM cond AS t1 LEFT OUTER JOIN anatEntityTaxonConstraint AS t2 ON t1.anatEntityId = t2.anatEntityId AND (t2.speciesId IS NULL OR t1.speciesId = t2.speciesId) LEFT OUTER JOIN anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId WHERE t2.anatEntityId IS NULL ORDER BY t1.speciesId" >> $@.tmp 2> $@.warn
	@echo >> $@.tmp
	@echo "Conditions with dev. stage not existing in related species:" >> $@.tmp
	@$(MYSQL) -e "SELECT DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain FROM cond AS t1 LEFT OUTER JOIN stageTaxonConstraint AS t2 ON t1.stageId = t2.stageId AND (t2.speciesId IS NULL OR t1.speciesId = t2.speciesId) LEFT OUTER JOIN anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId WHERE t2.stageId IS NULL ORDER BY t1.speciesId" >> $@.tmp 2>> $@.warn
	@echo "Conditions with cell type not existing in related species:" >> $@.tmp
	@$(MYSQL) -e "SELECT DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.cellTypeId, t5.anatEntityName, t1.sex, t1.strain FROM cond AS t1 LEFT OUTER JOIN stageTaxonConstraint AS t2 ON t1.stageId = t2.stageId AND (t2.speciesId IS NULL OR t1.speciesId = t2.speciesId) LEFT OUTER JOIN anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId LEFT OUTER JOIN anatEntity as t5 ON t5.anatEntityId = t1.cellTypeId WHERE t5.anatEntityId IS NULL and t5.anatEntityId IS NOT NULL ORDER BY t1.speciesId" >> $@.tmp 2>> $@.warn
	@$(MV) $@.tmp $@

insert_expression: 
	# Insert the expression summaries
	@perl 3Insertion/insert_scrna_seq_expression.pl -bgee=$(BGEECMD) > $@.tmp 2>warnings.$@
	@echo
	@echo "Searching for incorrect updates of scRnaSeqFullLengthResult/expression tables (there should be none)" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM scRnaSeqFullLengthResult WHERE reasonForExclusion = 'not excluded' AND expressionId IS NULL" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM scRnaSeqFullLengthResult WHERE reasonForExclusion != 'not excluded' AND expressionId IS NOT NULL" >> $@.tmp
	@echo >> $@.tmp
	@echo "Statistics for rscRnaSeqFullLengthResult" >> $@.tmp
	@$(MYSQL) -e "SELECT speciesId, detectionFlag, COUNT(*) AS scRnaSeqFullLengthResultCount FROM scRnaSeqFullLengthResult AS t1 INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId GROUP BY speciesId detectionFlag ORDER BY speciesId, detectionFlag" >> $@.tmp
	@echo "Statistics for expression table" >> $@.tmp
	@$(MYSQL) -e "SELECT speciesId, COUNT(*), COUNT(DISTINCT t1.bgeeGeneId), COUNT(DISTINCT t1.conditionId), COUNT(DISTINCT t3.anatEntityId, t3.stageId), COUNT(DISTINCT t3.anatEntityId, t3.stageId, t3.sex) FROM expression AS t1 INNER JOIN cond AS t3 ON t1.conditionId = t3.conditionId GROUP BY speciesId ORDER BY speciesId" >> $@.tmp
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@
