PIPELINEROOT := ../../
DIR_NAME := scRNA_Seq/Full_length/

include $(PIPELINEROOT)Makefile.common

################## STEPS TO RUN ON CURNAGL SERVER ##################

## Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations submodule and store them in the source_files directory
get_annot:
	@$(GIT) submodule init
	@$(GIT) submodule update
	@$(CP) $(ANNOT_SUBMODULE_DIR)/Strains/StrainMapping.tsv             $(STRAIN_MAPPING_FILE)
	@$(CP) $(ANNOT_SUBMODULE_DIR)/scRNA_Seq/scRNASeqLibrary_merged.tsv  $(SC_RNASEQ_LIB_FL_FILEPATH)
	@$(CP) $(ANNOT_SUBMODULE_DIR)/scRNA_Seq/scRNASeqExperiment.tsv      $(SC_RNASEQ_EXP_FL_FILEPATH)
	@touch $@

## Control the annotation (run rule in front)
filter_annotation: get_annot 0Preparation/filter_annotation.R
	@echo --- Filter FL libraries from annotation files ---
	@$(CLUSTER_LOAD_APPTAINER) $(CLUSTER_CONTAINER_CMD_WITH_BIND_CURNAGL) Rscript 0Preparation/filter_annotation.R scRNASeqLibrary=\"$(SC_RNASEQ_LIB_FL_FILEPATH)\" scRNASeqExperiment=\"$(SC_RNASEQ_EXP_FL_FILEPATH)\" filteredLibrariesFile=\"$(SC_RNASEQ_FL_LIB_FILTERED_FILEPATH)\" filteredExperimentsFile=\"$(SC_RNASEQ_FL_EXP_FILTERED_FILEPATH)\" > $@.tmp 2> $@.err
	@echo --- Map strain names ---
	@$(CLUSTER_SOFTWARE_BIN) perl -e 'use lib "../.."; use Utils; Utils::map_strain_names("$(SC_RNASEQ_FL_LIB_FILTERED_FILEPATH)", "$(STRAIN_MAPPING_FILE)")' > $(SC_RNASEQ_FL_LIB_FILTERED_FILEPATH).ori
	@$(MV) $(SC_RNASEQ_FL_LIB_FILTERED_FILEPATH).ori $(SC_RNASEQ_FL_LIB_FILTERED_FILEPATH)
	@$(MV) $@.tmp $@

## Retrieve metadata (run rule with sbatch)
## NOTE: metadata here are information allowing to download fastq libraries a posteriori
THREADS_RETRIEVE_METADATA := 16
retrieve_metadata: filter_annotation 0Preparation/retrieve_metadata.R
	@srun --output=${PWD}/$@.out --error=${PWD}/$@.err --partition=${CLUSTER_PARTITION} --account=${CLUSTER_ACCOUNT} --job-name=retrieve_metadata_fulllength --mem=2G --time=1-00:00:00 --nodes=1 --ntasks=1 --export=ALL --cpus-per-task=$(THREADS_RETRIEVE_METADATA) bash -c "$(CLUSTER_SOFT_UNIL) $(CLUSTER_SOFT_APPTAINER) $(CLUSTER_CONTAINER_CMD_WITH_BIND_CURNAGL) Rscript ${PWD}/0Preparation/retrieve_metadata.R filteredLibraryAnnotation=\'$(SC_RNASEQ_FL_LIB_FILTERED_FILEPATH)\' metadataFile=\'$(SC_RNASEQ_METADATA_FL_FILEPATH)\' threads=\'$(THREADS_RETRIEVE_METADATA)\'"
	$(MV) $@.out $@

# Parallelize download of libraries and do not download again already downloaded libraries
parallelized_download: retrieve_metadata
	@$(CLUSTER_SOFTWARE_BIN) perl 0Preparation/parallelized_download_SRA.pl -metadataFile=$(SC_RNASEQ_METADATA_FL_FILEPATH) -parallelJobs=50 -downloadedLibraries=$(SC_RNASEQ_DOWNLOADED_LIB_FILEPATH) -excludedLibraries=$(SC_RNASEQ_EXCLUDED_LIB_FILEPATH) -outputDir=$(SC_RNASEQ_DOWNLOAD_PATH_FULL_LENGTH) --loadApptainerModule="$(CLUSTER_LOAD_APPTAINER)" --containerCmd="$(CLUSTER_CONTAINER_CMD_WITH_BIND_CURNAGL)" -queue=$(CLUSTER_PARTITION) -account=$(CLUSTER_ACCOUNT) >$@.tmp 2> $@.err
	@$(MV) $@.tmp $@

commit_annotation_and_metadata: list_new_downloads
	# Commit the library information file that will be used for the rest of the pipeline
	@$(GIT) add $(SC_RNASEQ_LIB_FL_FILEPATH)
	@$(GIT) add $(SC_RNASEQ_LIB_PASS_FILEPATH)
	@$(GIT) add $(SC_RNASEQ_LIB_NOT_PASS_FILEPATH)
	@$(GIT) add $(SC_RNASEQ_METADATA_FILEPATH)
	@$(GIT) add $(SC_RNASEQ_METADATA_NOT_MATCH_FILEPATH)
	@$(GIT) add $(SC_RNASEQ_DOWNLOADED_LIB_FILEPATH)
	@$(GIT) commit -m 'Update metadata for scRNASeq for $(DBNAME)' $(SC_RNASEQ_LIB_FL_FILEPATH) $(SC_RNASEQ_LIB_PASS_FILEPATH) $(SC_RNASEQ_LIB_NOT_PASS_FILEPATH) $(SC_RNASEQ_METADATA_FILEPATH) $(SC_RNASEQ_METADATA_NOT_MATCH_FILEPATH) $(SC_RNASEQ_DOWNLOADED_LIB_FILEPATH) || true
	@$(GIT) push
	@echo -e "All information is ready, you can go to the sensitive cluster to prepare the scRNA-Seq information file and continue the pipeline \n"
	@touch $@


## NOTE: Copy all git repository to the sensitive cluster
################## STEPS TO RUN ON SENSITIVE CLUSTER ##################

clusterJURA:
	@echo -e "\tBe sure everything is up-to-date before running full-length single cell RNASeq pipeline\n"
	@touch $@

check_tools: clusterJURA
	@echo -e "\tGo to 'cd pipeline/scRNA_Seq/Full_Length_Protocols/' and be prepared to work\n"
	# Check if all required tools/libs are available
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) which R                          > $@.tmp
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("Biostrings")'    >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("HelpersMG")'     >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("LaplacesDemon")' >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("data.table")'    >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("dplyr")'         >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("forcats")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("gghighlight")'   >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("ggplot2")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("gridExtra")'     >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("plyr")'          >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("reshape2")'      >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("rjson")'         >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("stringr")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) R -e 'library("tools")'         >> $@.tmp  2>/dev/null
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) which xz                        >> $@.tmp
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) which sbatch                    >> $@.tmp
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) which kallisto                  >> $@.tmp
	@$(CLUSTER_CONTAINER_CMD_WITH_BIND) which fastp                     >> $@.tmp
	@$(MV) $@.tmp $@

## Prepare single cell info file
prepare_singlecell_info: $(SC_RNASEQ_LIB_PASS_FILEPATH) 0Preparation/prepare_scrna_seq_sample_info.R
	@echo --- prepare single cell information file ---
	@srun --export=ALL --partition=${SENSITIVE_CLUSTER_PARTITION} --account=${SENSITIVE_CLUSTER_ACCOUNT} --job-name=$@ --output=${PWD}/$@.out --error=${PWD}/$@.err  --nodes=1 --time=12:00:00 --ntasks=1 --cpus-per-task=1 --mem=10G bash -c "$(CLUSTER_CONTAINER_CMD_WITH_BIND) Rscript $(PWD)/0Preparation/prepare_scrna_seq_sample_info.R filteredLibraryAnnotation=\'$(PWD)/$(SC_RNASEQ_FL_LIB_FILTERED_FILEPATH)\' metadataFile=\'$(PWD)/$(SC_RNASEQ_METADATA_FL_FILEPATH)\' fastqDir=\'$(SC_RNASEQ_FASTQ_FULL_LENGTH)\' outputDir=\'$(SC_RNASEQ_FL_FOLDER_GENERATED)\'"
	@$(MV) $@.out $@

# generate kallisto indexes with BgeeCall. Indexes already generated during the bulk RNA-Seq pipeline can be reused.
# run script from bulk RNASeq pipeline
indexes_bgeecall: prepare_singlecell_info
	# generate BgeeCall input file from the scrna_seq_sample_info.txt file
	# Run in front because it only parse a small file in order to create a new one.
	$(CLUSTER_SOFTWARE_BIN) $(CLUSTER_CONTAINER_CMD_WITH_BIND) perl 1Run/sc_sample_info_to_bgeecall.pl -sc_sample_info_file=$(SC_RNASEQ_SAMPINFO_FL_FILEPATH) -transcriptome_dir=$(RNASEQ_CLUSTER_GTF) -annotation_dir=$(RNASEQ_CLUSTER_GTF) -fastq_dir=$(SC_RNASEQ_FASTQ_FULL_LENGTH) -output_dir=$(SC_RNASEQ_CLUSTER_BGEECALL) -bgeecall_file=$(SC_RNASEQ_BGEECALL_FILE_FULL_LENGTH) -ref_intergenic_dir=$(CLUSTER_REF_INTERGENIC_FOLDER) > $@.tmp 2> $@.err
	# generate kallisto indexes with BgeeCall. Use same working directory than for bulk RNASeq
	$(CLUSTER_CONTAINER_CMD_WITH_BIND) Rscript $(PIPELINEROOT)$(RNASEQPATH)1Run/bgeecall_index.R  bgeecall_input_file=\"$(SC_RNASEQ_BGEECALL_FILE_FULL_LENGTH)\" account=\"$(SENSITIVE_CLUSTER_ACCOUNT)\" time=\"10:00:00\" partition=\"$(SENSITIVE_CLUSTER_PARTITION)\" working_path=\"$(RNASEQ_CLUSTER_BGEECALL_OUTPUT)\" container_cmd=\"$(CLUSTER_CONTAINER_CMD_WITHOUT_BIND)\" >> $@.tmp 2>> $@.err
	@$(MV) $@.tmp $@

# Use BgeeCall R package to generate present/absent calls.
# run script from bulk RNASeq pipeline
presence_absence_bgeecall: indexes_bgeecall
	# generate again the BgeeCall input file containing libraries to process. It allows to filter libraries already processed before running BgeeCall
	$(CLUSTER_SOFTWARE_BIN) $(CLUSTER_CONTAINER_CMD_WITH_BIND) perl 1Run/sc_sample_info_to_bgeecall.pl -sc_sample_info_file=$(SC_RNASEQ_SAMPINFO_FL_FILEPATH) -transcriptome_dir=$(RNASEQ_CLUSTER_GTF) -annotation_dir=$(RNASEQ_CLUSTER_GTF) -fastq_dir=$(SC_RNASEQ_FASTQ_FULL_LENGTH) -output_dir=$(SC_RNASEQ_CLUSTER_BGEECALL) -bgeecall_file=$(SC_RNASEQ_BGEECALL_FILE_FULL_LENGTH) -ref_intergenic_dir=$(CLUSTER_REF_INTERGENIC_FOLDER) -filter_processed > $@.tmp 2> $@.err
	$(CLUSTER_CONTAINER_CMD_WITH_BIND) Rscript $(PIPELINEROOT)$(RNASEQPATH)1Run/bgeecall_calls.R bgeecall_input_file=\"$(SC_RNASEQ_BGEECALL_FILE_FULL_LENGTH)\" account=\"$(SENSITIVE_CLUSTER_ACCOUNT)\" container_cmd=\"$(CLUSTER_CONTAINER_CMD_WITHOUT_BIND)\" time=\"2-00:00:00\" partition=\"$(SENSITIVE_CLUSTER_PARTITION)\" working_path=\"$(RNASEQ_CLUSTER_BGEECALL_OUTPUT)\" decrypt_file_path=\"\" >> $@.tmp 2>> $@.err
	@$(MV) $@.tmp $@

## Quality control by cell population to perform the calls (not pass if not bimodal)
QC_cellPopulation: presence_absence_bgeecall 1Run/QC_cellPopulation.R
	@echo --- Starting the quality control per cell population ---
	@sed -i 's@--output=.*@--output=${PWD}/QC_cellPopulation.out@'                                                     1Run/QC_cellPopulation.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/QC_cellPopulation.err@'                                                       1Run/QC_cellPopulation.sbatch
	@sed -i 's@--partition=.*@--partition=${SENSITIVE_CLUSTER_PARTITION}@'                                             1Run/QC_cellPopulation.sbatch
	@sed -i 's@--account=.*@--account=${SENSITIVE_CLUSTER_ACCOUNT}@'                                                   1Run/QC_cellPopulation.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                                      1Run/QC_cellPopulation.sbatch
	@sed -i 's@export scrna_seq_sample_info=.*@export scrna_seq_sample_info=$(SC_RNASEQ_SAMPINFO_FILEPATH)@'           1Run/QC_cellPopulation.sbatch
	@sed -i 's@export cells_folder=.*@export cells_folder=$(SC_RNASEQ_CLUSTER_BGEECALL)@'                              1Run/QC_cellPopulation.sbatch
	@sed -i 's@export sample_info_pass=.*@export sample_info_pass=$(SC_RNASEQ_SAMPINFO_PASS_FILEPATH)@'                1Run/QC_cellPopulation.sbatch
	@sed -i 's@export sample_info_discarded=.*@export sample_info_discarded=$(SC_RNASEQ_SAMPINFO_NOT_PASS_FILEPATH)@'  1Run/QC_cellPopulation.sbatch
	@sed -i 's@export modality_info=.*@export modality_info=$(SC_RNASEQ_MODALITY_CELL_POP)@'                           1Run/QC_cellPopulation.sbatch
	@sed -i 's@export calls_file_name=.*@export calls_file_name=$(ABUNDANCEFILE)@'                                     1Run/QC_cellPopulation.sbatch
	@sed -i 's@export plot=.*@export plot=yes@'                                                                        1Run/QC_cellPopulation.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_QC_CELLPOP_RES)@'                                        1Run/QC_cellPopulation.sbatch
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
	# Touch all files so that they are not removed from $(RNASEQ_CLUSTER_DATA)
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
	@perl 3Insertion/insert_scrna_seq.pl -bgee=$(BGEECMD) -scRnaSeqExperiment=$(SC_RNASEQ_EXPERIMENT_FILEPATH) -library_info=$(SC_RNASEQ_SAMPINFO_PASS_FILEPATH) -excluded_libraries=$(SC_RNASEQ_SAMPINFO_NOT_PASS_FILEPATH) -library_stats=$(SC_RNASEQ_SAMP_STATS_FL) -report_info=$(SC_RNASEQ_REPORT_INFO_FL) -all_results=$(SC_RNASEQ_ALL_RES_FL) -sex_info=$(UBERON_SEX_INFO_FILE_PATH) -extraMapping=$(EXTRAMAPPING_FILEPATH) -Aport=$(IDMAPPINGPORT) -Sport=$(STGMAPPINGPORT) > $@.tmp 2>warnings.$@
	@$(IDMAPPING)  $(IDMAPPINGPORT) &
	@$(STGMAPPING) $(STGMAPPINGPORT) &
	@sleep 50 # sleep because mappers need time to load Uberon
	# Insert RNA Seq data
	@perl 3Insertion/insert_scrna_seq.pl -bgee=$(BGEECMD) -scRnaSeqExperiment=$(SC_RNASEQ_EXPERIMENT_FILEPATH) -library_info=$(SC_RNASEQ_SAMPINFO_PASS_FILEPATH) -excluded_libraries=$(SC_RNASEQ_SAMPINFO_NOT_PASS_FILEPATH) -library_stats=$(SCRNASEQSAMPSTATS) -report_info=$(SCRNASEQREPORTINFO) -all_results=$(SCRNASEQALLRES) -sex_info=$(UBERON_SEX_INFO_FILE_PATH) -extraMapping=$(EXTRAMAPPING_FILEPATH) -Aport=$(IDMAPPINGPORT) -Sport=$(STGMAPPINGPORT) > $@.tmp 2>warnings.$@
	@echo "Delete RNA-Seq experiments for which no RNA-Seq libraries have been inserted" >> $@.tmp
	@$(MYSQL) -e "DELETE t1 FROM rnaSeqExperiment AS t1 WHERE NOT EXISTS (SELECT 1 FROM rnaSeqLibrary AS t2 WHERE t1.rnaSeqExperimentId = t2.rnaSeqExperimentId)" >> $@.tmp
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

