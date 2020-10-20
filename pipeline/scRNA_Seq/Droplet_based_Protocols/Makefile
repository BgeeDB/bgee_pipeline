PIPELINEROOT := ../
DIR_NAME := scRNA_Seq/

include $(PIPELINEROOT)Makefile.common

################## STEPS TO RUN ON AXIOM SERVER ##################

## Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations and store them in the source_files directory
## NOTE: this files should be already in the source_files (see whats happen with download)
get_annot:
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeqLibrary.tsv               && $(MV) scRNASeqLibrary.tsv                $(INPUT_DIR)/scRNASeqLibrary.tsv                2>/dev/null  || rm -f scRNASeqLibrary.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeqExperiment.tsv            && $(MV) scRNASeqExperiment.tsv             $(INPUT_DIR)/scRNASeqExperiment.tsv             2>/dev/null  || rm -f scRNASeqExperiment.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeq_barcode.tsv              && $(MV) scRNASeq_barcode.tsv               $(INPUT_DIR)/scRNASeq_barcode.tsv               2>/dev/null  || rm -f scRNASeq_barcode.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeq_markers_10X.tsv          && $(MV) scRNASeq_markers_10X.tsv           $(INPUT_DIR)/scRNASeq_markers_10X.tsv           2>/dev/null  || rm -f scRNASeq_markers_10X.tsv
	@touch $@

## Retrieve metadata (run rule with sbatch)
retrieve_metadata: get_annot 0Preparation/retrieve_metadata.R
	@echo --- Retrieve metadata information ---
	@sed -i 's@\(BASE *= \).*@\1"$(SC_RNASEQ_DOWNLOAD_PATH)";@'                                             0Preparation/retrieve_metadata.R
	@sed -i 's@--output=.*@--output=${PWD}/retrieve_metadata.out@'                                          0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/retrieve_metadata.err@'                                            0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                            0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export scRNASeqExperiment=.*@export scRNASeqExperiment=$(SC_RNASEQ_EXPERIMENT_FILEPATH)@'    0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export scRNASeqLibrary=.*@export scRNASeqLibrary=$(SC_RNASEQ_LIB_FILEPATH)@'                 0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_PATH_GENERATED)@'                                 0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                              0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export R_LIBS_USER=.*@export R_LIBS_USER=$(R_LIBS_PATH_AXIOM)@'                              0Preparation/retrieve_metadata.sbatch
	@sbatch 0Preparation/retrieve_metadata.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@


## Download cleaning data that still is not downloaded/present in Jura server (run rule with sbatch)
## NOTE: Here in the Make file we have one rule after the another, but the download from the 3 repositories can be done in parallel

## Download from HCA
download_HCA: retrieve_metadata 0Preparation/download_HCA.sbatch
	@echo --- Starting downloading the data from Human Cell Atlas ---
	@sed -i 's@\(BASE *= \).*@\1"$(SC_RNASEQ_DOWNLOAD_PATH)";@'                                                       0Preparation/download_HCA.sbatch
	@sed -i 's@--output=.*@--output=${PWD}/download_HCA.out@'                                                         0Preparation/download_HCA.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/download_HCA.err@'                                                           0Preparation/download_HCA.sbatch
	@sed -i 's@export manifest_file=.*@export manifest_file=${SC_RNASEQ_MANIFEST_FILE_FILEPATH}@'                     0Preparation/download_HCA.sbatch
	@sed -i 's@export tmp_folder_Download_data=.*@export tmp_folder_Download_data=/tmp/DOWNLOAD_HCA_DATA.$RANDOM@'    0Preparation/download_HCA.sbatch
	@sed -i 's@export final_destination=.*@export final_destination=${SC_RNASEQ_DOWNLOAD_PATH_DROPLET}@'              0Preparation/download_HCA.sbatch
	@sbatch 0Preparation/download_HCA.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Download from EBI arrayExpress
download_EBI: retrieve_metadata 0Preparation/download_reformat_EBI.R
	@echo --- Starting downloading the data from EBI ArrayExpress ---
	@sed -i 's@\(BASE *= \).*@\1"$(SC_RNASEQ_DOWNLOAD_PATH)";@'                                                   0Preparation/download_reformat_EBI.R
	@sed -i 's@--output=.*@--output=${PWD}/download_reformat_EBI.out@'                                            0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/download_reformat_EBI.err@'                                              0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                  0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export metadata_info=.*@export metadata_info=$(SC_RNASEQ_METADATA_10X_FILEPATH)@'                  0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export librariesDownloadedJura=.*@export librariesDownloadedJura=$(SC_RNASEQ_JURA_LIB_FILEPATH)@'  0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export bamtofastq=.*@export bamtofastq=$(SCRNASEQ_SOFTWARE_BAMTOFASTQ)@'                           0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                    0Preparation/download_reformat_EBI.sbatch
	@sbatch 0Preparation/download_reformat_EBI.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Download from SRA
download_SRA: retrieve_metadata 0Preparation/download_SRA.R
	@echo --- Starting downloading the data from SRA ---
	@sed -i 's@\(BASE *= \).*@\1"$(SC_RNASEQ_DOWNLOAD_PATH)";@'                                                   0Preparation/download_SRA.R
	@sed -i 's@--output=.*@--output=${PWD}/download_SRA.out@'                                                     0Preparation/download_SRA.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/download_SRA.err@'                                                       0Preparation/download_SRA.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                  0Preparation/download_SRA.sbatch
	@sed -i 's@export metadata_info=.*@export metadata_info=$(SC_RNASEQ_METADATA_10X_FILEPATH)@'                  0Preparation/download_SRA.sbatch
	@sed -i 's@export librariesDownloadedJura=.*@export librariesDownloadedJura=$(SC_RNASEQ_JURA_LIB_FILEPATH)@'  0Preparation/download_SRA.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                0Preparation/download_SRA.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                    0Preparation/download_SRA.sbatch
	@sbatch 0Preparation/download_SRA.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Add new libraries download to the jura file list (run rule in front)
#NOTE *make -j3 list_new_downloads*  should run the 3 downloads in parallel
list_new_downloads: download_SRA download_EBI download_HCA
	@find $(SC_RNASEQ_DOWNLOAD_PATH_DROPLET) -type f -name \*.fastq.gz\* | xargs -r dirname | sed -e 's@^.*/@@' | sort -u > /tmp/new_downloads
	@cat $(SC_RNASEQ_JURA_LIB_FILEPATH) >>/tmp/new_downloads
	@sort -u /tmp/new_downloads >$(SC_RNASEQ_JURA_LIB_FILEPATH)
	@rm -f /tmp/new_downloads
	@$(GIT) add $(SC_RNASEQ_JURA_LIB_FILEPATH)
	@$(GIT) commit -m 'Add new downloaded libraries' $(SC_RNASEQ_JURA_LIB_FILEPATH) || true
	@touch $@

# Commit the library information file that will be used for the rest of the pipeline
commit_annotation_and_metadata: list_new_downloads
	@$(GIT) add $(INPUT_DIR)scRNASeqLibrary.tsv
	@$(GIT) add $(INPUT_DIR)scRNASeqExperiment.tsv
	@$(GIT) add $(INPUT_DIR)scRNASeq_barcode.tsv
	@$(GIT) add $(INPUT_DIR)scRNASeq_markers_10X.tsv
	@$(GIT) add $(OUTPUT_DIR)metadata_info_10X.txt
	@$(GIT) add $(OUTPUT_DIR)metadata_notMatch_10X.txt
	@$(GIT) commit -m 'Update metadata for droplet-based protocols for scRNASeq for $(DBNAME)' $(INPUT_DIR)scRNASeqLibrary.tsv $(INPUT_DIR)scRNASeqExperiment.tsv $(INPUT_DIR)scRNASeq_barcode.tsv $(INPUT_DIR)scRNASeq_markers_10X.tsv $(OUTPUT_DIR)metadata_info_10X.txt $(OUTPUT_DIR)metadata_notMatch_10X.txt || true
	@$(GIT) push
	@echo -e "All information is ready, you can go to JURA cluster to continue the pipeline of the droplet protocols\n"
	@touch $@


## NOTE: Copy all git repository to JURA!
################## STEPS TO RUN ON JURA SERVER ##################

clusterJURA:
	@echo -e "\tBe sure everything is up-to-date before running single cell RNASeq pipeline for Droplet protocols\n"
	@touch $@


check_tools: clusterJURA
	@echo -e "\tGo to 'cd pipeline/scRNA_Seq/Droplet_based_Protocols/' and be prepared to work\n"
	@echo -e "\n\tRun this command to give access to all modules installed on vital-it\n\tmodule add Bioinformatics/Software/vital-it\n"
	# Check if all required tools/libs are available
	@module add Bioinformatics/Software/vital-it || true
	@$(CLUSTER_R_CMD) which R                          > $@.tmp
	@$(CLUSTER_R_CMD) R -e 'library("data.table")'    >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("stringr")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("dplyr")'         >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("devtools")'      >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("ggplot2")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("magrittr")'      >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("gridExtra")'     >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("biomaRt")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("UniprotR")'      >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("lattice")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("Matrix")'        >> $@.tmp  2>/dev/null
	### Note this 3 packages maybe are not installed in Vital-it module
	@$(CLUSTER_R_CMD) R -e 'library("BUSpaRse")'      >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("DropletUtils")'  >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("Seurat")'        >> $@.tmp  2>/dev/null
	@which xz                                         >> $@.tmp
	@which sbatch                                     >> $@.tmp
	@$(CLUSTER_KALLISTO_CMD)   which kallisto         >> $@.tmp
	@$(CLUSTER_BUSTOOLS_CMD)   which bustools         >> $@.tmp
	@$(MV) $@.tmp $@

## generate informative files: transcript_to_gene_with_intergenic + gene_to_biotype_with_intergenic for each species
generate_info: check_tools 0Run/generateInfo.R
	@echo --- Start generating the informative files ---
	@sed -i 's@\(BASE *= \).*@\1"$(SENSITIVE_GENERAL_PATH)";@'                     1Run/generateInfo.R
	@sed -i 's@--output=.*@--output=${PWD}/generateInfo.out@'                      1Run/generateInfo.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/generateInfo.err@'                        1Run/generateInfo.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                   1Run/generateInfo.sbatch
	@sed -i 's@export folder_gtf=.*@export folder_gtf=$(RNASEQ_CLUSTER_GTF)@'      1Run/generateInfo.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'   1Run/generateInfo.sbatch
	@sbatch 1Run/generateInfo.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@


## Run Kallisto bus per library
kallisto_bus: generate_info 1Run/Kallisto_bus.R
	@echo --- Start running kallisto bus for all libraries ---
	@sed -i 's@\(BASE *= \).*@\1"$(SENSITIVE_GENERAL_PATH)";@'                                    1Run/Kallisto_bus.R
	@sed -i 's@--output=.*@--output=${PWD}/Kallisto_bus.out@'                                     1Run/Kallisto_bus.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/Kallisto_bus.err@'                                       1Run/Kallisto_bus.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                  1Run/Kallisto_bus.sbatch
	@sed -i 's@export metadata_file=.*@export metadata_file=$(SC_RNASEQ_METADATA_10X_FILEPATH)@'  1Run/Kallisto_bus.sbatch
	@sed -i 's@export annotation_file=.*@export annotation_file=$(SC_RNASEQ_LIB_FILEPATH)@'       1Run/Kallisto_bus.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_FASTQ_DROPLET)@'              1Run/Kallisto_bus.sbatch
	@sed -i 's@export folderSupport=.*@export folderSupport=$(RNASEQ_CLUSTER_GTF)@'               1Run/Kallisto_bus.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_PATH_GENERATED)@'                       1Run/Kallisto_bus.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                  1Run/Kallisto_bus.sbatch
	@sbatch 1Run/Kallisto_bus.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

# NOTE: Before run the process_busFile rule (uncompress the barcodes files, because in github we cannot have files with more than 100Mb size)
uncompress_barcodes: kallisto_bus $SC_RNASEQ_PATH_SOURCE/barcode_whitelist_10X_v2.txt.zip $SC_RNASEQ_PATH_SOURCE/barcode_whitelist_10X_v3.txt.zip
	cd $(SC_RNASEQ_PATH_SOURCE)
	unzip barcode_whitelist_10X_v2.txt.zip
	unzip barcode_whitelist_10X_v3.txt.zip

## Process bus files
process_busFile: uncompress_barcodes 1Run/process_busFile.R
	@echo --- Starting the process of bus files ---
	@sed -i 's@\(BASE *= \).*@\1"$(SENSITIVE_GENERAL_PATH)";@'                                    1Run/process_busFile.R
	@sed -i 's@--output=.*@--output=${PWD}/process_busFile.out@'                                  1Run/process_busFile.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/process_busFile.err@'                                    1Run/process_busFile.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                  1Run/process_busFile.sbatch
	@sed -i 's@export scRNASeq_Info=.*@export scRNASeq_Info=$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)@'  1Run/process_busFile.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_FASTQ_DROPLET)@'              1Run/process_busFile.sbatch
	@sed -i 's@export folderSupport=.*@export folderSupport=$(RNASEQ_CLUSTER_GTF)@'               1Run/process_busFile.sbatch
	@sed -i 's@export whiteList_Path=.*@export whiteList_Path=$(SC_RNASEQ_PATH_SOURCE)@'          1Run/process_busFile.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                  1Run/process_busFile.sbatch
	@sbatch 1Run/process_busFile.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Quality control + cell type identification
qc_cellType: process_busFile 1Run/QC_CellTypeIdentification.R
	@echo --- Starting the quality control and the cell type identification ---
	@sed -i 's@\(BASE *= \).*@\1"$(SENSITIVE_GENERAL_PATH)";@'                                     1Run/QC_CellTypeIdentification.R
	@sed -i 's@--output=.*@--output=${PWD}/QC_CellTypeIdentification.out@'                         1Run/QC_CellTypeIdentification.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/QC_CellTypeIdentification.err@'                           1Run/QC_CellTypeIdentification.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                   1Run/QC_CellTypeIdentification.sbatch
	@sed -i 's@export scRNASeq_Info=.*@export scRNASeq_Info=$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)@'   1Run/QC_CellTypeIdentification.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_FASTQ_DROPLET)@'               1Run/QC_CellTypeIdentification.sbatch
	@sed -i 's@export infoFolder=.*@export infoFolder=$(SC_RNASEQ_PATH_SOURCE)@'                   1Run/QC_CellTypeIdentification.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_CLUSTER_QC_CELLTYPE)@'                   1Run/QC_CellTypeIdentification.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                   1Run/QC_CellTypeIdentification.sbatch
	@sbatch 1Run/QC_CellTypeIdentification.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Sum raw UMI per cell population and the correspondent normalized values
sum_raw_UMI: qc_cellType 1Run/sum_raw_UMI.R
	@echo --- Starting the summing of the raw UMI counts ---
	@sed -i 's@\(BASE *= \).*@\1"$(SENSITIVE_GENERAL_PATH)";@'                                                         1Run/sum_raw_UMI.R
	@sed -i 's@--output=.*@--output=${PWD}/sum_raw_UMI.out@'                                                           1Run/sum_raw_UMI.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/sum_raw_UMI.err@'                                                             1Run/sum_raw_UMI.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                       1Run/sum_raw_UMI.sbatch
	@sed -i 's@export scRNASeq_Info=.*@export scRNASeq_Info=$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)@'                       1Run/sum_raw_UMI.sbatch
	@sed -i 's@export InformationAllLibraries=.*@export InformationAllLibraries=$(SC_RNASEQ_INFO_ALL_LIBRARIES_10X)@'  1Run/sum_raw_UMI.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_FASTQ_DROPLET)@'                                   1Run/sum_raw_UMI.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_CLUSTER_SUM_RAW_UMI)@'                                       1Run/sum_raw_UMI.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                                       1Run/sum_raw_UMI.sbatch
	@sbatch 1Run/sum_raw_UMI.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@


## calls per individual cell and cell population
call_cell_cellPop: sum_raw_UMI 1Run/Call_PresentGenes_indivCell_cellPop.R
	@echo --- Starting the calls per individual cell and cell population ---
	@sed -i 's@\(BASE *= \).*@\1"$(SENSITIVE_GENERAL_PATH)";@'                                                         1Run/Call_PresentGenes_indivCell_cellPop.R
	@sed -i 's@--output=.*@--output=${PWD}/Call_PresentGenes_indivCell_cellPop.out@'                                   1Run/Call_PresentGenes_indivCell_cellPop.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/Call_PresentGenes_indivCell_cellPop.err@'                                     1Run/Call_PresentGenes_indivCell_cellPop.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                       1Run/Call_PresentGenes_indivCell_cellPop.sbatch
	@sed -i 's@export scRNASeq_Info=.*@export scRNASeq_Info=$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)@'                       1Run/Call_PresentGenes_indivCell_cellPop.sbatch
	@sed -i 's@export InformationAllLibraries=.*@export InformationAllLibraries=$(SC_RNASEQ_INFO_ALL_LIBRARIES_10X)@'  1Run/Call_PresentGenes_indivCell_cellPop.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_FASTQ_DROPLET)@'                                   1Run/Call_PresentGenes_indivCell_cellPop.sbatch
	@sed -i 's@export desired_r_cutoff=.*@export desired_r_cutoff=$(SC_RNASEQ_RATIO_PROPORTION)@'                      1Run/Call_PresentGenes_indivCell_cellPop.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_CLUSTER_CALLS_10X)@'                                         1Run/Call_PresentGenes_indivCell_cellPop.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                                       1Run/Call_PresentGenes_indivCell_cellPop.sbatch
	@sbatch 1Run/Call_PresentGenes_indivCell_cellPop.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

final_status: call_cell_cellPop
	@$(GIT) status

tar_all: final_status
	cd $(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)
	tar -cvfSp tarball_scRNASeq_droplet_10X.tar .
	gzip -9 tarball_scRNASeq_droplet_10X.tar
	# TODO: cp tarball_scRNASeq_droplet_10X.tar.gz to archive

