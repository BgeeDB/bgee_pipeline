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
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeq_barcode.tsv              && $(MV) scRNASeq_barcode.tsv               $(INPUT_DIR)/scRNASeq_barcode.tsv               2>/dev/null  || rm -f scRNASeq_barcode.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeq_markers_10X.tsv          && $(MV) scRNASeq_markers_10X.tsv           $(INPUT_DIR)/scRNASeq_markers_10X.tsv           2>/dev/null  || rm -f scRNASeq_markers_10X.tsv
	@touch $@

## Retrieve metadata (run rule with sbatch)
retrieve_metadata: get_annot 0Preparation/retrieve_metadata.R
	@echo --- Retrieve metadata information ---
	@sed -i 's@--output=.*@--output=${PWD}/retrieve_metadata.out@'                                                      0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/retrieve_metadata.err@'                                                        0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION}@'                                                    	0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT}@'                                                          	0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                        0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export scRNASeqExperiment=.*@export scRNASeqExperiment=../$(SC_RNASEQ_EXPERIMENT_FILEPATH)@'             0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export scRNASeqLibrary=.*@export scRNASeqLibrary=../$(SC_RNASEQ_LIB_FILEPATH)@'                          0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export metadata_file_path=.*@export output=../$(SC_RNASEQ_METADATA_10X_FILEPATH)@'                       0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export metadata_notmatch_file_path=.*@export output=../$(SC_RNASEQ_METADATA_10X_NOT_MATCH_FILEPATH)@'    0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                          0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export R_LIBS_USER=.*@export R_LIBS_USER=$(R_LIBS_PATH_AXIOM)@'                                          0Preparation/retrieve_metadata.sbatch
	@sbatch 0Preparation/retrieve_metadata.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@


## Download cleaning data that still is not downloaded/present in Jura server (run rule with sbatch)
## NOTE: Here in the Make file we have one rule after the another, but the download from the 3 repositories can be done in parallel

## Download from HCA
#XXX: Why this rule is the only one already using the dcsr environment?
download_HCA: retrieve_metadata 0Preparation/download_HCA.sbatch
	@echo --- Starting downloading the data from Human Cell Atlas ---
	@sed -i 's@--output=.*@--output=${PWD}/download_HCA.out@'                                                         0Preparation/download_HCA.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/download_HCA.err@'                                                           0Preparation/download_HCA.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION}@'                                                      0Preparation/download_HCA.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT}@'                                                            0Preparation/download_HCA.sbatch
	@sed -i 's@export manifest_file=.*@export manifest_file=../${SC_RNASEQ_MANIFEST_FILE_FILEPATH}@'                  0Preparation/download_HCA.sbatch
	@sed -i 's@export tmp_folder_Download_data=.*@export tmp_folder_Download_data=/tmp/DOWNLOAD_HCA_DATA.$$RANDOM@'   0Preparation/download_HCA.sbatch
	@sed -i 's@export final_destination=.*@export final_destination=${SC_RNASEQ_DOWNLOAD_PATH_DROPLET}@'              0Preparation/download_HCA.sbatch
	@sbatch 0Preparation/download_HCA.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Download from EBI arrayExpress
download_EBI: retrieve_metadata 0Preparation/download_reformat_EBI.R
	@echo --- Starting downloading the data from EBI ArrayExpress ---
	@sed -i 's@--output=.*@--output=${PWD}/download_reformat_EBI.out@'                                               0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/download_reformat_EBI.err@'                                                 0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION}@'                                                     0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT}@'                                                           0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                     0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export metadata_info=.*@export metadata_info=../$(SC_RNASEQ_METADATA_10X_FILEPATH)@'                  0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export librariesDownloadedJura=.*@export librariesDownloadedJura=../$(SC_RNASEQ_JURA_LIB_FILEPATH)@'  0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                   0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export bamtofastq=.*@export bamtofastq=$(SCRNASEQ_SOFTWARE_BAMTOFASTQ)@'                              0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                       0Preparation/download_reformat_EBI.sbatch
	@sbatch 0Preparation/download_reformat_EBI.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Download from SRA
download_SRA: retrieve_metadata 0Preparation/download_SRA.R
	@echo --- Starting downloading the data from SRA ---
	@sed -i 's@--output=.*@--output=${PWD}/download_SRA.out@'                                                        0Preparation/download_SRA.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/download_SRA.err@'                                                          0Preparation/download_SRA.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION}@'                                                     0Preparation/download_SRA.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT}@'                                                           0Preparation/download_SRA.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                     0Preparation/download_SRA.sbatch
	@sed -i 's@export metadata_info=.*@export metadata_info=../$(SC_RNASEQ_METADATA_10X_FILEPATH)@'                  0Preparation/download_SRA.sbatch
	@sed -i 's@export librariesDownloadedJura=.*@export librariesDownloadedJura=../$(SC_RNASEQ_JURA_LIB_FILEPATH)@'  0Preparation/download_SRA.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                   0Preparation/download_SRA.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                       0Preparation/download_SRA.sbatch
	@sbatch 0Preparation/download_SRA.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Add new libraries download to the jura file list (run rule in front)
#NOTE *make -j3 list_new_downloads*  should run the 3 downloads in parallel
list_new_downloads: download_SRA download_EBI download_HCA
	@find $(SC_RNASEQ_DOWNLOAD_PATH_DROPLET) -type f -name \*.fastq.gz\* | xargs -r dirname | sed -e 's@^.*/@@' | sort -u > /tmp/new_downloads
	@cat ../$(SC_RNASEQ_JURA_LIB_FILEPATH) >>/tmp/new_downloads
	@sort -u /tmp/new_downloads >../$(SC_RNASEQ_JURA_LIB_FILEPATH)
	@rm -f /tmp/new_downloads
	@$(GIT) add ../$(SC_RNASEQ_JURA_LIB_FILEPATH)
	@$(GIT) commit -m 'Add new downloaded libraries' ../$(SC_RNASEQ_JURA_LIB_FILEPATH) || true
	@touch $@

# Commit the library information file that will be used for the rest of the pipeline
commit_annotation_and_metadata: list_new_downloads
	@$(GIT) add ../$(INPUT_DIR)scRNASeqLibrary.tsv
	@$(GIT) add ../$(INPUT_DIR)scRNASeqExperiment.tsv
	@$(GIT) add ../$(INPUT_DIR)scRNASeq_barcode.tsv
	@$(GIT) add ../$(INPUT_DIR)scRNASeq_markers_10X.tsv
	@$(GIT) add ../$(OUTPUT_DIR)metadata_info_10X.txt
	@$(GIT) add ../$(OUTPUT_DIR)metadata_notMatch_10X.txt
	@$(GIT) commit -m 'Update metadata for droplet-based protocols for scRNASeq for $(DBNAME)' ../$(INPUT_DIR)scRNASeqLibrary.tsv ../$(INPUT_DIR)scRNASeqExperiment.tsv ../$(INPUT_DIR)scRNASeq_barcode.tsv ../$(INPUT_DIR)scRNASeq_markers_10X.tsv ../$(OUTPUT_DIR)metadata_info_10X.txt ../$(OUTPUT_DIR)metadata_notMatch_10X.txt || true
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
	@export R_LIBS_USER=$(R_LIBS_PATH_JURA); $(CLUSTER_R_CMD) R -e 'library("UniprotR")'      >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("lattice")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("Matrix")'        >> $@.tmp  2>/dev/null
	@export R_LIBS_USER=$(R_LIBS_PATH_JURA); $(CLUSTER_R_CMD) R -e 'library("BUSpaRse")'      >> $@.tmp  2>/dev/null
	@export R_LIBS_USER=$(R_LIBS_PATH_JURA); $(CLUSTER_R_CMD) R -e 'library("DropletUtils")'  >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("Seurat")'        >> $@.tmp  2>/dev/null
	@which xz                                         >> $@.tmp
	@which sbatch                                     >> $@.tmp
	@$(CLUSTER_KALLISTO_CMD)   which kallisto         >> $@.tmp
	@$(CLUSTER_BUSTOOLS_CMD)   which bustools         >> $@.tmp
	@$(MV) $@.tmp $@

## generate informative files: transcript_to_gene_with_intergenic + gene_to_biotype_with_intergenic for each species
## not the informative files contain intergenic regions
##XXX : tx2gene and gene2biotype files are generated by BgeeCall during the RNA-Seq pipeline. It is propably better to reuse them rather than create new on. Are geneid2geneName mandatory?
generate_info: check_tools 0Preparation/generateInfo.R
	@echo --- Start generating the informative files ---
	@sed -i 's@--output=.*@--output=${PWD}/generateInfo.out@'                                               0Preparation/generateInfo.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/generateInfo.err@'                                                 0Preparation/generateInfo.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION_SENSITIVE}@'                                  0Preparation/generateInfo.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT_SENSITIVE}@'                                        0Preparation/generateInfo.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                           0Preparation/generateInfo.sbatch
	@sed -i 's@export folder_gtf=.*@export folder_gtf=$(RNASEQ_CLUSTER_GTF)@'                               0Preparation/generateInfo.sbatch
	@sed -i 's@export metadata_info_file=.*@export metadata_info_file=${PWD}/$(SC_RNASEQ_METADATA_10X_FILEPATH)@'  0Preparation/generateInfo.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                            0Preparation/generateInfo.sbatch
	@sbatch 0Preparation/generateInfo.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@


## Run Kallisto bus per library
##TODO: since Bgee15 kallisto indexes containing transcripts and ref intergenic sequences are generated using BgeeCall. They will not be anymore in the $(RNASEQ_CLUSTER_GTF) folder. It has to be updated
kallisto_bus: generate_info 1Run/Kallisto_bus.R
	@echo --- Start running kallisto bus for all libraries ---
	@sed -i 's@--output=.*@--output=${PWD}/Kallisto_bus.out@'                                         1Run/Kallisto_bus.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/Kallisto_bus.err@'                                           1Run/Kallisto_bus.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION_SENSITIVE}@'                            1Run/Kallisto_bus.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT_SENSITIVE}@'                                  1Run/Kallisto_bus.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                     1Run/Kallisto_bus.sbatch
	@sed -i 's@export metadata_file=.*@export metadata_file=$(SC_RNASEQ_METADATA_10X_FILEPATH)@'      1Run/Kallisto_bus.sbatch
	@sed -i 's@export annotation_file=.*@export annotation_file=$(SC_RNASEQ_LIB_FILEPATH)@'           1Run/Kallisto_bus.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_FASTQ_DROPLET)@'                  1Run/Kallisto_bus.sbatch
	@sed -i 's@export folderSupport=.*@export folderSupport=$(RNASEQ_CLUSTER_GTF)@'                   1Run/Kallisto_bus.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_PATH_GENERATED)@'                           1Run/Kallisto_bus.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                      1Run/Kallisto_bus.sbatch
	@sbatch 1Run/Kallisto_bus.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

# NOTE: Before run the process_busFile rule (uncompress the barcodes files, because in github we cannot have files with more than 100Mb size)
uncompress_barcodes: kallisto_bus ../$(SC_RNASEQ_PATH_SOURCE)/barcode_whitelist_10X_v2.txt.zip ../$(SC_RNASEQ_PATH_SOURCE)/barcode_whitelist_10X_v3.txt.zip
	cd ../$(SC_RNASEQ_PATH_SOURCE) && \
	unzip barcode_whitelist_10X_v2.txt.zip && \
	unzip barcode_whitelist_10X_v3.txt.zip
	@touch $@

## Process bus files
process_busFile: uncompress_barcodes 1Run/process_busFile.R
	@echo --- Starting the process of bus files ---
	@sed -i 's@--output=.*@--output=${PWD}/process_busFile.out@'                                     1Run/process_busFile.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/process_busFile.err@'                                       1Run/process_busFile.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION_SENSITIVE}@'                           1Run/process_busFile.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT_SENSITIVE}@'                                 1Run/process_busFile.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                    1Run/process_busFile.sbatch
	@sed -i 's@export scRNASeq_Info=.*@export scRNASeq_Info=../$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)@'  1Run/process_busFile.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_FASTQ_DROPLET)@'                 1Run/process_busFile.sbatch
	@sed -i 's@export folderSupport=.*@export folderSupport=$(RNASEQ_CLUSTER_GTF)@'                  1Run/process_busFile.sbatch
	@sed -i 's@export whiteList_Path=.*@export whiteList_Path=../$(SC_RNASEQ_PATH_SOURCE)@'          1Run/process_busFile.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                     1Run/process_busFile.sbatch
	@sbatch 1Run/process_busFile.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Quality control (Knee + PCA + UMAP for each library) + cell type identification
qc_cellType: process_busFile 1Run/QC_CelltypeIdentification.R
	@echo --- Starting the quality control and the cell type identification ---
	@sed -i 's@--output=.*@--output=${PWD}/QC_CelltypeIdentification.out@'                            1Run/QC_CelltypeIdentification.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/QC_CelltypeIdentification.err@'                              1Run/QC_CelltypeIdentification.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION_SENSITIVE}@'                            1Run/QC_CelltypeIdentification.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT_SENSITIVE}@'                                  1Run/QC_CelltypeIdentification.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                     1Run/QC_CelltypeIdentification.sbatch
	@sed -i 's@export scRNASeq_Info=.*@export scRNASeq_Info=../$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)@'   1Run/QC_CelltypeIdentification.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_FASTQ_DROPLET)@'                  1Run/QC_CelltypeIdentification.sbatch
	@sed -i 's@export infoFolder=.*@export infoFolder=../$(SC_RNASEQ_PATH_SOURCE)@'                   1Run/QC_CelltypeIdentification.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_CLUSTER_QC_CELLTYPE)@'                      1Run/QC_CelltypeIdentification.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                      1Run/QC_CelltypeIdentification.sbatch
	@sed -i 's@export R_LIBS_USER=.*@export R_LIBS_USER=$(R_LIBS_PATH_JURA)@'                         1Run/QC_CelltypeIdentification.sbatch
	@sbatch 1Run/QC_CelltypeIdentification.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@
	
## Quality control based on the bimodality of the cell-population per library/cell-type pop
bimodality_cellPop: qc_cellType 1Run/QC_cellPop_bimodality.R
	@echo --- Starting the quality control for the cell population based on the bimodality ---
	@sed -i 's@--output=.*@--output=${PWD}/QC_cellPop_bimodality_TB.out@'                                 1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/QC_cellPop_bimodality_TB.err@'                                   1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION_SENSITIVE}@'                                1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT_SENSITIVE}@'                                      1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                         1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@export scRNASeq_Info=.*@export scRNASeq_Info=../$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)@'       1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_CLUSTER_QC_CELLTYPE)@'                1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@export output_folder=.*@export output_folder=$(SC_RNASEQ_CLUSTER_QC_CELLTYPE)@'            1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                          1Run/QC_cellPop_bimodality.sbatch
	@sbatch 1Run/QC_cellPop_bimodality.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@
	
## Call present and absent genes per library/cell population that pass the bimodality test
calls: bimodality_cellPop 1Run/Calls_cell_pop_per_library.R
	@echo --- Starting the calls of present and absent genes ---
	@sed -i 's@--output=.*@--output=${PWD}/Calls_library_cellPop.out@'                                                             1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/Calls_library_cellPop.err@'                                                               1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION_SENSITIVE}@'                                                         1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT_SENSITIVE}@'                                                               1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                                                  1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export scRNASeq_Info=.*@export scRNASeq_Info=../$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)@'                                1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export bimodalityFile=.*@export bimodalityFile=$(SC_RNASEQ_CLUSTER_QC_CELLTYPE)/bimodality_targetBased.txt@'        1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_CLUSTER_QC_CELLTYPE)@'                                         1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export folder_refIntergenic=.*@export folder_refIntergenic=../$(CLUSTER_REF_INTERGENIC_FOLDER)@'                    1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export desired_pValue_cutoff=.*@export desired_pValue_cutoff=../$(PVALUE_CUTOFF)@'                                  1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export output_folder=.*@export output_folder=$(SC_RNASEQ_CLUSTER_CALLS_10X)@'                                       1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                                                   1Run/Calls_cell_pop_per_library.sbatch
	@sbatch 1Run/Calls_cell_pop_per_library.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

final_status: calls
	@$(GIT) status
	
final_commit: final_status
	# Commit the scRNASeq info files after QC and information file about the stats of all libraries
	@$(GIT) add ../$(OUTPUT_DIR)InformationAllLibraries.txt
	@$(GIT) add ../$(OUTPUT_DIR)markerGenes_Validated.txt
	@$(GIT) add ../$(OUTPUT_DIR)variabilityClusters.txt
	@$(GIT) add ../$(OUTPUT_DIR)bimodality_targetBased.txt
	@$(GIT) add ../$(OUTPUT_DIR)All_cellPopulation_stats_10X.tsv
	@$(GIT) commit -m 'Update info files for scRNASeq target-based for $(DBNAME)' ../$(OUTPUT_DIR)InformationAllLibraries.txt ../$(OUTPUT_DIR)markerGenes_Validated.txt ../$(OUTPUT_DIR)variabilityClusters.txt ../$(OUTPUT_DIR)bimodality_targetBased.txt ../$(OUTPUT_DIR)All_cellPopulation_stats_10X.tsv || true
	@$(GIT) push
	@echo -e "All information is ready, you can make a tar of the results. \n"
	@touch $@

tar_all: final_commit
	cd $(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)
	tar -cvfSp tarball_scRNASeq_droplet_10X.tar .
	gzip -9 tarball_scRNASeq_droplet_10X.tar
	# TODO: cp tarball_scRNASeq_droplet_10X.tar.gz to archive

