PIPELINEROOT := ../../
DIR_NAME := scRNA_Seq/Droplet_based_Protocols

include $(PIPELINEROOT)Makefile.common

################## STEPS TO RUN ON THE CURNAGL SERVER ##################

## Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations submodule and store them in the source_files directory
get_annot:
	@$(GIT) submodule init
	@$(GIT) submodule update
	@$(CP) $(ANNOT_SUBMODULE_DIR)/Strains/StrainMapping.tsv                 $(STRAIN_MAPPING_FILE)
	@$(CP) $(ANNOT_SUBMODULE_DIR)/scRNA_Seq/scRNASeqTBLibrary.tsv           $(INPUT_DIR)/Target_based/scRNASeqTBLibrary.tsv
	@$(CP) $(ANNOT_SUBMODULE_DIR)/scRNA_Seq/scRNASeqExperiment.tsv          $(INPUT_DIR)/scRNASeqExperiment.tsv
	@$(CP) $(ANNOT_SUBMODULE_DIR)/scRNA_Seq/scRNASeq_barcode_*.tsv          $(INPUT_DIR)/Target_based/
	@touch $@

## Run in the front as it does not require lots of CPU and GPU
## Filter experiment/libraries to consider for the target-based pipeline
## This filtering :
##   - removes full-length experiment/libraries
##   - removes protocols not accepted in Bgee
##   - filter on speciesIds AND libraryIds if such arguments are provided (e.g speciesIds="7227,9606")
filter_annot:
	@$(CLUSTER_R_CMD) R CMD BATCH --no-save --no-restore '--args scRNASeqExperiment="$(SC_RNASEQ_EXPERIMENT_FILEPATH)" scRNASeqTBLibrary="$(SC_RNASEQ_LIB_TB_FILEPATH)" strainMapping="$(STRAIN_MAPPING_FILE)" acceptedProtocols="$(SC_RNASEQ_ACCEPTED_PROTOCOLS)" outputDir="$(SC_RNASEQ_TB_FOLDER_GENERATED)"' 0Preparation/filter_annotation.R $(SC_RNASEQ_TB_FOLDER_GENERATED)filter_annotation.Rout
	@touch $@

## Retrieve metadata (run rule with sbatch)
retrieve_metadata: get_annot 0Preparation/retrieve_metadata.R
	@echo --- Retrieve metadata information ---
	@sed -i 's@--output=.*@--output=${PWD}/retrieve_metadata.out@'                                                 0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/retrieve_metadata.err@'                                                   0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION}@'                                                   0Preparation/retrieve_metadata.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT}@'                                                         0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}@'                                                   0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export scRNASeqExperiment=.*@export scRNASeqExperiment=$(SC_RNASEQ_EXPERIMENT_FILEPATH_FILTERED)@'  0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export scRNASeqTBLibrary=.*@export scRNASeqTBLibrary=$(SC_RNASEQ_LIB_TB_FILEPATH_FILTERED)@'        0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export metadata_file=.*@export metadata_file=$(SC_RNASEQ_METADATA_10X_FILEPATH)@'                   0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export information_file=.*@export information_file=$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)@'             0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                     0Preparation/retrieve_metadata.sbatch
	@sed -i 's@export R_LIBS_USER=.*@export R_LIBS_USER=$(R_LIBS_PATH_CURNAGL)@'                                   0Preparation/retrieve_metadata.sbatch
	@sbatch 0Preparation/retrieve_metadata.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@


## Download cleaning data that still is not downloaded/present in Jura server (run rule with sbatch)
## NOTE: Here in the Make file we have one rule after the another, but the download from the 3 repositories can be done in parallel

## Download from HCA
#XXX: Why this rule is the only one already using the dcsr environment?
download_HCA: retrieve_metadata 0Preparation/download_HCA.R
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
	@sed -i 's@export librariesDownloadedJura=.*@export librariesDownloadedJura=../$(SC_RNASEQ_SENSITIVE_LIB_FILEPATH)@'  0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export output=.*@export output=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                   0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export bamtofastq=.*@export bamtofastq=$(SCRNASEQ_SOFTWARE_BAMTOFASTQ)@'                              0Preparation/download_reformat_EBI.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_DOWNLOAD_PATH_DROPLET)@'                                       0Preparation/download_reformat_EBI.sbatch
	@sbatch 0Preparation/download_reformat_EBI.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

# Rule used to parallelize download of libraries and to homogenize fastq file names
# it replaces the old rules not parallelized and downloading files depending on the sources
# (see scripts 0Preparation/download_SRA.R, 0Preparation/download_HCA.R and 0Preparation/download_reformat_EBI.R)
parallelized_download:
	@perl 0Preparation/parallelized_download_SRA.pl -metadataFile=$(SC_RNASEQ_METADATA_10X_FILEPATH) -parallelJobs=50 -outputDir=$(CLUSTER_CURNAGL_PATH)scRNA_Seq_All/scRNASeq_libraries_Droplet_10X -bamtofastq=$(SCRNASEQ_SOFTWARE_BAMTOFASTQ) -queue=$(CLUSTER_PARTITION) -account=$(CLUSTER_ACCOUNT) >$@.tmp 2> $@.err
	@$(MV) $@.tmp $@

## Add new libraries download to the jura file list (run rule in front)
#NOTE *make -j3 list_new_downloads*  should run the 3 downloads in parallel
list_new_downloads: download_SRA download_EBI download_HCA
	@find $(SC_RNASEQ_DOWNLOAD_PATH_DROPLET) -type f -name \*.fastq.gz\* | xargs -r dirname | sed -e 's@/FASTQ@@' -e 's@^.*/@@' | sort -u > /tmp/new_downloads
	@cat $(SC_RNASEQ_SENSITIVE_LIB_FILEPATH) >>/tmp/new_downloads
	@sort -u /tmp/new_downloads >$(SC_RNASEQ_SENSITIVE_LIB_FILEPATH)
	@rm -f /tmp/new_downloads
	@$(GIT) add $(SC_RNASEQ_SENSITIVE_LIB_FILEPATH)
	@$(GIT) commit -m 'Add new downloaded libraries' $(SC_RNASEQ_SENSITIVE_LIB_FILEPATH) || true
	@touch $@

# Commit the library information file that will be used for the rest of the pipeline
commit_annotation_and_metadata: list_new_downloads
	@$(GIT) add $(INPUT_DIR)/Target_based/scRNASeqTBLibrary.tsv
	@$(GIT) add $(INPUT_DIR)scRNASeqExperiment.tsv
	@$(GIT) add $(INPUT_DIR)/Target_based/scRNASeq_barcode_*
	@$(GIT) add $(SC_RNASEQ_TB_FOLDER_GENERATED)metadata_info_10X.txt
	@$(GIT) add $(SC_RNASEQ_TB_FOLDER_GENERATED)metadata_notMatch_10X.txt
	@$(GIT) commit -m 'Update metadata for droplet-based protocols for scRNASeq for $(DBNAME)' || true
	@$(GIT) push
	@echo -e "All information is ready, you can go to SENSITIVE cluster to continue the pipeline of the droplet protocols\n"
	@touch $@

## NOTE: Copy all git repository to sensitive cluster!
################## STEPS TO RUN ON SENSITIVE SERVER ##################

clusterSensitive:
	@echo -e "\tBe sure everything is up-to-date before running single cell RNASeq pipeline for Droplet protocols\n"
	@touch $@


check_tools: clusterSensitive
	@echo -e "\tGo to 'cd pipeline/scRNA_Seq/Droplet_based_Protocols/' and be prepared to work\n"
	@echo -e "\n\tRun this command to give access to all modules installed on vital-it\n\tmodule use /software/module/\n"
	# Check if all required tools/libs are available
	@module use /software/module/ || true
	@$(CLUSTER_R_CMD) which R                          > $@.tmp
	@$(CLUSTER_R_CMD) R -e 'library("biomaRt")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("Biostrings")'    >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("data.table")'    >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("devtools")'      >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("dplyr")'         >> $@.tmp  2>/dev/null
	# this package does not exist. Is it still used?
	#@$(CLUSTER_R_CMD) R -e 'library("HCAExplorer")'   >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("HelpersMG")'     >> $@.tmp  2>/dev/null
	#@$(CLUSTER_R_CMD) R -e 'library("ggExtra")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("ggplot2")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("gridExtra")'     >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("LaplacesDemon")' >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("lattice")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("magrittr")'      >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("Matrix")'        >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("mclust")'        >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("Seurat")'        >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("stringr")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("tibble")'        >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("tools")'         >> $@.tmp  2>/dev/null
	@export R_LIBS_USER=$(R_LIBS_PATH_SENSITIVE); $(CLUSTER_R_CMD) R -e 'library("BUSpaRse")'      >> $@.tmp  2>/dev/null
	@export R_LIBS_USER=$(R_LIBS_PATH_SENSITIVE); $(CLUSTER_R_CMD) R -e 'library("DropletUtils")'  >> $@.tmp  2>/dev/null
	@export R_LIBS_USER=$(R_LIBS_PATH_SENSITIVE); $(CLUSTER_R_CMD) R -e 'library("UniprotR")'      >> $@.tmp  2>/dev/null
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
	@sed -i 's@--output=.*@--output=${PWD}/generateInfo.out@'                                                      0Preparation/generateInfo.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/generateInfo.err@'                                                        0Preparation/generateInfo.sbatch
	@sed -i 's@--partition=.*@--partition=${SENSITIVE_CLUSTER_PARTITION}@'                                         0Preparation/generateInfo.sbatch
	@sed -i 's@--account=.*@--account=${SENSITIVE_CLUSTER_ACCOUNT}@'                                               0Preparation/generateInfo.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                                  0Preparation/generateInfo.sbatch
	@sed -i 's@export folder_gtf=.*@export folder_gtf=$(RNASEQ_CLUSTER_GTF)@'                                      0Preparation/generateInfo.sbatch
	@sed -i 's@export metadata_info_file=.*@export metadata_info_file=${PWD}/$(SC_RNASEQ_METADATA_10X_FILEPATH)@'  0Preparation/generateInfo.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                                   0Preparation/generateInfo.sbatch
	@sbatch 0Preparation/generateInfo.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Cleaning duplicates in the barcodes annotation files before mapping and analyis
cleaning_barcodes: generate_info 0Preparation/cleaning_barcodes.R
	@echo --- Cleaning duplicated barcodes per experimentID/libraryID ---
	@sed -i 's@--output=.*@--output=${PWD}/cleaning_barcodes.out@'                                               0Preparation/cleaning_barcodes.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/cleaning_barcodes.err@'                                                 0Preparation/cleaning_barcodes.sbatch
	@sed -i 's@--partition=.*@--partition=${SENSITIVE_CLUSTER_PARTITION}@'                                       0Preparation/cleaning_barcodes.sbatch
	@sed -i 's@--account=.*@--account=${SENSITIVE_CLUSTER_ACCOUNT}@'                                             0Preparation/cleaning_barcodes.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                                0Preparation/cleaning_barcodes.sbatch
	@sed -i 's@export barcodesFolder=.*@export barcodesFolder=$(SC_RNASEQ_TB_FOLDER_SOURCE)@'                    0Preparation/cleaning_barcodes.sbatch
	@sed -i 's@export output=.*@export output=${PWD}/$(SC_RNASEQ_CLEANED_BARCODES)@'                             0Preparation/cleaning_barcodes.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                                 0Preparation/cleaning_barcodes.sbatch
	@sbatch 0Preparation/cleaning_barcodes.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## generate a gtf file containing non matured transcripts. This file will be used only
## for single nucleus libraries
generate_gtf_nascent: cleaning_barcodes 0Preparation/generate_GTF_nascent.R
	@$(CLUSTER_R_CMD) R CMD BATCH --no-save --no-restore '--args gtf_dir="$(RNASEQ_CLUSTER_GTF)" metadata_file="$(SC_RNASEQ_METADATA_10X_FILEPATH)" scRNASeqInfoFile="$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)"' 0Preparation/generate_GTF_nascent.R $(SC_RNASEQ_TB_FOLDER_GENERATED)generate_GTF_nascent.Rout
	@touch $@

prepare_single_nucleus_index: generate_gtf_nascent 0Preparation/slurm_single_nucleus_index.pl
	#Preparing single nucleus indexed transcriptome for every species
	@perl 0Preparation/slurm_single_nucleus_index.pl -transcriptome_folder=$(RNASEQ_CLUSTER_GTF) -output_log_folder=$(SC_RNASEQ_TB_FOLDER_GENERATED) -account=$(SENSITIVE_CLUSTER_ACCOUNT) -partition=$(SENSITIVE_CLUSTER_PARTITION) -cluster_kallisto_cmd="$(CLUSTER_KALLISTO_CMD)" -cluster_tophat_cmd="$(CLUSTER_TOPHAT_CMD)" >$@.tmp 2>&1
	@$(MV) $@.tmp $@

## Run Kallisto bus per taxon_iD/library
kallisto_bus: prepare_single_nucleus_index 1Run/kallisto_bus_one_lib.R
	@echo --- Start running kallisto bus for all libraries ---
	@$(SENSITIVE_PERL_CMD) perl 1Run/parallelized_kallisto_bus.pl -metadataFile=$(SC_RNASEQ_METADATA_10X_FILEPATH) -parallelJobs=100 -fastqDir=$(SC_RNASEQ_FASTQ_DROPLET) -gtfDir=$(RNASEQ_CLUSTER_GTF) -scRNASeqInfoFile=$(SC_RNASEQ_SAMPINFO_10X_FILEPATH) -kallistoResults=$(SC_RNASEQ_CLUSTER_RES_KALLISTO_DROPLET) -queue=$(SENSITIVE_CLUSTER_PARTITION) -account=$(SENSITIVE_CLUSTER_ACCOUNT) -pathToScript=${PWD}/1Run/kallisto_bus_one_lib.R
	@touch $@

# generate stats from fastp and kallisto (e.g. read length, reads mapped, ...)
generate_stats: kallisto_bus 1Run/generate_stats.R
	@$(CLUSTER_R_CMD) R CMD BATCH --no-save --no-restore '--args metadata_file="$(SC_RNASEQ_METADATA_10X_FILEPATH)" kallisto_dir="$(SC_RNASEQ_CLUSTER_RES_KALLISTO_DROPLET)" fastq_dir="$(SC_RNASEQ_FASTQ_DROPLET)" output_file="$(SC_RNASEQ_STATS_LIBRARIES_10X)"' 1Run/generate_stats.R $(SC_RNASEQ_TB_FOLDER_GENERATED)generate_stats.Rout
	@touch $@

# NOTE: Before run the process_busFile rule (uncompress the barcodes files, because in github we cannot have files with more than 100Mb size)
uncompress_barcodes: generate_stats $(SC_RNASEQ_TB_FOLDER_SOURCE)barcode_whitelist_10X_v2.txt.zip $(SC_RNASEQ_TB_FOLDER_SOURCE)barcode_whitelist_10X_v3.txt.zip
	cd $(SC_RNASEQ_TB_FOLDER_SOURCE) && \
	unzip -o barcode_whitelist_10X_v2.txt.zip && \
	unzip -o barcode_whitelist_10X_v3.txt.zip
	@touch $@

## Process bus files
process_busFile: uncompress_barcodes 1Run/process_busfile_one_lib.R
	@echo --- Starting the process of bus files ---
	@perl 1Run/parallelized_process_busfile.pl -metadataFile=$(SC_RNASEQ_METADATA_10X_FILEPATH) -parallelJobs=100 -fastqDir=$(SC_RNASEQ_FASTQ_DROPLET) -gtfDir=$(RNASEQ_CLUSTER_GTF) -scRNASeqInfoFile=$(SC_RNASEQ_SAMPINFO_10X_FILEPATH) -whiteListPath=$(SC_RNASEQ_TB_FOLDER_SOURCE) -kallistoResults=$(SC_RNASEQ_CLUSTER_RES_KALLISTO_DROPLET) -queue=$(SENSITIVE_CLUSTER_PARTITION) -account=$(SENSITIVE_CLUSTER_ACCOUNT) -pathToScript=${PWD}/1Run/process_busfile_one_lib.R > $@.tmp 2> $@.err
	@mv $@.tmp $@


## Quality control (Knee + PCA + UMAP for each library) + cell type identification
qc_cellType: process_busFile 1Run/parallelized_qcCellIdentification.pl
	@echo --- Starting the quality control and the cell type identification ---
	@perl 1Run/parallelized_qcCellIdentification.pl -metadataFile=$(SC_RNASEQ_METADATA_10X_FILEPATH)  -parallelJobs=100 -gtfDir=$(RNASEQ_CLUSTER_GTF) -barcodeFolder=$(SC_RNASEQ_CLEANED_BARCODES) -kallistoResults=$(SC_RNASEQ_CLUSTER_RES_KALLISTO_DROPLET) -outputDir=$(SC_RNASEQ_CLUSTER_QC_CELLTYPE) -queue=$(SENSITIVE_CLUSTER_PARTITION) -account=$(SENSITIVE_CLUSTER_ACCOUNT) -pathToScript=${PWD}/1Run/celltype_one_lib.R -rLibs=$(R_LIBS_PATH_SENSITIVE) > $@.tmp 2>> $@.err
	@mv $@.tmp $@

## Quality control based on the bimodality of the cell-population per library/cell-type pop
## This QC is NOT USED ANYMORE to filter libraries because cell type annotation are not always
## very granular. It means a celltype can correspond to different cells. Then gene expression
## does not always follow a bimodal distribution.
## Kept this rule for sake of visualisation of the plots
bimodality_cellPop: qc_cellType 1Run/QC_cellPop_bimodality.R
	@echo --- Starting the quality control for the cell population based on the bimodality ---
	@sed -i 's@--output=.*@--output=${PWD}/QC_cellPop_bimodality_TB.out@'                                 1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/QC_cellPop_bimodality_TB.err@'                                   1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@--partition=.*@--partition=${SENSITIVE_CLUSTER_PARTITION}@'                                1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@--account=.*@--account=${SENSITIVE_CLUSTER_ACCOUNT}@'                                      1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                         1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@export scRNASeq_Info=.*@export scRNASeq_Info=$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)@'          1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_CLUSTER_QC_CELLTYPE)@'                1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@export output_folder=.*@export output_folder=$(SC_RNASEQ_CLUSTER_QC_CELLTYPE)@'            1Run/QC_cellPop_bimodality.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                          1Run/QC_cellPop_bimodality.sbatch
	@sbatch 1Run/QC_cellPop_bimodality.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

## Call present and absent genes per library/cell population that pass the bimodality test
## Not yet parallelized
calls: qc_cellType 1Run/Calls_cell_pop_per_library.R
	@echo --- Starting the calls of present and absent genes ---
	@sed -i 's@--output=.*@--output=${PWD}/Calls_library_cellPop.out@'                                                       1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/Calls_library_cellPop.err@'                                                         2Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@--partition=.*@--partition=${SENSITIVE_CLUSTER_PARTITION}@'                                                   1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@--account=.*@--account=${SENSITIVE_CLUSTER_ACCOUNT}@'                                                         1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export SCRIPT_PATH=.*@export SCRIPT_PATH=${PWD}/@'                                                            1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export scRNASeq_Info=.*@export scRNASeq_Info=$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)@'                             1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export bimodalityFile=.*@export bimodalityFile=$(SC_RNASEQ_CLUSTER_QC_CELLTYPE)/bimodality_targetBased.txt@'  1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export folder_data=.*@export folder_data=$(SC_RNASEQ_CLUSTER_QC_CELLTYPE)@'                                   1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export folder_gtf=.*@export folder_gtf=$(RNASEQ_CLUSTER_GTF)@'  1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export desired_pValue_cutoff=.*@export desired_pValue_cutoff=$(PVALUE_CUTOFF)@'                               1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export output_folder=.*@export output_folder=$(SC_RNASEQ_CLUSTER_CALLS_10X)@'                                 1Run/Calls_cell_pop_per_library.sbatch
	@sed -i 's@export ROUT=.*@export ROUT=$(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)@'                                             1Run/Calls_cell_pop_per_library.sbatch
	@sbatch 1Run/Calls_cell_pop_per_library.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo --- DONE ---
	@touch $@

final_status: calls
	@$(GIT) status

final_commit: final_status
	# Commit the scRNASeq info files after QC and information file about the stats of all libraries
	@$(GIT) add ../$(OUTPUT_DIR)$(SC_RNASEQ_SAMPINFO_10X_FILEPATH)
	@$(GIT) add ../$(OUTPUT_DIR)InformationAllLibraries.txt
	@$(GIT) add ../$(OUTPUT_DIR)bimodality_targetBased.txt
	@$(GIT) add ../$(OUTPUT_DIR)All_cellPopulation_stats_10X.tsv
	@$(GIT) commit -m 'Update info files for scRNASeq target-based for $(DBNAME)' ../$(OUTPUT_DIR)$(SC_RNASEQ_SAMPINFO_10X_FILEPATH) ../$(OUTPUT_DIR)InformationAllLibraries.txt ../$(OUTPUT_DIR)bimodality_targetBased.txt ../$(OUTPUT_DIR)All_cellPopulation_stats_10X.tsv || true
	@$(GIT) push
	@echo -e "All information is ready, you can make a tar of the results. \n"
	@touch $@

## Currently parallelize job using different thread to write in the database
## (argument coreNumber). Could optimize by also creating one Job per experiment
## (on unil cluster using slurm) if this rule takes too much time.
## Then, a perl script should take care of managing parallelization. It could also
## be possible to use the rslurm R package to manage parallelization
insert_cell_abundance:
	@perl 3Insertion/insert_scrna_seq.pl -bgee=$(BGEECMD) -targetBaseLibrary=$(SC_RNASEQ_LIB_TB_FILEPATH_FILTERED) -singleCellExperiment=$(SC_RNASEQ_EXPERIMENT_FILEPATH_FILTERED) -bgeeLibraryInfo=$(SC_RNASEQ_METADATA_10X_FILEPATH) -filteredBarcodeDir=$(SC_RNASEQ_CLEANED_BARCODES) -pipelineCallsSummary=$(SC_RNASEQ_SERVER_INFO_CALLS_10X) -pipelineReportFile=$(SC_RNASEQ_STATS_LIBRARIES_10X)  -kallistoResults=$(SC_RNASEQ_SERVER_RES_KALLISTO_DROPLET) -callsResults=$(SC_RNASEQ_SERVER_CALLS_10X) -sourceDir=$(SC_RNASEQ_TB_FOLDER_SOURCE) -sexInfo=$(UBERON_SEX_INFO_FILE_PATH) -numberCore=20 > $@.tmp 2> $@.err
	@if [[ ! -s $@.err ]]; then $(RM) $@.err; fi
	@$(MV) $@.tmp $@

insert_expression:
	@perl 3Insertion/insert_scrna_seq_expression.pl -bgee=$(BGEECMD) -number_threads=40> $@.tmp 2> $@.err
	@if [[ ! -s $@.err ]]; then $(RM) $@.err; fi
	@$(MV) $@.tmp $@

tar_all: final_commit
	cd $(SC_RNASEQ_CLUSTER_ALL_RES_DROPLET)
	tar -cvfSp tarball_scRNASeq_droplet_10X.tar .
	gzip -9 tarball_scRNASeq_droplet_10X.tar
	# TODO cp tarball_scRNASeq_droplet_10X.tar.gz to archive

