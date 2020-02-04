PIPELINEROOT := ../
DIR_NAME := scRNA_Seq/

include $(PIPELINEROOT)Makefile.common

################## STEPS RUN ON BGEE SERVER ##################

get_annot:
	# Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations. Store them in the source_files directory
	@$(WGET) $(ANNOTATION_GIT_URL)/scRNA_Seq/scRNASeqLibrary.tsv               && $(MV) scRNASeqLibrary.tsv                $(INPUT_DIR)/scRNASeqLibrary.tsv                2>/dev/null  || rm -f scRNASeqLibrary.tsv
	@touch $@

## Control the annotation 
control_annotation: get_annot 0Preparation/pre_process_control_annotation.R
	@echo --- Starting the control of the annotation  ---
	R CMD BATCH --no-save --no-restore '--args scRNASeqLibrary="${PWD}/$(INPUT_DIR)/scRNASeqLibrary.tsv" output_folder="${PWD}/$(OUTPUT_DIR)"' 0Preparation/pre_process_control_annotation.R 0Preparation/pre_process_control_annotation.Rout > $@.tmp
	@echo --- DONE ---
	@mv $@.tmp $@

## Retrieve metadata
retrieve_metadata: control_annotation 0Preparation/retrieve_metadata.R
	@echo --- Retrieve metadata information ---
	R CMD BATCH --no-save --no-restore '--args NEW_scRNASeqLibrary="${PWD}/$(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv" output_folder="${PWD}/$(OUTPUT_DIR)"' 0Preparation/retrieve_metadata.R 0Preparation/retrieve_metadata.Rout > $@.tmp
	@echo --- DONE ---
	@mv $@.tmp $@

commit_annotation_and_metadata: retrieve_metadata
	# Commit the library information file thta will be used for the rest of the pipeline
	@$(GIT) add $(INPUT_DIR)scRNASeqLibrary.tsv
	@$(GIT) add $(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv
	@$(GIT) add $(OUTPUT_DIR)metadata_info.txt
	@$(GIT) commit -m 'Update metadata for scRNASeq for $(DBNAME)' $(INPUT_DIR)scRNASeqLibrary.tsv $(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv $(OUTPUT_DIR)metadata_info.txt
	@$(GIT) push
	@echo -e "\t$(RNASEQ_SAMPINFO_FILEPATH) is ready, you can go to regular cluster to download new SRA files with get_SRA.pl as *admin* user\n\tDo a 'git pull' before starting\n"
	@$(MV) $@.tmp $@
	
################## STEPS RUN ON AXIOM ##################

pull_annotation_and_metadata: commit_annotation_and_metadata
	@$(GIT) pull
	@touch $@

## Download cleaning data
download_cleaning_data : pull_annotation_and_metadata 0Preparation/download_cleaning_data.R
	@echo --- Starting downloading the data ---
	R CMD BATCH --no-save --no-restore '--args metadata_info="${PWD}/$(OUTPUT_DIR)metadata_info.txt" output_folder="$(SC_RNASEQ_DOWNLOAD_LIB_DIR)"' 0Preparation/download_cleaning_data.R 0Preparation/download_cleaning_data.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@
	
## Prepare single cell info file
prepare_singlecell_info: download_cleaning_data $(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv 0Preparation/prepare_scrna_seq_sample_info.R
	@echo --- prepare single cell information file ---
	R CMD BATCH --no-save --no-restore '--args NEW_scRNASeqLibrary="$(OUTPUT_DIR)NEW_scRNASeqLibrary.tsv" cells_folder="$(SC_RNASEQ_DOWNLOAD_LIB_DIR)" output_folder="${PWD}/$(OUTPUT_DIR)"' 0Preparation/prepare_scrna_seq_sample_info.R 0Preparation/prepare_scrna_seq_sample_info.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@

################## STEPS RUN ON JURA ##################

# Start here on cluster
cluster1:

	@echo -e "\tBe sure everything is up-to-date before running RNASeq pipeline\n"
	@touch $@

# MUST be on a machine with read/write access to /data/ul/dee/bgee (e.g. rserv01 or dev), as bbgee user
# Ask for a cluster installation if a tool/library is missing
check_tools: cluster1
	@echo -e "\n\tFirst of all, go to '$(RNASEQ_CLUSTER_READONLY)GIT/' and do a 'git pull'"
	@echo -e "\tThen 'cd pipeline/RNA_Seq/' and be prepared to work\n"
	@echo -e "\n\tRun this command to give access to all modules installed on vital-it\n\tmodule add Bioinformatics/Software/vital-it\n"
	@echo -e "\n\tRun this command to prevent errors with Utils.pm\n\tmodule add Development/Ensembl_API/$(ENSRELEASE);\n"
	# Check if logged on cluster
	@if [[ `hostname -d` != 'axiom.unil.ch' ]]; then false; fi
	# Check if all required tools/libs are available
	@which perl                                        > $@.tmp
	@perl -Mdiagnostics           -e 1                >> $@.tmp
	@perl -MDBI                   -e 1                >> $@.tmp
	@perl -MDigest::SHA           -e 1                >> $@.tmp
	@perl -MFile::Basename        -e 1                >> $@.tmp
	@perl -MFile::Path            -e 1                >> $@.tmp
	@perl -MFile::Slurp           -e 1                >> $@.tmp
	@perl -MFile::Spec            -e 1                >> $@.tmp
	@perl -MFindBin               -e 1                >> $@.tmp
	@perl -MGetopt::Long          -e 1                >> $@.tmp
	@perl -MIO::Socket            -e 1                >> $@.tmp
	@perl -Mlib                   -e 1                >> $@.tmp
	@perl -MList::MoreUtils       -e 1                >> $@.tmp
	@perl -MLWP::Simple           -e 1                >> $@.tmp
	@perl -MSpreadsheet::Read     -e 1                >> $@.tmp
	@perl -Mstrict                -e 1                >> $@.tmp
	@perl -MTime::localtime       -e 1                >> $@.tmp
	@perl -Mwarnings              -e 1                >> $@.tmp
	@perl -MCpanel::JSON::XS      -e 1                >> $@.tmp
	@module add Bioinformatics/Software/vital-it
	@$(CLUSTER_R_CMD) which R                         >> $@.tmp
	@$(CLUSTER_R_CMD) R -e 'library("chipseq")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("edgeR")'         >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("GenomicFeatures")' >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("limma")'         >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("mclust")'        >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("rhdf5")'         >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("tximport")'      >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("rtracklayer")'   >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("biomaRt")'       >> $@.tmp  2>/dev/null
	@$(CLUSTER_R_CMD) R -e 'library("jsonlite")'      >> $@.tmp  2>/dev/null
	@which xz                                         >> $@.tmp
	@which sbatch                                     >> $@.tmp
	@$(CLUSTER_ASPERA_CMD)     which ascp             >> $@.tmp
	@$(CLUSTER_TOPHAT_CMD)     which gtf_to_fasta     >> $@.tmp
	@$(CLUSTER_SRATOOLKIT_CMD) which fastq-dump       >> $@.tmp
	@$(CLUSTER_KALLISTO_CMD)   which kallisto         >> $@.tmp
	@mkdir -p $(RNASEQ_CLUSTER_GTF)
	@if [[ -d $(RNASEQ_CLUSTER_GTF) ]]; then echo -n; else echo "Directory [$(RNASEQ_CLUSTER_GTF)] does not exist"  >> $@.tmp; false; fi
	@if [[ -w $(RNASEQ_CLUSTER_GTF) ]]; then echo -n; else echo "Directory [$(RNASEQ_CLUSTER_GTF)] is not writable" >> $@.tmp; false; fi
	@mkdir -p $(RNASEQ_CLUSTER_SCRATCH)
	@if [[ -d $(RNASEQ_CLUSTER_SCRATCH) ]]; then echo -n; else echo "Directory [$(RNASEQ_CLUSTER_SCRATCH)] does not exist"  >> $@.tmp; false; fi
	@if [[ -w $(RNASEQ_CLUSTER_SCRATCH) ]]; then echo -n; else echo "Directory [$(RNASEQ_CLUSTER_SCRATCH)] is not writable" >> $@.tmp; false; fi
	@$(MV) $@.tmp $@
















