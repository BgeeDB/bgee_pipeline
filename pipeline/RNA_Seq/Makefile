PIPELINEROOT := ../
DIR_NAME := RNA_Seq/
include $(PIPELINEROOT)Makefile.common
#TODO Add steps to recover file with sensitive infos (logins and passwords)

all: $(VERIFICATIONFILE)


################################ RETRIEVE BGEE ANNOTATIONS ################################

#FIXME git submodule should allow the get_annot step to work while protecting the annotation git from unallowed users!
get_annot:
	# Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations. Store them in the source_files directory
	@$(WGET) $(ANNOTATION_GIT_URL)/Strains/StrainMapping.tsv                && $(MV) StrainMapping.tsv                $(STRAIN_MAPPING_FILE)
	@$(WGET) $(ANNOTATION_GIT_URL)/RNA_Seq/RNASeqLibrary.tsv                && $(MV) RNASeqLibrary.tsv                $(SOURCE_FILES_DIR)/$(DIR_NAME)/RNASeqLibrary.tsv                2>/dev/null  || rm -f RNASeqLibrary.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/RNA_Seq/RNASeqExperiment.tsv             && $(MV) RNASeqExperiment.tsv             $(SOURCE_FILES_DIR)/$(DIR_NAME)/RNASeqExperiment.tsv             2>/dev/null  || rm -f RNASeqExperiment.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/RNA_Seq/RNASeqLibraryPlatformChecks.tsv  && $(MV) RNASeqLibraryPlatformChecks.tsv  $(SOURCE_FILES_DIR)/$(DIR_NAME)/RNASeqLibraryPlatformChecks.tsv  2>/dev/null  || rm -f RNASeqLibraryPlatformChecks.tsv
# Download wormbase annotation as well. Will be merged with our annotation in create_rna_seq_sample_info step
	@$(WGET) $(ANNOTATION_GIT_URL)/RNA_Seq/RNASeqLibrary_worm.tsv           && $(MV) RNASeqLibrary_worm.tsv           $(SOURCE_FILES_DIR)/$(DIR_NAME)/RNASeqLibrary_worm.tsv           2>/dev/null  || rm -f RNASeqLibrary_worm.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/RNA_Seq/RNASeqLibrary_worm_exclusion.tsv && $(MV) RNASeqLibrary_worm_exclusion.tsv $(SOURCE_FILES_DIR)/$(DIR_NAME)/RNASeqLibrary_worm_exclusion.tsv 2>/dev/null  || rm -f RNASeqLibrary_worm_exclusion.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/RNA_Seq/RNASeqExperiment_worm.tsv        && $(MV) RNASeqExperiment_worm.tsv        $(SOURCE_FILES_DIR)/$(DIR_NAME)/RNASeqExperiment_worm.tsv        2>/dev/null  || rm -f RNASeqExperiment_worm.tsv
#TODO clarify the steps from WormBase raw files to formatted file: $(ANNOTATION_GIT_URL)/RNA_Seq/RNASeqLibrary_worm.tsv
	@touch $@

check_annot: get_annot $(RNASEQ_EXPERIMENT_FILEPATH) $(RNASEQ_EXPERIMENT_FILEPATH_WORM) $(RNASEQ_LIB_FILEPATH) $(RNASEQ_LIBRARY_FILEPATH_WORM)
	# Concatenate our annotation file with wormbase annotations (without header, and sorted by experiments and libraries)
	@tail -n+2 $(RNASEQ_EXPERIMENT_FILEPATH_WORM) | sort -k1,1       | cat $(RNASEQ_EXPERIMENT_FILEPATH) - > $(RNASEQ_EXPERIMENT_FILEPATH_FULL)
	@tail -n+2 $(RNASEQ_LIB_FILEPATH_WORM)        | sort -k2,2 -k1,1 | cat $(RNASEQ_LIB_FILEPATH)        - > $(RNASEQ_LIB_FILEPATH_FULL).ori
	# Map strain names
	@perl -e 'use lib ".."; use Utils; Utils::map_strain_names("$(RNASEQ_LIB_FILEPATH_FULL).ori", "$(STRAIN_MAPPING_FILE)")' >$(RNASEQ_LIB_FILEPATH_FULL)
	@$(RM) $(RNASEQ_LIB_FILEPATH_FULL).ori
	# First check of annotations
	perl 0Before/check_rna_seq_curation.pl -bgee=$(BGEECMD) -RNAseqExperiment=$(RNASEQ_EXPERIMENT_FILEPATH_FULL) -RNAseqLib=$(RNASEQ_LIB_FILEPATH_FULL) -allRes=$(RNASEQALLRES) before >$@.tmp 2>&1
	@echo -e "Check file \"check_annot\" for the output of the script 0Before/check_rna_seq_curation.pl, which indicates potential errors to correct in the annotation files.\n"
	@$(MV) $@.tmp $@
#TODO the script needs better check of annotation files for leading and trailing spaces (only done for experimentId now, but needs to be done for other fields as well)

create_rna_seq_sample_info: check_annot $(EXTRAMAPPING_FILEPATH)
	# Generate rna_seq_sample_info.txt from RNASeq lib annotation file
	@perl 0Before/create_rna_seq_sample_info.pl -bgee=$(BGEECMD) -RNAseqLib=$(RNASEQ_LIB_FILEPATH_FULL) -RNAseqLibChecks=$(RNASEQ_LIB_CHECKS_FILEPATH) -RNAseqLibWormExclusion=$(RNASEQ_LIB_EXCLUSION_FILEPATH_WORM) -extraMapping=$(EXTRAMAPPING_FILEPATH) -outFile=$(RNASEQ_SAMPINFO_FILEPATH) >$@.tmp 2>&1
	@$(MV) $@.tmp $@
	@echo -e "Check file \"create_rna_seq_sample_info\" for the output of the script 0Before/create_rna_seq_sample_info.pl, which indicates potential errors to correct in the annotation files.\n"

commit_rna_seq_sample_info: create_rna_seq_sample_info $(RNASEQ_SAMPINFO_FILEPATH) $(STRAIN_MAPPING_FILE)
	# Commit the library information file thta will be used for the rest of the pipeline
	@$(GIT) add $(RNASEQ_SAMPINFO_FILEPATH) $(STRAIN_MAPPING_FILE) $(RNASEQ_EXPERIMENT_FILEPATH_FULL) $(RNASEQ_LIB_FILEPATH_FULL)
	@$(GIT) commit -m 'Update $(RNASEQ_SAMPINFO_FILEPATH) and $(STRAIN_MAPPING_FILE) for $(DBNAME)' $(RNASEQ_SAMPINFO_FILEPATH) || true
	@$(GIT) push
	@echo -e "\t$(RNASEQ_SAMPINFO_FILEPATH) is ready, you can go to regular cluster to download new SRA files with get_SRA.pl as *admin* user\n\tDo a 'git pull' before starting.\n"
	@touch $@


################################ SEND FILES TO CLUSTER ################################

send_files_to_cluster: commit_rna_seq_sample_info $(RNASEQ_LIB_FILEPATH_FULL) $(RNASEQ_EXPERIMENT_FILEPATH_FULL) $(RNASEQ_LIB_CHECKS_FILEPATH)
	# Send not-versionned annotation files to cluster  with login/password $(CLUSTERLOGIN)/$(CLUSTERPASSW)
	@scp get_annot \
		 check_annot \
		 create_rna_seq_sample_info \
		 $(RNASEQ_EXPERIMENT_FILEPATH_FULL) \
		 $(RNASEQ_LIB_FILEPATH_FULL) \
		 $(CLUSTERLOGIN)@$(CLUSTERHOST):$(RNASEQ_CLUSTER_READONLY)GIT/pipeline/RNA_Seq/
	@echo -e "\tNext steps have to be done on cluster\n\tGo to '$(RNASEQ_CLUSTER_READONLY)GIT/', do a 'git pull' and copy missing files in '$(RNASEQ_CLUSTER_READONLY)GIT/pipeline/RNA_Seq/'\n\tThen restart with 'make cluster1' to start the pipeline run step\n"
	@touch $@


################################ DOWNLOAD RNA-SEQ LIBRARIES ################################

##############################################################################################
# RNASeq library download
#NOTE to be done on regular cluster, with network access
get_sra: $(RNASEQ_SAMPINFO_FILEPATH) $(RNASEQ_ALREADY_DOWNLOADED)
	@$(GIT) pull
	@sed -i 's@\(BASE *= \).*@\1"$(RNASEQ_DOWNLOAD_LIB_DIR)";@'                      0Before/get_SRA.pl
	@sed -i 's@--output=.*@--output=${PWD}/get_sra.out@'                             0Before/download_lib.sbatch
	@sed -i 's@--partition=.*@--partition=${CLUSTER_PARTITION}@'                     0Before/download_lib.sbatch
	@sed -i 's@--account=.*@--account=${CLUSTER_ACCOUNT}@'                           0Before/download_lib.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/get_sra.err@'                               0Before/download_lib.sbatch
	@sed -i 's@SCRIPT_PATH=.*@SCRIPT_PATH=${PWD}@'                                   0Before/download_lib.sbatch
	@sed -i 's@ANNOTATION_FILE=.*@ANNOTATION_FILE=$(RNASEQ_SAMPINFO_FILEPATH)@'      0Before/download_lib.sbatch
	@sed -i 's@DONE_FILE=.*@DONE_FILE=$(RNASEQ_ALREADY_DOWNLOADED)@'                 0Before/download_lib.sbatch
	@sed -i 's@^.*module add UHTS/Analysis/sratoolkit/.*@$(CLUSTER_SRATOOLKIT_CMD)@' 0Before/download_lib.sbatch
	@sed -i 's@^.*module add Utility/aspera_connect/.*@#$(CLUSTER_ASPERA_CMD)@'      0Before/download_lib.sbatch
	@sed -i 's@^.*module add UHTS/Quality_control/fastp/.*@$(CLUSTER_FASTP_CMD)@'    0Before/download_lib.sbatch
	@sed -i 's@^.*module add R/.*@$(CLUSTER_R_CMD)@'                                 0Before/download_lib.sbatch
	@sbatch 0Before/download_lib.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo '!!! Rerun this step several times to complete downloads !!!'
	@touch $@

check_new_downloads: get_sra
	@echo 'If the following commands return something, check FASTQ files and/or rerun *get_sra* with more memory!'
	@find $(RNASEQ_DOWNLOAD_LIB_DIR)/FASTQ/RNAseq/ -type f -name \*.R.stat -exec wc -l {} \; | grep -v '^2 ' | cat
	@find $(RNASEQ_DOWNLOAD_LIB_DIR)/FASTQ/RNAseq/ -type f -name \*.R.stat | xargs grep -Hv '^#' | cut -f1 | grep ':0$$' | cat
	@touch $@

list_new_downloads: check_new_downloads
	@find $(RNASEQ_DOWNLOAD_LIB_DIR)/FASTQ/RNAseq/ -type f -name \*.fastq.gz\* | xargs -r dirname | sed -e 's@^.*/@@' | sort -u > /tmp/new_downloads
	@cat $(RNASEQ_ALREADY_DOWNLOADED) >>/tmp/new_downloads
	@sort -u /tmp/new_downloads >$(RNASEQ_ALREADY_DOWNLOADED)
	@rm -f /tmp/new_downloads
	@$(GIT) add $(RNASEQ_ALREADY_DOWNLOADED)
	@$(GIT) commit -m 'Add new downloaded libraries' $(RNASEQ_ALREADY_DOWNLOADED) || true
	@$(GIT) push
	@touch $@

#TODO Add a step to tranfer folders to sensitive cluster!
check_sra: list_new_downloads
	@echo -e "\tRe-run 0Before/get_SRA.pl to be sure all SRA are downloaded and FASTQ prepared\n"
	@scp  $<  $(CLUSTERLOGIN)@$(CLUSTERHOST):$(RNASEQ_CLUSTER_READONLY)GIT/pipeline/RNA_Seq/
	@touch $@


################################ DOWNLOAD GENOMES AND ANNOTATIONS ################################

#NOTE to be done on regular cluster, with network access
get_GTF: clean_cluster_folders create_rna_seq_sample_info
	# Get GTF files from Ensembl and NCBI FTP
	@perl 0Before/get_GTF_files.pl  -RNAseqSample=$(RNASEQ_SAMPINFO_FILEPATH) -ensRelease=$(ENSRELEASE) -ensMetazoaRelease=$(ENSMETAZOARELEASE) -outDir=$(RNASEQ_DOWNLOAD_GTF)  >$@.tmp 2>&1
	@$(MV) $@.tmp $@

# this rule is created to remove some elements from NCBI/RefSeq gtf files in order to make them compatible with our pipeline
update_GTF: get_GTF
	perl 0Before/update_GTF.pl -path_to_gtf_folder=$(RNASEQ_CLUSTER_GTF) -sample_info_file=$(RNASEQ_SAMPINFO_FILEPATH)>$@.tmp 2>&1
	@$(MV) $@.tmp $@

get_genome: update_GTF
	# Download genomes from Ensembl and NCBI databases
	@perl 0Before/get_genome_files.pl  -GTF_dir=$(RNASEQ_DOWNLOAD_GTF) -ensRelease=$(ENSRELEASE) -ensMetazoaRelease=$(ENSMETAZOARELEASE) -outDir=$(RNASEQ_DOWNLOAD_GTF) >$@.tmp 2>&1
	@echo "You can transfer GTF and genome files from UNIL cluster [$(RNASEQ_DOWNLOAD_GTF)] to sensitive cluster [$(RNASEQ_CLUSTER_GTF)]"
	@$(MV) $@.tmp $@


################################ PREPARE CLUSTER ENVIRONMENT ################################

# Start here on sensitive cluster
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
	@if [[ `hostname -d` != 'chuv.vital-it.ch' ]]; then false; fi
	# Check if all required tools/libs are available
	@module add Bioinformatics/Software/vital-it
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
	@$(CLUSTER_TOPHAT_CMD)     which gtf_to_fasta     >> $@.tmp
	@$(CLUSTER_SRATOOLKIT_CMD) which fastq-dump       >> $@.tmp
	@$(CLUSTER_FASTP_CMD)      which fastp            >> $@.tmp
	@$(CLUSTER_KALLISTO_CMD)   which kallisto         >> $@.tmp
	@mkdir -p $(RNASEQ_CLUSTER_GTF)
	@if [[ -d $(RNASEQ_CLUSTER_GTF) ]]; then echo -n; else echo "Directory [$(RNASEQ_CLUSTER_GTF)] does not exist"  >> $@.tmp; false; fi
	@if [[ -w $(RNASEQ_CLUSTER_GTF) ]]; then echo -n; else echo "Directory [$(RNASEQ_CLUSTER_GTF)] is not writable" >> $@.tmp; false; fi
	@mkdir -p $(RNASEQ_CLUSTER_SCRATCH)
	@if [[ -d $(RNASEQ_CLUSTER_SCRATCH) ]]; then echo -n; else echo "Directory [$(RNASEQ_CLUSTER_SCRATCH)] does not exist"  >> $@.tmp; false; fi
	@if [[ -w $(RNASEQ_CLUSTER_SCRATCH) ]]; then echo -n; else echo "Directory [$(RNASEQ_CLUSTER_SCRATCH)] is not writable" >> $@.tmp; false; fi
	@$(MV) $@.tmp $@

clean_cluster_folders: check_tools
	# Clean folders on cluster
	@echo "rm -Rf $(RNASEQ_CLUSTER_GTF)*.genome.* $(RNASEQ_CLUSTER_GTF)*.gtf.gz"
	@echo "xz -9 $(RNASEQ_CLUSTER_GTF)*.gtf_all $(RNASEQ_CLUSTER_GTF)*.transcriptome.* $(RNASEQ_CLUSTER_GTF)*.gene2transcript $(RNASEQ_CLUSTER_GTF)*.gene2biotype"
	@echo "rm -Rf $(RNASEQ_CLUSTER_SCRATCH)*"
	@echo
	@echo "Those scripts have to be run from $(RNASEQ_CLUSTER_SCRIPTS)$(RNASEQPATH) on the frontal (can write on /data/)!"
	@touch $@


################################ GENERATE INTERGENIC SEQUENCES FOR EACH SPECIES ################################
prepare_GTF: update_GTF
	# Prepare GTF files : $(CLUSTER_R_CMD)
	@perl 0Before/slurm_prepare_GTF.pl -gtf_dir=$(RNASEQ_CLUSTER_GTF) -block_size_N=31 -account=$(CLUSTER_ACCOUNT_SENSITIVE) -partition=$(CLUSTER_PARTITION_SENSITIVE) -proportion_N=0.05 -output_gtf_path=$(RNASEQ_CLUSTER_GTF) -output_log_folder=$(OUTPUT_DIR) -cluster_R_cmd="$(CLUSTER_R_CMD)" >$@.tmp 2>&1
	@echo rm -f $(RNASEQ_CLUSTER_GTF)/*.gtf.gz
	@$(MV) $@.tmp $@


prepare_indexed_transcriptome: prepare_GTF
	#Preparing indexed transcriptome for every species
	#Extract transcriptome.fa from gtf_all and genome.fa files
	#perl one liner to remove arbitrary numbering in fasta header from gtf_to_fasta
	# Prepare indexes for kallisto: one with default k-mer size $(RNASEQ_KALLISTO_KMER_DEFAULT), one with short k-mer size $(RNASEQ_KALLISTO_KMER_SHORT)
	@perl 1Run/slurm_index_creation.pl -transcriptome_folder=$(RNASEQ_CLUSTER_GTF) -output_log_folder=$(OUTPUT_DIR) -account=$(CLUSTER_ACCOUNT_SENSITIVE) -partition=$(CLUSTER_PARTITION_SENSITIVE) -short_index_length=$(RNASEQ_KALLISTO_KMER_SHORT)  -cluster_kallisto_cmd="$(CLUSTER_KALLISTO_CMD)" -cluster_tophat_cmd="$(CLUSTER_TOPHAT_CMD)" >$@.tmp 2>&1
	#TODO Use short k-mer size = 21nt instead of 15
	#TODO The gffread utility in cufflinks package seems more flexible and reliable
	#TODO add step to copy .passw file to bbgee's home on cluster / or maybe just echo a message to tell the user scp it? Path and name of this file on cluster is stored in ENCRYPT_PASSWD_FILE. We should store this file somewhere (devbioinfo?) because it cannot be added to the gitlab project
	#Or just put it in /home/bbgee/? in read only mode, only for bbgee user
	@$(MV) $@.tmp $@

#NOTE Better to run this step in screen or with nohup as submission is done X jobs per X jobs to not overload the system (and not decrease our user priority)!
abundance_all_intergenic: prepare_indexed_transcriptome $(RNASEQ_SAMPINFO_FILEPATH)
	# Running the pipeline
	@perl 1Run/slurm_scheduler.pl -sample_info_file=$(RNASEQ_SAMPINFO_FILEPATH) -exclude_sample_file=$(RNASEQ_SAMPEXCLUDED_FILEPATH) -output_log_folder=$(RNASEQ_CLUSTER_LOG) -account=$(CLUSTER_ACCOUNT_SENSITIVE) -partition=$(CLUSTER_PARTITION_SENSITIVE) -index_folder=$(RNASEQ_CLUSTER_GTF) -fastq_folder=$(RNASEQ_SENSITIVE_FASTQ) -kallisto_out_folder=$(RNASEQ_CLUSTER_ALL_RES) -enc_passwd_file=$(ENCRYPT_PASSWD_FILE) -cluster_kallisto_cmd='$(CLUSTER_KALLISTO_CMD)' -cluster_R_cmd='$(CLUSTER_R_CMD)' >$@.tmp 2>&1
	@echo "TODO: At the end it is a good idea to relaunch the abundance_all_intergenic step to be sure everything was run!"
	@$(MV) $@.tmp $@

check_abundance_all_intergenic: abundance_all_intergenic
	# Check SLURM stderr files
	@cat $(RNASEQ_CLUSTER_ALL_RES)/*/*.err > $@.tmp
	# Check kallisto runs without any mappings
	#@grep 'Kallisto results include only NAs'    $(RNASEQ_CLUSTER_ALL_RES)/*/*.Rout  >>$@.tmp
	# Check kallisto runs with low mappings, with few reads, ...
	# Check length of reads, ...
	@grep 'Problem:'                             $(RNASEQ_CLUSTER_ALL_RES)/*/*.txt   >>$@.tmp
	@grep 'Warning:'                             $(RNASEQ_CLUSTER_ALL_RES)/*/*.txt   >>$@.tmp
	@echo
	@echo "It is probably easier to relaunch the problematic samples manually, notably those requiring memory or runtime extreme limits"
	@echo
	@echo -e "TODO: Flagged / excluded samples with low % genes mapped / low number of reads mapped / for which mapping failed\nAdd them manually to file: $(RNASEQ_SAMPEXCLUDED_FILEPATH)"
	@$(MV) $@.tmp $@

#This rule was not run for Bgee 15.0
export_length: check_abundance_all_intergenic
	#Export transcript length to a file to export
	$(eval TIME_LENGTH := $(shell date +'%Y%m%d-%H%M%S'))
	@echo "#!/bin/sh" > $@_$(TIME_LENGTH).sh
	@echo "perl 3Insertion/export_feature_length.pl -library_info=$(RNASEQ_SAMPINFO_FILEPATH) -excluded_libraries=$(RNASEQ_SAMPEXCLUDED_FILEPATH) -tx2gene_dir=$(RNASEQ_CLUSTER_GTF) -all_results=$(RNASEQ_CLUSTER_ALL_RES) -length_info=$(RNASEQ_LENGTH_INFO_FILEPATH)" >> $@_$(TIME_LENGTH).sh
	@sbatch --account=$(CLUSTER_ACCOUNT_SENSITIVE) --time=2:00:00 --mem=2G --partition=$(CLUSTER_PARTITION_SENSITIVE) --error=$@_$(TIME_LENGTH).err --output=$@_$(TIME_LENGTH).out --nodes=1 --ntasks=1 --cpus-per-task=1 --job-name=$@_$(TIME_LENGTH) $@_$(TIME_LENGTH).sh > $@.tmp
	NUMBER_JOBS=$$(squeue --name=$@_$(TIME_LENGTH) --noheader | wc -l) ; \
	while [ $$NUMBER_JOBS -gt 0 ] ; do \
		sleep 30 ; \
		NUMBER_JOBS=$$(squeue --name=$@_$(TIME_LENGTH) --noheader | wc -l) ; \
	done ; \

	#@rm -r $@_$(TIME_LENGTH).*
	@echo -e "\tTranscript length are now exported into $(RNASEQ_LENGTH_INFO_FILEPATH) file\n"
	@$(MV) $@.tmp $@

finalize_pipeline: check_abundance_all_intergenic
	# Collect infos from .report files
	$(eval TIME_FINALIZE := $(shell date +'%Y%m%d-%H%M%S'))
	@echo "#!/bin/sh" > $@_$(TIME_FINALIZE).sh
	@echo "perl 3Insertion/create_rna_seq_report_info.pl -library_info=$(RNASEQ_SAMPINFO_FILEPATH) -excluded_libraries=$(RNASEQ_SAMPEXCLUDED_FILEPATH) -report_info=$(RNASEQ_CLUSTER_REPORTINFO) -all_results=$(RNASEQ_CLUSTER_ALL_RES)" >> $@_$(TIME_FINALIZE).sh
	@sbatch --account=$(CLUSTER_ACCOUNT_SENSITIVE) --time=2:00:00 --mem=2G --partition=$(CLUSTER_PARTITION_SENSITIVE) --error=$@_$(TIME_FINALIZE) --output=$@_$(TIME_FINALIZE).out --nodes=1 --ntasks=1 --cpus-per-task=1 --job-name=$@_$(TIME_FINALIZE) $@_$(TIME_FINALIZE).sh > $@.tmp
	# Wait end of the job
	NUMBER_JOBS=$$(squeue --name=$@_$(TIME_FINALIZE) --noheader | wc -l) ; \
	while [ $$NUMBER_JOBS -gt 0 ] ; do \
		sleep 30 ; \
		NUMBER_JOBS=$$(squeue --name=$@_$(TIME_FINALIZE) --noheader | wc -l) ; \
	done ; \

	# Touch all files so that they are not removed from /scratch/temporary
	@find $(RNASEQ_CLUSTER_SCRATCH) -exec touch {} \;
	# We use tail -n+1 instead of cat because it writes the name of the file in the concatenated file
	@tail -n+1 $(RNASEQ_CLUSTER_LOG)/*/*.err    > $(RNASEQ_CLUSTER_LOG)/all_std_err.txt 2>>warnings.$@
	@tail -n+1 $(RNASEQ_CLUSTER_LOG)/*/*.out    > $(RNASEQ_CLUSTER_LOG)/all_std_out.txt 2>>warnings.$@
	@tail -n+1 $(RNASEQ_CLUSTER_LOG)/*/*.report > $(RNASEQ_CLUSTER_LOG)/all_reports.txt 2>>warnings.$@
	# Back-up all data (not final but it's worth doing an intermediate backup here)
	@$(RM)    $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar*
	@tar -C $(RNASEQ_CLUSTER_SCRATCH)/ -cf $(RNASEQ_CLUSTER_SCRATCH)/all_results_$(DBNAME).tar  all_results_$(DBNAME)/  2>>warnings.$@
	@gzip -9  $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar                                                        2>>warnings.$@
	# Move back-up data to /data/
	@mkdir -p $(RNASEQ_CLUSTER_ALL_RES_BACKUP)/
	@$(MV) $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar.gz $(RNASEQ_CLUSTER_ALL_RES_BACKUP) 2>>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

# Infer blood libraries: provide info for each library (rule can be run in front)
infer_blood_samples: check_abundance_all_intergenic
	@echo --- start infering the blood libraries ---
	@$(CLUSTER_R_CMD) R CMD BATCH --no-save --no-restore '--args RNASeqLibrary="$(RNASEQ_LIB_FILEPATH)" globin_file="$(RNASEQ_GLOBIN_FILEPATH)" kallisto_count_folder="$(RNASEQ_CLUSTER_ALL_RES)" output="$(GENERATED_FILES_DIR)$(RNASEQPATH)"' 1Run/blood_protocols_inference.R 1Run/blood_protocols_inference.Rout > $@.tmp
	@echo --- DONE ---
	@$(MV) $@.tmp $@


#NOTE For iterative updates the  sum_by_species  step has to be skipped if you want to reuse the previous gaussian curves as they are!
#     So use the previous release $(RNASEQ_CLUSTER_GAUSSIAN_CHOICE) file and $(RNASEQ_CLUSTER_SUM_RES) result folder (Just update the database version if required)
#     for the  presence_absence  step!
sum_by_species: check_abundance_all_intergenic $(RNASEQ_SAMPINFO_FILEPATH) $(RNASEQ_SAMPEXCLUDED_FILEPATH)
	# Script using all data from each species to deconvolute the coding genes and intergenic regions underlying distributions
	@mkdir -p $(RNASEQ_CLUSTER_SUM_RES)

	$(eval TIME_SUM := $(shell date +'%Y%m%d-%H%M%S'))
	@echo "#!/bin/sh" > $@_$(TIME_SUM).sh
	@echo "$(CLUSTER_R_CMD) R CMD BATCH --vanilla --slave '--args rna_seq_sample_info=\"$(RNASEQ_SAMPINFO_FILEPATH)\" rna_seq_sample_excluded=\"$(RNASEQ_SAMPEXCLUDED_FILEPATH)\" kallisto_count_folder=\"$(RNASEQ_CLUSTER_ALL_RES)\" tx2gene_folder=\"$(RNASEQ_CLUSTER_GTF)\" sum_by_species_folder=\"$(RNASEQ_CLUSTER_SUM_RES)\"' 1Run/rna_seq_sum_by_species.R $(OUTPUT_DIR)rna_seq_sum_by_species.Rout" >> $@_$(TIME_SUM).sh
	@sbatch --account=$(CLUSTER_ACCOUNT_SENSITIVE) --time=2-00:00:00 --mem=30G --partition=$(CLUSTER_PARTITION_SENSITIVE) --output=$@_$(TIME_SUM).out --error=$@_$(TIME_SUM).err  --nodes=1 --ntasks=1 --cpus-per-task=1 --job-name=$@_$(TIME_SUM) $@_$(TIME_SUM).sh > $@.tmp
	NUMBER_JOBS=$$(squeue --name=$@_$(TIME_SUM) --noheader | wc -l) ; \
	echo $$NUMBER_JOBS ; \
	while [ $$NUMBER_JOBS -gt 0 ] ; do \
		sleep 30 ; \
		NUMBER_JOBS=$$(squeue --name=$@_$(TIME_SUM) --noheader | wc -l) ; \
	done ; \

	@echo -e "You HAVE TO create/update MANUALLY a file [$(RNASEQ_CLUSTER_GAUSSIAN_CHOICE)] with selected gaussians for coding and intergenic regions!\n"
	@$(MV) $@.tmp $@

## variables used only to generate intergenic sequences.
## should potentialy be moved at the top of the file or in the Makefile.common file
TRANSCRIPTOME_COMPRESSION_EXT	:= .xz
SUM_ABUNDANCE_FILE_PATH			:= $(RNASEQ_CLUSTER_SUM_RES)sum_abundance_gene_level+fpkm+intergenic+classification_SPECIES_ID.tsv
## Generate 2 intergenic fasta files for each species. One for reference intergenic sequences and one for other reference intergenic sequences
##
generate_intergenic_sequences: sum_by_species
	#generate intergenic dirs with all rights for user and group members
	@mkdir -p -m 0770 $(CLUSTER_REF_INTERGENIC_FOLDER)
	@mkdir -p -m 0770 $(CLUSTER_OTHER_INTERGENIC_FOLDER)

	$(eval TIME_INTERGENIC := $(shell date +'%Y%m%d-%H%M%S'))
	@echo "#!/bin/sh" > $@_$(TIME_INTERGENIC).sh
	@echo "perl 1Run/create_intergenic_fasta.pl -sample_info_path=$(RNASEQ_SAMPINFO_FILEPATH) -transcriptomes_folder=$(RNASEQ_CLUSTER_GTF) -transcriptome_compression_ext=$(TRANSCRIPTOME_COMPRESSION_EXT) -sum_abundance_file_path=$(SUM_ABUNDANCE_FILE_PATH) -gaussian_file_path=$(RNASEQ_CLUSTER_GAUSSIAN_CHOICE) -ref_intergenic_dir=$(CLUSTER_REF_INTERGENIC_FOLDER) -other_intergenic_dir=$(CLUSTER_OTHER_INTERGENIC_FOLDER) > $@.tmp 2> $@.warn" >> $@_$(TIME_INTERGENIC).sh
	@sbatch --account=$(CLUSTER_ACCOUNT_SENSITIVE) --time=2-00:00:00 --mem=5G --partition=$(CLUSTER_PARTITION_SENSITIVE) --error=$@_$(TIME_INTERGENIC).err --output=$@_$(TIME_INTERGENIC).out --nodes=1 --ntasks=1 --cpus-per-task=1 --job-name=$@_$(TIME_INTERGENIC) $@_$(TIME_INTERGENIC).sh > $@.tmp
	NUMBER_JOBS=$$(squeue --name=$@_$(TIME_INTERGENIC) --noheader | wc -l) ; \
	while [ $$NUMBER_JOBS -gt 0 ] ; do \
		sleep 30 ; \
		NUMBER_JOBS=$$(squeue --name=$@_$(TIME_INTERGENIC) --noheader | wc -l) ; \
	done ; \

	# compress intergenic fasta files
	@find $(CLUSTER_REF_INTERGENIC_FOLDER) -type f -name '*_intergenic.fa' -exec gzip --verbose --best {} \; >> $@.tmp 2> $@.warn
	@find $(CLUSTER_OTHER_INTERGENIC_FOLDER) -type f -name '*_intergenic.fa' -exec gzip --verbose --best {} \; >> $@.tmp 2> $@.warn
	@$(MV) $@.tmp $@

################################ GENERATE PRESENT/ABSENT CALLS WITH BgeeCall ################################

indexes_bgeecall: 
	# generate BgeeCall input file from the rna_seq_sample_info.txt file
	@perl 1Run/sample_info_to_bgeecall.pl -sample_info_file=$(RNASEQ_SAMPINFO_FILEPATH) -sample_excluded=$(RNASEQ_SAMPEXCLUDED_FILEPATH) -output_dir=$(RNASEQ_CLUSTER_BGEECALL_OUTPUT) -transcriptome_dir=$(RNASEQ_CLUSTER_GTF) -annotation_dir=$(RNASEQ_CLUSTER_GTF) -fastq_dir=$(RNASEQ_SENSITIVE_FASTQ) -bgeecall_file=$(RNASEQ_BGEECALL_FILE) -ref_intergenic_dir=$(CLUSTER_REF_INTERGENIC_FOLDER) >>$@.tmp 2>&1
	# generate kallisto indexes with BgeeCall
	@$(CLUSTER_R_CMD) R CMD BATCH --vanilla --slave '--args bgeecall_input_file="$(RNASEQ_BGEECALL_FILE)" account="$(CLUSTER_ACCOUNT_SENSITIVE)" time="24:00:00" partition="$(CLUSTER_PARTITION_SENSITIVE)" working_path="$(RNASEQ_CLUSTER_BGEECALL_OUTPUT)"' 1Run/bgeecall_index.R $(OUPUT_DIR)bgeecall_index.Rout
	@$(MV) $@.tmp $@

presence_absence_bgeecall: indexes_bgeecall
	# generate present/absent expression calls with BgeeCall
	@$(CLUSTER_R_CMD) R CMD BATCH --vanilla --slave '--args bgeecall_input_file="$(RNASEQ_BGEECALL_FILE)" account="$(CLUSTER_ACCOUNT_SENSITIVE)" time="2-00:00:00" partition="$(CLUSTER_PARTITION_SENSITIVE)" working_path="$(RNASEQ_CLUSTER_BGEECALL_OUTPUT)" decrypt_file_path="$(ENCRYPT_PASSWD_FILE)"' 1Run/bgeecall_calls.R $(OUTPUT_DIR)bgeecall_calls.Rout
	# TODO generate summary plot once all calls have been generated
	@touch $@

check_presence_absence: presence_absence_bgeecall 
	# check that presence absence calls have been generated for all libraries present in the rna_seq_sample_info file
	# generate file containing calls info for all libraries and use this file to generate different plots
	@$(CLUSTER_R_CMD) R CMD BATCH --vanilla --slave '--args bgeecall_sample_info="$(RNASEQ_BGEECALL_FILE)" calls_dir="$(RNASEQ_CLUSTER_BGEECALL_OUTPUT)all_results"' 1Run/rna_seq_calls_plot.R $(OUTPUT_DIR)rna_seq_calls_plot.Rout
	@touch $@

#NOTE if output directory is not $(RNASEQ_CLUSTER_ALL_RES): copy all files with calls to this directory, as well as summary stats and plot files
#for folder in *; do echo $folder; /bin/cp $folder/* ../all_results_bgee_v15/$folder/; done
# /bin/cp used because cp is an alias to cp -i
# Also copy presence_absence_all_samples.txt, presence_absence_all_samples.RDa and presence_absence_boxplots.pdf

save_and_send_results_back: finalize_pipeline check_presence_absence export_length
	# Touch all files so that they are not removed from $(RNASEQ_CLUSTER_SCRATCH)
	@find $(RNASEQ_CLUSTER_SCRATCH) -exec touch {} \; 2>warnings.$@
	# Back-up all data
	@$(RM) $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar* $(RNASEQ_CLUSTER_SCRATCH)presence_absence_$(DBNAME).tar* $(RNASEQ_CLUSTER_SCRATCH)sum_by_species_$(MAJOR_RELEASE).tar*
	@tar -C $(RNASEQ_CLUSTER_SCRATCH) -cf $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar            all_results_$(DBNAME)/            2>>warnings.$@
	@tar -C $(RNASEQ_CLUSTER_SCRATCH) -cf $(RNASEQ_CLUSTER_SCRATCH)presence_absence_$(DBNAME).tar       presence_absence_$(DBNAME)/       2>>warnings.$@
	@tar -C $(RNASEQ_CLUSTER_SCRATCH) -cf $(RNASEQ_CLUSTER_SCRATCH)sum_by_species_$(MAJOR_RELEASE).tar  sum_by_species_$(MAJOR_RELEASE)/  2>>warnings.$@
	@gzip -9 $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar            2>>warnings.$@
	@gzip -9 $(RNASEQ_CLUSTER_SCRATCH)presence_absence_$(DBNAME).tar       2>>warnings.$@
	@gzip -9 $(RNASEQ_CLUSTER_SCRATCH)sum_by_species_$(MAJOR_RELEASE).tar  2>>warnings.$@
	# Move back-up data to $(RNASEQ_CLUSTER_ALL_RES_BACKUP)
	@$(MV) $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar.gz $(RNASEQ_CLUSTER_SCRATCH)presence_absence_$(DBNAME).tar.gz $(RNASEQ_CLUSTER_SCRATCH)sum_by_species_$(MAJOR_RELEASE).tar.gz $(RNASEQ_CLUSTER_ALL_RES_BACKUP) 2>>warnings.$@
	# Whole archive is probably too big to be copied to our servers (>100Gb). Commit only the final gene-level expression + calls files:
	@$(CP) $(RNASEQ_CLUSTER_REPORTINFO)                                    $(RNASEQREPORTINFO) 2>>warnings.$@
	@$(CP) $(RNASEQ_CLUSTER_PRESENCE_RES)/presence_absence_all_samples.txt $(RNASEQSAMPSTATS)  2>>warnings.$@
	#To do to complete
	@echo -e "\tTODO: Commit/Push $(RNASEQREPORTINFO) and $(RNASEQSAMPSTATS) files\n"
	@echo -e "\tTODO: Save tarballs to $(CLUSTER_ARCHIVE_PATH)/rna_seq/all_results_$(DBNAME)/\n"
	@echo -e "\tTODO: Copy tarballs to development server ($(PIPEHOST)) for db insertion in [$(RNASEQALLRES)]"
	@echo -e "\t      For insertion only  '*abundance_gene_level+new_tpm+new_fpkm+calls.tsv' $(RNASEQ_CLUSTER_REPORTINFO) $(RNASEQ_CLUSTER_PRESENCE_RES)/presence_absence_all_samples.txt $(RNASEQ_CLUSTER_PRESENCE_RES)/presence_absence_all_samples.RDa $(RNASEQ_CLUSTER_PRESENCE_RES)/presence_absence_boxplots.pdf  look to be required\n"
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@echo -e "\tYou can go out of the sensitive cluster now\n"
	@touch $@


################################ INSERT CALLS IN DATABASE ################################

#NOTE For iterative updates die commands after insert/update/delete statement failure have to be changed to warn
#     because with iterative updates the database may contain duplicates in regard to what you are inserting/updating/deleting!
insert_RNA_seq: $(RNASEQALLRES) $(RNASEQSAMPSTATS) $(RNASEQREPORTINFO) $(RNASEQ_LIB_FILEPATH_FULL) $(RNASEQ_EXPERIMENT_FILEPATH_FULL) $(UBERON_SEX_INFO_FILE_PATH) $(CUSTOM_UBERON_FILE_PATH) $(DEV_STAGE_ONT_FILE_PATH) $(RNASEQ_SAMPINFO_FILEPATH) $(RNASEQ_SAMPEXCLUDED_FILEPATH) $(EXTRAMAPPING_FILEPATH)
	# Launch the organ stage mapping tool (using $(CUSTOM_UBERON_FILE_PATH) and $(DEV_STAGE_ONT_FILE_PATH))
	@$(IDMAPPING)  $(IDMAPPINGPORT) &
	@$(STGMAPPING) $(STGMAPPINGPORT) &
	@sleep 50 # sleep because mappers need time to load Uberon
	# Insert RNA Seq data
	@perl 3Insertion/insert_rna_seq.pl -bgee=$(BGEECMD) -rnaSeqLibrary=$(RNASEQ_LIB_FILEPATH_FULL) -rnaSeqExperiment=$(RNASEQ_EXPERIMENT_FILEPATH_FULL) -library_info=$(RNASEQ_SAMPINFO_FILEPATH) -excluded_libraries=$(RNASEQ_SAMPEXCLUDED_FILEPATH) -library_stats=$(RNASEQSAMPSTATS) -report_info=$(RNASEQREPORTINFO) -all_results=$(RNASEQALLRES) -sex_info=$(UBERON_SEX_INFO_FILE_PATH) -extraMapping=$(EXTRAMAPPING_FILEPATH) -Aport=$(IDMAPPINGPORT) -Sport=$(STGMAPPINGPORT) > $@.tmp 2>warnings.$@
	@echo "Check inconsistencies between condition species and gene species (there should be none): " >> $@.tmp
	@$(MYSQL) -e "SELECT t1.* FROM rnaSeqResult AS t1 INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId INNER JOIN rnaSeqLibrary AS t3 ON t1.rnaSeqLibraryId = t3.rnaSeqLibraryId INNER JOIN cond AS t4 on t3.conditionId = t4.conditionId WHERE t2.speciesId != t4.speciesId" >> $@.tmp
	@echo >> $@.tmp
	@echo "Distinct strains in RNA-Seq conditions, check that they are correct (e.g., no 'wild type' instead of 'wild-type')" >> $@.tmp
	@$(MYSQL) -e "SELECT DISTINCT t1.strain, t1.speciesId FROM cond AS t1 INNER JOIN rnaSeqLibrary AS t2 ON t1.conditionId = t2.conditionId ORDER BY t1.strain" >> $@.tmp
	@echo >> $@.tmp
	@echo "Statistics on libraries" >> $@.tmp
	@$(MYSQL) -e "SELECT t1.speciesId, COUNT(DISTINCT t2.rnaSeqLibraryId), MIN(t2.tmmFactor), MAX(t2.tmmFactor), AVG(t2.tmmFactor), MIN(t2.fpkmThreshold), MAX(t2.fpkmThreshold), AVG(t2.fpkmThreshold), MIN(t2.tpmThreshold), MAX(t2.tpmThreshold), AVG(t2.tpmThreshold), MIN(t2.allGenesPercentPresent), MAX(t2.allGenesPercentPresent), AVG(t2.allGenesPercentPresent), MIN(t2.proteinCodingGenesPercentPresent), MAX(t2.proteinCodingGenesPercentPresent), AVG(t2.proteinCodingGenesPercentPresent), MIN(t2.mappedReadsCount), MAX(t2.mappedReadsCount), AVG(t2.mappedReadsCount) FROM cond AS t1 INNER JOIN rnaSeqLibrary AS t2 ON t1.conditionId = t2.conditionId GROUP BY t1.speciesId ORDER BY COUNT(DISTINCT t2.rnaSeqLibraryId) DESC" >> $@.tmp
	@echo >> $@.tmp
	@echo "Statistics on conditions" >> $@.tmp
	@$(MYSQL) -e "SELECT t1.speciesId, COUNT(DISTINCT t2.rnaSeqLibraryId), COUNT(DISTINCT t2.conditionId), COUNT(DISTINCT t1.exprMappedConditionId), COUNT(DISTINCT t1.anatEntityId), COUNT(DISTINCT t1.stageId), COUNT(DISTINCT t1.anatEntityId, t1.stageId, t1.sex), GROUP_CONCAT(DISTINCT t1.sex ORDER BY t1.sex SEPARATOR ', '), GROUP_CONCAT(DISTINCT t1.strain ORDER BY t1.strain SEPARATOR ', ') FROM cond AS t1 INNER JOIN rnaSeqLibrary AS t2 ON t1.conditionId = t2.conditionId GROUP BY t1.speciesId ORDER BY COUNT(DISTINCT t2.rnaSeqLibraryId) DESC" >> $@.tmp
	@echo >> $@.tmp
	@echo "Same condition information, but for mapped conditions of expression tables" >> $@.tmp
	@$(MYSQL) -e "SELECT t3.speciesId, COUNT(DISTINCT t2.rnaSeqLibraryId), COUNT(DISTINCT t3.anatEntityId), COUNT(DISTINCT t3.stageId), COUNT(DISTINCT t3.anatEntityId, t3.stageId, t3.sex), GROUP_CONCAT(DISTINCT t3.sex ORDER BY t3.sex SEPARATOR ', '), GROUP_CONCAT(DISTINCT t3.strain ORDER BY t3.strain SEPARATOR ', ') FROM cond AS t1 INNER JOIN rnaSeqLibrary AS t2 ON t1.conditionId = t2.conditionId INNER JOIN cond AS t3 ON t1.exprMappedConditionId = t3.conditionId GROUP BY t3.speciesId ORDER BY COUNT(DISTINCT t2.rnaSeqLibraryId) DESC" >> $@.tmp
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@
#TODO
#keep files in all_results (rsync)
#rsync all rna_seq folder as done for Affymetrix

check_conditions: insert_RNA_seq
	@echo "Conditions with anat. entity not existing in related species:" > $@.tmp
	@$(MYSQL) -e "SELECT DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain FROM cond AS t1 LEFT OUTER JOIN anatEntityTaxonConstraint AS t2 ON t1.anatEntityId = t2.anatEntityId AND (t2.speciesId IS NULL OR t1.speciesId = t2.speciesId) LEFT OUTER JOIN anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId LEFT OUTER JOIN rnaSeqLibrary AS t10 ON t1.conditionId = t10.conditionId WHERE t2.anatEntityId IS NULL AND t10.conditionId IS NOT NULL ORDER BY t1.speciesId" >> $@.tmp 2> $@.warn
	@echo >> $@.tmp
	@echo "Conditions with dev. stage not existing in related species:" >> $@.tmp
	@$(MYSQL) -e "SELECT DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain FROM cond AS t1 LEFT OUTER JOIN stageTaxonConstraint AS t2 ON t1.stageId = t2.stageId AND (t2.speciesId IS NULL OR t1.speciesId = t2.speciesId) LEFT OUTER JOIN anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId LEFT OUTER JOIN rnaSeqLibrary AS t10 ON t1.conditionId = t10.conditionId WHERE t2.stageId IS NULL AND t10.conditionId IS NOT NULL ORDER BY t1.speciesId" >> $@.tmp 2>> $@.warn
	@$(MV) $@.tmp $@

TMP_TABLE_CLAUSE   := SELECT t10.expressionId, IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'present' AND t11.rnaSeqExperimentCallQuality = 'high quality') >= 2, 'pstGold', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'present' AND t11.rnaSeqExperimentCallQuality = 'high quality') = 1, 'pstSilverFromHigh', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'present' AND t11.rnaSeqExperimentCallQuality = 'poor quality') >= 2, 'pstSilverFromLow', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'present' AND t11.rnaSeqExperimentCallQuality = 'poor quality') = 1, 'pstBronze', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'absent' AND t11.rnaSeqExperimentCallQuality = 'high quality') >= 2, 'absGold', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'absent' AND t11.rnaSeqExperimentCallQuality = 'high quality') = 1, 'absSilverFromHigh', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'absent' AND t11.rnaSeqExperimentCallQuality = 'poor quality') >= 2, 'absSilverFromLow', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'absent' AND t11.rnaSeqExperimentCallQuality = 'poor quality') = 1, 'absBronze', 'IMPOSSIBLE')))))))) AS expressionCall FROM rnaSeqExperimentExpression AS t10

insert_expression: check_conditions insert_RNA_seq
	# Insert the expression summaries
	@perl 3Insertion/insert_rna_seq_expression.pl -bgee=$(BGEECMD) > $@.tmp 2>warnings.$@
	@echo
	@echo "Searching for incorrect updates of rnaSeqResult/expression tables (there should be none)" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM rnaSeqResult WHERE reasonForExclusion = 'not excluded' AND expressionId IS NULL" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM rnaSeqResult WHERE reasonForExclusion != 'not excluded' AND expressionId IS NOT NULL" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM rnaSeqResult AS t1 WHERE expressionId IS NOT NULL AND NOT EXISTS(SELECT 1 FROM rnaSeqExperimentExpression AS t2 WHERE t2.expressionId = t1.expressionId)" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM rnaSeqExperimentExpression AS t1 WHERE NOT EXISTS (SELECT 1 FROM rnaSeqResult AS t2 WHERE t2.expressionId = t1.expressionId)" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM rnaSeqExperimentExpression WHERE (presentHighRNASeqLibraryCount + presentLowRNASeqLibraryCount + absentHighRNASeqLibraryCount + absentLowRNASeqLibraryCount) = 0" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM rnaSeqExperimentExpression WHERE (presentHighRNASeqLibraryCount > 0 AND (rnaSeqExperimentCallDirection != 'present' OR rnaSeqExperimentCallQuality != 'high quality')) OR (presentHighRNASeqLibraryCount = 0 AND presentLowRNASeqLibraryCount > 0 AND (rnaSeqExperimentCallDirection != 'present' OR rnaSeqExperimentCallQuality != 'poor quality')) OR (presentHighRNASeqLibraryCount = 0 AND presentLowRNASeqLibraryCount = 0 AND absentHighRNASeqLibraryCount > 0 AND (rnaSeqExperimentCallDirection != 'absent' OR rnaSeqExperimentCallQuality != 'high quality')) OR (presentHighRNASeqLibraryCount = 0 AND presentLowRNASeqLibraryCount = 0 AND absentHighRNASeqLibraryCount = 0 AND (rnaSeqExperimentCallDirection != 'absent' OR rnaSeqExperimentCallQuality != 'poor quality'))" >> $@.tmp
	@$(MYSQL) -e "SELECT t1.* FROM rnaSeqExperiment AS t1 WHERE NOT EXISTS (SELECT 1 FROM rnaSeqExperimentExpression AS t2 WHERE t2.rnaSeqExperimentId = t1.rnaSeqExperimentId)" >> $@.tmp
	@echo >> $@.tmp
	@echo "Statistics for rnaSeqResult" >> $@.tmp
	@$(MYSQL) -e "SELECT speciesId, expressionCall, detectionFlag, rnaSeqData, COUNT(*) AS rnaSeqResultCount FROM rnaSeqResult AS t1 INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId INNER JOIN ($(TMP_TABLE_CLAUSE)) AS tmpTable ON t1.expressionId = tmpTable.expressionId GROUP BY speciesId, expressionCall, detectionFlag, rnaSeqData ORDER BY speciesId, expressionCall, detectionFlag, rnaSeqData" >> $@.tmp
	@echo "Statistics for expression table" >> $@.tmp
	@$(MYSQL) -e "SELECT speciesId, expressionCall, COUNT(*), COUNT(DISTINCT t1.bgeeGeneId), COUNT(DISTINCT t1.conditionId), COUNT(DISTINCT t3.anatEntityId, t3.stageId), COUNT(DISTINCT t3.anatEntityId, t3.stageId, t3.sex) FROM expression AS t1 INNER JOIN cond AS t3 ON t1.conditionId = t3.conditionId INNER JOIN ($(TMP_TABLE_CLAUSE)) AS tmpTable ON t1.expressionId = tmpTable.expressionId GROUP BY speciesId, expressionCall ORDER BY speciesId, expressionCall" >> $@.tmp
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

#NOTE For iterative updates die commands after insert/update/delete statement failure have to be changed to warn
#     because with iterative updates the database may contain duplicates in regard to what you are inserting/updating/deleting!
insert_feature_length: $(RNASEQ_LENGTH_INFO_FILEPATH)
	# Insert the feature length information
	@perl 3Insertion/insert_feature_length.pl -bgee=$(BGEECMD) -length_info=$(RNASEQ_LENGTH_INFO_FILEPATH) > $@.tmp 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

#NOTE For iterative updates die commands to warn
#     because with iterative updates tries to calculate TMM for already processed RNA-Seq libraries not available here!
launch_calculate_TMM_factors: insert_RNA_seq $(RNASEQALLRES)
	# Launch calculation of TMM factors for RNA-seq
	@mkdir -p $(RNASEQTMMTARG) $(RNASEQTMMPATH)
	@perl 3Insertion/launch_calculate_TMM_factors.pl -bgee=$(BGEECMD) -path_generes=$(RNASEQALLRES) -path_target=$(RNASEQTMMTARG) -path_processed=$(RNASEQTMMPATH) >$@.tmp 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

check_TMM_factors: launch_calculate_TMM_factors $(RNASEQTMMPATH)
	# Check results: all were calculated (none :0$ expected)
	@grep -H -c 'proc.time' $(RNASEQTMMPATH)/*.log | grep ':0$$'  >$@.tmp   || true
	# Check warnings: no warnings
	@grep 'Warning' -A 5 $(RNASEQTMMPATH)/*.log                   >>$@.tmp  || true
	@$(MV) $@.tmp $@
#TODO check that all libraries in database have a TMM factor calculated:
# cat $(RNASEQTMMPATH)/*.tsv | grep -v rnaSeqExperimentId | wc -l
# This should correspond to: SELECT count(*) FROM rnaSeqLibrary;

insert_TMM_factors: check_TMM_factors $(RNASEQTMMPATH)
	# Insert TMM factors into rnaSeqLibrary table
	@perl 3Insertion/insert_TMM_factors.pl -bgee=$(BGEECMD) -tmm_results=$(RNASEQTMMPATH) >$@.tmp 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@


#TODO Keep a way to use an external annotation file
#FIXME Useful? As Wormbase annotations are merged within RNASeq*_full.tsv files used in the main part!
wormbase:
	# Get Annotation file from WormBase
	@$(WGET) 'http://athena.caltech.edu/MrExpTable.csv' && $(MV) MrExpTable.csv 4External/$@.tsv 2>/dev/null  || rm -f MrExpTable.csv
	# Filter useful experiments, species/strains, conditions, ...
	@perl 4External/extract_rnaseq_info.pl 4External/$@.tsv > 4External/$@.map 2> warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	# Make them compatible with the pipeline
	@echo TODO

external: wormbase


#TODO Add 1st steps
#TODO add a verification that each transcriptId is uniquely used in a given species
#(see comment for key 'unique(transcriptId, bgeeGeneId)' in bgeeConstraint.sql)
$(VERIFICATIONFILE): send_files_to_cluster check_conditions insert_expression
	# Copy run_info_file.txt for next release comparison (check_runs step)
	@$(CP) $(RNASEQ_RUNINFO_FILEPATH) $(RNASEQ_RUNINFO_FILEPATH)_v$(RELEASE)
#TODO this needs to be updated

	# Check RNA Seq data insertion & expression
	# TODO: to update for Bgee 14. See existing queries in specific targets of this Makefile
#	@echo 'Get number of RNA-Seq libraries and conditions studied per species' > $@.tmp
#	@$(MYSQL) -e "SELECT (SELECT t3.speciesId FROM gene AS t3 INNER JOIN rnaSeqResult AS t2 ON t3.geneId = t2.geneId WHERE t2.rnaSeqLibraryId = t1.rnaSeqLibraryId LIMIT 1) AS speciesId, (SELECT t4.speciesCommonName FROM species AS t4 INNER JOIN gene AS t3 ON t3.speciesId = t4.speciesId INNER JOIN rnaSeqResult AS t2 ON t3.geneId = t2.geneId WHERE t2.rnaSeqLibraryId = t1.rnaSeqLibraryId limit 1) AS speciesName, COUNT(DISTINCT t1.rnaSeqLibraryId) AS libraryCount, COUNT(DISTINCT t1.anatEntityId, t1.stageId) AS numberOfConditions, COUNT(DISTINCT t1.anatEntityId) AS organCount, COUNT(DISTINCT t1.stageId) AS stageCount FROM rnaSeqLibrary AS t1 GROUP BY speciesId ORDER BY speciesId;" >> $@.tmp
#	@echo 'Get info for the table "rnaSeqResult"' >> $@.tmp
#	@$(MYSQL) -e "SELECT t2.speciesId, t3.speciesCommonName, detectionFlag, rnaSeqData, COUNT(*) AS resultCount FROM rnaSeqResult AS t1 INNER JOIN gene AS t2 ON t1.geneId = t2.geneId INNER JOIN species AS t3 ON t2.speciesId = t3.speciesId GROUP BY t2.speciesId, detectionFlag, rnaSeqData;" >> $@.tmp
#	@echo 'Check expression insertion' >> $@.tmp
#	@$(MYSQL) -e "SELECT t2.speciesId, t3.speciesCommonName, t1.rnaSeqData,             COUNT(DISTINCT t1.expressionId),   COUNT(DISTINCT t1.geneId) FROM expression   AS t1 INNER JOIN gene AS t2 ON t1.geneId = t2.geneId INNER JOIN species AS t3 ON t2.speciesId = t3.speciesId WHERE rnaSeqData             != 'no data' GROUP BY t2.speciesId, rnaSeqData;" >> $@.tmp
	@$(MV) $@.tmp $@
# add the verification file to git
	-@$(GIT) add $@
	-@$(GIT) commit $@ -m "Commit verification file for RNA Seq" || true
	-@$(GIT) push

#TODO Clean added/removed steps
#TODO Remove unused scripts/files in folders and sub-folders
clean:
	-@$(RM) -R $(VERIFICATIONFILE)  check_annot check_tools create_rna_seq_sample_info get_annot get_GTF prepare_GTF get_genome prepare_indexed_transcriptome generate_genome_info prepare_run_info_final create_config_file cluster get_sra check_new_downloads list_new_downloads check_sra check_curation check_runs insert_RNA_seq insert_expression wormbase external $(RNASEQGTFDATAPATH)/*.gtf $(RNASEQGTFDATAPATH)/*.gtf.gz *.tmp check_conditions insert_expression
	-@$(GIT) rm $(VERIFICATIONFILE)
	-@$(GIT) commit $(VERIFICATIONFILE) -m "Removing verification file for RNA Seq" || true
	-@$(GIT) push

# special target to remove RNA_seq from database. We make the name long on purpose, and does not
# incude it in the clean target, to not wipe the database accidentally
.PHONY := $(.PHONY) deleteRNASeq
deleteRNASeq:
	@$(MYSQL) -e "DELETE FROM rnaSeqResult; DELETE FROM rnaSeqRun; DELETE FROM rnaSeqLibrary; DELETE FROM rnaSeqExperimentToKeyword; DELETE FROM rnaSeqExperiment; DELETE FROM rnaSeqPlatform;"
#TODO check if new tables were added
