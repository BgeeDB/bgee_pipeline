PIPELINEROOT := ../
DIR_NAME := RNA_Seq/
include $(PIPELINEROOT)Makefile.common
#TODO Add steps to recover file with sensitive infos (logins and passwords)

all: $(VERIFICATIONFILE)


################################ RETRIEVE BGEE ANNOTATIONS ################################

#FIXME Fix strain names to be more consistent with Anne's strain file!
get_annot:
	# Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations. Store them in the source_files directory
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
	@tail -n+2 $(RNASEQ_LIB_FILEPATH_WORM)        | sort -k2,2 -k1,1 | cat $(RNASEQ_LIB_FILEPATH)        - > $(RNASEQ_LIB_FILEPATH_FULL)
	# First check of annotations
	perl 0Before/check_rna_seq_curation.pl -bgee=$(BGEECMD) -RNAseqExperiment=$(RNASEQ_EXPERIMENT_FILEPATH_FULL) -RNAseqLib=$(RNASEQ_LIB_FILEPATH_FULL) -allRes=$(RNASEQALLRES) before >$@.tmp 2>&1
	@echo "Check file \"check_annot\" for the output of the script 0Before/check_rna_seq_curation.pl, which indicates potential errors to correct in the annotation files.\n"
	@$(MV) $@.tmp $@
#TODO the script needs better check of annotation files for leading and trailing spaces (only done for experimentId now, but needs to be done for other fields as well)

create_rna_seq_sample_info: check_annot $(EXTRAMAPPING_FILEPATH)
	# Generate rna_seq_sample_info.txt from RNASeq lib annotation file
	@perl 0Before/create_rna_seq_sample_info.pl -bgee=$(BGEECMD) -RNAseqLib=$(RNASEQ_LIB_FILEPATH_FULL) -RNAseqLibChecks=$(RNASEQ_LIB_CHECKS_FILEPATH) -RNAseqLibWormExclusion=$(RNASEQ_LIB_EXCLUSION_FILEPATH_WORM) -extraMapping=$(EXTRAMAPPING_FILEPATH) -outFile=$(RNASEQ_SAMPINFO_FILEPATH) >$@.tmp 2>&1
	@$(MV) $@.tmp $@
	@echo "Check file \"create_rna_seq_sample_info\" for the output of the script 0Before/create_rna_seq_sample_info.pl, which indicates potential errors to correct in the annotation files.\n"

commit_rna_seq_sample_info: create_rna_seq_sample_info $(RNASEQ_SAMPINFO_FILEPATH)
	# Commit the library information file thta will be used for the rest of the pipeline
	@$(GIT) add $(RNASEQ_SAMPINFO_FILEPATH)
	@$(GIT) commit -m 'Update $(RNASEQ_SAMPINFO_FILEPATH) for $(DBNAME)' $(RNASEQ_SAMPINFO_FILEPATH)
	@$(GIT) push
	@echo -e "\t$(RNASEQ_SAMPINFO_FILEPATH) is ready, you can go to regular cluster to download new SRA files with get_SRA.pl as *admin* user\n\tDo a 'git pull' before starting\n"
	@$(MV) $@.tmp $@


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
	@sed -i 's@\(BASE *= \).*@\1"$(RNASEQ_DOWNLOAD_LIB_DIR)";@'                 0Before/get_SRA.pl
	@sed -i 's@--output=.*@--output=${PWD}/get_sra.out@'                        0Before/download_lib.sbatch
	@sed -i 's@--error=.*@--error=${PWD}/get_sra.err@'                          0Before/download_lib.sbatch
	@sed -i 's@SCRIPT_PATH=.*@SCRIPT_PATH=${PWD}@'                              0Before/download_lib.sbatch
	@sed -i 's@ANNOTATION_FILE=.*@ANNOTATION_FILE=$(RNASEQ_SAMPINFO_FILEPATH)@' 0Before/download_lib.sbatch
	@sed -i 's@DONE_FILE=.*@DONE_FILE=$(RNASEQ_ALREADY_DOWNLOADED)@'            0Before/download_lib.sbatch
	@sbatch 0Before/download_lib.sbatch
	@echo 'Check with  squeue/sacct -j <JOB_ID>  the job status'
	@echo '!!! Rerun this step several times to complete downloads !!!'
	@touch $@

list_new_downloads: get_sra
	@find $(RNASEQ_DOWNLOAD_LIB_DIR) -type f -name \*.fastq.gz\* | xargs dirname | sed -e 's@^.*/@@' | sort -u > /tmp/new_downloads
	@cat $(RNASEQ_ALREADY_DOWNLOADED) >>/tmp/new_downloads
	@sort -u /tmp/new_downloads >$(RNASEQ_ALREADY_DOWNLOADED)
	@rm -f /tmp/new_downloads
	@$(GIT) add $(RNASEQ_ALREADY_DOWNLOADED)
	@$(GIT) commit -m 'Add new downloaded libraries' $(RNASEQ_ALREADY_DOWNLOADED)
	@$(GIT) push
	@touch $@

#TODO Add a step to tranfer folders to sensitive cluster!
check_sra: list_new_downloads
	@echo -e "\tRe-run 0Before/get_SRA.pl to be sure all SRA are downloaded and FASTQ prepared\n"
	@scp  $<  $(CLUSTERLOGIN)@$(CLUSTERHOST):$(RNASEQ_CLUSTER_READONLY)GIT/pipeline/RNA_Seq/
	@touch $@


################################ PREPARE CLUSTER ENVIRONMENT ################################

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

clean_cluster_folders: check_tools
	# Clean folders on cluster
	@echo "rm -Rf $(RNASEQ_CLUSTER_GTF)*.genome.* $(RNASEQ_CLUSTER_GTF)*.gtf.gz"
	@echo "xz -9 $(RNASEQ_CLUSTER_GTF)*.gtf_all $(RNASEQ_CLUSTER_GTF)*.transcriptome.* $(RNASEQ_CLUSTER_GTF)*.gene2transcript $(RNASEQ_CLUSTER_GTF)*.gene2biotype"
	@echo "rm -Rf $(RNASEQ_CLUSTER_SCRATCH)*"
	@echo
	@echo "Those scripts have to be run from $(RNASEQ_CLUSTER_SCRIPTS)$(RNASEQPATH) on dev/frt/rserv01 (can write on /data/)!"
	@touch $@


################################ DOWNLOAD GENOMES AND ANNOTATIONS ################################

get_GTF: clean_cluster_folders create_rna_seq_sample_info
	# Get GTF files from Ensembl FTP
	@perl 0Before/get_GTF_files.pl  -RNAseqSample=$(RNASEQ_SAMPINFO_FILEPATH) -ensRelease=$(ENSRELEASE) -ensMetazoaRelease=$(ENSMETAZOARELEASE) -outDir=$(RNASEQ_CLUSTER_GTF)  >$@.tmp 2>&1
	@$(MV) $@.tmp $@
#TODO implement use of new "database" field in rna_seq_sample_info.txt (instead of trying both Ensembl and Ensembl Metazoa URLs)

get_genome: get_GTF
	# Download genome from Ensembl databases
	@perl 0Before/get_genome_files.pl  -GTF_dir=$(RNASEQ_VITALIT_GTF) -ensRelease=$(ENSRELEASE) -ensMetazoaRelease=$(ENSMETAZOARELEASE) -outDir=$(RNASEQ_VITALIT_GTF) >$@.tmp 2>&1
	@$(MV) $@.tmp $@
#TODO implement use of new "database" field in rna_seq_sample_info.txt (instead of trying both Ensembl and Ensembl Metazoa URLs)


################################ GENERATE REFERENCE INTERGENIC SEQUENCES ################################

prepare_GTF: get_genome
	# Prepare GTF files
	@for gtf in $(RNASEQ_CLUSTER_GTF)/*.gtf.gz; do \
		echo "#!/bin/sh" > run.sh; \
		echo "$(CLUSTER_R_CMD) R CMD BATCH --vanilla --slave '--args gene_gtf_path=\"$$gtf\" output_gtf_path=\"$${gtf/.gtf.gz/}\"'  0Before/prepare_GTF.R  $@.tmp" >> run.sh; \
		sbatch  --account=mrobinso_bgee --time=12:00:00 --mem=20G --partition=ax-normal --nodes=1 --ntasks=1 --cpus-per-task=1 run.sh ; \
	done
	@$(RM) run.sh
	@echo rm -f $(RNASEQ_CLUSTER_GTF)/*.gtf.gz
	@touch $@

prepare_indexed_transcriptome: get_genome prepare_GTF
	#Preparing indexed transcriptome for every species
	#Extract transcriptome.fa from gtf_all and genome.fa files
	#perl one liner to remove arbitrary numbering in fasta header from gtf_to_fasta
    # Prepare indexes for kallisto: one with default k-mer size $(RNASEQ_KALLISTO_KMER_DEFAULT), one with short k-mer size $(RNASEQ_KALLISTO_KMER_SHORT)
	@perl 1Run/slurm_index_creation.pl -transcriptome_folder=$(RNASEQ_CLUSTER_GTF) -output_log_folder=$(OUTPUT_DIR) -ens_release=$(ENSRELEASE) -ens_metazoa_release=$(ENSMETAZOARELEASE) -short_index_length=$(KALLISTO_SHORT_KMER_SIZE)  -cluster_kallisto_cmd="$(CLUSTER_KALLISTO_CMD)" -cluster_tophat_cmd="$(CLUSTER_TOPHAT_CMD)" >$@.tmp 2>&1
	#TODO Use short k-mer size = 21nt instead of 15
	#TODO The gffread utility in cufflinks package seems more flexible and reliable
	#TODO add step to copy .passw file to bbgee's home on cluster / or maybe just echo a message to tell the user scp it? Path and name of this file on cluster is stored in ENCRYPT_PASSWD_FILE. We should store this file somewhere (devbioinfo?) because it cannot be added to the gitlab project
	#Or just put it in /home/bbgee/? in read only mode, only for bbgee user
	@$(MV) $@.tmp $@
#TODO Use short k-mer size = 21nt instead of 15
#TODO The gffread utility in cufflinks package seems more flexible and reliable
#TODO add step to copy .passw file to bbgee's home on cluster / or maybe just echo a message to tell the user scp it? Path and name of this file on cluster is stored in ENCRYPT_PASSWD_FILE. We should store this file somewhere (devbioinfo?) because it cannot be added to the gitlab project

#NOTE Better to run this step in screen or with nohup as submission is done X jobs per X jobs to not overload the system (and not decrease our user priority)!
abundance_all_intergenic: prepare_indexed_transcriptome $(RNASEQ_SAMPINFO_FILEPATH)
	# Running the pipeline
	@perl 1Run/slurm_scheduler.pl -sample_info_file=$(RNASEQ_SAMPINFO_FILEPATH) -exclude_sample_file=$(RNASEQ_SAMPEXCLUDED_FILEPATH) -output_log_folder=$(RNASEQ_CLUSTER_LOG) -index_folder=$(RNASEQ_CLUSTER_GTF) -fastq_folder=$(RNASEQ_BIGBGEE_FASTQ) -kallisto_out_folder=$(RNASEQ_CLUSTER_ALL_RES) -ens_release=$(ENSRELEASE) -ens_metazoa_release=$(ENSMETAZOARELEASE) -enc_passwd_file=$(ENCRYPT_PASSWD_FILE) -cluster_kallisto_cmd="$(CLUSTER_KALLISTO_CMD)" -cluster_R_cmd="$(CLUSTER_R_CMD)" >$@.tmp 2>&1
	@echo "TODO: At the end it is a good idea to relaunch the abundance_all_intergenic step to be sure everything was run!"
	@$(MV) $@.tmp $@

check_abundance_all_intergenic: abundance_all_intergenic
#TODO - Add a steps to automatically look for problems? See README.md in pipeline/RNA_Seq/ folder
#     - It is probably easier to relaunch the problematic samples manually.
#     - Flagged / excluded samples with low % genes mapped / low number of reads mapped / for which mapping failed: added manually to file $(RNASEQ_SAMPEXCLUDED_FILEPATH).

export_length: check_abundance_all_intergenic
	#Export transcript length to a file to export
	@perl 3Insertion/export_feature_length.pl -library_info=$(RNASEQ_SAMPINFO_FILEPATH) -excluded_libraries=$(RNASEQ_SAMPEXCLUDED_FILEPATH) -tx2gene_dir=$(OUTPUT_DIR) -all_results=$(RNASEQ_CLUSTER_ALL_RES) -length_info=$(RNASEQ_LENGTH_INFO_FILEPATH) > $@.tmp 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@echo -e "\tTranscript length are now exported into $(RNASEQ_LENGTH_INFO_FILEPATH) file\n"
	@$(MV) $@.tmp $@

finalize_pipeline: check_abundance_all_intergenic
	# Collect infos from .report files
	@perl 3Insertion/create_rna_seq_report_info.pl -library_info=$(RNASEQ_SAMPINFO_FILEPATH) -excluded_libraries=$(RNASEQ_SAMPEXCLUDED_FILEPATH) -report_info=$(RNASEQ_CLUSTER_REPORTINFO) -all_results=$(RNASEQ_CLUSTER_ALL_RES) > $@.tmp 2>warnings.$@
	# Touch all files so that they are not removed from /scratch/temporary
	@find $(RNASEQ_CLUSTER_SCRATCH) -exec touch {} \;
# We use tail -n+1 instead of cat because it writes the name of the file in the concatenated file
	@tail -n+1 $(RNASEQ_CLUSTER_LOG)/*/*.err    > $(RNASEQ_CLUSTER_LOG)/all_std_err.txt 2>>warnings.$@
	@tail -n+1 $(RNASEQ_CLUSTER_LOG)/*/*.out    > $(RNASEQ_CLUSTER_LOG)/all_std_out.txt 2>>warnings.$@
	@tail -n+1 $(RNASEQ_CLUSTER_LOG)/*/*.report > $(RNASEQ_CLUSTER_LOG)/all_reports.txt 2>>warnings.$@
	# Back-up all data (not final but it's worth doing an intermediate backup here)
	@tar -cf $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar $(RNASEQ_CLUSTER_ALL_RES) 2>>warnings.$@
	@gzip -9 $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar                           2>>warnings.$@
	# Move back-up data to /data/
	@$(MV) $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar.gz $(RNASEQ_CLUSTER_ALL_RES_BACKUP) 2>>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@


#NOTE For iterative updates the  sum_by_species  step has to be skipped if you want to reuse the previous gaussian curves as they are!
#     So use the previous release $(RNASEQ_CLUSTER_GAUSSIAN_CHOICE) file and $(RNASEQ_CLUSTER_SUM_RES) result folder (Just update the database version if required)
#     for the  presence_absence  step!
sum_by_species: check_pipeline $(RNASEQ_SAMPINFO_FILEPATH) $(RNASEQ_SAMPEXCLUDED_FILEPATH)
# Script using all data from each species to deconvolute the coding genes and intergenic regions underlying distributions
	@mkdir -p $(RNASEQ_CLUSTER_SUM_RES)
	@$(CLUSTER_R_CMD) R CMD BATCH --vanilla --slave '--args rna_seq_sample_info="$(RNASEQ_SAMPINFO_FILEPATH)" rna_seq_sample_excluded="$(RNASEQ_SAMPEXCLUDED_FILEPATH)" kallisto_count_folder="$(RNASEQ_CLUSTER_ALL_RES)" sum_by_species_folder="$(RNASEQ_CLUSTER_SUM_RES)"' $(RNASEQ_CLUSTER_SCRIPTS)$(RNASEQPATH)1Run/rna_seq_sum_by_species.R $(RNASEQ_CLUSTER_SUM_RES)rna_seq_sum_by_species.Rout
	@echo -e "You HAVE TO create/update MANUALLY a file [$(RNASEQ_CLUSTER_GAUSSIAN_CHOICE)] with selected gaussians for coding and intergenic regions!\n"
	@touch $@

## variables used only to generate intergenic sequences.
## should potentialy be moved at the top of the file or in the Makefile.common file
TRANSCRIPTOME_COMPRESSION_EXT	:= .xz
SUM_ABUNDANCE_FILE_PATH			:= $(RNASEQ_VITALIT_SUM_RES)sum_abundance_gene_level+fpkm+intergenic+classification_SPECIES_ID.tsv
## Generate 2 intergenic fasta files for each species. One for reference intergenic sequences and one for other reference intergenic sequences
##
generate_intergenic_sequences: sum_by_species
	@perl create_intergenic_fasta.pl -bgee=$(BGEECMD) -ens_release=$(ENSRELEASE) -transcriptomes_folder=$(RNASEQ_VITALIT_GTF) -transcriptome_compression_ext=$(TRANSCRIPTOME_FILE_EXT) -sum_abundance_file_path=$(SUM_ABUNDANCE_FILE_PATH) -gaussian_file_path=$(RNASEQ_VITALIT_GAUSSIAN_CHOICE) -ref_intergenic_dir=$(VITALIT_REF_INTERGENIC_FOLDER) -other_intergenic_dir=$(VITALIT_OTHER_INTERGENIC_FOLDER) >> $@.tmp 2> $@.warn
	@perl create_intergenic_fasta.pl -bgee=$(BGEECMD) -ens_release=$(ENSMETAZOARELEASE) -transcriptomes_folder=$(RNASEQ_VITALIT_GTF) -transcriptome_compression_ext=$(TRANSCRIPTOME_FILE_EXT) -sum_abundance_file_path=$(SUM_ABUNDANCE_FILE_PATH) -gaussian_file_path=$(RNASEQ_VITALIT_GAUSSIAN_CHOICE) -ref_intergenic_dir=$(VITALIT_REF_INTERGENIC_FOLDER) -other_intergenic_dir=$(VITALIT_OTHER_INTERGENIC_FOLDER) >> $@.tmp 2> $@.warn

	# compress intergenic fasta files
	@find $(VITALIT_REF_INTERGENIC_FOLDER) -type f -name '*_intergenic.fa' -exec gzip --verbose --best {} \; >> $@.tmp 2> $@.warn
	@find $(VITALIT_OTHER_INTERGENIC_FOLDER) -type f -name '*_intergenic.fa' -exec gzip --verbose --best {} \; >> $@.tmp 2> $@.warn


################################ GENERATE PRESENT/ABSENT CALLS ################################

# Run BgeeCall R package on all libraries
# generate index for all species
#      	- based on sample info file outside of BgeeCall?
#	   	- based on sample info file using one function present in BgeeCall (to develop) ?
# reuse script generating abundance with all intergenic and adapt it to generate abundance per library using BgeeCall


presence_absence: sum_by_species $(RNASEQ_SAMPINFO_FILEPATH) $(RNASEQ_SAMPEXCLUDED_FILEPATH) $(RNASEQ_VITALIT_GAUSSIAN_CHOICE)
#	Script to calculate presence/absence threshold for all samples and call genes present/absent + generate density plots + generate statistics on % present (cutoff info files)
	@mkdir -p $(RNASEQ_CLUSTER_PRESENCE_RES)
	@$(CLUSTER_R_CMD) R CMD BATCH --vanilla --slave '--args rna_seq_sample_info="$(RNASEQ_SAMPINFO_FILEPATH)" rna_seq_sample_excluded="$(RNASEQ_SAMPEXCLUDED_FILEPATH)" kallisto_count_folder="$(RNASEQ_CLUSTER_ALL_RES)" sum_by_species_folder="$(RNASEQ_CLUSTER_SUM_RES)" gaussian_choice="$(RNASEQ_CLUSTER_GAUSSIAN_CHOICE)" out_folder="$(RNASEQ_CLUSTER_PRESENCE_RES)" plot_only=FALSE desired_r_cutoff="0.05"' $(RNASEQ_CLUSTER_SCRIPTS)$(RNASEQPATH)1Run/rna_seq_presence_absence.R $(RNASEQ_CLUSTER_PRESENCE_RES)rna_seq_presence_absence.Rout
	@touch $@
#NOTE if output directory is not $(RNASEQ_CLUSTER_ALL_RES): copy all files with calls to this directory, as well as summary stats and plot files
#for folder in *; do echo $folder; /bin/cp $folder/* ../all_results_bgee_v15/$folder/; done
# /bin/cp used because cp is an alias to cp -i
# Also copy presence_absence_all_samples.txt, presence_absence_all_samples.RDa and presence_absence_boxplots.pdf

save_and_send_results_back: finalize_pipeline presence_absence export_length
	# Touch all files so that they are not removed from /scratch/temporary
	@find $(RNASEQ_CLUSTER_SCRATCH) -exec touch {} \; 2>warnings.$@
	# Back-up all data
	@tar -cf $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar --absolute-names $(RNASEQ_CLUSTER_ALL_RES) $(RNASEQ_CLUSTER_SUM_RES) $(RNASEQ_CLUSTER_PRESENCE_RES) 2>>warnings.$@
	@gzip -9 $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar 2>>warnings.$@
	# Move back-up data to /data/
	@$(MV) $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar.gz $(RNASEQ_CLUSTER_ALL_RES_BACKUP) 2>>warnings.$@
	# Send results back to our server: $(PIPELOGIN)@$(PIPEHOST)  password: $(PIPEPASSWD)
	@ssh $(PIPELOGIN)@$(PIPEHOST) mkdir -p $(RNASEQALLRES)
	# @scp $(RNASEQ_CLUSTER_SCRATCH)all_results_$(DBNAME).tar.gz  $(PIPELOGIN)@$(PIPEHOST):$(RNASEQALLRES) >$@.tmp 2>&1
	# Whole archive is probably too big to be copied to our servers (>100Gb). Copy only the final gene-level expression + calls files:
	@rsync -rav -e ssh --include='*/' --include='*abundance_gene_level+new_tpm+new_fpkm+calls.tsv' --exclude='*' $(RNASEQ_CLUSTER_PRESENCE_RES) $(PIPELOGIN)@$(PIPEHOST):$(RNASEQALLRES) >>$@.tmp 2>&1
	@$(CP) $(RNASEQ_CLUSTER_REPORTINFO) $(RNASEQREPORTINFO)
	@$(CP) $(RNASEQ_CLUSTER_PRESENCE_RES)/presence_absence_all_samples.txt $(RNASEQSAMPSTATS)
	@rsync -av -e ssh $(RNASEQ_CLUSTER_REPORTINFO) $(RNASEQ_CLUSTER_PRESENCE_RES)/presence_absence_all_samples.txt $(RNASEQ_CLUSTER_PRESENCE_RES)/presence_absence_all_samples.RDa $(RNASEQ_CLUSTER_PRESENCE_RES)/presence_absence_boxplots.pdf $(PIPELOGIN)@$(PIPEHOST):$(RNASEQALLRES)/ >>$@.tmp 2>&1
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

cluster2: save_and_send_results_back
	@echo -e "\tYou can go out of cluster now\n"
	@touch $@


################################ INSERT CALLS IN DATABASE ################################

backup: $(RNASEQALLRES)
	# Rsync files for backup by '$(DATALOGIN)' user on '$(DATAHOST)' (see wiki for password)
	@echo -e "\t\e[1m\e[37m\e[41mRsync files in '$(RNASEQALLRES)' to backup (NAS) storage!\e[m\n"
	@ssh $(DATALOGIN)@$(DATAHOST) mkdir -p $(DATAPATH)/rna_seq/all_results_$(DBNAME)/
	@rsync -av --del -f "- /lost+found" --exclude=.svn --exclude=.git $(RNASEQALLRES)  $(DATALOGIN)@$(DATAHOST):$(DATAPATH)/rna_seq/all_results_$(DBNAME)/  >   $@.tmp
	@$(MV) $@.tmp $@
#TODO backup log files too + backup other files? We should keep as much as we can otherwise it's impossible to trace back things after a run

#NOTE For iterative updates die commands after insert/update/delete statement failure have to be changed to warn
#     because with iterative updates the database may contain duplicates in regard to what you are inserting/updating/deleting!
insert_RNA_seq: backup $(RNASEQSAMPSTATS) $(RNASEQREPORTINFO) $(RNASEQ_LIB_FILEPATH_FULL) $(RNASEQ_EXPERIMENT_FILEPATH_FULL) $(UBERON_SEX_INFO_FILE_PATH) $(CUSTOM_UBERON_FILE_PATH) $(DEV_STAGE_ONT_FILE_PATH) $(RNASEQ_SAMPINFO_FILEPATH) $(RNASEQ_SAMPEXCLUDED_FILEPATH) $(EXTRAMAPPING_FILEPATH)
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
	@$(MYSQL) -e "SELECT DISTINCT t1.strain FROM cond AS t1 INNER JOIN rnaSeqLibrary AS t2 ON t1.conditionId = t2.conditionId ORDER BY t1.strain" >> $@.tmp
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

TMP_TABLE_CLAUSE   := SELECT t10.expressionId, IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'present' AND t11.rnaSeqExperimentCallQuality = 'high quality') >= 2, 'pstGold', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'present' AND t11.rnaSeqExperimentCallQuality = 'high quality') = 1, 'pstSilverFromHigh', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'present' AND t11.rnaSeqExperimentCallQuality = 'poor quality') >= 2, 'pstSilverFromLow', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'present' AND t11.rnaSeqExperimentCallQuality = 'poor quality') = 1, 'pstBronze', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'absent' AND t11.rnaSeqExperimentCallQuality = 'high quality') >= 2, 'absGold', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'absent' AND t11.rnaSeqExperimentCallQuality = 'high quality') = 1, 'absSilverFromHigh', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'absent' AND t11.rnaSeqExperimentCallQuality = 'poor quality') >= 2, 'absSilverFromLow', IF((SELECT COUNT(rnaSeqExperimentId) FROM rnaSeqExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.rnaSeqExperimentCallDirection = 'absent' AND t11.rnaSeqExperimentCallQuality = 'poor quality') = 1, 'absBronze', 'IMPOSSIBLE')))))))) AS expressionCall FROM rnaSeqExperimentExpression AS t10

check_conditions: insert_RNA_seq
	@echo "Conditions with anat. entity not existing in related species: " > $@.tmp
	@$(MYSQL) -e "SELECT DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain FROM cond AS t1 LEFT OUTER JOIN anatEntityTaxonConstraint AS t2 ON t1.anatEntityId = t2.anatEntityId AND (t2.speciesId IS NULL OR t1.speciesId = t2.speciesId) LEFT OUTER JOIN anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId LEFT OUTER JOIN rnaSeqLibrary AS t10 ON t1.conditionId = t10.conditionId WHERE t2.anatEntityId IS NULL AND t10.conditionId IS NOT NULL ORDER BY t1.speciesId" >> $@.tmp 2> $@.warn
	@echo "Conditions with dev. stage not existing in related species: " >> $@.tmp
	@$(MYSQL) -e "SELECT DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain FROM cond AS t1 LEFT OUTER JOIN stageTaxonConstraint AS t2 ON t1.stageId = t2.stageId AND (t2.speciesId IS NULL OR t1.speciesId = t2.speciesId) LEFT OUTER JOIN anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId LEFT OUTER JOIN rnaSeqLibrary AS t10 ON t1.conditionId = t10.conditionId WHERE t2.stageId IS NULL AND t10.conditionId IS NOT NULL ORDER BY t1.speciesId" >> $@.tmp 2>> $@.warn
	@$(MV) $@.tmp $@

insert_expression: check_conditions insert_RNA_seq
	# Insert the expression summaries
	@perl 3Insertion/insert_rna_seq_expression.pl -bgee=$(BGEECMD) > $@.tmp 2>warnings.$@
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

insert_feature_length: $(RNASEQ_LENGTH_INFO_FILEPATH)
	# Insert the feature length information
	@perl 3Insertion/insert_feature_length.pl -bgee=$(BGEECMD) -length_info=$(RNASEQ_LENGTH_INFO_FILEPATH) > $@.tmp 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

launch_calculate_TMM_factors: insert_RNA_seq $(RNASEQALLRES)
	# Launch calculation of TMM factors for RNA-seq
	@perl 3Insertion/launch_calculate_TMM_factors.pl -bgee=$(BGEECMD) -path_generes=$(RNASEQALLRES) -path_target=$(RNASEQTMMTARG) -path_processed=$(RNASEQTMMPATH) >$@.tmp 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

check_TMM_factors: launch_calculate_TMM_factors $(RNASEQTMMPATH)
	# Check results
	@grep -H -c 'proc.time' $(RNASEQTMMPATH)/*.log | grep ':0$$'  >$@.tmp
	# Check warnings
	@grep 'Warning' -A 5 $(RNASEQTMMPATH)/*.log                   >>$@.tmp
	@$(MV) $@.tmp $@
#TODO check that all libraries in database have a TMM factor calculated:
# cat /var/bgee/extra/pipeline/rna_seq/processed_TMM_$(DBNAME)/*.tsv | grep -v rnaSeqExperimentId | wc -l
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
	-@$(GIT) commit $@ -m "Commit verification file for RNA Seq"
	-@$(GIT) push

#TODO Clean added/removed steps
#TODO Remove unused scripts/files in folders and sub-folders
clean:
	-@$(RM) -R $(VERIFICATIONFILE)  check_annot check_tools create_rna_seq_sample_info get_annot get_GTF prepare_GTF get_genome prepare_indexed_transcriptome generate_genome_info prepare_run_info_final create_config_file cluster backup check_curation check_runs insert_RNA_seq insert_expression wormbase external $(RNASEQGTFDATAPATH)/*.gtf $(RNASEQGTFDATAPATH)/*.gtf.gz *.tmp check_conditions insert_expression
	-@$(GIT) rm $(VERIFICATIONFILE)
	-@$(GIT) commit $(VERIFICATIONFILE) -m "Removing verification file for RNA Seq"
	-@$(GIT) push

# special target to remove RNA_seq from database. We make the name long on purpose, and does not
# incude it in the clean target, to not wipe the database accidentally
.PHONY := $(.PHONY) deleteRNASeq
deleteRNASeq:
	@$(MYSQL) -e "DELETE FROM rnaSeqResult; DELETE FROM rnaSeqRun; DELETE FROM rnaSeqLibrary; DELETE FROM rnaSeqExperimentToKeyword; DELETE FROM rnaSeqExperiment; DELETE FROM rnaSeqPlatform;"
#TODO check if new tables were added
