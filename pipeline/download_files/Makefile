PIPELINEROOT := ../
DIR_NAME := download_files/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)

############
# Variables to create specific files
###########
EXPR_SIMPLE_CAT := expr_simple
EXPR_COMPLETE_CAT := expr_advanced
ALL_CALLS_CATS := $(EXPR_SIMPLE_CAT)$(VAL_SEP)$(EXPR_COMPLETE_CAT)

RNA_SEQ_ANNOT_CAT := rnaseq_annot
RNA_SEQ_DATA_CAT := rnaseq_data
ALL_RAW_DATA_CATS := $(RNA_SEQ_ANNOT_CAT)$(VAL_SEP)$(RNA_SEQ_DATA_CAT)

ORGAN_PARAM := ANAT_ENTITY_ID
CELLTYPE_PARAM := CELL_TYPE_ID
STAGE_PARAM := DEV_STAGE_ID
SEX_PARAM := SEX_ID
STRAIN_PARAM := STRAIN_ID
ALL_COND_PARAM := $(ORGAN_PARAM)$(ENTRY_SEP)$(STAGE_PARAM)$(ENTRY_SEP)$(SEX_PARAM)$(ENTRY_SEP)$(STRAIN_PARAM)

######################################################
# VARIABLES CREATED, INSERT
######################################################

# =========================
# Processed data download files
# =========================
$(OUTPUT_DIR)generate_processed_files:
	@perl generate_ref_expr_files.pl -bgee=$(BGEECMD) -speciesArg=- -parallelJobs=10 -rnaSeqDir=$(DOWNLOAD_FILE_ROOT_PATH)$(PROC_EXPR_VALUES_PATH)rna_seq_processed/ -flScRnaSeqDir=$(DOWNLOAD_FILE_ROOT_PATH)$(PROC_EXPR_VALUES_PATH)sc_rna_seq_fl_processed/ -dbScRnaSeqDir=$(DOWNLOAD_FILE_ROOT_PATH)$(PROC_EXPR_VALUES_PATH)sc_rna_seq_db_processed/ -bgeeVersion=$(DBNAME) > $@.tmp 2>$@.warnings
	# generate one tarball per species for bulk RNA-Seq. This tarball is used by the BgeeDB R package.
	cd $(DOWNLOAD_FILE_ROOT_PATH)$(PROC_EXPR_VALUES_PATH)rna_seq_processed/ && \
	for d in */ ; do \
	  species="$$d"; species="$${species%/}"; \
	  if compgen -G "$${species}/$${species}_RNA-Seq_read_counts_TPM_*.tsv.gz" > /dev/null; then \
	    files=("$${species}/$${species}_RNA-Seq_read_counts_TPM_"*.tsv.gz); \
	    tar -czvf "$${species}_RNA-Seq_read_counts_TPM.tar.gz" "$${files[@]}" >> $@.tmp 2>> $@.warnings; \
	  fi; \
	done
	# generate one tarball per species for full-length single-cell RNA-Seq. This tarball is used by the BgeeDB R package.
	@cd $(DOWNLOAD_FILE_ROOT_PATH)$(PROC_EXPR_VALUES_PATH)sc_rna_seq_fl_processed/ && \
	for d in */ ; do \
	  species="$$d"; species="$${species%/}"; \
	  if compgen -G "$${species}/$${species}_Full-Length_SC_RNA-Seq_read_counts_TPM_*.tsv.gz" > /dev/null; then \
	    files=("$${species}/$${species}_Full-Length_SC_RNA-Seq_read_counts_TPM_"*.tsv.gz); \
	    tar -czvf "$${species}_Full-Length_SC_RNA-Seq_read_counts_TPM.tar.gz" "$${files[@]}" >> $@.tmp 2>> $@.warnings; \
	  fi; \
	done
	# generate one tarball per species for droplet-based single-cell RNA-Seq. This tarball is used by the BgeeDB R package.
	@cd $(DOWNLOAD_FILE_ROOT_PATH)$(PROC_EXPR_VALUES_PATH)sc_rna_seq_db_processed/ && \
	for d in */ ; do \
	  species="$$d"; species="$${species%/}"; \
	  if compgen -G "$${species}/$${species}_Droplet-Based_SC_RNA-Seq_read_counts_CPM_*.tsv.gz" > /dev/null; then \
	    files=("$${species}/$${species}_Droplet-Based_SC_RNA-Seq_read_counts_CPM_"*.tsv.gz); \
	    tar -czvf "$${species}_Droplet-Based_SC_RNA-Seq_read_counts_CPM.tar.gz" "$${files[@]}" >> $@.tmp 2>> $@.warnings; \
	  fi; \
	done

	@$(MV) $@.tmp $@

# =========================
# full-length and droplet-based h5ad download files
# =========================
# Path to the uv project containing the Python environment for h5ad generation
H5AD_PROJECT_DIR := $(TMPDIR)/h5ad_env
H5AD_PYTHON := $(H5AD_PROJECT_DIR)/.venv/bin/python

# Create the virtual environment from the uv project (installs dependencies from pyproject.toml)
$(H5AD_PYTHON): $(OUTPUT_DIR)generate_processed_files
	cd $(H5AD_PROJECT_DIR) && uv sync

#TODO: generate files with directory structure and file names as in the FTP server
$(OUTPUT_DIR)generate_h5ad_files: $(H5AD_PYTHON)
	@$(H5AD_PYTHON) generate_h5ad_download_files.py --output_dir $(DOWNLOAD_FILE_ROOT_PATH)$(H5AD_PATH) --server $(DBHOST) --db $(DBNAME) --usr $(DBUSER) --pwd $(DBPASS) --result_dir $(SC_RNASEQ_SERVER_RES_KALLISTO_DROPLET) --intergenic_prefixes $(INTERGENIC_PREFIXES) > $@.tmp 2>$@.warnings
	@$(MV) $@.tmp $@

# =========================
# Gene expression call download files
# =========================
# With the new OTF approach of Bgee 16.0 all observed data are used to generate these files but 2 fitlering are used before writing the calls:
#     * keep only calls with a meta stages (starting with the UBERON namespace)
#     * keep only calls with an observed anat. entity.
#XXX generation of all_cond files are memory intensive especially for the M. musculus one. Options to avoid OOM are:
#     * increase Xmx in Makefile.common
#     * decrease threads used
#     * decrease number of genes in a batch (last attribut of GenerateBasicExprFile)
# The generation of the anatomy call files took less than 1 day and the generation of the all cond. calls files took ~ 3 days for Bgee 16.0 (with 20 threads and 10 genes per batch).
# The generation of the M. musculus all_cond $(EXPR_COMPLETE_CAT) file alone took more than 1.5 days. That file was 98Gb big before compression.
$(OUTPUT_DIR)generate_anatomy_call_files: $(OUTPUT_DIR)generate_h5ad_files
	@$(JAVA_NO_JAR) -Djava.util.concurrent.ForkJoinPool.common.parallelism=10 -jar $(PIPELINEROOT)$(BGEEPIPELINEJAR)  GenerateBasicExprFile $(EMPTY_LIST) $(EXPR_SIMPLE_CAT)$(ENTRY_SEP)$(EXPR_COMPLETE_CAT) $(DOWNLOAD_FILE_ROOT_PATH)$(EXPR_CALLS_PATH) $(ORGAN_PARAM) 100 > $@.tmp 2>$@.warnings
	@$(MV) $@.tmp $@

$(OUTPUT_DIR)generate_all_cond_params_call_files: $(OUTPUT_DIR)generate_anatomy_call_files
	@$(JAVA_NO_JAR) -Djava.util.concurrent.ForkJoinPool.common.parallelism=20 -jar $(PIPELINEROOT)$(BGEEPIPELINEJAR)  GenerateBasicExprFile - $(EXPR_SIMPLE_CAT) $(DOWNLOAD_FILE_ROOT_PATH)$(EXPR_CALLS_PATH) $(ALL_COND_PARAM) 10 > $@.tmp 2>$@.warnings
	@$(MV) $@.tmp $@

# =========================
# topAnat expression call files
# =========================
# Six files per species, read by topAnat: two call files, gzip files with one member per stage,
# with their index, and two files describing the conditions. topAnat cannot propagate a whole
# species inside a web request. The web application only generates missing files at start-up when
# its property org.bgee.webapp.generateTopAnatCallFiles is true, one gene batch at a time; this
# step produces them all, several batches at a time.
# The parallelism is passed to the jar rather than set through the parallelism of the common
# ForkJoinPool, the generation using a pool of its own.
# Bgee 16, parallelism of 10: about 16 h for mouse, the largest species, and 6 h 40 for human.
$(OUTPUT_DIR)generate_topanat_call_files: $(OUTPUT_DIR)generate_all_cond_params_call_files
	@mkdir -p $(TOPANAT_CALL_FILES_PATH)
	@$(JAVA) GenerateTopAnatFiles $(TOPANAT_CALL_FILES_PATH) $(EMPTY_LIST) false $(TOPANAT_PARALLELISM) > $@.temp 2>$@.warnings
	@echo "Content of directory $(TOPANAT_CALL_FILES_PATH)topanat_calls/" >> $@.temp
	@ls -l $(TOPANAT_CALL_FILES_PATH)topanat_calls/ >> $@.temp
	@$(MV) $@.temp $@

# The file is generated to retrieve species information from the database. Until Bgee 15.2 the R package was querying the Java API to retrieve that info.
# Since Bgee 15.2 we added droplet based data. As in the Java API we do not distinguish between full length and droplet based single cell data (they are
# both considered single cell), we had either to create a new datatype EnumSet specific to the download file (called DownloadFileDatatypes) or to generate
# a file available on the FTP containing all the information required (list of species, datatypes available for each species). We decided to go to the direction
# of a file because it was less time consuming and we wanted the R package to be able to download single cell processed expression values and h5ad as soon as
# possible.
#TODO: In the future it would be cleaner to avoid generation of files (species file generated here and release_vX.tsv file) and always query the Java API.
# Since Bgee 16 there are no Affymetrix and EST data any more, and no table to query for them: their columns
# are kept, always false, as BgeeDB reads them.
$(OUTPUT_DIR)bgeeDBRpackageSpeciesFile:
	@$(MYSQL) -e "select t1.speciesId AS \`ID\`, t1.genus AS \`GENUS\`, t1.species AS \`SPECIES_NAME\`, t1.speciesCommonName AS \`COMMON_NAME\`, 'F' AS \`AFFYMETRIX\`, 'F' AS \`EST\`, IF( EXISTS( select 1 from inSituSpot as t2 inner join cond as t3 on t2.conditionId = t3.conditionId where t3.speciesId = t1.speciesId limit 1), 'T', 'F') AS \`IN_SITU\`, IF( EXISTS( select 1 from rnaSeqLibraryAnnotatedSample as t2 inner join cond as t3 on t2.conditionId = t3.conditionId inner join rnaSeqLibrary AS t4 ON t2.rnaSeqLibraryId = t4.rnaSeqLibraryId where t3.speciesId = t1.speciesId AND t4.rnaSeqTechnologyIsSingleCell = 0 limit 1), 'T', 'F') AS \`RNA_SEQ\`, IF( EXISTS( select 1 from rnaSeqLibraryAnnotatedSample as t2 inner join cond as t3 on t2.conditionId = t3.conditionId inner join rnaSeqLibrary AS t4 ON t2.rnaSeqLibraryId = t4.rnaSeqLibraryId where t3.speciesId = t1.speciesId AND t4.rnaSeqTechnologyIsSingleCell = 1 AND t2.multipleLibraryIndividualSample = 0 limit 1), 'T', 'F') AS \`SC_FULL_LENGTH\`, IF( EXISTS( select 1 from rnaSeqLibraryAnnotatedSample as t2 inner join cond as t3 on t2.conditionId = t3.conditionId inner join rnaSeqLibrary AS t4 ON t2.rnaSeqLibraryId = t4.rnaSeqLibraryId where t3.speciesId = t1.speciesId AND t4.rnaSeqTechnologyIsSingleCell = 1 AND t2.multipleLibraryIndividualSample = 1 limit 1), 'T', 'F') AS \`SC_DROPLET_BASED\` from species AS t1 group by ID" > $(R_SPECIES_INFO) 2> $@.warnings
	@echo -e "\n\n\n\n\n" >> $(R_SPECIES_INFO)
	@touch $@

#TODO check it, not used because of usage of different servers
$(OUTPUT_DIR)gzip_call_files: $(OUTPUT_DIR)bgeeDBRpackageSpeciesFile
	@find $(DOWNLOAD_FILE_ROOT_PATH)$(EXPR_CALLS_PATH) -type f ! -name '*.gz' ! -name '*.tmp' -exec gzip --verbose --best {} \;

# =========================
# Insert info about download files
# =========================
#TODO instead of using insert_download_file_info.sql to insert data into downloadFile table, fix InsertSpeciesDataGroups
$(OUTPUT_DIR)insert_download_files_info: $(OUTPUT_DIR)gzip_call_files $(OUTPUT_DIR)generate_processed_files $(OUTPUT_DIR)generate_anatomy_call_files $(OUTPUT_DIR)generate_all_cond_params_call_files
	#insert information about species data groups
	@$(MYSQL) < insert_download_file_info.sql > $@.tmp 2>$@.warnings
	#insert information about download files
	@perl insert_download_file_info.pl -bgee=$(BGEECMD) -processed_dir=$(DOWNLOAD_FILE_ROOT_PATH)$(PROC_EXPR_VALUES_PATH) -calls_dir=$(DOWNLOAD_FILE_ROOT_PATH)$(CALLS_PATH) --h5ad_dir=$(DOWNLOAD_FILE_ROOT_PATH)$(H5AD_PATH)> $@.tmp 2>$@.warnings
	@$(MYSQL) -e "SELECT (SELECT COUNT(*) FROM speciesDataGroup) as \`data group count\`, (SELECT COUNT(distinct speciesId) FROM speciesToDataGroup) as \`species in data group count\`, (SELECT COUNT(*) FROM downloadFile) as \`download file count\`" >> $@.tmp 2>>$@.warnings
	@$(MV) $@.tmp $@

######################################################
######################################################

# =========================
# XRefs with expression information
# =========================
# IDs of data source Uniprot/SPTREMBL and Uniprot/SWISSPROT in our database 
# (see pipeline/db_creation/insert_data_sources.sql) 
# Very important that these are correct.
TREMBL_ID := 4
SWISSPROT_ID := 5

# XRef file: 
XREF_FILE := $(OUTPUT_DIR)XRefBgee.txt

# The SQL query run to generate the XRefFile, in readable format: 
# 
# SELECT DISTINCT t1.XRefId as uniprotId, uniqueMapping.uniqueGeneId as ensemblId, uniqueMapping.speciesId as speciesId
# FROM geneXRef AS t1 
# INNER JOIN (
# -- Some xRefIds are mapped to multiple genes, we retrieve a unique gene ID with expression data
#     SELECT t1.XRefId, (SELECT t4.geneId FROM geneXRef AS t2 
#                        INNER JOIN gene AS t4 ON t2.bgeeGeneId = t4.bgeeGeneId
#                        WHERE t2.XRefId = t1.XRefId 
#                        AND EXISTS (SELECT 1 FROM expression as t3 WHERE t3.bgeeGeneId = t2.bgeeGeneId)
#                        ORDER BY t4.geneId LIMIT 1) AS uniqueGeneId,
#						 (SELECT t4.geneId FROM geneXRef AS t2 
#                        INNER JOIN gene AS t4 ON t2.bgeeGeneId = t4.bgeeGeneId
#                        WHERE t2.XRefId = t1.XRefId 
#                        AND EXISTS (SELECT 1 FROM expression as t3 WHERE t3.bgeeGeneId = t2.bgeeGeneId)
#                        ORDER BY t4.geneId LIMIT 1) AS speciesId 
#     FROM geneXRef AS t1 
# ) AS uniqueMapping ON t1.XRefId = uniqueMapping.XRefId 
# WHERE t1.dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND t1.XRefId IS NOT NULL AND t1.XRefId != '' AND t1.XRefId NOT LIKE '%\_%' 
# AND uniqueMapping.uniqueGeneId IS NOT NULL 
# ORDER BY uniqueMapping.uniqueGeneId, t1.XRefId;

XREF_SQL_QUERY := SELECT DISTINCT t1.XRefId as uniprotId, uniqueMapping.uniqueGeneId as geneId, uniqueMapping.speciesId as speciesId FROM geneXRef AS t1 INNER JOIN (SELECT t1.XRefId, (SELECT t4.geneId FROM geneXRef AS t2 INNER JOIN gene AS t4 ON t2.bgeeGeneId = t4.bgeeGeneId WHERE t2.XRefId = t1.XRefId AND EXISTS (SELECT 1 FROM expression as t3 WHERE t3.bgeeGeneId = t2.bgeeGeneId) ORDER BY t4.geneId LIMIT 1) AS uniqueGeneId, (SELECT t4.speciesId FROM geneXRef AS t2 INNER JOIN gene AS t4 ON t2.bgeeGeneId = t4.bgeeGeneId WHERE t2.XRefId = t1.XRefId AND EXISTS (SELECT 1 FROM expression as t3 WHERE t3.bgeeGeneId = t2.bgeeGeneId) ORDER BY t4.geneId LIMIT 1) AS speciesId  FROM geneXRef AS t1) AS uniqueMapping ON t1.XRefId = uniqueMapping.XRefId WHERE t1.dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND t1.XRefId IS NOT NULL AND t1.XRefId != '' AND t1.XRefId NOT LIKE '%\_%' AND uniqueMapping.uniqueGeneId IS NOT NULL ORDER BY uniqueMapping.uniqueGeneId, t1.XRefId;

# Query Wikidata directly (federated SERVICE calls are refused by the Bgee endpoint).
# WDQS ignores the format parameter, hence the Accept header in the wget call.
WIKIDATA_UBERON_IDS := "https://query.wikidata.org/sparql?query=SELECT%20DISTINCT%20%3Fuberon_id_prefixed%20WHERE%20%7B%0A%20%20%7B%20%3Fitem%20%3Chttp%3A%2F%2Fwww.wikidata.org%2Fprop%2Fdirect%2FP1554%3E%20%3Fuberon_id%20.%0A%20%20%20%20BIND%28CONCAT%28%22UBERON%3A%22%2C%20STR%28%3Fuberon_id%29%29%20AS%20%3Fuberon_id_prefixed%29%20%7D%0A%20%20UNION%0A%20%20%7B%20%3Fitem%20%3Chttp%3A%2F%2Fwww.wikidata.org%2Fprop%2Fdirect%2FP7963%3E%20%3Fonto_cl_id%20.%0A%20%20%20%20BIND%28REPLACE%28STR%28%3Fonto_cl_id%29%2C%20%22_%22%2C%20%22%3A%22%29%20AS%20%3Fuberon_id_prefixed%29%20%7D%0A%7D"

UBERON_WIKIDATA_FILE_PATH := $(OUTPUT_DIR)uberon_wikidata.txt

XREFS_TO_GENERATE := WIKIDATA,UNIPROT,GENE_CARDS,BGEE_GENE_SUMMARY

##TODO Remove the XREF_SQL_QUERY query and implement all steps of uniprot XRef generation in Java
#Generate XRefs for different datasources. As for Bgee 15.0 these datasources are uniprot, genecards and wikidata
$(OUTPUT_DIR)generate_xrefs: $(OUTPUT_DIR)insert_download_files_info
	@$(MYSQL) -e "$(XREF_SQL_QUERY)" > $(XREF_FILE).tmp 2>$@.warnings
	@$(WGET) --header="Accept: text/csv" -O $(UBERON_WIKIDATA_FILE_PATH) $(WIKIDATA_UBERON_IDS) > $@.temp 2> $@.warnings
	@$(JAVA_NO_JAR) -Djava.util.concurrent.ForkJoinPool.common.parallelism=2 -Dorg.bgee.core.version.minor=$(MINOR_RELEASE) -Dorg.bgee.core.version.major=$(MAJOR_RELEASE) -jar $(PIPELINEROOT)$(BGEEPIPELINEJAR) GenerateXRefs $(XREF_FILE).tmp $(XREFS_PATH) $(UBERON_WIKIDATA_FILE_PATH) $(XREFS_TO_GENERATE) > $@.temp 2> $@.warnings
	# create a verification file
	# First, we display the targeted data sources, so that we can detect 
	# incorrect data source IDs used (TREMBL_ID and SWISSPROT_ID)
	@echo "Data sources used for mapping \nIt is very important that it corresponds to TREMBL and SWISSPROT: " >> $@.temp
	@$(MYSQL) -e "SELECT * FROM dataSource WHERE dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID))" >> $@.temp  2>>$@.warnings
	@echo "\n\nNumber of expected mappings in XRef file: " >> $@.temp
	@$(MYSQL) -e "SELECT COUNT(DISTINCT t1.XRefId) FROM geneXRef AS t1 WHERE dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND EXISTS (SELECT 1 FROM expression AS t2 WHERE t2.bgeeGeneId = t1.bgeeGeneId)" >> $@.temp 2>>$@.warnings
	@echo "\nNumber of lines in UniProt XRef file: " >> $@.temp
	@wc -l $(XREF_FILE).tmp >> $@.temp
	@echo "\nNumber of mappings per species: " >> $@.temp
	@$(MYSQL) -e "SELECT t3.speciesId, t4.speciesCommonName, CONCAT(t4.genus, ' ', t4.species) AS species, COUNT(DISTINCT t1.XRefId) FROM geneXRef AS t1 INNER JOIN gene AS t3 ON t1.bgeeGeneId = t3.bgeeGeneId INNER JOIN species AS t4 ON t3.speciesId = t4.speciesId WHERE t1.dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND EXISTS (SELECT 1 FROM expression AS t2 WHERE t2.bgeeGeneId = t1.bgeeGeneId) GROUP BY t3.speciesId" >> $@.temp 2>>$@.warnings
	@echo "\n\nNumber of ambiguous mappings per species (we forced the choice of one mapping in these cases): " >> $@.temp
	@$(MYSQL) -e "SELECT multipleGeneCount.speciesId, t4.speciesCommonName, CONCAT(t4.genus, ' ', t4.species) AS species, COUNT(*) FROM (SELECT t1.XRefId, t2.speciesId FROM geneXRef AS t1 INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId WHERE t1.dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND EXISTS (SELECT 1 FROM expression as t3 WHERE t3.bgeeGeneId = t1.bgeeGeneId) GROUP BY t1.XRefId HAVING COUNT(DISTINCT t1.bgeeGeneId) > 1) AS multipleGeneCount INNER JOIN species AS t4 ON multipleGeneCount.speciesId = t4.speciesId  GROUP BY multipleGeneCount.speciesId" >> $@.temp 2>>$@.warnings
	@echo "\nTotal number of ambiguous mappings: " >> $@.temp
	@$(MYSQL) -e "SELECT COUNT(*) FROM (SELECT t1.XRefId FROM geneXRef AS t1 WHERE t1.dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND EXISTS (SELECT 1 FROM expression as t3 WHERE t3.bgeeGeneId = t1.bgeeGeneId) GROUP BY t1.XRefId HAVING COUNT(DISTINCT t1.bgeeGeneId) > 1) AS multipleGeneCount" >> $@.temp 2>>$@.warnings
	@$(MV) $@.temp $@

$(OUTPUT_DIR)insert_expression_summary:
	@perl insert_expression_summary.pl -bgee=$(BGEECMD) -expression_summary_file=$(GENE_SUMMARY_FILE_PATH) > $@.tmp 2> $@.err
	@$(MV) $@.tmp $@

# ======================
# Verification file
# ======================
$(VERIFICATIONFILE): $(OUTPUT_DIR)insert_download_files_info $(OUTPUT_DIR)xref_uniprot_info $(OUTPUT_DIR)rank_files_info $(OUTPUT_DIR)generate_processed_files $(OUTPUT_DIR)generate_anatomy_call_files $(OUTPUT_DIR)generate_anatomy_development_call_files $(OUTPUT_DIR)gzip_call_files $(OUTPUT_DIR)generate_topanat_call_files

clean:
	-@$(RM) -R $(VERIFICATIONFILE) $(OUTPUT_DIR)xref_uniprot_info.tmp $(OUTPUT_DIR)xref_uniprot_info $(XREF_FILE).tmp.noheader $(XREF_FILE).tmp $(XREF_FILE) $(OUTPUT_DIR)insert_download_files_info $(OUTPUT_DIR)rank_files_info $(OUTPUT_DIR)generate_processed_files $(OUTPUT_DIR)generate_anatomy_call_files $(OUTPUT_DIR)generate_anatomy_development_call_files $(OUTPUT_DIR)gzip_call_files $(OUTPUT_DIR)generate_topanat_call_files
