########
# Run this makefile for each mapping file, see commented and uncommented $(GENEX_MAPPING_FILE) variables at $(PIPELINEROOT)Makefile.common
# The mappings are split to avoid query processing overload with the MySQL server during tiple materialization.
########
PIPELINEROOT := ../
DIR_NAME := SPARQL-RDF/
include $(PIPELINEROOT)Makefile.common

FTP_RELEASE   := $(subst .,_,$(RELEASE))
# Number of slices of genes materialized in parallel by generate_triples
RDF_SLICES    := 12

all: $(VERIFICATIONFILE)


download_mapping_file:
	# Download mapping file
	@$(WGET) $(GENEX_MAPPING_URL) -O $(GENEX_MAPPING_FILE) && $(MV) $(GENEX_MAPPING_FILE) $(INPUT_DIR) > $@.tmp
	@$(RM) $(GENEX_MAPPING_FILE)
	@touch $@
	$(MV) $@.tmp $@

download_ontology_file: download_mapping_file
	# Download OWL ontology file
	@$(WGET) $(GENEX_ONTOLOGY_URL) -O $(GENEX_ONTOLOGY_FILE) && $(MV) $(GENEX_ONTOLOGY_FILE) $(INPUT_DIR) > $@.tmp
	@$(RM) $(GENEX_ONTOLOGY_FILE)
	$(MV) $@.tmp $@

download_and_install_ontop: download_ontology_file
	# Download ontop
	@wget $(ONTOP_DOWNLOAD_URL) -O $(ONTOP_DOWNLOAD_FILE) && $(MV) $(ONTOP_DOWNLOAD_FILE) $(INPUT_DIR) 1>$@.tmp 2>$@.err
	#check ontop dir already exists
	@rm -rf $(ONTOP_INSTALL_DIR)
	@mkdir -p $(ONTOP_INSTALL_DIR)
	# Uncompress ontop
	@unzip $(INPUT_DIR)$(ONTOP_DOWNLOAD_FILE) -d $(ONTOP_INSTALL_DIR) 1>> $@.tmp 2>> $@.err
	#add MySQL connector lib for Ontop
	@wget $(ONTOP_MYSQL_JDBC_URL) -O $(ONTOP_INSTALL_DIR)lib/mysql_connector.jar 1>> $@.tmp 2>> $@.err
	$(RM) $(INPUT_DIR)$(ONTOP_DOWNLOAD_FILE)
	$(MV) $@.tmp $@

create_property_file: $(ONTOP_BGEE_PROP_TEMPLATE)
	# Create a property file for server connection to database
	# useCursorFetch and defaultFetchSize make the MySQL driver read the rows of a query 10,000 at a time,
	# instead of loading the whole result into memory, which exceeds the Java heap for the largest tables
	@$(CP) $(ONTOP_BGEE_PROP_TEMPLATE) $(OUTPUT_DIR)bgee_connection.properties
	@echo -e 'jdbc.url=jdbc:mysql://$(DBHOST_EASYBGEE):$(DBPORT)/$(DBNAME_EASYBGEE)?useCursorFetch=true&defaultFetchSize=10000\njdbc.driver=com.mysql.jdbc.Driver\njdbc.user=$(DBUSER_EASYBGEE)\njdbc.name=$(DBNAME_EASYBGEE)\njdbc.password=$(DBPASS_EASYBGEE)' >> $(OUTPUT_DIR)bgee_connection.properties
	@touch $@

# One Ontop process works on a single thread: for Bgee 16.0 (1.4 billion rows in globalExpression) the sequential
# bgee_rdf.sh would have taken days. bgee_rdf_parallel.sh starts one Ontop process per slice of genes, plus one
# for the small tables. Each process holds only a few GB, -Xmx is a limit: RDF_SLICES + 1 processes x 12GB.
# The query gives the first gene of each slice, so that the slices hold about the same number of rows of
# globalExpression: k is the slice of a gene, from the number of rows of the genes before it.
# To check the split on a few genes first, see the test mode (-l) of bgee_rdf_parallel.sh.
generate_triples: create_property_file download_and_install_ontop download_ontology_file download_mapping_file
	@bounds=$$($(MYSQL_EASYBGEE) -N -e "SELECT GROUP_CONCAT(first_gene ORDER BY k) FROM (SELECT k, MIN(bgeeGeneId) AS first_gene FROM (SELECT bgeeGeneId, FLOOR((SUM(n) OVER (ORDER BY bgeeGeneId) - n) * $(RDF_SLICES) / SUM(n) OVER ()) AS k FROM (SELECT bgeeGeneId, COUNT(*) AS n FROM globalExpression GROUP BY bgeeGeneId) AS t) AS s GROUP BY k) AS b") && \
	./bgee_rdf_parallel.sh -m $(INPUT_DIR)$(GENEX_MAPPING_FILE) -o $(TTL_DIRECTORY_FULL_PATH) -p $(OUTPUT_DIR)bgee_connection.properties -x $(ONTOP_INSTALL_DIR) -t $(INPUT_DIR)$(GENEX_ONTOLOGY_FILE) -c $(INPUT_DIR)catalog-v001.xml -g "$$bounds" -j '-Xmx12G' > $@.tmp 2> $@.err
	$(MV) $@.tmp $@


create_isql_load_script: generate_triples
	@echo "Generating iSQL script to load RDF data into a virtuoso data store..."
	#delete $(BGEE_VIRTUOSO_FILE) if already exists
	-f $(BGEE_VIRTUOSO_FILE) && $(RM) $(BGEE_VIRTUOSO_FILE)
	#Uncompress ontop
	# find rather than a glob: with the ~36,000 files of Bgee 16.0 the glob exceeds the maximum length of a command
	@find $(TTL_DIRECTORY_FULL_PATH) -maxdepth 1 -name '*.ttl' -printf '%f\n' | sort > $(TTL_FILE_LIST)
	@touch $(BGEE_VIRTUOSO_FILE)
	while IFS= read -r line ; do \
		echo "DB.DBA.TTLP_MT(file_to_string_output('$(TTL_DIRECTORY_FULL_PATH)$$line'), '', 'https://www.bgee.org/rdf_v$(FTP_RELEASE)');" >>  $(BGEE_VIRTUOSO_FILE) ; \
		echo "SPARQL LOAD <http://purl.org/lscr#> INTO <https://www.bgee.org/rdf_v$(FTP_RELEASE)>;" >>  $(BGEE_VIRTUOSO_FILE)  ; \
	done < $(TTL_FILE_LIST)

	@echo "SPARQL LOAD <http://purl.org/genex#> INTO <https://www.bgee.org/rdf_v$(FTP_RELEASE)>;" >>  $(BGEE_VIRTUOSO_FILE)
	@echo "SPARQL LOAD <http://purl.org/lscr#> INTO <https://www.bgee.org/rdf_v$(FTP_RELEASE)>;" >>  $(BGEE_VIRTUOSO_FILE)
	@$(RM) $(TTL_FILE_LIST)
	@touch $@

	#@echo "Please execute manually iSQL script to load RDF data into a virtuoso data store..."
	#@$(ISQL_PATH) $(VIRTUOSO_HOST) $(VIRTUSER) $(VIRTPASS) "EXEC=LOAD $(BGEE_VIRTUOSO_FILE)" > $@.tmp 2> $@.err

zip_and_store: create_isql_load_script
	#zip TTL directory and copy resulting archive to Bgee FTP
	@zip -9 -r  $(ZIP_TTL_FILE) $(TTL_DIRECTORY_FULL_PATH) > $@.tmp $@.err && find $(TTL_DIRECTORY_FULL_PATH) -maxdepth 1 -name '*.ttl' -delete && $(RM) $(BGEE_VIRTUOSO_FILE)
	@scp $(ZIP_TTL_FILE) $(FTPUSER)@$(FTPHOST):$(FTP_SAVE_DUMP_PATH)
	@echo "RDF file was stored on FTP. Please manually remove local version with the command : rm -f $(ZIP_TTL_FILE)"
	@$(MV) $@.tmp $@

$(VERIFICATIONFILE): zip_and_store
	@touch $@

clean:
	-@$(RM) -R $(VERIFICATIONFILE) insert_to_triplestore generate_triples create_property_file download_mapping_file download_ontology_file download_and_install_ontop

