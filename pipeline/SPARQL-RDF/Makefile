########
# Run this makefile for each mapping file, see commented and uncommented $(GENEX_MAPPING_FILE) variables at $(PIPELINEROOT)Makefile.common
# The mappings are split to avoid query processing overload with the MySQL server during tiple materialization.
########
PIPELINEROOT := ../
DIR_NAME := SPARQL-RDF/
include $(PIPELINEROOT)Makefile.common

FTP_RELEASE   := $(subst .,_,$(RELEASE))

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
	[ -d "$(ONTOP_INSTALL_DIR)" ] && rm -rf $(ONTOP_INSTALL_DIR)
	@mkdir $(ONTOP_INSTALL_DIR)
	# Uncompress ontop
	@unzip $(INPUT_DIR)$(ONTOP_DOWNLOAD_FILE) -d $(ONTOP_INSTALL_DIR) 1>> $@.tmp 2>> $@.err
	#add MySQL connector lib for Ontop
	@wget $(ONTOP_MYSQL_JDBC_URL) -O $(ONTOP_INSTALL_DIR)lib/mysql_connector.jar 1>> $@.tmp 2>> $@.err
	$(RM) $(INPUT_DIR)$(ONTOP_DOWNLOAD_FILE)
	$(MV) $@.tmp $@

create_property_file: $(ONTOP_BGEE_PROP_TEMPLATE)
	# Create a property file for server connection to database
	@$(CP) $(ONTOP_BGEE_PROP_TEMPLATE) $(OUTPUT_DIR)bgee_connection.properties
	@echo -e 'jdbc.url=jdbc:mysql://$(DBHOST_EASYBGEE):$(DBPORT)/$(DBNAME_EASYBGEE)\njdbc.driver=com.mysql.jdbc.Driver\njdbc.user=$(DBUSER_EASYBGEE)\njdbc.name=$(DBNAME_EASYBGEE)\njdbc.password=$(DBPASS_EASYBGEE)' >> $(OUTPUT_DIR)bgee_connection.properties
	@touch $@

generate_triples: create_property_file download_and_install_ontop download_ontology_file download_mapping_file
	@./bgee_rdf.sh -m $(INPUT_DIR)$(GENEX_MAPPING_FILE) -o $(TTL_DIRECTORY_FULL_PATH) -p $(OUTPUT_DIR)bgee_connection.properties -x $(ONTOP_INSTALL_DIR) -t $(INPUT_DIR)$(GENEX_ONTOLOGY_FILE) -v $(FTP_RELEASE) -j '-Xmx128G' > $@.tmp 2> $@.err
	$(MV) $@.tmp $@


create_isql_load_script: generate_triples
	@echo "Generating iSQL script to load RDF data into a virtuoso data store..."
	#delete $(BGEE_VIRTUOSO_FILE) if already exists
	-f $(BGEE_VIRTUOSO_FILE) && $(RM) $(BGEE_VIRTUOSO_FILE)
	#Uncompress ontop
	@ls $(TTL_DIRECTORY_FULL_PATH)*.ttl > $(TTL_FILE_LIST)
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
	@zip -r  $(ZIP_TTL_FILE) $(TTL_DIRECTORY_FULL_PATH) > $@.tmp $@.err && $(RM) $(TTL_DIRECTORY_FULL_PATH)*.ttl $(BGEE_VIRTUOSO_FILE)
	@scp $(ZIP_TTL_FILE) $(FTPUSER)@$(FTPHOST):$(FTP_SAVE_DUMP_PATH)
	@echo "RDF file was stored on FTP. Please manually remove local version with the command : rm -f $(ZIP_TTL_FILE)"
	@$(MV) $@.tmp $@

$(VERIFICATIONFILE): zip_and_store
	@touch $@

clean:
	-@$(RM) -R $(VERIFICATIONFILE) insert_to_triplestore generate_triples create_property_file download_mapping_file download_ontology_file download_and_install_ontop

