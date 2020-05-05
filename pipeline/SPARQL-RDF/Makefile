PIPELINEROOT := ../
DIR_NAME := SPARQL-RDF/
include $(PIPELINEROOT)Makefile.common

FTP_RELEASE   := $(subst .,_,$(RELEASE))

all: $(VERIFICATIONFILE)


download_mapping_file:
	# Download mapping file
	@$(WGET) $(GENEX_MAPPING_URL) -O $(GENEX_MAPPING_FILE) && $(MV) $(GENEX_MAPPING_FILE) $(INPUT_DIR)/
	@$(RM) $(GENEX_MAPPING_FILE)
	@touch $@

download_ontology_file: 
	# Download OWL ontology file
	@$(WGET) $(GENEX_ONTOLOGY_URL) -O $(GENEX_ONTOLOGY_FILE) && $(MV) $(GENEX_ONTOLOGY_FILE) $(INPUT_DIR)/
	@$(RM) $(GENEX_ONTOLOGY_FILE)
	@touch $@

download_and_install_ontop: 
	# Download ontop
	@$(WGET) $(ONTOP_DOWNLOAD_URL) -O $(ONTOP_DOWNLOAD_FILE) && $(MV) $(ONTOP_DOWNLOAD_FILE) $(INPUT_DIR)/
	#check ontop dir already exists
	! -d $(ONTOP_INSTALL_DIR) && mkdir -p $(ONTOP_INSTALL_DIR)
	# Uncompress ontop
	@unzip $(INPUT_DIR)$(ONTOP_DOWNLOAD_FILE) -d $(ONTOP_INSTALL_DIR)
	$(RM) $(INPUT_DIR)$(ONTOP_DOWNLOAD_FILE)
	@touch $@

create_property_file: $(ONTOP_BGEE_PROP_TEMPLATE)
	# Create a property file for server connection to database
	@$(CP) $(ONTOP_BGEE_PROP_TEMPLATE) $(OUTPUT_DIR)bgee_connection.properties
	@echo -e 'jdbc.url=jdbc:mysql://$(DBHOST_EASYBGEE):$(DBPORT)/$(DBNAME_EASYBGEE)\njdbc.driver=com.mysql.jdbc.Driver\njdbc.user=$(DBUSER_EASYBGEE)\njdbc.name=$(DBNAME_EASYBGEE)\njdbc.password=$(DBPASS_EASYBGEE)' >> $(OUTPUT_DIR)bgee_connection.properties
	@touch $@

generate_triples: create_property_file download_and_install_ontop download_ontology_file download_mapping_file
	@bgee_rdf.sh -m $(INPUT_DIR)$(GENEX_MAPPING_FILE) -o $(TTL_DIRECTORY) -p $(OUTPUT_DIR)bgee_connection.properties -x $(ONTOP_INSTALL_DIR) -t $(INPUT_DIR)$(GENEX_ONTOLOGY_FILE) -v '$(FTP_RELEASE)' -j '-Xmx128G'" > $@.tmp 2> $@.err
	$(MV) $@.tmp $@


TTL_FILE_LIST      := $(TTL_DIRECTORY)ttl_files.txt
BGEE_VIRTUOSO_FILE := $(TTL_DIRECTORY)add_bgee_virtuoso.sql"

insert_to_triplestore: generate_triples
	@echo "Generating iSQL script to load RDF data into a virtuoso data store..."
	
	ls $(TTL_DIRECTORY)*.ttl > $(TTL_FILE_LIST) ;\
	while IFS= read -r line ; do\
		echo "DB.DBA.TTLP_MT(file_to_string_output('$(TTL_DIRECTORY)$line'), '', 'https://bgee.org/rdf_v$(FTP_RELEASE)');" >>  $(BGEE_VIRTUOSO_FILE) ;\
	done < $(TTL_FILE_LIST) ;\

	@echo "SPARQL LOAD <http://purl.org/genex#> INTO <https://bgee.org/rdf_v$(FTP_RELEASE)>;" >>  $(BGEE_VIRTUOSO_FILE)
	@$(RM) $(TTL_FILE_LIST)

	@echo "Executing iSQL script to load RDF data into a virtuoso data store..."
	@$(ISQL_PATH) $(VIRTUOSO_HOST) $(VIRTUSER) $(VIRTPASS) "EXEC=LOAD $(BGEE_VIRTUOSO_FILE)" > $@.tmp 2> $@.err

	@echo "Finished. All files were loaded into a virtuoso data store."
	$(MV) $@.tmp $@

zip_and_store: insert_to_triplestore
	@zip $(BGEE_VIRTUOSO_FILE) $(BGEE_VIRTUOSO_FILE).zip > $@.tmp $@.err
	@scp $(BGEE_VIRTUOSO_FILE).zip $(FTPUSER)@$(FTPHOST):$(FTP_SAVE_DUMP_PATH) > $@.tmp 2> $@.err
	@echo "RDF file was stored on FTP. Please manually remove local version with the command : \nrm -f $(BGEE_VIRTUOSO_FILE) $(BGEE_VIRTUOSO_FILE).zip"
	$(MV) $@.tmp $@

$(VERIFICATIONFILE): zip_and_store
	@touch $@

clean:
	-@$(RM) -R $(VERIFICATIONFILE) insert_to_triplestore generate_triples create_property_file download_mapping_file download_ontology_file download_and_install_ontop
