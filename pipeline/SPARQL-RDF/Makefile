PIPELINEROOT := ../
DIR_NAME := SPARQL-RDF/
include $(PIPELINEROOT)Makefile.common

FTP_RELEASE   := $(subst .,_,$(RELEASE))
ONTOP_VERSION := '3.0.0-beta-2'
TOMCAT_FOLDER := '/var/lib/tomcat8'

all: $(VERIFICATIONFILE)

get_bgeeLite_dump:
	# Download bgee lite dump
	@$(WGET) "ftp://ftp.bgee.org/bgee_v${FTP_RELEASE}/sql_lite_dump.tar.gz"  && $(MV) sql_lite_dump.tar.gz /tmp/  || rm -f sql_lite_dump.tar.gz
#Do not store sql_lite_dump.tar.gz in source_files/ because big and in our FTP
	@tar xvfz /tmp/sql_lite_dump.tar.gz && rm -f /tmp/sql_lite_dump.tar.gz
	@touch $@

create_bgeeLite_db: get_bgeeLite_dump
	# Create a database with bgee lite dump file
	@$(MYSQLNODBNAME_BGEELITE) -e "CREATE DATABASE $(DBNAME_BGEELITE)" > $@.tmp 2> $@.warnings
	@$(MYSQL_BGEELITE) < sql_lite_dump.sql >> $@.tmp 2>> $@.warnings
	@if [[ ! -s $@.warnings ]]; then $(RM) $@.warnings; fi
	# Remove dump file
	@$(RM) sql_lite_dump.sql
	@$(MV) $@.tmp $@


create_property_file:
	# Create a property file for server connection to database
	@echo -e "jdbc.url=jdbc:mysql://$(DBHOST_BGEELITE):$(DBPORT)/$(DBNAME_BGEELITE)\njdbc.driver=com.mysql.jdbc.Driver\njdbc.user=$(DBUSER_BGEELITE)\njdbc.name=$(DBNAME_BGEELITE)\njdbc.password=$(DBPASS_BGEELITE)" > bgee_connection.properties
	@touch $@

download_mapping_file:
	# Download mapping file
	@$(WGET) 'https://raw.githubusercontent.com/biosoda/bioquery/master/Bgee_OBDA_mappings/bgee_genex.obda'  && $(MV) bgee_genex.obda $(SOURCE_FILES_DIR)/$(DIR_NAME)/  || rm -f bgee_genex.obda
	@touch $@

download_ontology_file:
	# Download OWL ontology file
	@$(WGET) 'https://raw.githubusercontent.com/biosoda/genex/master/genex.owl'  && $(MV) genex.owl $(SOURCE_FILES_DIR)/$(DIR_NAME)/  || rm -f genex.owl
	@touch $@

download_ontop_webapps:
	# Download ontop webapps
	@$(WGET) 'https://sourceforge.net/projects/ontop4obda/files/ontop-$(ONTOP_VERSION)/ontop-webapps-$(ONTOP_VERSION).zip/download' -O ontop-webapps-$(ONTOP_VERSION).zip  && $(MV) ontop-webapps-$(ONTOP_VERSION).zip $(SOURCE_FILES_DIR)/$(DIR_NAME)/  || rm -f ontop-webapps-$(ONTOP_VERSION).zip
	@touch $@

install_RDF_webapp: create_property_file download_mapping_file download_ontology_file download_ontop_webapps
	@echo -e 'Copy extracted files from $(SOURCE_FILES_DIR)/$(DIR_NAME)ontop-webapps-$(ONTOP_VERSION).zip\n\tin the tomcat webapp folder: **$(TOMCAT_FOLDER)/webapps/**\n'
	@echo -e 'Configure tomcat proxy through apache accordingly\n\t# Ontop settings\n\tJkMount /rdf4j-server      default\n\tJkMount /rdf4j-server/*    default\n\tJkMount /rdf4j-workbench   default\n\tJkMount /rdf4j-workbench/* default\n'
	@echo -e 'Restart apache service\n'
	@echo -e 'Copy bgee_connection.properties, $(SOURCE_FILES_DIR)/$(DIR_NAME)genex.owl and $(SOURCE_FILES_DIR)/$(DIR_NAME)bgee_genex.obda files\n\tin **$(TOMCAT_FOLDER)/bgee_webapps/$(DBNAME_BGEELITE)_conf/**\n'
	@echo -e 'Browse to <YOUR SERVER URL>/rdf4j-workbench\n'
	@echo -e 'On the <YOUR SERVER URL>/rdf4j-workbench web interface click on\n\t**SYSTEM** link, then on\n\t**Repositories -> New repository**\n'
	@echo -e 'IF IT DOES NOT WORK!\nCheck in tomcat log for permission issues regarding\n\tlog folder, .RDF4J/ folder, ...\n'
	@echo -e 'Create a new repository\n\tType: \tOntop Virtual RDF Store\n\tID: \t$(DBNAME_BGEELITE)\n\tTitle: \t$(DBNAME_BGEELITE)\n'
	@echo -e 'Click on **Next**\n'
	@echo -e 'Fill the input text fields in the form according to the directory path $(TOMCAT_FOLDER)/bgee_webapps/$(DBNAME_BGEELITE)_conf/\n\tOWL ontology file: \t$(TOMCAT_FOLDER)/bgee_webapps/$(DBNAME_BGEELITE)_conf/genex.owl\n\tMapping file: \t\t$(TOMCAT_FOLDER)/bgee_webapps/$(DBNAME_BGEELITE)_conf/bgee_genex.obda\n\tOntop properties file: \t$(TOMCAT_FOLDER)/bgee_webapps/$(DBNAME_BGEELITE)_conf/bgee_connection.properties\n'
	@echo -e 'Click on **Create**\n'
	@echo -e '... \n'
	@echo -e 'rdf4j-workbench MUST be accessible ONLY internally for admin tasks!\nrdf4j-server    MUST be accessible to all RDF/SPARQL users!\n'


drop_bgeeLite_db:
	# Drop $(DBNAME_BGEELITE) database in $(DBHOST_BGEELITE)
	@$(MYSQLNODBNAME_BGEELITE) -e "DROP DATABASE $(DBNAME_BGEELITE)" > $@.tmp 2> $@.warnings
	@if [[ ! -s $@.warnings ]]; then $(RM) $@.warnings; fi
	@$(MV) $@.tmp $@


$(VERIFICATIONFILE): create_bgeeLite_db install_RDF_webapp
	@touch $@

clean:
	-@$(RM) -R $(VERIFICATIONFILE) get_bgeeLite_dump create_bgeeLite_db create_property_file download_mapping_file download_ontology_file download_ontop_webapps install_RDF_webapp drop_bgeeLite_db
