PIPELINEROOT := ../
DIR_NAME := SPARQL-RDF/
include $(PIPELINEROOT)Makefile.common

FTP_RELEASE := `echo $(RELEASE)| sed -e 's@\.@_@'`


all: $(VERIFICATIONFILE)

get_bgeeLite_dump:
	# Download bgee lite dump
	@$(WGET) "ftp://ftp.bgee.org/bgee_v${FTP_RELEASE}/sql_lite_dump.tar.gz" && $(MV) sql_lite_dump.tar.gz /tmp/ || rm -f sql_lite_dump.tar.gz
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


drop_bgeeLite_db:
	# Drop $(DBNAME_BGEELITE) database in $(DBHOST_BGEELITE)
	@$(MYSQLNODBNAME_BGEELITE) -e "DROP DATABASE $(DBNAME_BGEELITE)" > $@.tmp 2> $@.warnings
	@if [[ ! -s $@.warnings ]]; then $(RM) $@.warnings; fi
	@$(MV) $@.tmp $@


$(VERIFICATIONFILE): create_bgeeLite_db
	@touch $@

clean:
	-@$(RM) -R $(VERIFICATIONFILE) get_bgeeLite_dump create_bgeeLite_db drop_bgeeLite_db
