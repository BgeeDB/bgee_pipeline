PIPELINEROOT := ../
DIR_NAME := easybgee_creation/
include $(PIPELINEROOT)Makefile.common



all: $(OUTPUT_DIR)extract_data_from_bgee $(OUTPUT_DIR)import_to_easybgee

$(OUTPUT_DIR)create_schema: easyBgeeSchema.sql
	@echo "create the database $(DBNAME_EASYBGEE) in $(DBHOST_EASYBGEE)" > $@.tmp 2> $@.warnings
	@$(MYSQLNODBNAME_EASYBGEE) -e "CREATE DATABASE $(DBNAME_EASYBGEE)" >> $@.tmp 2>> $@.warnings
	@echo "Load the database schema for $(DBNAME_EASYBGEE) in $(DBHOST_EASYBGEE)" >> $@.tmp 2>> $@.warnings
	@$(MYSQL_EASYBGEE) < easyBgeeSchema.sql >> $@.tmp 2>> $@.warnings
	@echo "get the number of tables as a verification" >> $@.tmp 2>> $@.warnings
	@$(MYSQL_EASYBGEE) -e "SELECT (SELECT COUNT(*) FROM information_schema.tables WHERE table_schema = '$(DBNAME_EASYBGEE)' and TABLE_TYPE='BASE TABLE') as \`Table count\`" >> $@.tmp 2>> $@.warnings
	@$(MV) $@.tmp $@

$(OUTPUT_DIR)extract_data_from_bgee:
	# extract data from the $(DBNAME) database in $(DBHOST)
	@$(JAVA_NO_JAR) -Djava.util.concurrent.ForkJoinPool.common.parallelism=4 -jar $(PIPELINEROOT)$(BGEEPIPELINEJAR) GenerateEasyBgee extractFromBgee $(EASYBGEE_PATH)  - > $@.tmp 2>> $@.warnings
	@$(MV) $@.tmp $@

$(OUTPUT_DIR)import_to_easybgee: $(OUTPUT_DIR)create_schema
	# import data to the $(DBNAME_EASYBGEE) database in $(DBHOST_EASYBGEE)
	@$(JAVA_EASYBGEE) GenerateEasyBgee tsvToEasyBgee $(EASYBGEE_PATH) > $@.tmp 2> $@.warnings
	@$(MV) $@.tmp $@

clean:
	-@$(RM) $(OUTPUT_DIR)create_schema $(OUTPUT_DIR)create_schema.tmp $(OUTPUT_DIR)extract_data_from_bgee $(OUTPUT_DIR)extract_data_from_bgee.tmp $(OUTPUT_DIR)import_to_easybgee $(OUTPUT_DIR)import_to_easybgee.tmp drop_easybgee drop_easybgee.tmp

drop_easybgee:
	# drop $(DBNAME_EASYBGEE) database in $(DBHOST_EASYBGEE)
	@$(MYSQLNODBNAME_EASYBGEE) -e "DROP DATABASE $(DBNAME_EASYBGEE)" > $@.tmp 2> $@.warnings
	@$(MV) $@.tmp $@

.PHONY := $(.PHONY) dropDatabaseBgee$(RELEASE)
