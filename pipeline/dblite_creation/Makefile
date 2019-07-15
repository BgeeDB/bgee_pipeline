PIPELINEROOT := ../
DIR_NAME := dblite_creation/
include $(PIPELINEROOT)Makefile.common



all: $(OUTPUT_DIR)extract_data_from_bgee $(OUTPUT_DIR)import_to_bgeelite

$(OUTPUT_DIR)create_schema: bgeeLiteSchema.sql
	# create the database $(DBNAME_BGEELITE) in $(DBHOST_BGEELITE)
	@$(MYSQLNODBNAME_BGEELITE) -e "CREATE DATABASE $(DBNAME_BGEELITE)" > $@.tmp 2> $@.warnings
	@$(MYSQL_BGEELITE) < bgeeLiteSchema.sql >> $@.tmp 2>> $@.warnings
	# get the number of tables as a verification
	@$(MYSQL_BGEELITE) -e "SELECT (SELECT COUNT(*) FROM information_schema.tables WHERE table_schema = '$(DBNAME_BGEELITE)' and TABLE_TYPE='BASE TABLE') as \`Table count\`" >> $@.tmp 2>> $@.warnings
	@$(MV) $@.tmp $@

$(OUTPUT_DIR)extract_data_from_bgee:
	# extract data from the $(DBNAME) database in $(DBHOST)
	@$(JAVA) GenerateBgeeLite extractFromBgee $(OUTPUT_DIR) - > $@.tmp 2> $@.warnings
	@$(MV) $@.tmp $@

$(OUTPUT_DIR)import_to_bgeelite: $(OUTPUT_DIR)create_schema
	# import data to the $(DBNAME_BGEELITE) database in $(DBHOST_BGEELITE)
	@$(JAVA_BGEE_LITE) GenerateBgeeLite tsvToBgeeLite $(OUTPUT_DIR) > $@.tmp 2> $@.warnings
	@$(MV) $@.tmp $@

clean:
	-@$(RM) $(OUTPUT_DIR)create_schema $(OUTPUT_DIR)create_schema.tmp $(OUTPUT_DIR)extract_data_from_bgee $(OUTPUT_DIR)extract_data_from_bgee.tmp $(OUTPUT_DIR)import_to_bgeelite $(OUTPUT_DIR)import_to_bgeelite.tmp drop_bgee_lite drop_bgee_lite.tmp

drop_bgee_lite:
	# drop $(DBNAME_BGEELITE) database in $(DBHOST_BGEELITE)
	@$(MYSQLNODBNAME_BGEELITE) -e "DROP DATABASE $(DBNAME_BGEELITE)" > $@.tmp 2> $@.warnings
	@$(MV) $@.tmp $@

.PHONY := $(.PHONY) dropDatabaseBgee$(RELEASE)
