PIPELINEROOT := ../
DIR_NAME := dblight_creation/
include $(PIPELINEROOT)Makefile.common



all: $(OUTPUT_DIR)import_to_bgeelight

$(OUTPUT_DIR)create_schema: bgeeLightSchema.sql
	# create the database $(DBNAME_BGEELIGHT) in $(DBHOST_BGEELIGHT)
	@$(MYSQLNODBNAME_BGEELIGHT) -e "CREATE DATABASE $(DBNAME_BGEELIGHT)" > $@.tmp 2>> $@.warnings
	@$(MYSQL_BGEELIGHT) < bgeeLightSchema.sql > $@.tmp 2>> $@.warnings
	# get the number of tables as a verification
	@$(MYSQL_BGEELIGHT) -e "SELECT (SELECT COUNT(*) FROM information_schema.tables WHERE table_schema = '$(DBNAME_BGEELIGHT)' and TABLE_TYPE='BASE TABLE') as \`Table count\`" > $@.tmp 2>> $@.warnings
	@$(MV) $@.tmp $@

# As there is really a file named update_data_sources.sql in this directory, this target will never be executed,
# unless explicitly called. This is the expected behavior, as this file should be used at the very end of the pipeline run.
# this step also insert data into table dataSourceToSpecies
$(OUTPUT_DIR)extract_data_from_bgee:
	# extract data from the $(DBNAME) database in $(DBHOST)
	@$(JAVA) GenerateBgeeLight extractFromBgee $(OUTPUT_DIR) - > $@.tmp 2>> $@.warnings
	@$(MV) $@.tmp $@
	
$(OUTPUT_DIR)import_to_bgeelight: $(OUTPUT_DIR)create_schema
	# import data to the $(DBNAME_BGEELIGHT) database in $(DBHOST_BGEELIGHT)
	@$(JAVA_BGEE_LIGHT) GenerateBgeeLight tsvToBgeeLight $(OUTPUT_DIR) > $@.tmp 2>> $@.warnings
	@$(MV) $@.tmp $@

clean:
	-@$(RM) $(OUTPUT_DIR)create_schema $(OUTPUT_DIR)create_schema.tmp $(OUTPUT_DIR)extract_data_from_bgee $(OUTPUT_DIR)extract_data_from_bgee.tmp $(OUTPUT_DIR)import_to_bgeelight $(OUTPUT_DIR)import_to_bgeelight.tmp drop_bgee_light drop_bgee_light.tmp
	
drop_bgee_light:
	# drop $(DBNAME_BGEELIGHT) database in $(DBHOST_BGEELIGHT)
	@$(MYSQLNODBNAME_BGEELIGHT) -e "DROP DATABASE $(DBNAME_BGEELIGHT)" > $@.tmp 2>> $@.warnings
	@$(MV) $@.tmp $@

.PHONY := $(.PHONY) dropDatabaseBgee$(RELEASE)
dropDatabaseBgeeLight$(RELEASE):
	@$(MYSQLLIGHT) -e "DROP DATABASE $(DBLIGHTNAME)"
