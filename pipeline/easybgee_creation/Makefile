PIPELINEROOT := ../
DIR_NAME := easybgee_creation/
include $(PIPELINEROOT)Makefile.common

# compressed dump of the Easy Bgee database
EASYBGEE_DUMP_FILE := $(OUTPUT_DIR)$(DBNAME_EASYBGEE).sql.gz

all: $(OUTPUT_DIR)extract_data_from_bgee $(OUTPUT_DIR)import_to_easybgee $(OUTPUT_DIR)dump_easybgee

$(OUTPUT_DIR)create_schema: easyBgeeSchema.sql
	@echo "create the database $(DBNAME_EASYBGEE) in $(DBHOST_EASYBGEE)" > $@.tmp 2> $@.warnings
	@$(MYSQLNODBNAME_EASYBGEE) -e "CREATE DATABASE IF NOT EXISTS $(DBNAME_EASYBGEE)" >> $@.tmp 2>> $@.warnings
	@echo "Load the database schema for $(DBNAME_EASYBGEE) in $(DBHOST_EASYBGEE)" >> $@.tmp 2>> $@.warnings
	@$(MYSQL_EASYBGEE) < easyBgeeSchema.sql >> $@.tmp 2>> $@.warnings
	@echo "get the number of tables as a verification" >> $@.tmp 2>> $@.warnings
	@$(MYSQL_EASYBGEE) -e "SELECT (SELECT COUNT(*) FROM information_schema.tables WHERE table_schema = '$(DBNAME_EASYBGEE)' and TABLE_TYPE='BASE TABLE') as \`Table count\`" >> $@.tmp 2>> $@.warnings
	@$(MV) $@.tmp $@

$(OUTPUT_DIR)extract_data_from_bgee:
	@echo "extract data from the $(DBNAME) database in $(DBHOST)"
	@$(JAVA_NO_JAR) -Djava.util.concurrent.ForkJoinPool.common.parallelism=10 -jar $(PIPELINEROOT)$(BGEEPIPELINEJAR) GenerateEasyBgee extractFromBgee $(EASYBGEE_PATH) - true > $@.tmp 2>> $@.warnings
	@$(MV) $@.tmp $@

$(OUTPUT_DIR)import_to_easybgee: $(OUTPUT_DIR)create_schema $(OUTPUT_DIR)extract_data_from_bgee
	@echo "import data to the $(DBNAME_EASYBGEE) database in $(DBHOST_EASYBGEE)"
	@$(JAVA_EASYBGEE) GenerateEasyBgee tsvToEasyBgee $(EASYBGEE_PATH) 6 > $@.tmp 2> $@.warnings
	@$(MV) $@.tmp $@

# dump the Easy Bgee database and compress it on the fly.
# pipefail: make the rule fail if mysqldump fails, not only if gzip fails.
$(OUTPUT_DIR)dump_easybgee: $(OUTPUT_DIR)import_to_easybgee
	@echo "dump $(DBNAME_EASYBGEE) from $(DBHOST_EASYBGEE) to $(EASYBGEE_DUMP_FILE)" > $@.tmp 2> $@.warnings
	@set -o pipefail; mysqldump -u $(DBUSER_EASYBGEE) -p$(DBPASS_EASYBGEE) -h $(DBHOST_EASYBGEE) -P $(DBPORT) --single-transaction --routines --triggers $(DBNAME_EASYBGEE) 2>> $@.warnings | gzip -9 > $(EASYBGEE_DUMP_FILE)
	@echo "check the integrity of the compressed dump" >> $@.tmp 2>> $@.warnings
	@gzip -t $(EASYBGEE_DUMP_FILE) >> $@.tmp 2>> $@.warnings
	@ls -lh $(EASYBGEE_DUMP_FILE) >> $@.tmp 2>> $@.warnings
	@$(MV) $@.tmp $@

clean:
	-@$(RM) $(OUTPUT_DIR)create_schema $(OUTPUT_DIR)create_schema.tmp $(OUTPUT_DIR)create_schema.warnings $(OUTPUT_DIR)extract_data_from_bgee $(OUTPUT_DIR)extract_data_from_bgee.tmp $(OUTPUT_DIR)extract_data_from_bgee.warnings $(OUTPUT_DIR)import_to_easybgee $(OUTPUT_DIR)import_to_easybgee.tmp $(OUTPUT_DIR)import_to_easybgee.warnings $(OUTPUT_DIR)dump_easybgee $(OUTPUT_DIR)dump_easybgee.tmp $(OUTPUT_DIR)dump_easybgee.warnings $(EASYBGEE_DUMP_FILE)

# also remove the marker files of the steps that fill the database, otherwise make
# would consider them up to date after the database has been dropped
drop_easybgee:
	@echo "drop $(DBNAME_EASYBGEE) database in $(DBHOST_EASYBGEE)"
	@$(MYSQLNODBNAME_EASYBGEE) -e "DROP DATABASE IF EXISTS $(DBNAME_EASYBGEE)"
	@$(RM) $(OUTPUT_DIR)create_schema $(OUTPUT_DIR)import_to_easybgee $(OUTPUT_DIR)dump_easybgee

.PHONY := $(.PHONY) drop_easybgee
