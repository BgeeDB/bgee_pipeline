PIPELINEROOT := ../
DIR_NAME := db_creation/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)

$(VERIFICATIONFILE): bgeeSchema.sql bgeeConstraint.sql bgeeIndex.sql bgeeForeignKey.sql insert_data_sources.sql
# create the database and insert data sources
	@$(MYSQLNODBNAME) -e "CREATE DATABASE $(DBNAME)"
	@$(MYSQL) < bgeeSchema.sql
	@$(MYSQL) < insert_data_sources.sql
	@$(MYSQL) < bgeeConstraint.sql
	@$(MYSQL) < bgeeIndex.sql
	@$(MYSQL) < bgeeForeignKey.sql
# get the number of tables and of data sources inserted as a verification
	@$(MYSQL) -e "SELECT (SELECT COUNT(*) FROM information_schema.tables WHERE table_schema = '$(DBNAME)' and TABLE_TYPE='BASE TABLE') as \`Table count\`, (SELECT COUNT(*) FROM dataSource) as \`Data source count\`" > $@.temp
	@$(MV) $@.temp $@

# As there is really a file named update_data_sources.sql in this directory, this target will never be executed,
# unless explicitly called. This is the expected behavior, as this file should be used at the very end of the pipeline run.
# this step also insert data into table dataSourceToSpecies
update_data_sources.sql: $(VERIFICATIONFILE)
	@$(MYSQL) < $@
	@$(MYSQL) < insert_data_sources_to_species.sql
	@echo "----------Looking for species with no data sources: " > $<.temp
	@$(MYSQL) -e "SELECT t1.speciesId AS 'ERROR: SPECIES IDS WITH NO DATA SOURCES' FROM species AS t1 LEFT OUTER JOIN dataSourceToSpecies AS t2 ON t1.speciesId = t2.speciesId WHERE t2.speciesId IS NULL" >> $<.temp
	@echo "----------Looking for species with no data types defined through data sources: " >> $<.temp
	@$(MYSQL) -e "SELECT t1.speciesId AS 'ERROR: SPECIES IDS WITH NO DATA TYPES' FROM species AS t1 LEFT OUTER JOIN dataSourceToSpecies AS t2 ON t1.speciesId = t2.speciesId AND t2.infoType = 'data' WHERE t2.speciesId IS NULL" >> $<.temp
	@$(MYSQL) -e "SELECT * FROM dataSource ORDER BY displayOrder" >> $<.temp
	@$(MYSQL) -e "SELECT * FROM dataSourceToSpecies" >> $<.temp
	@$(CAT) $<.temp >> $<
	@$(RM) $<.temp

clean:
	-@$(RM) $(VERIFICATIONFILE).temp $(VERIFICATIONFILE)

.PHONY := $(.PHONY) dropDatabaseBgee$(RELEASE)
dropDatabaseBgee$(RELEASE):
	@$(MYSQL) -e "DROP DATABASE $(DBNAME)"
