PIPELINEROOT := ../
DIR_NAME := annotations/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)
HOMONTOLOGY := homology_ontology.obo
ECOONTOLOGY := eco.obo
CONFONTOLOGY := conf_information.obo

# Condition remapping
$(OUTPUT_DIR)remap_cond: $(INPUT_DIR)condition_remapping.tsv $(UBERON_SEX_INFO_FILE_PATH)
# First, we create or identify the remapped conditions
# thanks to the script remap_conditions.pl
	@$(MYSQL) -e "DELETE FROM remapCond" > $@.tmp 2>$@.warnings
	@perl remap_conditions.pl -bgee=$(BGEECMD) -remapping_file=$(INPUT_DIR)condition_remapping.tsv -sex_info=$(UBERON_SEX_INFO_FILE_PATH) > $@.tmp 2>$@.warnings
# Now, we update the tables using the SQL script remap_condition.sql
# WARNING: check the script remap_conditions.sql before using it,
# in case new tables using an expressionId or a conditinId have been added
	@$(MYSQL) < remap_conditions.sql > $@.tmp 2>$@.warnings
# Then it is your responsibility to delete rows in the table remapCond,
# no big deal if you don't. We don't do it here in case you wanna check things.
# WARNING: then it is necessary to remove the unused expression calls,
# then the unused conditions (see below)
	@$(MV) $@.tmp $@

# WARNING: check the script delete_unused_expressions.sql before using it,
# in case new tables using an expressionId have been added
$(OUTPUT_DIR)delete_unused_expressions:
	@$(MYSQL) < delete_unused_expressions.sql > $@.tmp 2>$@.warnings
	@$(MV) $@.tmp $@

# WARNING: check the script delete_unused_conditions.sql before using it,
# in case new tables using a conditionId have been added
$(OUTPUT_DIR)delete_unused_conditions: $(OUTPUT_DIR)delete_unused_expressions
	@$(MYSQL) < delete_unused_conditions.sql > $@.tmp 2>$@.warnings
	@$(MV) $@.tmp $@

$(SIMILARITYCLEANFILE): $(SIMILARITYFILE) $(TAXONCONSTRAINTSFILEPATH) $(UBERONFILEPATH) $(TAXONOMYFILEPATH) $(HOMONTOLOGY) $(ECOONTOLOGY) $(CONFONTOLOGY)
	@$(JAVA) SimilarityAnnotation generateReleaseFile $^ $@.tmp
	@$(MV) $@.tmp $@

clean:
	-@$(RM) $(SIMILARITYCLEANFILE).tmp $(SIMILARITYCLEANFILE) $(OUTPUT_DIR)remap_cond.tmp $(OUTPUT_DIR)remap_cond $(OUTPUT_DIR)remap_cond.warnings
