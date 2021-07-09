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

# Actually this should be part of the similarity annotation repository,
# not here.
$(SIMILARITYCLEANFILE): $(SIMILARITYFILE) $(TAXONCONSTRAINTSFILEPATH) $(UBERONFILEPATH) $(TAXONOMYFILEPATH) $(HOMONTOLOGY) $(ECOONTOLOGY) $(CONFONTOLOGY)
	@$(JAVA) SimilarityAnnotation generateReleaseFile $^ $@.tmp
	@$(MV) $@.tmp $@

SIM_ANNOT_DOWNLOAD_URL_START := https://raw.githubusercontent.com/BgeeDB/anatomical-similarity-annotations/master/release/

$(RAW_SIMILARITY_FILE):
	@$(CURL) $(SIM_ANNOT_DOWNLOAD_URL_START)raw_similarity_annotations.tsv $(APPEND_CURL_COMMAND)
$(SUMMARY_SIMILARITY_FILE):
	@$(CURL) $(SIM_ANNOT_DOWNLOAD_URL_START)summary_similarity_annotations.tsv $(APPEND_CURL_COMMAND)

# TODO: manage insertion of ECO and CIO ontologies before this step
$(OUTPUT_DIR)insert_sim_annots: $(RAW_SIMILARITY_FILE) $(SUMMARY_SIMILARITY_FILE)
	@$(JAVA) InsertSimilarityAnnotation $(RAW_SIMILARITY_FILE) $(SUMMARY_SIMILARITY_FILE) > $@.tmp 2> $@.warn
# We need to add the root of all cell types to the similarity mapping
	@$(MYSQL) -e "INSERT INTO summarySimilarityAnnotation (taxonId, negated, CIOId) VALUES (1, 0, 'CIO:0000003'); INSERT INTO similarityAnnotationToAnatEntityId (summarySimilarityAnnotationId, anatEntityId) VALUES ((SELECT MAX(summarySimilarityAnnotationId) FROM summarySimilarityAnnotation), '$(CELL_TYPE_ROOT_ID)'); INSERT INTO rawSimilarityAnnotation (summarySimilarityAnnotationId, negated, ECOId, CIOId, referenceId, referenceTitle, supportingText, assignedBy, curator, annotationDate) VALUES ((SELECT MAX(summarySimilarityAnnotationId) FROM summarySimilarityAnnotation), 0, 'ECO:0000306', 'CIO:0000003', '', '', '', 'Bgee', 'FBB', '2021-06-14')" >> $@.tmp 2>> $@.warn
	@$(MV) $@.tmp $@

clean:
	-@$(RM) $(SIMILARITYCLEANFILE).tmp $(SIMILARITYCLEANFILE) $(OUTPUT_DIR)remap_cond.tmp $(OUTPUT_DIR)remap_cond $(OUTPUT_DIR)remap_cond.warnings
