PIPELINEROOT := ../
DIR_NAME := annotations/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)
HOMONTOLOGY := homology_ontology.obo
ECOONTOLOGY := eco.obo
CONFONTOLOGY := conf_information.obo

$(OUTPUT_DIR)remap_cond: $(INPUT_DIR)condition_remapping.tsv $(UBERON_SEX_INFO_FILE_PATH)
	@perl remap_conditions.pl -bgee=$(BGEECMD) -remapping_file=$(INPUT_DIR)condition_remapping.tsv -sex_info=$(UBERON_SEX_INFO_FILE_PATH) > $@.tmp 2>$@.warnings
	@$(MV) $@.tmp $@

$(SIMILARITYCLEANFILE): $(SIMILARITYFILE) $(TAXONCONSTRAINTSFILEPATH) $(UBERONFILEPATH) $(TAXONOMYFILEPATH) $(HOMONTOLOGY) $(ECOONTOLOGY) $(CONFONTOLOGY)
	@$(JAVA) SimilarityAnnotation generateReleaseFile $^ $@.tmp
	@$(MV) $@.tmp $@

clean:
	-@$(RM) $(SIMILARITYCLEANFILE).tmp $(SIMILARITYCLEANFILE) $(OUTPUT_DIR)remap_cond.tmp $(OUTPUT_DIR)remap_cond $(OUTPUT_DIR)remap_cond.warnings
