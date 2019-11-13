PIPELINEROOT := ../../
DIR_NAME := collaboration/oncoMX/
include $(PIPELINEROOT)Makefile.common

all : $(VERIFICATIONFILE)

$(OUTPUT_DIR)generateFiles :
	# Launch file generation
	@$(JAVA) GenerateOncoMXFile $(INPUT_DIR)human_doid_slim_uberon_mapping.csv $(OUTPUT_DIR) 9606//UBERON:0000113,10090//UBERON:0000113 > $@.tmp
	# compress generated files
	@find $(OUTPUT_DIR) -type f -not -name '*.gz' -not -name '*.tmp' -name '*.tsv' -execdir gzip --verbose --best {} \;
	@$(MV) $@.tmp $@

# ======================
# Verification file
# ======================
$(VERIFICATIONFILE): $(OUTPUT_DIR)generateFiles
	@touch $(VERIFICATIONFILE)

clean :
	-@$(RM) $(OUTPUT_DIR)generateFiles $(VERIFICATIONFILE)
