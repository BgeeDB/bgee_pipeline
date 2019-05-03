PIPELINEROOT := ../../
DIR_NAME := collaboration/genofish/
include $(PIPELINEROOT)Makefile.common

all : $(VERIFICATIONFILE)

$(OUTPUT_DIR)homology :
	# Launch file generation
	@$(JAVA) GenoFishProject homology 7955,9606,10090,9913,10116,9544,9593,9598 UBERON:0000955,UBERON:0000922,UBERON:0011150,UBERON:0000948,UBERON:0001474,UBERON:0000160,UBERON:0002113,UBERON:0002107,UBERON:0001630,UBERON:0000992,UBERON:0000473 $(OUTPUT_DIR) > $@.tmp
	@$(MV) $@.tmp $@

# ======================
# Verification file
# ======================
$(VERIFICATIONFILE): $(OUTPUT_DIR)homology
	@touch $(VERIFICATIONFILE)

clean :
	-@$(RM) $(OUTPUT_DIR)homology $(VERIFICATIONFILE)