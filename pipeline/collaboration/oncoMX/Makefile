PIPELINEROOT := ../../
DIR_NAME := collaboration/oncoMX/
include $(PIPELINEROOT)Makefile.common

all : $(VERIFICATIONFILE)

oncoMX :
	perl generate_oncoMX_expression_files.pl -bgee=$(BGEECMD) -speciesArg=10090 -devStageArg=UBERON:0000113 -outputDir=$(OUTPUT_DIR) -bgeeVersion=$(DBNAME) 1 > $@.tmp
	@$(MV) $@.tmp $@



# ======================
# Verification file
# ======================
$(VERIFICATIONFILE): oncoMX
	@touch $(VERIFICATIONFILE)

clean :
	-@$(RM) oncoMX $(VERIFICATIONFILE)
