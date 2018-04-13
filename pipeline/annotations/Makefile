PIPELINEROOT := ../
DIR_NAME := annotations/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)
HOMONTOLOGY := homology_ontology.obo
ECOONTOLOGY := eco.obo
CONFONTOLOGY := conf_information.obo


$(SIMILARITYCLEANFILE): $(SIMILARITYFILE) $(TAXONCONSTRAINTSFILEPATH) $(UBERONFILEPATH) $(TAXONOMYFILEPATH) $(HOMONTOLOGY) $(ECOONTOLOGY) $(CONFONTOLOGY)
	@$(JAVA) SimilarityAnnotation generateReleaseFile $^ $@.tmp
	@$(MV) $@.tmp $@

clean:
	-@$(RM) $(SIMILARITYCLEANFILE).tmp $(SIMILARITYCLEANFILE)
