PIPELINEROOT := ../
DIR_NAME := species/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)

# Generate a file with only the taxon IDs present in the similarity annotation file.
$(ANNOTTAXIDSFILEPATH):
# download similarity file in TMP dir
	@$(CURL) $(RAW_SIM_ANNOT_DOWNLOAD_URL) -z $(OUTPUT_DIR)similarity_annotations.tsv -o $(OUTPUT_DIR)similarity_annotations.tsv.tmp && $(MV) $(OUTPUT_DIR)similarity_annotations.tsv.tmp $(OUTPUT_DIR)similarity_annotations.tsv 2>/dev/null || rm -f $(OUTPUT_DIR)similarity_annotations.tsv.tmp
	@$(JAVA) SimilarityAnnotation extractTaxonIds $(OUTPUT_DIR)similarity_annotations.tsv $@.tmp
# need to add a header
	@echo "taxon ID" | cat - $@.tmp > $@
	@$(RM) $@.tmp

#Generate a file containing all species IDs from Ensembl and ensmebl Metazoa and some RefSeq species.
#This file will also contain information about species from NCBI listed in the
# input file $(NCBI_SPECIES_FILEPATH) that has to be manually updated before each
# major Bgee release
$(OUTPUT_DIR)generate_species_file: $(ANNOTTAXIDSFILEPATH) $(REFSEQ_SPECIES_FILEPATH)
	@perl get_ensembl_species_info.pl -ensembl_species_info_file=$(BGEE_SPECIES_FILEPATH).tmp > $@.tmp 2> $@.warnings
# remove first line (header) of the file listing refseq species
	@sed 1d $(REFSEQ_SPECIES_FILEPATH) > $(REFSEQ_SPECIES_FILEPATH).noheader
# merge ensembl and refseq species info files
	@cat $(BGEE_SPECIES_FILEPATH).tmp $(REFSEQ_SPECIES_FILEPATH).noheader > $(BGEE_SPECIES_FILEPATH)
	@rm -f $(BGEE_SPECIES_FILEPATH).tmp $(BGEE_SPECIES_FILEPATH).merged $(REFSEQ_SPECIES_FILEPATH).noheader
# delete warnings file if empty
	@[ -s $@.warnings ] || rm -f $@.warnings
	@mv $@.tmp $@


# generate a file containing the IDs of all species used in Bgee, all taxa
# used in the similarity annotations, all taxa used in Uberon
$(ALLTAXIDSFILEPATH): $(OUTPUT_DIR)generate_species_file $(BGEE_SPECIES_FILEPATH)
# download the composite-metazoan uberon ontology
	@$(CURL) $(COMPOSITE_DOWNLOAD_URL) -z $(GENERATED_FILES_DIR)composite-metazoan.owl -o $(GENERATED_FILES_DIR)composite-metazoan.owl.tmp && $(MV) $(GENERATED_FILES_DIR)composite-metazoan.owl.tmp $(GENERATED_FILES_DIR)composite-metazoan.owl 2>/dev/null || rm -f $(GENERATED_FILES_DIR)composite-metazoan.owl.tmp
	@$(JAVA) Uberon extractTaxonIds $(GENERATED_FILES_DIR)composite-metazoan.owl $(GENERATED_FILES_DIR)uberonTaxIds_composite.tsv.tmp
# get the annotation taxon IDs without the first header line
	@sed 1d $(ANNOTTAXIDSFILEPATH) > $(GENERATED_FILES_DIR)$(ANNOTTAXIDSFILE).noheader
# concatenate these annotation taxa with the Uberon taxa, no header in this file
	@cat $(GENERATED_FILES_DIR)$(ANNOTTAXIDSFILE).noheader $(GENERATED_FILES_DIR)uberonTaxIds_composite.tsv.tmp | sort | uniq > $@.tmp
# now concatenate with the file containing the IDs of species used in Bgee, and add its header
	@cat $(BGEE_SPECIES_FILEPATH) $@.tmp > $@
	@$(RM) $(GENERATED_FILES_DIR)$(ANNOTTAXIDSFILE).noheader $(GENERATED_FILES_DIR)uberonTaxIds_composite.tsv.tmp $@.tmp

# update display order. First all model organism and then order the remaining ones per
# english common name
$(OUTPUT_DIR)update_display_order: $(OUTPUT_DIR)generate_species_file
	@Rscript update_display_order.R bgeeSpeciesFile=\"$(BGEE_SPECIES_FILEPATH)\" outputFile=\"$(BGEE_SPECIES_FILEPATH)\" > $@.tmp 2> $@.err
	@[ -s $@.err ] || rm -f $@.err
	@mv $@.tmp $@


$(OUTPUT_DIR)installRobot: $(OUTPUT_DIR)update_display_order
	@if [ ! -d $(ROBOT_DIR_PATH) ]; then mkdir $(ROBOT_DIR_PATH); fi
	@$(CURL) $(ROBOT_JAR_URL) -z $(ROBOT_JAR_PATH) -o $(ROBOT_JAR_PATH).tmp && $(MV) $(ROBOT_JAR_PATH).tmp $(ROBOT_JAR_PATH) > $@.tmp 2> $@.err || rm -f $(ROBOT_JAR_PATH).tmp
	@$(CURL) $(ROBOT_SH_URL) -z $(ROBOT_SH_PATH) -o $(ROBOT_SH_PATH).tmp && $(MV) $(ROBOT_SH_PATH).tmp $(ROBOT_SH_PATH) > $@.tmp 2> $@.err || rm -f $(ROBOT_SH_PATH).tmp
	@chmod u+x $(ROBOT_SH_PATH)
# delete error file if empty
	@[ -s $@.err ] || rm -f $@.err
	@mv $@.tmp $@

# generate the .owl ncbi taxon ontology using ROBOT and the ncbitaxon github repository
NCBITAXON_OWL_FILE := ncbitaxon.owl

$(OUTPUT_DIR)download_ncbi_taxonomy: $(OUTPUT_DIR)installRobot
	@$(CURL) $(NCBITAXON_OBOFOUNDRY_URL) -z $(OUTPUT_DIR)$(NCBITAXON_OWL_FILE) -o $(OUTPUT_DIR)$(NCBITAXON_OWL_FILE).tmp && $(MV) $(OUTPUT_DIR)$(NCBITAXON_OWL_FILE).tmp $(OUTPUT_DIR)$(NCBITAXON_OWL_FILE) 2>/dev/null || rm -f $(OUTPUT_DIR)$(NCBITAXON_OWL_FILE).tmp
	@touch $@

# THIS RULE IS NOT USED IN THE PIPELINE. IT ALLOWS TO GENERATE THE NCBI TAXONOMY ONTOLOGY.
# IT CAN BE USED IN CASE OF ISSUE WITH THE OWL FILE AVAILABLE ON OBOFOUNDRY AND WILL THEN RELACE
# THE RULE $(OUTPUT_DIR)download_ncbi_taxonomy
$(OUTPUT_DIR)generate_ncbi_taxonomy: $(OUTPUT_DIR)installRobot
# First clone the ncbitaxon github repository
	@if [ ! -d $(NCBITAXON_DIR_PATH) ]; then mkdir -p $(NCBITAXON_DIR_PATH); fi
	@if [ -d $(NCBITAXON_DIR_PATH) ] && [ -z "$$(ls -A $(NCBITAXON_DIR_PATH))" ]; then git clone $(NCBITAXON_GITHUB_URL) $(NCBITAXON_DIR_PATH); fi
# Then update the ROBOT variable from the makefile to the ROBOT_SH_PATH
	@sed -i "s|ROBOT=robot|ROBOT=$(ROBOT_SH_PATH)|" $(NCBITAXON_DIR_PATH)/Makefile
# Generate the owl ontology from the ncbitaxon.owl rule of the ncbitaxon makefile
	@make $(NCBITAXON_OWL_FILE) -C $(NCBITAXON_DIR_PATH)
	@mv $(NCBITAXON_DIR_PATH)/$(NCBITAXON_OWL_FILE) $(OUTPUT_DIR)$(NCBITAXON_OWL_FILE)
	@rm -rf $(NCBITAXON_DIR_PATH)
	@touch $@

$(TAXONOMY_FILE_PATH): $(OUTPUT_DIR)download_ncbi_taxonomy
# keep only first column
	@cut -f1 $(ALLTAXIDSFILEPATH) > $(OUTPUT_DIR)extractTaxonIds.tmp
# remove header
	@sed 1d $(OUTPUT_DIR)extractTaxonIds.tmp > $(OUTPUT_DIR)extractTaxonIds.noheader && rm -f $(OUTPUT_DIR)extractTaxonIds.tmp
# add prefix NCBITaxon:
	@sed 's/^/NCBITaxon:/' $(OUTPUT_DIR)extractTaxonIds.noheader > $(OUTPUT_DIR)extractTaxonIds.txt && rm -f $(OUTPUT_DIR)extractTaxonIds.noheader
# Filter ontology by extracting only terms that are present in the file containing all IDs of all taxa and species used in Bgee or
# in the homology annotations
	@$(ROBOT_SH_PATH) extract --input $(PWD)/$(OUTPUT_DIR)$(NCBITAXON_OWL_FILE) --method BOT --term-file $(OUTPUT_DIR)extractTaxonIds.txt --annotate-with-source true --output $@_wo_disjoint_classes.owl
# Generate the final Bgee taxonomy ontology by generating disjoint classes axioms for sibling taxons
	@$(JAVA) GenerateTaxonOntology createTaxonDisjointAxioms $@_wo_disjoint_classes.owl $@.tmp.owl
# Remove temporary files
	@rm	-rf $(OUTPUT_DIR)$(NCBITAXON_OWL_FILE) $@_wo_disjoint_classes.owl
# move the Bgee taxonomy ontology to its final destination
	@$(MV) $@.tmp.owl $@

# Generate a file containing all IDs of all taxa and species used in Bgee, all taxa
# used in our annotations, all taxa used in Uberon,
# *along with their ancestor IDs*. These IDs will be used by the Uberon team to generate the taxon constraints.
ALL_TAX_IDS_FOR_UBERON := $(OUTPUT_DIR)taxon_ids_for_constraints.txt
$(ALL_TAX_IDS_FOR_UBERON): $(TAXONOMY_FILE_PATH)
	@$(JAVA) OntologyTools ExtractAllIds $< $@

$(VERIFICATIONFILE): $(ALL_TAX_IDS_FOR_UBERON)
	@$(JAVA) InsertTaxa $(BGEE_SPECIES_FILEPATH) $(ALLTAXIDSFILEPATH) $(TAXONOMY_FILE_PATH)
# insert sex-related info
	@$(MYSQL) < $(INPUT_DIR)insert_species_sex_info.sql
# get the number of species and of taxa inserted as a verification
	@$(MYSQL) -e "SELECT (SELECT COUNT(*) FROM species) as \`species count\`, (SELECT COUNT(*) FROM taxon) as \`all taxon count\`, (SELECT COUNT(*) FROM taxon WHERE bgeeSpeciesLCA = TRUE) as \`LCA taxon count\`" > $@.temp
# Display species with no sex info
	@echo "----------Looking for species with no sex information: " >> $@.temp
	@$(MYSQL) -e "SELECT t1.speciesId AS \`ERROR: SPECIES IDS WITH NO SEX INFO\` FROM species AS t1 LEFT OUTER JOIN speciesToSex AS t2 ON t1.speciesId = t2.speciesId WHERE t2.speciesId IS NULL" >> $@.temp
# Display all the species
	@echo "----------All species inserted: " >> $@.temp
	@$(MYSQL) -e "SELECT * FROM species order by speciesId" >> $@.temp
# Display the least common ancestors of species stored in Bgee
	@echo "----------Species least common ancestors: " >> $@.temp
	@$(MYSQL) -e "SELECT * FROM taxon where bgeeSpeciesLCA = TRUE order by taxonLeftBound" >> $@.temp
# Display all the taxa in taxonomy order
	@echo "----------All taxa inserted: " >> $@.temp
	@$(MYSQL) -e "SELECT * FROM taxon order by taxonLeftBound" >> $@.temp
	@$(MV) $@.temp $@

clean:
	-@$(RM) $(ANNOTTAXIDSFILEPATH) $(OUTPUT_DIR)generate_species_file $(BGEE_SPECIES_FILEPATH) $(ALLTAXIDSFILEPATH) $(OUTPUT_DIR)update_display_order $(OUTPUT_DIR)installRobot $(ROBOT_DIR_PATH) $(OUTPUT_DIR)download_ncbi_taxonomy $(OUTPUT_DIR)$(NCBITAXON_OWL_FILE) $(TAXONOMY_FILE_PATH) $(ALL_TAX_IDS_FOR_UBERON) $(VERIFICATIONFILE)

# special target to remove species and taxa from database. We make the name long on purpose, and does not
# incude it in the clean target, to not wipe the database accidentally
.PHONY := $(.PHONY) deleteSpeciesAndTaxa
deleteSpeciesAndTaxa:
	@$(MYSQL) -e "DELETE FROM species; DELETE FROM taxon;"

