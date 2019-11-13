PIPELINEROOT := ../
DIR_NAME := species/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)

# Generate a file with only the taxon IDs present in the similarity annotation file.
$(ANNOTTAXIDSFILEPATH):
# download similarity file in TMP dir
	@$(CURL) $(RAW_SIM_ANNOT_DOWNLOAD_URL) -z $(TMPDIR)similarity_annotations.tsv -o $(TMPDIR)similarity_annotations.tsv.tmp && $(MV) $(TMPDIR)similarity_annotations.tsv.tmp $(TMPDIR)similarity_annotations.tsv 2>/dev/null || rm -f $(TMPDIR)similarity_annotations.tsv.tmp
	@$(JAVA) SimilarityAnnotation extractTaxonIds $(TMPDIR)similarity_annotations.tsv $@.tmp
# need to add a header
	@echo "taxon ID" | cat - $@.tmp > $@
	@$(RM) $@.tmp

# generate a file containing the IDs of all species used in Bgee, all taxa
# used in the similarity annotations, all taxa used in Uberon
$(ALLTAXIDSFILEPATH): $(ANNOTTAXIDSFILEPATH) $(SPECIESFILEPATH)
# download all Uberon versions with different taxon constraints in tmp dir
	@$(CURL) $(EXT_DOWNLOAD_URL) -z $(TMPDIR)ext.owl -o $(TMPDIR)ext.owl.tmp && $(MV) $(TMPDIR)ext.owl.tmp $(TMPDIR)ext.owl 2>/dev/null || rm -f $(TMPDIR)ext.owl.tmp
	@$(CURL) $(COMPOSITE_DOWNLOAD_URL) -z $(TMPDIR)composite-metazoan.owl -o $(TMPDIR)composite-metazoan.owl.tmp && $(MV) $(TMPDIR)composite-metazoan.owl.tmp $(TMPDIR)composite-metazoan.owl 2>/dev/null || rm -f $(TMPDIR)composite-metazoan.owl.tmp
	@$(JAVA) Uberon extractTaxonIds $(TMPDIR)ext.owl $(TMPDIR)uberonTaxIds.tsv.tmp
	@$(JAVA) Uberon extractTaxonIds $(TMPDIR)composite-metazoan.owl $(TMPDIR)uberonTaxIds_composite.tsv.tmp
# get the annotation taxon IDs without the first header line
	@sed 1d $(ANNOTTAXIDSFILEPATH) > $(TMPDIR)$(ANNOTTAXIDSFILE).noheader
# concatenate these annotation taxa with the Uberon taxa, no header in this file
	@cat $(TMPDIR)$(ANNOTTAXIDSFILE).noheader $(TMPDIR)uberonTaxIds.tsv.tmp $(TMPDIR)uberonTaxIds_composite.tsv.tmp | sort | uniq > $@.tmp
# now concatenate with the file containing the IDs of species used in Bgee, and add its header
	@cat $(SPECIESFILEPATH) $@.tmp > $@
	@$(RM) $(TMPDIR)$(ANNOTTAXIDSFILE).noheader $(TMPDIR)uberonTaxIds.tsv.tmp $(TMPDIR)uberonTaxIds_composite.tsv.tmp $@.tmp

# name of the source file for the taxonomy
TAXDAT := $(INPUT_DIR)taxonomy.dat
$(TAXONOMYFILEPATH): $(ALLTAXIDSFILEPATH)
	@$(CURL) $(TAX_DAT_DOWNLOAD_URL) -z $(TAXDAT) -o $(TAXDAT).tmp && $(MV) $(TAXDAT).tmp $(TAXDAT) 2>/dev/null || rm -f $(TAXDAT).tmp
# the extension of the ontology file must be .owl. That is why we add the ".tmp.owl" suffix rather than the usual .tmp
	@$(JAVA) GenerateTaxonOntology $(TAXDAT) $^ $@.tmp.owl
	@xz -9 $(TAXDAT)
	@$(MV) $@.tmp.owl $@
	

# Generate a file containing all IDs of all taxa and species used in Bgee, all taxa
# used in our annotations, all taxa used in Uberon,
# *along with their ancestor IDs*. These IDs will be used by the Uberon team to generate the taxon constraints.
ALL_TAX_IDS_FOR_UBERON := $(OUTPUT_DIR)taxon_ids_for_constraints.txt
$(ALL_TAX_IDS_FOR_UBERON): $(TAXONOMYFILEPATH)
	@$(JAVA) OntologyTools ExtractAllIds $< $@

$(VERIFICATIONFILE): $(SPECIESFILEPATH) $(ANNOTTAXIDSFILEPATH) $(ALLTAXIDSFILEPATH) $(TAXONOMYFILEPATH) $(ALL_TAX_IDS_FOR_UBERON)
	@$(JAVA) InsertTaxa $(SPECIESFILEPATH) $(ALLTAXIDSFILEPATH) $(TAXONOMYFILEPATH)
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
	-@$(RM) $(TAXDAT).tmp $(TAXDAT) $(TAXDAT).xz $(ALLTAXIDSFILEPATH).tmp $(ALLTAXIDSFILEPATH) $(ANNOTTAXIDSFILEPATH).tmp $(ANNOTTAXIDSFILEPATH) $(TAXONOMYFILEPATH) $(VERIFICATIONFILE) $(ALL_TAX_IDS_FOR_UBERON)

# special target to remove species and taxa from database. We make the name long on purpose, and does not
# incude it in the clean target, to not wipe the database accidentally
.PHONY := $(.PHONY) deleteSpeciesAndTaxa
deleteSpeciesAndTaxa:
	@$(MYSQL) -e "DELETE FROM species; DELETE FROM taxon;"

