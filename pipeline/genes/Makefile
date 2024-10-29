PIPELINEROOT := ../
DIR_NAME := genes/
include $(PIPELINEROOT)Makefile.common


all: $(VERIFICATIONFILE)

#FIXME is go still useful ???
$(OUTPUT_DIR)go:
	# Download Gene Ontology (obo format-version: 1.2)
	@$(WGET) 'http://purl.obolibrary.org/obo/go/go-basic.obo' -O $(INPUT_DIR)go.tmp && $(MV) $(INPUT_DIR)go.tmp $(INPUT_DIR)go 2>/dev/null  || rm -f $(INPUT_DIR)go.tmp
	# Insert Gene Ontology
	@$(JAVA) InsertGO $(INPUT_DIR)go >$@.tmp 2>&1
	@$(MV) $@.tmp $@

#TODO Need a way to distinct species "hidden" in Bgee, and species with data in Bgee
#TODO SPECIESFILEPATH name has changed!
genes: $(SPECIESFILEPATH)
	@$(MYSQL) -e "SELECT CONCAT(s.speciesId, '__', s.genomeSpeciesId, '__', s.genus, '_', REPLACE(s.species, ' ', '_'), '__', d.dataSourceName) FROM species s, dataSource d WHERE s.dataSourceId=d.dataSourceId ORDER BY s.speciesId" | grep -v 'speciesId' >$@.species
	@if [[ -s $@.species ]]; then echo -n; else echo "NO species in file"; exit 2; fi
	# Insert gene information for Ensembl species (may be long: between 1 and 10 hours per species)
	@for sp in `cat $@.species | grep -i __Ensembl`; do \
		date; \
		echo "species: $$sp"; \
		export spe=`echo $$sp | grep -i __Ensembl | sed -e 's|__Ensembl.*||; s|^.*__||' | tr 'A-Z' 'a-z'`; \
		if [[ $$spe == 'heterocephalus_glaber' ]]; then export spe='heterocephalus_glaber_female'; fi; \
		export url='https://ftp.ensembl.org/pub/release-$(ENSRELEASE)'; \
		if [[ `echo $$sp | grep 'EnsemblMetazoa$$'` ]]; then export url='http://ftp.ensemblgenomes.org/pub/metazoa/release-$(ENSMETAZOARELEASE)'; fi; \
		$(WGET) "$$url/json/$$spe/$$spe.json" -O $(INPUT_DIR)$$spe.tmp && $(MV) $(INPUT_DIR)$$spe.tmp $(INPUT_DIR)$$spe.json 2>/dev/null || $(RM) $(INPUT_DIR)$$spe.tmp; \
		perl insert_ensembl_genes.pl -species=$$sp -bgee=$(BGEECMD) $(INPUT_DIR)$$spe.json >> $@.tmp 2>>$@.err; \
		$(RM) $(INPUT_DIR)$$spe.json \
	done
	# Insert gene information for non-Ensembl species
	@for sp in `cat $@.species | grep -vi __Ensembl | sed -e 's| .*||'`; do \
		date; \
		echo "species: $$sp"; \
		perl insert_nonensembl_genes.pl -species=$$sp -bgee=$(BGEECMD) -bgeeSpecies=$(SPECIESFILEPATH) >> $@.tmp 2>>$@.err; \
	done
	@if [[ -s $@.err ]]; then echo -n; else $(RM) $@.err; fi
	@$(MV) $@.tmp $@

# Please check that GeneCards datasource ID did not change in the database before running updateXrefs
GENECARDS_ID := 40

#FIXME still needed with the Ensembl JSON parsing? Maybe only for non-Ensembl species
updateXrefs : genes
	# UPDATE UNIPROT XREFS
	# first remove already existing xrefs_from_uniprot file
	@rm -f $(XREFS_FROM_UNIPROT_FILE)
	# generate tsv file containing all information to update uniprot xrefs in the database. Also retrieve gene name.
	#XXX for now it uses the parallel R package and can be run on any server. Attribute "number_cores"
	#    has to be updated depending on the server the script is run on. It corresponds to the number of species queried
	#    in parallel
	#XXX Could be updated to directly run on UNIL cluster
	#XXX Could be updated to retrieve other information from Uniprot (e.g GeneCards XRefs,...)
	@R CMD BATCH --no-save --no-restore '--args species_file="$(SPECIESFILEPATH)" number_cores=28 bgee="$(BGEECMD)" output_file="$(XREFS_FROM_UNIPROT_FILE)"' update_genes_from_UniProt.R > $@.tmp
	# Now update/insert XRefs in the database
	@perl insert_uniprot_xrefs.pl -bgee=$(BGEECMD) -uniprot_xrefs=$(XREFS_FROM_UNIPROT_FILE) > $@.tmp 2>>$@.err
	# ADD GENECARDS XREFS
	@$(MYSQL) -e "INSERT INTO geneXRef (bgeeGeneId, XRefId, dataSourceId) SELECT bgeeGeneId, geneId, $(GENECARDS_ID) from gene where speciesId = 9606;" > $@.tmp 2> $@.warning
	@$(MV) $@.tmp $@

sameIdGeneCount: genes
	@$(MYSQL) -e "CREATE TEMPORARY TABLE geneCount (PRIMARY KEY(geneId)) SELECT geneId, COUNT(*) AS geneCount FROM gene GROUP BY geneId; UPDATE gene INNER JOIN geneCount ON gene.geneId = geneCount.geneId SET gene.geneMappedToGeneIdCount = geneCount.geneCount;" > $@.tmp 2> $@.warning
	@$(MV) $@.tmp $@

#Name of the file on OMA server
ORTHOXMLFILE := hogs_for_bgee.orthoxml
hogs: genes
# Download OMA OrthoXML file
	@$(WGET) http://omabrowser.org/All/$(ORTHOXMLFILE).gz -O $(ORTHOXMLFILE).gz.tmp && $(MV) $(ORTHOXMLFILE).gz.tmp $(INPUT_DIR)$(ORTHOXMLFILE).gz 2>/dev/null  || rm -f $(ORTHOXMLFILE).gz.tmp
# Gunzip it
	@$(GUNZIP) $(INPUT_DIR)$(ORTHOXMLFILE).gz
# Insert hierarchical groups and update genes
	@$(JAVA) ParseOrthoXML $(INPUT_DIR)$(ORTHOXMLFILE)
# Get the number of hierarchical groups inserted as a verification
	@$(MYSQL) -e "SELECT (SELECT COUNT(*) FROM OMAHierarchicalGroup) AS 'Total HOGs count', (SELECT COUNT(*) FROM OMAHierarchicalGroup where taxonId is not null) AS 'Number of HOGs with a taxonId not null', (SELECT COUNT(distinct taxonId) FROM OMAHierarchicalGroup where taxonId is not null) AS 'number of distinct taxonIds used', (SELECT COUNT(*) FROM gene WHERE OMAParentNodeId is not null) AS 'Number of genes with OMAParentNodeId not null', (SELECT COUNT(distinct OMAParentNodeId) FROM gene WHERE OMAParentNodeId is not null) AS 'Number of distinct OMAParentNodeId not null'"> $@.temp
	@xz -9 $(INPUT_DIR)$(ORTHOXMLFILE)
	@$(MV) $@.temp $@


$(VERIFICATIONFILE): $(OUTPUT_DIR)go genes sameIdGeneCount hogs $(OUTPUT_DIR)insert_oncoMX_XRefs
	# Check overview for all species
	@$(MYSQL) -e "SELECT t1.speciesId, COUNT(DISTINCT t1.geneId) AS geneCount, COUNT(t2.geneId) AS synonymCount, COUNT(t2.geneId)/COUNT(DISTINCT t1.geneId) AS synonymsByGene FROM gene AS t1 LEFT OUTER JOIN geneNameSynonym AS t2 ON t1.geneId = t2.geneId GROUP BY speciesId" > $@.temp
	# Check insertion of genes
	@$(MYSQL) -e "SELECT speciesId FROM species ORDER BY speciesId" | grep -v 'speciesId' >species
	# Check insertion of oncoMX XRefs
	@$(MYSQL) -e "select COUNT(DISTINCT t1.bgeeGeneId, t1.XRefId) AS "number of oncoMX XRefs" from geneXRef AS t1 INNER JOIN dataSource AS t2 ON t1.datasourceId = t2.datasourceId WHERE t2.dataSourceName = 'OncoMX' > $@.temp
	@for sp in `cat species`; do \
		$(MYSQL) -e "SELECT t1.speciesId, t2.geneBioTypeName, COUNT(DISTINCT t1.geneId) AS geneCount, COUNT(t3.geneId) AS xrefCount, COUNT(t3.geneId)/COUNT(DISTINCT t1.geneId) AS XRefsByGene FROM gene AS t1 LEFT OUTER JOIN geneBioType AS t2 ON t1.geneBioTypeId = t2.geneBioTypeId LEFT OUTER JOIN geneXRef AS t3 ON t1.geneId = t3.geneId WHERE t1.speciesId = $$sp GROUP BY t1.speciesId, t2.geneBioTypeName ORDER BY geneCount DESC, xrefCount DESC" >> $@.temp; \
		$(MYSQL) -e "SELECT t1.speciesId, t2.dataSourceName, COUNT(DISTINCT t1.geneId) AS geneCount, COUNT(t3.geneId) AS xrefCount, COUNT(t3.geneId)/COUNT(DISTINCT t1.geneId) AS XRefsByGene FROM gene AS t1 LEFT OUTER JOIN geneXRef AS t3 ON t1.geneId = t3.geneId LEFT OUTER JOIN dataSource AS t2 ON t2.dataSourceId = t3.dataSourceId WHERE t1.speciesId = $$sp GROUP BY t1.speciesId, t2.dataSourceName ORDER BY geneCount DESC, xrefCount DESC" >> $@.temp; \
	done
	@$(RM) species
	# Add output of check of insertion of hierarchical groups
	@$(CAT) hogs >> $@.temp
	@$(MV) $@.temp $@
	# Add the verification file to GIT
	-@$(GIT) add $@
	-@$(GIT) commit $@ -m "Commit verification file for gene and hierarchical group insertion"
	@$(GIT) push


#number of processes, default is 3 maximum is 4
N_PROCESSES = 4

#Log mode (debug or info)
LOG_HOMOLOGS = debug

#Skip species pairs stated in the {orthologs|paralogs}_pairs_with_results.tmp (see README.md)
#The tmp files should be placed in the TEMPORARY_CONFIG_FILES_DIRECTORY. Possible values: True and False
REMOVE_SPECIES_PAIRS = False

#For more information about the config template, please see the README.md file of the orthologs/paralogs app
define CONFIG_TEMPLATE
[DEFAULT]
bgee_sparql_endpoint = https://www.bgee.org/sparql/
oma_sparql_endpoint = https://sparql.omabrowser.org/sparql/
ncbi_gene2ensembl_URL = $(NCBI_GENE2ENSEMBL_URL)
log = $(LOG_HOMOLOGS)
remove_species_pairs = $(REMOVE_SPECIES_PAIRS)

[OUTPUT]
file_directory_output = $(HOMOLOGS_DIR)

[RDF_VOCABULARY]
ncbi_gene_crossref_property = lscr:xrefNCBIGene
ensembl_gene_crossref_property = lscr:xrefEnsemblGene

[SPECIES]
species_file = $(ALL_SPECIES_FILE)
ncbi_species_file = $(ONLY_NCBI_SPECIES_FILE)
flybase_species_file = $(ONLY_FLYBASE_SPECIES_FILE)
taxon_namespace = http://purl.uniprot.org/taxonomy/
species_to_gene_id_prefix = {6239: 'WBGene', 7740: 'BL'}
species_ncbi_to_ensembl = [105023]

[SPECIES_MAPPING]
7237 = 46245
9593 = 9595
endef

export CONFIG_TEMPLATE

$(OUTPUT_DIR)config_file:
	@echo "$$CONFIG_TEMPLATE" > $(GENERATE_HOMOLOGS_APP_DIR)/template.properties
	@touch $@

$(OUTPUT_DIR)generate_homologs_app: config_file
	@cd $(GENERATE_HOMOLOGS_APP_DIR)
	@$(PIPENV) install 1>$@.tmp 2>$@.err
	@$(PIPENV) run $(PYTHON3) $(GENERATE_HOMOLOGS_APP_FILE) -c $(TEMPORARY_CONFIG_FILES_DIR) -n $(N_PROCESSES) >> $@.tmp 2>> $@.err
	@$(MV) $(HOMOLOGS_DIR)/orthologs* $(ORTHOLOGS_DIR)
	@$(MV) $(HOMOLOGS_DIR)/paralogs* $(PARALOGS_DIR)
	@$(RM) -r $(HOMOLOGS_DIR)
	@$(MV) $@.tmp $@

$(OUTPUT_DIR)insert_gene_homology:
	# Insert orthology/paralogy  into the database
	@perl insert_homology.pl -bgee=$(BGEECMD) -paralogs_dir_path=$(PARALOGS_DIR) -orthologs_dir_path=$(ORTHOLOGS_DIR) >>$@.tmp 2>>$@.err
	@if [[ -s $@.err ]]; then echo -n; else $(RM) $@.err; fi
	@$(MV) $@.tmp $@

# Init variables used to insert OncoMX XRefs in the database
ONCOMX_DATASOURCE_NAME := OncoMX

$(OUTPUT_DIR)insert_oncoMX_XRefs:
	# Download oncoMX ids (uniprot ID and gene name)
	# TODO : OncoMX is supposed to update the tsv file containing their IDs. Please be sure that the URL used to download this dataset is still up to date
	@$(WGET) --no-check-certificate 'https://data.oncomx.org/ln2wwwdata/reviewed/biomarker_masterlist.csv' -O $(INPUT_DIR)oncoMX.csv.tmp >>$@.tmp 2>>$@.err && $(MV) $(INPUT_DIR)oncoMX.csv.tmp $(INPUT_DIR)oncoMX.csv 2>/dev/null  || rm -f $(INPUT_DIR)oncoMX.csv.tmp
	# Insert OncoMX XRefs into the database
	@perl insert_oncoMX_xrefs.pl -bgee=$(BGEECMD) -oncoMX_ids_file=$(INPUT_DIR)oncoMX.csv -oncoMX_datasource_name=$(ONCOMX_DATASOURCE_NAME) >>$@.tmp 2>>$@.err
	@if [[ -s $@.err ]]; then echo -n; else $(RM) $@.err; fi
	@$(MV) $@.tmp $@

# Rule to run if oncoMX ids changed. Will remove all oncoMX XRefs already present in the Bgee database
$(OUTPUT_DIR)delete_oncoMX_XRefs:
	@$(MYSQL) -e "DELETE FROM geneXRef INNER JOIN dataSource ON dataSource.dataSourceId = geneXRef.dataSourceId WHERE dataSource.dataSourceName = $(ONCOMX_DATASOURCE_NAME);" > $@.tmp
	@$(MV) $@.tmp $@

clean:
	-@$(RM) $(VERIFICATIONFILE) $(OUTPUT_DIR)go genes genes.err species *.tmp $(INPUT_DIR)$(ORTHOXMLFILE).gz $(INPUT_DIR)$(ORTHOXMLFILE) hogs
	-@$(GIT) rm $(VERIFICATIONFILE)
	-@$(GIT) commit $(VERIFICATIONFILE) -m "Removing verification file for gene (and related) and hierarchical group insertion"
	@$(GIT) push

.PHONY := $(.PHONY) deleteGene, deleteOMAHierarchicalGroup
# special target to remove genes, geneBioTypes, geneNameSynonyms, geneOntologyTerms, geneToGeneOntolongyTerms from database.
# It is not included in the clean target, to not wipe the database accidentally
deleteGene:
	@$(MYSQL) -e "DELETE FROM geneNameSynonym; DELETE FROM geneXRef; \
		          DELETE FROM geneOntologyTermAltId; DELETE FROM geneOntologyTerm; DELETE FROM geneBioType; ALTER TABLE geneBioType AUTO_INCREMENT=1; DELETE FROM gene;"
deleteOMAHierarchicalGroup:
	@$(MYSQL) -e "UPDATE gene SET OMAParentNodeId=null; DELETE FROM OMAHierarchicalGroup;"
