PIPELINEROOT := ../
DIR_NAME := genes/
include $(PIPELINEROOT)Makefile.common


all: $(VERIFICATIONFILE)

$(OUTPUT_DIR)go:
	# Download Gene Ontology
	@$(WGET) 'ftp://ftp.geneontology.org/pub/go/ontology/obo_format_1_2/gene_ontology.1_2.obo' -O $(INPUT_DIR)go.tmp && $(MV) $(INPUT_DIR)go.tmp $(INPUT_DIR)go 2>/dev/null  || rm -f $(INPUT_DIR)go.tmp
	# Insert Gene Ontology
	@$(JAVA) InsertGO $(INPUT_DIR)go >$@.tmp 2>&1
	@$(MV) $@.tmp $@

go.obsolete: $(OUTPUT_DIR)go
	# Get obsolete GO term list
	@$(JAVA) OntologyTools ExtractObsoleteIds $(INPUT_DIR)go $@.tmp
	@$(MV) $@.tmp $@

# Deal with different Ensembl API versions for Ensembl and Ensembl Metazoa/Genomes
# Keep default for Ensembl
genes: $(SPECIESFILEPATH) $(OUTPUT_DIR)go go.obsolete
	@$(MYSQL) -e "SELECT CONCAT(s.speciesId, '__', s.genomeSpeciesId, '__', s.genus, '_', REPLACE(s.species, ' ', '_'), '__', d.dataSourceName) FROM species s, dataSource d WHERE s.dataSourceId=d.dataSourceId ORDER BY s.speciesId" | grep -v 'speciesId' >$@.tmp
	@if [[ -s $@.tmp ]]; then echo -n; else echo "NO species in file"; exit 2; fi
	# Insert gene information for all species (may be long: between 1 and 10 hours per species)
	@for sp in `cat $@.tmp`; do \
		export PERL5LIBORI=$$PERL5LIB; \
		echo "species: $$sp"; \
		if [[ `echo $$sp | grep 'Ensembl$$'` ]]; then \
			export PERL5LIB=`echo $$PERL5LIB | perl -ne 'print join ":", map { s/\/ensembl\//\/ensembl_$(ENS_API_RELEASE)\//; $$_ } grep { /\w/ } split(/:/, $$_)'`; \
		elif [[ `echo $$sp | grep 'EnsemblMetazoa$$'` ]]; then \
			export PERL5LIB=`echo $$PERL5LIB | perl -ne 'print join ":", map { s/\/ensembl\//\/ensembl_$(ENSMETAZOA_API_RELEASE)\//; $$_ } grep { /\w/ } split(/:/, $$_)'`; \
		fi; \
		perl insert_genes.pl -species=$$sp -bgee=$(BGEECMD) -ensembl=$(ENSCMD) -obsGO=go.obsolete >> $@.tmp 2>>$@.err; \
		export PERL5LIB=$$PERL5LIBORI; \
	done
	@if [[ -s $@.err ]]; then echo -n; else $(RM) $@.err; fi
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


$(VERIFICATIONFILE): $(OUTPUT_DIR)go go.obsolete genes sameIdGeneCount hogs $(OUTPUT_DIR)insert_oncoMX_XRefs
	# Check overview for all species
	@$(MYSQL) -e "SELECT t1.speciesId, COUNT(DISTINCT t1.geneId) AS geneCount, COUNT(t2.geneId) AS synonymCount, COUNT(t2.geneId)/COUNT(DISTINCT t1.geneId) AS synonymsByGene FROM gene AS t1 LEFT OUTER JOIN geneNameSynonym AS t2 ON t1.geneId = t2.geneId GROUP BY speciesId" > $@.temp
	# Check insertion of genes & generation of geneToTerm
	@$(MYSQL) -e "SELECT speciesId FROM species ORDER BY speciesId" | grep -v 'speciesId' >species
	# Check insertion of oncoMX XRefs
	@$(MYSQL) -e "select COUNT(DISTINCT t1.bgeeGeneId, t1.XRefId) AS "number of oncoMX XRefs" from geneXRef AS t1 INNER JOIN dataSource AS t2 ON t1.datasourceId = t2.datasourceId WHERE t2.dataSourceName = 'OncoMX' > $@.temp
	@for sp in `cat species`; do \
		$(MYSQL) -e "SELECT t1.speciesId, t2.geneBioTypeName, COUNT(DISTINCT t1.geneId) AS geneCount, COUNT(t3.geneId) AS xrefCount, COUNT(t3.geneId)/COUNT(DISTINCT t1.geneId) AS XRefsByGene FROM gene AS t1 LEFT OUTER JOIN geneBioType AS t2 ON t1.geneBioTypeId = t2.geneBioTypeId LEFT OUTER JOIN geneXRef AS t3 ON t1.geneId = t3.geneId WHERE t1.speciesId = $$sp GROUP BY t1.speciesId, t2.geneBioTypeName ORDER BY geneCount DESC, xrefCount DESC" >> $@.temp; \
		$(MYSQL) -e "SELECT t1.speciesId, t2.dataSourceName, COUNT(DISTINCT t1.geneId) AS geneCount, COUNT(t3.geneId) AS xrefCount, COUNT(t3.geneId)/COUNT(DISTINCT t1.geneId) AS XRefsByGene FROM gene AS t1 LEFT OUTER JOIN geneXRef AS t3 ON t1.geneId = t3.geneId LEFT OUTER JOIN dataSource AS t2 ON t2.dataSourceId = t3.dataSourceId WHERE t1.speciesId = $$sp GROUP BY t1.speciesId, t2.dataSourceName ORDER BY geneCount DESC, xrefCount DESC" >> $@.temp; \
		$(MYSQL) -e "SELECT t1.speciesId, t2.geneBioTypeName, COUNT(t3.geneId) AS termCount, COUNT(t3.geneId)/COUNT(DISTINCT t1.geneId) AS termsByGene FROM gene AS t1 LEFT OUTER JOIN geneBioType AS t2 ON t1.geneBioTypeId = t2.geneBioTypeId LEFT OUTER JOIN geneToTerm AS t3 ON t1.geneId = t3.geneId WHERE t1.speciesId = $$sp GROUP BY t1.speciesId, t2.geneBioTypeName ORDER BY termCount DESC" >> $@.temp; \
	done
	@$(RM) species
	# Add output of check of insertion of hierarchical groups
	@$(CAT) hogs >> $@.temp
	@$(MV) $@.temp $@
	# Add the verification file to GIT
	-@$(GIT) add $@
	-@$(GIT) commit $@ -m "Commit verification file for gene and hierarchical group insertion"
	@$(GIT) push

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
	-@$(RM) $(VERIFICATIONFILE) $(OUTPUT_DIR)go go.obsolete genes genes.err species *.tmp $(INPUT_DIR)$(ORTHOXMLFILE).gz $(INPUT_DIR)$(ORTHOXMLFILE) hogs
	-@$(GIT) rm $(VERIFICATIONFILE)
	-@$(GIT) commit $(VERIFICATIONFILE) -m "Removing verification file for gene (and related) and hierarchical group insertion"
	@$(GIT) push

.PHONY := $(.PHONY) deleteGene, deleteOMAHierarchicalGroup
# special target to remove genes, geneBioTypes, geneNameSynonyms, geneOntologyTerms, geneToGeneOntolongyTerms from database.
# It is not included in the clean target, to not wipe the database accidentally
deleteGene:
	@$(MYSQL) -e "DELETE FROM geneNameSynonym; DELETE FROM geneXRef; DELETE FROM geneToGeneOntologyTerm; DELETE FROM geneToTerm; \
		          DELETE FROM geneOntologyTermAltId; DELETE FROM geneOntologyTerm; DELETE FROM geneBioType; ALTER TABLE geneBioType AUTO_INCREMENT=1; DELETE FROM gene;"
deleteOMAHierarchicalGroup:
	@$(MYSQL) -e "UPDATE gene SET OMAParentNodeId=null; DELETE FROM OMAHierarchicalGroup;"
