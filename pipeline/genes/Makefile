PIPELINEROOT := ../
include $(PIPELINEROOT)Makefile.common


all: $(VERIFICATIONFILE)

go:
	# Download Gene Ontology
	@$(WGET) 'http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo' -O $@.tmp && $(MV) $@.tmp $(SOURCE_FILES_DIR)/$@ 2>/dev/null  || rm -f $@.tmp
	# Insert Gene Ontology
	@$(JAVA) InsertGO $(SOURCE_FILES_DIR)/$@  >$@ 2>&1

go.obsolete: go
	# Get obsolete GO term list
	@$(JAVA) OntologyTools ExtractObsoleteIds $(SOURCE_FILES_DIR)/$^ $@

# Deal with different Ensembl API versions for Ensembl and Ensembl Metazoa/Genomes
# Keep default for Ensembl
genes: $(SPECIESFILEPATH) go go.obsolete
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
	@$(WGET) http://omabrowser.org/All/$(ORTHOXMLFILE).gz -O $(ORTHOXMLFILE).gz.tmp && mv $(ORTHOXMLFILE).gz.tmp $(ORTHOXMLFILE).gz 2>/dev/null  || rm -f $(ORTHOXMLFILE).gz.tmp
# Gunzip it
	@$(GUNZIP) $(ORTHOXMLFILE).gz
# Insert hierarchical groups and update genes
	@$(JAVA) ParseOrthoXML $(ORTHOXMLFILE)
# Get the number of hierarchical groups inserted as a verification
	@$(MYSQL) -e "SELECT (SELECT COUNT(*) FROM OMAHierarchicalGroup) AS 'Total HOGs count', (SELECT COUNT(*) FROM OMAHierarchicalGroup where taxonId is not null) AS 'Number of HOGs with a taxonId not null', (SELECT COUNT(distinct taxonId) FROM OMAHierarchicalGroup where taxonId is not null) AS 'number of distinct taxonIds used', (SELECT COUNT(*) FROM gene WHERE OMAParentNodeId is not null) AS 'Number of genes with OMAParentNodeId not null', (SELECT COUNT(distinct OMAParentNodeId) FROM gene WHERE OMAParentNodeId is not null) AS 'Number of distinct OMAParentNodeId not null'"> $@.temp
	@$(MV) $@.temp $@


$(VERIFICATIONFILE): go go.obsolete genes sameIdGeneCount hogs
	# Check overview for all species
	@$(MYSQL) -e "SELECT t1.speciesId, COUNT(DISTINCT t1.geneId) AS geneCount, COUNT(t2.geneId) AS synonymCount, COUNT(t2.geneId)/COUNT(DISTINCT t1.geneId) AS synonymsByGene FROM gene AS t1 LEFT OUTER JOIN geneNameSynonym AS t2 ON t1.geneId = t2.geneId GROUP BY speciesId" > $@.temp
	# Check insertion of genes & generation of geneToTerm
	@$(MYSQL) -e "SELECT speciesId FROM species ORDER BY speciesId" | grep -v 'speciesId' >species
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

clean:
	-@$(RM) $(VERIFICATIONFILE) go go.obsolete genes genes.err species *.tmp $(ORTHOXMLFILE).gz $(ORTHOXMLFILE) hogs
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
