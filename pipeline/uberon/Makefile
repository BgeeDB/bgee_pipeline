PIPELINEROOT := ../
DIR_NAME := uberon/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)


$(UBERON_COMPOSITE_FILE_PATH):
	@$(CURL) $(COMPOSITE_DOWNLOAD_URL) $(APPEND_CURL_COMMAND)
$(UBERON_EXT_FILE_PATH):
	@$(CURL) $(EXT_DOWNLOAD_URL) $(APPEND_CURL_COMMAND)

# Note: this step requires lots of memory (was run with 512 Go)
$(OUTPUT_DIR)uberon_merged_with_taxonomy: $(UBERON_COMPOSITE_FILE_PATH) $(UBERON_EXT_FILE_PATH)
	@$(JAVA) TaxonConstraints mergeUberonAndTaxonomy $(UBERON_COMPOSITE_FILE_PATH) $(TAXONOMYFILEPATH) $(OUTPUT_DIR)merged_composite-metazoan_taxonomy > $@.tmp 2> $@.warn
	@$(JAVA) TaxonConstraints mergeUberonAndTaxonomy $(UBERON_EXT_FILE_PATH) $(TAXONOMYFILEPATH) $(OUTPUT_DIR)merged_ext_taxonomy >> $@.tmp 2>> $@.warn
	@$(MV) $@.tmp $@

# As of Bgee_v15, the composite-metazoan.owl does not contain "never_in_taxon" relations,
# see https://github.com/obophenotype/uberon/issues/705#issuecomment-233767158
# So only ext.owl is usable to generate taxon constraints. We will then miss constraints
# for all species-specific terms, and will need to enforce later them
# (we don't do it at this point, although it is possible, to generate the "original" constraints).
$(TAXONCONSTRAINTSFILEPATH): $(UBERON_EXT_FILE_PATH) $(TAXONOMYFILEPATH) $(SPECIESFILEPATH)
	@$(RM) $(OUTPUT_DIR)/intermediate_ontologies/*
	@$(JAVA) TaxonConstraints generateTaxonConstraints $^ $(SIMPLIFICATION_STEPS) - - $@.tmp $(OUTPUT_DIR)/intermediate_ontologies/
	@$(MV) $@.tmp $@

# # Uberon simplification
# ## Original command:
# java -Xmx2g -jar myJar Uberon simplifyUberon composite-metazoan.owl custom_composite simplification_composite.tsv UBERON:0001062,UBERON:0000465,UBERON:0000061,UBERON:0010000,UBERON:0008979 BFO:0000050,RO:0002202,RO:0002494 UBERON:0000922/UBERON:0002050,UBERON:0004716/UBERON:0000922,UBERON:0000467/UBERON:0000468,UBERON:0000475/UBERON:0000468,UBERON:0000479/UBERON:0000468,UBERON:0000480/UBERON:0000468,UBERON:0007688/UBERON:0000468,UBERON:0010707/UBERON:0000468,UBERON:0012641/UBERON:0000468,UBERON:0002199/UBERON:0000468,UBERON:0002416/UBERON:0000468,UBERON:0007376/UBERON:0000468,UBERON:0000463/UBERON:0000468,UBERON:0001048/UBERON:0000468,UBERON:0007567/UBERON:0000468,UBERON:0015119/UBERON:0000468 NBO:0000313,GO:0008150,ENVO:01000254,BFO:0000040,GO:0003674,PATO:0000001,CHEBI:24431,UBERON:0004458,UBERON:0000466,SO:0000704 UBERON:0013701,UBERON:0000026,UBERON:0000480,UBERON:0000479,UBERON:0000468,GO:0005575 grouping_class,non_informative,ubprop:upper_level,upper_level UBERON:0013701,UBERON:0000026,UBERON:0000480,UBERON:0000479,UBERON:0011676,GO:0005575
#
# See README.md in this folder for a up-to-date description of the options of this command
# roots of subgraph to filter after test: UBERON:0013701$(LIST_SEP)UBERON:0000026$(LIST_SEP)UBERON:0000480$(LIST_SEP)UBERON:0000479$(LIST_SEP)UBERON:0000468$(LIST_SEP)GO:0005575
# subsets for incoming edges filtering after test: grouping_class$(LIST_SEP)non_informative$(LIST_SEP)ubprop:upper_level$(LIST_SEP)upper_level
$(CUSTOM_UBERON_FILE_PATH): $(UBERON_COMPOSITE_FILE_PATH)
	@$(JAVA) Uberon simplifyUberon $(UBERON_COMPOSITE_FILE_PATH) $(UBERON_OUTPUT_PATH)$(CUSTOM_UBERON_PREFIX).tmp $(UBERON_OUTPUT_PATH)simplification_composite_info.tsv CARO:0000006$(LIST_SEP)UBERON:0010000$(LIST_SEP)UBERON:0008979$(LIST_SEP)GO:0016265$(LIST_SEP)BFO:0000030$(LIST_SEP)BFO:0000020$(LIST_SEP)BFO:0000004$(LIST_SEP)BFO:0000002$(LIST_SEP)UBERON:0001062$(LIST_SEP)BFO:0000040 BFO:0000050$(LIST_SEP)RO:0002202$(LIST_SEP)RO:0002494 UBERON:0010707$(KEY_VAL_SEP)UBERON:0000468$(ENTRY_SEP)UBERON:0000463$(KEY_VAL_SEP)UBERON:0000468$(ENTRY_SEP)UBERON:0000481$(KEY_VAL_SEP)NCBITaxon:6072$(ENTRY_SEP)UBERON:0005212$(KEY_VAL_SEP)NCBITaxon:7742 NBO:0000313$(LIST_SEP)CEPH:0000300$(LIST_SEP)CHEBI:24431$(LIST_SEP)PATO:0000001$(LIST_SEP)BFO:0000019$(LIST_SEP)ENVO:01000254$(LIST_SEP)GO:0003674$(LIST_SEP)CHEBI:23367$(LIST_SEP)CL:0000063$(LIST_SEP)UBERON:0000952$(LIST_SEP)WBls:0000075$(LIST_SEP)BFO:0000003$(LIST_SEP)FBdv:00007008$(LIST_SEP)CS:0$(LIST_SEP)ZFS:0100000$(LIST_SEP)OlatDv:0000010$(LIST_SEP)UBERON:0000466$(LIST_SEP)PORO:0000019$(LIST_SEP)BFO:0000141$(LIST_SEP)PORO:0000923$(LIST_SEP)PORO:0000963$(LIST_SEP)SO:0000704 UBERON:0000465$(LIST_SEP)GO:0005575 upper_level UBERON:0013701$(LIST_SEP)UBERON:0000026$(LIST_SEP)UBERON:0000480$(LIST_SEP)FBbt:00007276$(LIST_SEP)UBERON:0000479$(LIST_SEP)UBERON:0011676$(LIST_SEP)UBERON:0005423$(LIST_SEP)UBERON:0000463$(LIST_SEP)GO:0005575
	@$(MV) $(UBERON_OUTPUT_PATH)$(CUSTOM_UBERON_PREFIX).tmp.obo $(UBERON_OUTPUT_PATH)$(CUSTOM_UBERON_PREFIX).obo
	@$(MV) $(UBERON_OUTPUT_PATH)$(CUSTOM_UBERON_PREFIX).tmp.owl $@

#Extract sex-related info about anatomical terms from Uberon
$(UBERON_SEX_INFO_FILE_PATH): $(CUSTOM_UBERON_FILE_PATH)
	@$(JAVA) Uberon extractSexInfo $< $@.tmp
	@$(MV) $@.tmp $@

# The composite-metazoan dev. stage ontology is now built separately on the related tracker,
# see https://github.com/obophenotype/developmental-stage-ontologies/tree/master/external
# Note, to extract only the life cycle stages from Uberon, the command used was:
# @$(JAVA) UberonDevStage generateStageOntology $(UBERON_COMPOSITE_FILE_PATH) $(UBERON_OUTPUT_PATH)$(DEV_STAGE_ONT_PREFIX) UBERON:0000067$(LIST_SEP)UBERON:0000071$(LIST_SEP)UBERON:0000105$(LIST_SEP)UBERON:0000000$(LIST_SEP)FBdv:00007008$(LIST_SEP)MmusDv:0000041$(LIST_SEP)ZFS:0000000$(LIST_SEP)UBERON:0035944$(LIST_SEP)UBERON:0035945$(LIST_SEP)BFO:0000015 - UBERON:0000481$(KEY_VAL_SEP)NCBITaxon:6072 BFO:0000050$(LIST_SEP)BFO:0000062$(LIST_SEP)RO:0002087 UBERON:0000104$(LIST_SEP)WBls:0000075$(LIST_SEP)FBdv:00005259$(LIST_SEP)ZFS:0100000$(LIST_SEP)XAO:1000094$(LIST_SEP)NCBITaxon:1

#TODO See https://gitlab.sib.swiss/Bgee/bgee_apps/issues/192#note_8407
# Remove redundant species name in species specific stages
$(DEV_STAGE_ONT_FILE_PATH):
	@$(CURL) $(DEV_STAGE_DOWNLOAD_URL) $(APPEND_CURL_COMMAND)

$(OUTPUT_DIR)insert_stages: $(DEV_STAGE_ONT_FILE_PATH) $(TAXONOMYFILEPATH) $(TAXONCONSTRAINTSFILEPATH) $(SPECIESFILEPATH)
	@$(JAVA) InsertUberon insertStages $(DEV_STAGE_ONT_FILE_PATH) $(TAXONOMYFILEPATH) $(TAXONCONSTRAINTSFILEPATH) $(OVERRIDE_TAXON_CONSTRAINTS) $(ROOTS_OF_SUBGRAPHS_TO_IGNORE) $(SPECIESFILEPATH)
	@echo "----------Number of developmental stages inserted: " > $@.temp
	@$(MYSQL) -e "SELECT count(*) from stage" >> $@.temp
	@echo "----------Number of terms in the ontology: " >> $@.temp
	@grep "\\[Term\\]" $(DEV_STAGE_ONT_FILE_PATH) | wc >> $@.temp
	@echo "----------Stages existing in all species: " >> $@.temp
	@$(MYSQL) -e "SELECT t1.stageId, t1.stageName, t1.stageLeftBound, t1.stageRightBound, t1.stageLevel FROM stage AS t1 INNER JOIN stageTaxonConstraint AS t2 ON t1.stageId = t2.stageId WHERE t2.speciesId IS NULL ORDER BY t1.stageLeftBound" >> $@.temp
	@echo "----------UBERON stages with taxon constraints: " >> $@.temp
	@$(MYSQL) -e "SELECT t1.stageId, t1.stageName, t1.stageLeftBound, t1.stageRightBound, t1.stageLevel, GROUP_CONCAT(t2.speciesId ORDER BY t2.speciesId) AS speciesIdList FROM stage AS t1 INNER JOIN stageTaxonConstraint AS t2 ON t1.stageId = t2.stageId WHERE t1.stageId LIKE 'UBERON:%' GROUP BY t1.stageId ORDER BY t1.stageLeftBound" >> $@.temp
	@echo "----------Number of developmental stages per species: " >> $@.temp
	@$(MYSQL) -e "SELECT t1.speciesId, t1.genus, t1.species, COUNT(DISTINCT t2.stageId) AS allStageCount, COUNT(DISTINCT t3.stageId) AS stageWithTaxonConstraintsCount FROM species AS t1 INNER JOIN stageTaxonConstraint AS t2 ON t1.speciesId = t2.speciesId OR t2.speciesId IS NULL LEFT OUTER JOIN stageTaxonConstraint AS t3 ON t1.speciesId = t3.speciesId GROUP BY t1.speciesId ORDER BY t1.speciesId" >> $@.temp
	@echo "----------Developmental stages with a parent not existing in the species considered. There must be only 'UBERON:0000104 life cycle' (the parent of a term must always exist in the species where one of its children exist): " >> $@.temp
	@$(MYSQL) -e "select DISTINCT t1.stageId, t1.stageName, t1.stageLeftBound, t1.stageRightBound, t1.stageLevel, t2.speciesId from stage as t1 inner join stageTaxonConstraint as t2 on t1.stageId = t2.stageId where NOT EXISTS (SELECT 1 FROM stage AS t3 inner join stageTaxonConstraint as t4 on t3.stageId = t4.stageId where (t4.speciesId IS NULL OR t2.speciesId IS NULL OR t4.speciesId = t2.speciesId) AND t3.stageLeftBound <= t1.stageLeftBound AND t3.stageRightBound >= t1.stageRightBound AND t3.stageLevel = t1.stageLevel - 1) order by t2.speciesId, t1.stageLeftBound" >> $@.temp
# Remove the header line from bgeeSpecies.tsv, to iterate each species
	@sed 1d $(SPECIESFILEPATH) > $(TMPDIR)$(species).noheader
	@echo "----------Now, dev. stages per species in readable format: " >> $@.temp
	@while IFS=$$'\t' read -r speId genus species; do \
		echo "---Species: $$speId $$genus $$species" >> $@.temp; \
		$(MYSQL) -e "SELECT CONCAT(REPEAT(' ', (4*(t1.stageLevel-1))), t1.stageId, ' ', t1.stageName) from stage as t1 inner join stageTaxonConstraint as t2 on t1.stageId = t2.stageId where t2.speciesId IS NULL OR t2.speciesId = $$speId order by stageLeftBound" >> $@.temp; \
	done <$(TMPDIR)$(species).noheader
	@echo "----------Now, dev. stages per species: " >> $@.temp
	@while IFS=$$'\t' read -r speId genus species; do \
		echo "---Species: $$speId $$genus $$species" >> $@.temp; \
		$(MYSQL) -e "select t1.stageId, t1.stageName, t1.stageLeftBound, t1.stageRightBound, t1.stageLevel from stage as t1 inner join stageTaxonConstraint as t2 on t1.stageId = t2.stageId where t2.speciesId IS NULL OR t2.speciesId = $$speId order by stageLeftBound" >> $@.temp; \
	done <$(TMPDIR)$(species).noheader
	@$(RM) $(TMPDIR)$(species).noheader
	@$(MV) $@.temp $@

$(OUTPUT_DIR)insert_anatomy: $(CUSTOM_UBERON_FILE_PATH) $(TAXONCONSTRAINTSFILEPATH) $(SPECIESFILEPATH)
	@$(JAVA) InsertUberon insertAnatomy $(CUSTOM_UBERON_FILE_PATH) $(TAXONCONSTRAINTSFILEPATH) $(OVERRIDE_TAXON_CONSTRAINTS) $(ROOTS_OF_SUBGRAPHS_TO_IGNORE) $(SPECIESFILEPATH)
	@echo "Inconsistent taxon constraints between anatomical entities and anatomical entity relations (there should be none)" > $@.temp
	@$(MYSQL) -e "SELECT t1.*, GROUP_CONCAT(DISTINCT t2.speciesId ORDER BY t2.speciesId) AS relationTaxonConstraints, GROUP_CONCAT(DISTINCT t3.speciesId ORDER BY t3.speciesId) AS sourceTaxonContraints, GROUP_CONCAT(DISTINCT t4.speciesId ORDER BY t4.speciesId) AS targetTaxonConstraints FROM anatEntityRelation AS t1 INNER JOIN anatEntityRelationTaxonConstraint AS t2 ON t1.anatEntityRelationId = t2.anatEntityRelationId INNER JOIN anatEntityTaxonConstraint AS t3 ON t1.anatEntitySourceId = t3.anatEntityId INNER JOIN anatEntityTaxonConstraint AS t4 ON t1.anatEntityTargetId = t4.anatEntityId WHERE NOT EXISTS (SELECT 1 FROM anatEntityTaxonConstraint AS t5 WHERE t5.anatEntityId = t1.anatEntitySourceId AND (t5.speciesId IS NULL OR t2.speciesId IS NULL OR t5.speciesId = t2.speciesId)) OR NOT EXISTS (SELECT 1 FROM anatEntityTaxonConstraint AS t5 WHERE t5.anatEntityId = t1.anatEntityTargetId AND (t5.speciesId IS NULL OR t2.speciesId IS NULL OR t5.speciesId = t2.speciesId)) GROUP BY t1.anatEntityRelationId;" >> $@.temp 2>> $@.warn
	@echo "\n----------Number of anatomical entities: " >> $@.temp
	@$(MYSQL) -e 'SELECT (SELECT count(*) from anatEntity) AS `Number of anatomical entities`, (SELECT count(distinct t1.anatEntityId) from anatEntity as t1 inner join anatEntityTaxonConstraint as t2 on t1.anatEntityId = t2.anatEntityId) AS `With taxon constraints defined (should be equal)`, (SELECT count(*) from anatEntity where nonInformative != 0) AS `Number of non-informative anatomical entities`, (SELECT count(distinct anatEntityId) from anatEntityTaxonConstraint where speciesId is null) AS `Anatomical entities existing in all species`, (SELECT count(distinct anatEntityId) from anatEntityTaxonConstraint where speciesId is not null) AS `Anatomical entities with taxon-specific constraints`' >> $@.temp
	@echo "----------Number of anatomical entities per ontology: " >> $@.temp
	@$(MYSQL) -e 'SELECT SUBSTRING_INDEX(anatEntityId, ":", 1) AS ontology, count(DISTINCT anatEntityId) FROM anatEntity GROUP BY ontology ORDER BY ontology ASC' >> $@.temp
	@echo "----------Number of anatomical entities defined by taxon constraints: " >> $@.temp
	@$(MYSQL) -e "SET SESSION group_concat_max_len = 1024; SELECT (SELECT t6.taxonScientificName FROM taxon AS t6 where t6.taxonLeftBound <= temp1.minLeftBound AND t6.taxonRightBound >= temp1.maxRightBound ORDER BY t6.taxonLeftBound DESC LIMIT 1) AS commonAncestor, temp1.species_list, temp1.speciesCount, temp1.count, temp1.anatEntityExamples FROM (SELECT t3.species_list, t3.speciesCount, count(t3.anatEntityId) AS count, MIN(taxonLeftBound) as minLeftBound, MAX(taxonRightBound) AS maxRightBound, group_concat(t3.anatEntityName ORDER BY RAND()) as anatEntityExamples FROM (SELECT t1.anatEntityId, t4.anatEntityName, MIN(t5.taxonLeftBound) AS taxonLeftBound, MAX(t5.taxonRightBound) AS taxonRightBound, group_concat(IF(t1.speciesId IS NULL, 'All species', IF (t2.speciesCommonName != '', t2.speciesCommonName, CONCAT_WS(' ', t2.genus, t2.species))) ORDER BY IF(t1.speciesId IS NULL, 'All species', IF (t2.speciesCommonName != '', t2.speciesCommonName, CONCAT_WS(' ', t2.genus, t2.species)))) AS species_list, COUNT(DISTINCT t2.speciesId) AS speciesCount from anatEntityTaxonConstraint as t1 INNER JOIN anatEntity AS t4 ON t1.anatEntityId = t4.anatEntityId left outer join species as t2 on t1.speciesId = t2.speciesId LEFT OUTER JOIN taxon AS t5 ON t2.taxonId = t5.taxonId group by t1.anatEntityId) AS t3 group by t3.species_list) AS temp1 order by commonAncestor ASC, speciesCount DESC, count DESC" >> $@.temp
	@echo "----------Roots of the ontology by any relation: " >> $@.temp
	@$(MYSQL) -e 'select t1.* from anatEntity as t1 left outer join anatEntityRelation as t2 on t1.anatEntityId = t2.anatEntitySourceId and t2.relationStatus != "reflexive" where t2.anatEntitySourceId is null' >> $@.temp
	@echo "----------Roots of the ontology by is_a/part_of relations: " >> $@.temp
	@$(MYSQL) -e 'select t1.* from anatEntity as t1 left outer join anatEntityRelation as t2 on t1.anatEntityId = t2.anatEntitySourceId and t2.relationStatus != "reflexive" and t2.relationType = "is_a part_of" where t2.anatEntitySourceId is null' >> $@.temp
	@echo "----------Number of relations between anatomical entities: " >> $@.temp
	@$(MYSQL) -e 'SELECT (SELECT count(*) from anatEntityRelation) AS `Number of relations`, (SELECT count(distinct t1.anatEntityRelationId) from anatEntityRelation as t1 inner join anatEntityRelationTaxonConstraint as t2 on t1.anatEntityRelationId = t2.anatEntityRelationId) AS `With taxon constraints defined (should be equal)`, (SELECT count(distinct anatEntityRelationId) from anatEntityRelationTaxonConstraint where speciesId is null) AS `Relations existing in all species`, (SELECT count(distinct anatEntityRelationId) from anatEntityRelationTaxonConstraint where speciesId is not null) AS `Relations with taxon-specific constraints`' >> $@.temp
	@echo "----------Relations per status and type: " >> $@.temp
	@$(MYSQL) -e 'SELECT t1.relationStatus, t1.relationType, count(distinct t1.anatEntityRelationId) AS count from anatEntityRelation as t1 GROUP BY t1.relationStatus, t1.relationType ORDER BY count DESC' >> $@.temp
	@echo "----------Relations existing in all species per status and type: " >> $@.temp
	@$(MYSQL) -e 'SELECT t1.relationStatus, t1.relationType, count(distinct t1.anatEntityRelationId) AS count from anatEntityRelation as t1 inner join anatEntityRelationTaxonConstraint as t2 on t1.anatEntityRelationId = t2.anatEntityRelationId where t2.speciesId is null GROUP BY t1.relationStatus, t1.relationType ORDER BY count DESC' >> $@.temp
	@echo "----------Relations existing in some species per status and type: " >> $@.temp
	@$(MYSQL) -e 'SELECT t1.relationStatus, t1.relationType, count(distinct t1.anatEntityRelationId) AS count from anatEntityRelation as t1 inner join anatEntityRelationTaxonConstraint as t2 on t1.anatEntityRelationId = t2.anatEntityRelationId where t2.speciesId is not null GROUP BY t1.relationStatus, t1.relationType ORDER BY count DESC' >> $@.temp
	@echo "----------Number of relations defined by taxon constraints: " >> $@.temp
	@$(MYSQL) -e "SELECT (SELECT t6.taxonScientificName FROM taxon AS t6 where t6.taxonLeftBound <= temp1.minLeftBound AND t6.taxonRightBound >= temp1.maxRightBound ORDER BY t6.taxonLeftBound DESC LIMIT 1) AS commonAncestor, temp1.species_list, temp1.speciesCount, temp1.count AS relationCount FROM (SELECT t3.species_list, t3.speciesCount, count(t3.anatEntityRelationId) AS count, MIN(taxonLeftBound) as minLeftBound, MAX(taxonRightBound) AS maxRightBound FROM (SELECT t1.anatEntityRelationId, MIN(t5.taxonLeftBound) AS taxonLeftBound, MAX(t5.taxonRightBound) AS taxonRightBound, group_concat(IF(t1.speciesId IS NULL, 'All species', IF (t2.speciesCommonName != '', t2.speciesCommonName, CONCAT_WS(' ', t2.genus, t2.species))) ORDER BY IF(t1.speciesId IS NULL, 'All species', IF (t2.speciesCommonName != '', t2.speciesCommonName, CONCAT_WS(' ', t2.genus, t2.species)))) AS species_list, COUNT(DISTINCT t2.speciesId) AS speciesCount from anatEntityRelationTaxonConstraint as t1 left outer join species as t2 on t1.speciesId = t2.speciesId LEFT OUTER JOIN taxon AS t5 ON t2.taxonId = t5.taxonId group by t1.anatEntityRelationId) AS t3 group by t3.species_list) AS temp1 order by commonAncestor ASC, speciesCount DESC, count DESC" >> $@.temp
	@$(MV) $@.temp $@

$(OUTPUT_DIR)detect_cycles: $(OUTPUT_DIR)insert_anatomy
	@perl check_cycles.pl -bgee=$(BGEECMD) > $@.tmp
	@$(MV) $@.tmp $@

delete_anatomy_from_db:
	@$(MYSQL) -e "DELETE FROM anatEntityRelationTaxonConstraint; DELETE FROM anatEntityRelation; DELETE FROM anatEntityNameSynonym; DELETE FROM anatEntityXRef; DELETE FROM anatEntityTaxonConstraint; DELETE FROM anatEntity;"

delete_stage_from_db:
	@$(MYSQL) -e "DELETE FROM stageNameSynonym; DELETE FROM stageXRef; DELETE FROM stageTaxonConstraint; DELETE FROM stage;"

$(VERIFICATIONFILE): $(TAXONCONSTRAINTSFILEPATH) $(OUTPUT_DIR)insert_stages $(OUTPUT_DIR)insert_anatomy $(OUTPUT_DIR)detect_cycles $(UBERON_SEX_INFO_FILE_PATH)
# 	@$(MYSQL) -e "SELECT * FROM taxon order by taxonLeftBound" >> $@.temp
#	@$(MV) $@.temp $@

clean:
	-@$(RM) $(UBERONFILEPATH) $(DEV_STAGE_ONT_FILE_PATH) $(UBERON_COMPOSITE_FILE_PATH) $(UBERON_EXT_FILE_PATH) $(TAXONCONSTRAINTSFILEPATH).tmp $(TAXONCONSTRAINTSFILEPATH) $(OUTPUT_DIR)/intermediate_ontologies $(VERIFICATIONFILE).tmp $(VERIFICATIONFILE)
	-@$(RM) $(OUTPUT_DIR)insert_stages $(OUTPUT_DIR)insert_anatomy $(OUTPUT_DIR)detect_cycles $(OUTPUT_DIR)detect_cycles.tmp

