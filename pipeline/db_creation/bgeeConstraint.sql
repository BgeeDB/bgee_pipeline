-- This file contains the primary key definitions, and other constraints
-- such as unique indexes. Indexes used solely for performance issues are not
-- defined in this file, but in bgeeIndex.sql.

--  ****************************************************
--  GENERAL
--  ****************************************************
/*!40000 ALTER TABLE `author` DISABLE KEYS */;
alter table author
modify authorId smallInt unsigned not null auto_increment primary key;
/*!40000 ALTER TABLE `author` ENABLE KEYS */;

/*!40000 ALTER TABLE `dataSource` DISABLE KEYS */;
alter table dataSource
modify dataSourceId smallInt unsigned not null auto_increment primary key;
/*!40000 ALTER TABLE `dataSource` ENABLE KEYS */;

/*!40000 ALTER TABLE `dataSourceToSpecies` DISABLE KEYS */;
alter table dataSourceToSpecies 
add primary key(speciesId, dataSourceId, dataType, infoType);
/*!40000 ALTER TABLE `dataSourceToSpecies` ENABLE KEYS */;

/*!40000 ALTER TABLE `keyword` DISABLE KEYS */;
alter table keyword
modify keywordId int unsigned not null auto_increment primary key,
add unique(keyword);
/*!40000 ALTER TABLE `keyword` ENABLE KEYS */;

--  ****************************************************
--  TAXONOMY
--  ****************************************************
/*!40000 ALTER TABLE `taxon` DISABLE KEYS */;
alter table taxon
add primary key(taxonId),
add unique(taxonLeftBound),
add unique(taxonRightBound);
/*!40000 ALTER TABLE `taxon` ENABLE KEYS */;

/*!40000 ALTER TABLE `species` DISABLE KEYS */;
alter table species 
add primary key(speciesId), 
add unique(speciesDisplayOrder), 
add unique(species, genus);
/*!40000 ALTER TABLE `species` ENABLE KEYS */;

/*!40000 ALTER TABLE `speciesToSex` DISABLE KEYS */;
alter table speciesToSex
add primary key (speciesId, sex);
/*!40000 ALTER TABLE `speciesToSex` ENABLE KEYS */;

/*!40000 ALTER TABLE `speciesToKeyword` DISABLE KEYS */;
alter table speciesToKeyword
add primary key (speciesId, keywordId);
/*!40000 ALTER TABLE `speciesToKeyword` ENABLE KEYS */;

--  ****************************************************
--  CONFIDENCE AND EVIDENCE ONTOLOGIES
--  ****************************************************
/*!40000 ALTER TABLE `CIOStatement` DISABLE KEYS */;
alter table CIOStatement
add primary key(CIOId);
/*!40000 ALTER TABLE `CIOStatement` ENABLE KEYS */;

/*!40000 ALTER TABLE `evidenceOntology` DISABLE KEYS */;
alter table evidenceOntology
add primary key(ECOId);
/*!40000 ALTER TABLE `evidenceOntology` ENABLE KEYS */;

--  ****************************************************
--  ANATOMY AND DEVELOPMENT
--  ****************************************************
/*!40000 ALTER TABLE `stage` DISABLE KEYS */;
alter table stage
add primary key(stageId),
add unique(stageLeftBound),
add unique(stageRightBound);
/*!40000 ALTER TABLE `stage` ENABLE KEYS */;

/*!40000 ALTER TABLE `stageTaxonConstraint` DISABLE KEYS */;
alter table stageTaxonConstraint
-- not a primary key because speciesId can be null
add unique(stageId, speciesId);
/*!40000 ALTER TABLE `stageTaxonConstraint` ENABLE KEYS */;

/*!40000 ALTER TABLE `stageNameSynonym` DISABLE KEYS */;
alter table stageNameSynonym
add primary key (stageNameSynonym, stageId);
/*!40000 ALTER TABLE `stageNameSynonym` ENABLE KEYS */;

/*!40000 ALTER TABLE `stageXRef` DISABLE KEYS */;
alter table stageXRef
add primary key (stageXRefId, stageId);
/*!40000 ALTER TABLE `stageXRef` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntity` DISABLE KEYS */;
alter table anatEntity
add primary key(anatEntityId);
/*!40000 ALTER TABLE `anatEntity` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntityTaxonConstraint` DISABLE KEYS */;
alter table anatEntityTaxonConstraint
-- not a primary key because speciesId can be null
add unique(anatEntityId, speciesId);
/*!40000 ALTER TABLE `anatEntityTaxonConstraint` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntityNameSynonym` DISABLE KEYS */;
alter table anatEntityNameSynonym
add primary key (anatEntityNameSynonym, anatEntityId);
/*!40000 ALTER TABLE `anatEntityNameSynonym` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntityXRef` DISABLE KEYS */;
alter table anatEntityXRef
add primary key (anatEntityXRefId, anatEntityId);
/*!40000 ALTER TABLE `anatEntityXRef` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntityRelation` DISABLE KEYS */;
-- See https://stackoverflow.com/a/42822691/1768736 for motivation
-- for this design of PK and UNIQUE indexes
alter table anatEntityRelation
modify anatEntityRelationId int unsigned not null auto_increment,
-- we allow a same relation with same source and target, but different status (direct/indirect)
-- because they can have different taxon constraints.
add primary key (anatEntitySourceId, anatEntityTargetId, relationType, relationStatus),
add unique(anatEntityRelationId);
/*!40000 ALTER TABLE `anatEntityRelation` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntityRelationTaxonConstraint` DISABLE KEYS */;
alter table anatEntityRelationTaxonConstraint
-- not a primary key because speciesId can be null
add unique(anatEntityRelationId, speciesId);
/*!40000 ALTER TABLE `anatEntityRelationTaxonConstraint` ENABLE KEYS */;

--  ****************************************************
--  SIMILARITY ANNOTATIONS
--  ****************************************************
/*!40000 ALTER TABLE `summarySimilarityAnnotation` DISABLE KEYS */;
alter table summarySimilarityAnnotation
modify summarySimilarityAnnotationId mediumint unsigned not null auto_increment primary key;
/*!40000 ALTER TABLE `summarySimilarityAnnotation` ENABLE KEYS */;

/*!40000 ALTER TABLE `similarityAnnotationToAnatEntityId` DISABLE KEYS */;
alter table similarityAnnotationToAnatEntityId
add primary key (summarySimilarityAnnotationId, anatEntityId);
/*!40000 ALTER TABLE `similarityAnnotationToAnatEntityId` ENABLE KEYS */;

/*!40000 ALTER TABLE `rawSimilarityAnnotation` DISABLE KEYS */;
alter table rawSimilarityAnnotation
-- It is theoretically possible to have an annotation supported by several evidence lines
-- with same ECO ID, CIO ID, ref ID, and meaning (positive or negative), but it is
-- so unlikely that we make a constraint on these fields for now.
add primary key (summarySimilarityAnnotationId, negated, ECOId, CIOId, referenceId, supportingText(255));
/*!40000 ALTER TABLE `rawSimilarityAnnotation` ENABLE KEYS */;

-- ****************************************************
-- GENE AND TRANSCRIPT INFO
-- ****************************************************
/*!40000 ALTER TABLE `OMAHierarchicalGroup` DISABLE KEYS */;
alter table OMAHierarchicalGroup
-- The PK is only OMANodeId, but we also add OMAGroupId for queries using both fields,
-- it's better to use the clustered index in such cases.
-- TODO: to reevaluate when used, do we have any query requesting both?
-- If only OMAGroupId is mainly used, switch the PK to primary key (OMAGroupId, OMANodeId),
-- and remove index on OMAGroupId from bgeeIndex.sql.
-- If only OMANodeId is mainly used, switch the PK to primary key (OMANodeId),
-- and remove the unique index on OMANodeId?
add primary key (OMANodeId, OMAGroupId),
-- To check the PK, we also add a unique index on OMANodeId
add unique (OMANodeId),
add unique (OMANodeLeftBound),
add unique (OMANodeRightBound);
/*!40000 ALTER TABLE `OMAHierarchicalGroup` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneOntologyTerm` DISABLE KEYS */;
alter table geneOntologyTerm
add primary key (goId);
/*!40000 ALTER TABLE `geneOntologyTerm` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneOntologyTermAltId` DISABLE KEYS */;
alter table geneOntologyTermAltId
add primary key (goAltId, goId);
/*!40000 ALTER TABLE `geneOntologyTermAltId` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneOntologyRelation` DISABLE KEYS */;
alter table geneOntologyRelation
add primary key(goAllTargetId, goAllSourceId);
/*!40000 ALTER TABLE `geneOntologyRelation` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneBioType` DISABLE KEYS */;
alter table geneBioType
modify geneBioTypeId smallint unsigned not null auto_increment primary key,
add unique (geneBioTypeName);
/*!40000 ALTER TABLE `geneBioType` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneOrthologs` DISABLE KEYS */;
alter table geneOrthologs
add primary key(bgeeGeneId, targetGeneId);
/*!40000 ALTER TABLE `geneOrthologs` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneParalogs` DISABLE KEYS */;
alter table geneParalogs
add primary key(bgeeGeneId, targetGeneId);
/*!40000 ALTER TABLE `geneParalogs` ENABLE KEYS */;

/*!40000 ALTER TABLE `gene` DISABLE KEYS */;
alter table gene
modify bgeeGeneId mediumint unsigned not null auto_increment primary key,
add unique(geneId, speciesId);
/*!40000 ALTER TABLE `gene` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneNameSynonym` DISABLE KEYS */;
alter table geneNameSynonym
add primary key (geneNameSynonym, bgeeGeneId);
/*!40000 ALTER TABLE `geneNameSynonym` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneXRef` DISABLE KEYS */;
alter table geneXRef
add primary key (XRefId, bgeeGeneId, dataSourceId);
/*!40000 ALTER TABLE `geneXRef` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneToTerm` DISABLE KEYS */;
alter table geneToTerm
add primary key (term, bgeeGeneId);
/*!40000 ALTER TABLE `geneToTerm` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneToGeneOntologyTerm` DISABLE KEYS */;
alter table geneToGeneOntologyTerm
add primary key (bgeeGeneId, goId);
/*!40000 ALTER TABLE `geneToGeneOntologyTerm` ENABLE KEYS */;

/*!40000 ALTER TABLE `transcript` DISABLE KEYS */;
alter table transcript
modify bgeeTranscriptId int unsigned not null auto_increment primary key,
-- transcriptId is not unique because we use closely-related genomes for some species, 
-- e.g., chimp genome for bonobo. So a same transcript can be used in different species, 
-- associated with different genes.
-- Actually, a 'unique (transcriptId, speciesId)' would be more correct, 
-- but the speciesId is not stored in this table, using a trigger is too expensive performance-wise.
add unique(transcriptId, bgeeGeneId);
/*!40000 ALTER TABLE `transcript` ENABLE KEYS */;

-- ****************************************************
-- CONDITIONS
-- ****************************************************
/*!40000 ALTER TABLE `cond` DISABLE KEYS */;
-- See https://stackoverflow.com/a/42822691/1768736 for motivation
-- for this design of PK and UNIQUE indexes
alter table cond
modify conditionId mediumint unsigned not null auto_increment,
add primary key(conditionId),
add unique(anatEntityId, cellTypeId, stageId, speciesId, sex, sexInferred, strain);
/*!40000 ALTER TABLE `cond` ENABLE KEYS */;

/*!40000 ALTER TABLE `remapCond` DISABLE KEYS */;
alter table remapCond
add primary key(incorrectConditionId, remappedConditionId);
/*!40000 ALTER TABLE `remapCond` ENABLE KEYS */;

/*!40000 ALTER TABLE `remapExpression` DISABLE KEYS */;
alter table remapExpression
add primary key(incorrectExpressionId, remappedExpressionId);
/*!40000 ALTER TABLE `remapExpression` ENABLE KEYS */;

/*!40000 ALTER TABLE `globalCond` DISABLE KEYS */;
alter table globalCond
modify globalConditionId mediumint unsigned not null auto_increment primary key,
-- not a primary key as for table cond, because some field can be null
add unique(anatEntityId, cellTypeId, stageId, speciesId, sex, strain);
/*!40000 ALTER TABLE `globalCond` ENABLE KEYS */;

/*!40000 ALTER TABLE `globalCondToCond` DISABLE KEYS */;
alter table globalCondToCond
-- we set up this primary key using conditionRelationOrigin to benefit from the clustered index
add primary key (globalConditionId, conditionRelationOrigin, conditionId),
-- but actually the unique constraint is on globalConditionId and conditionId
add unique(globalConditionId, conditionId);
/*!40000 ALTER TABLE `globalCondToCond` ENABLE KEYS */;

-- ****************************************************
-- EXPRESSION DATA
-- ****************************************************
/*!40000 ALTER TABLE `expression` DISABLE KEYS */;
-- See https://stackoverflow.com/a/42822691/1768736 for motivation
-- for this design of PK and UNIQUE indexes
alter table expression
modify expressionId int unsigned not null auto_increment,
add primary key(bgeeGeneId, conditionId),
add unique(expressionId);
/*!40000 ALTER TABLE `expression` ENABLE KEYS */;

/*!40000 ALTER TABLE `globalExpression` DISABLE KEYS */;
-- See https://stackoverflow.com/a/42822691/1768736 for motivation
-- for this design of PK and UNIQUE indexes
alter table globalExpression
-- modify globalExpressionId int unsigned not null auto_increment,
add primary key(bgeeGeneId, globalConditionId);
-- add unique(globalExpressionId);
/*!40000 ALTER TABLE `globalExpression` ENABLE KEYS */;

-- ****************************************************
-- RAW EST DATA
-- ****************************************************
/*!40000 ALTER TABLE `estLibrary` DISABLE KEYS */;
alter table estLibrary
add primary key (estLibraryId);
/*!40000 ALTER TABLE `estLibrary` ENABLE KEYS */;

/*!40000 ALTER TABLE `estLibraryToKeyword` DISABLE KEYS */;
alter table estLibraryToKeyword
add primary key (estLibraryId, keywordId);
/*!40000 ALTER TABLE `estLibraryToKeyword` ENABLE KEYS */;

/*!40000 ALTER TABLE `expressedSequenceTag` DISABLE KEYS */;
alter table expressedSequenceTag
add primary key (estId);
/*!40000 ALTER TABLE `expressedSequenceTag` ENABLE KEYS */;

/*!40000 ALTER TABLE `estLibraryExpression` DISABLE KEYS */;
alter table estLibraryExpression
add primary key (expressionId, estLibraryId);
/*!40000 ALTER TABLE `estLibraryExpression` ENABLE KEYS */;

--  ****************************************************
--  RAW AFFYMETRIX DATA
--  ****************************************************
/*!40000 ALTER TABLE `microarrayExperiment` DISABLE KEYS */;
alter table microarrayExperiment
add primary key (microarrayExperimentId);
/*!40000 ALTER TABLE `microarrayExperiment` ENABLE KEYS */;

/*!40000 ALTER TABLE `microarrayExperimentToKeyword` DISABLE KEYS */;
alter table microarrayExperimentToKeyword
add primary key (microarrayExperimentId, keywordId);
/*!40000 ALTER TABLE `microarrayExperimentToKeyword` ENABLE KEYS */;

-- chipTypes in this table don't represent affymetrix chips only
-- (could for instance represent cDNA chip)
/*!40000 ALTER TABLE `chipType` DISABLE KEYS */;
alter table chipType
add unique(chipTypeName),
add unique(cdfName),
add primary key (chipTypeId);
/*!40000 ALTER TABLE `chipType` ENABLE KEYS */;

/*!40000 ALTER TABLE `affymetrixChip` DISABLE KEYS */;
alter table affymetrixChip
modify bgeeAffymetrixChipId mediumint unsigned not null auto_increment primary key,
add unique (affymetrixChipId, microarrayExperimentId);
/*!40000 ALTER TABLE `affymetrixChip` ENABLE KEYS */;

/*!40000 ALTER TABLE `affymetrixProbeset` DISABLE KEYS */;
alter table affymetrixProbeset
add primary key (bgeeAffymetrixChipId, affymetrixProbesetId);
/*!40000 ALTER TABLE `affymetrixProbeset` ENABLE KEYS */;

/*!40000 ALTER TABLE `microarrayExperimentExpression` DISABLE KEYS */;
alter table microarrayExperimentExpression
add primary key (expressionId, microarrayExperimentId);
/*!40000 ALTER TABLE `microarrayExperimentExpression` ENABLE KEYS */;

--  ****************************************************
--  RAW IN SITU DATA
--  ****************************************************
/*!40000 ALTER TABLE `inSituExperiment` DISABLE KEYS */;
alter table inSituExperiment
add primary key (inSituExperimentId);
/*!40000 ALTER TABLE `inSituExperiment` ENABLE KEYS */;

/*!40000 ALTER TABLE `inSituExperimentToKeyword` DISABLE KEYS */;
alter table inSituExperimentToKeyword
add primary key (inSituExperimentId, keywordId);
/*!40000 ALTER TABLE `inSituExperimentToKeyword` ENABLE KEYS */;

/*!40000 ALTER TABLE `inSituEvidence` DISABLE KEYS */;
alter table inSituEvidence
add primary key (inSituEvidenceId);
/*!40000 ALTER TABLE `inSituEvidence` ENABLE KEYS */;

/*!40000 ALTER TABLE `inSituSpot` DISABLE KEYS */;
alter table inSituSpot
add primary key (inSituSpotId);
/*!40000 ALTER TABLE `inSituSpot` ENABLE KEYS */;

/*!40000 ALTER TABLE `inSituExperimentExpression` DISABLE KEYS */;
alter table inSituExperimentExpression
add primary key (expressionId, inSituExperimentId);
/*!40000 ALTER TABLE `inSituExperimentExpression` ENABLE KEYS */;

--  ****************************************************
--  RAW RNA-SEQ DATA
--  ****************************************************
/*!40000 ALTER TABLE `rnaSeqExperiment` DISABLE KEYS */;
alter table rnaSeqExperiment
add primary key (rnaSeqExperimentId);
/*!40000 ALTER TABLE `rnaSeqExperiment` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqExperimentToKeyword` DISABLE KEYS */;
alter table rnaSeqExperimentToKeyword
add primary key (rnaSeqExperimentId, keywordId);
/*!40000 ALTER TABLE `rnaSeqExperimentToKeyword` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqLibrary` DISABLE KEYS */;
alter table rnaSeqLibrary
add primary key (rnaSeqLibraryId);
/*!40000 ALTER TABLE `rnaSeqLibrary` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqLibraryDiscarded` DISABLE KEYS */;
alter table rnaSeqLibraryDiscarded
add primary key (rnaSeqLibraryId);
/*!40000 ALTER TABLE `rnaSeqLibraryDiscarded` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqRun` DISABLE KEYS */;
alter table rnaSeqRun
add primary key (rnaSeqRunId);
/*!40000 ALTER TABLE `rnaSeqRun` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqLibraryAnnotatedSample` DISABLE KEYS */;
alter table rnaSeqLibraryAnnotatedSample
modify rnaSeqLibraryAnnotatedSampleId mediumint unsigned not null auto_increment primary key;
/*!40000 ALTER TABLE `rnaSeqLibraryAnnotatedSample` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqLibraryAnnotatedSampleGeneResult` DISABLE KEYS */;
alter table rnaSeqLibraryAnnotatedSampleGeneResult
add primary key (rnaSeqLibraryAnnotatedSampleId, bgeeGeneId);
/*!40000 ALTER TABLE `rnaSeqLibraryAnnotatedSampleGeneResult` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqLibraryIndividualSample` DISABLE KEYS */;
alter table rnaSeqLibraryIndividualSample
modify rnaSeqLibraryIndividualSampleId int unsigned not null auto_increment primary key;
/*!40000 ALTER TABLE `rnaSeqLibraryIndividualSample` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqLibraryIndividualSampleGeneResult` DISABLE KEYS */;
alter table rnaSeqLibraryIndividualSampleGeneResult
add primary key (rnaSeqLibraryIndividualSampleId, bgeeGeneId);
/*!40000 ALTER TABLE `rnaSeqLibraryIndividualSampleGeneResult` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqPopulationCapture` DISABLE KEYS */;
alter table rnaSeqPopulationCapture
add primary key (rnaSeqPopulationCaptureId);
/*!40000 ALTER TABLE `rnaSeqPopulationCapture` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls` DISABLE KEYS */;
alter table rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls
add primary key (rnaSeqPopulationCaptureId, geneBioTypeId);
/*!40000 ALTER TABLE `rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqPopulationCaptureSpeciesMaxRank` DISABLE KEYS */;
alter table rnaSeqPopulationCaptureSpeciesMaxRank
add primary key (speciesId, rnaSeqPopulationCaptureId);
/*!40000 ALTER TABLE `rnaSeqPopulationCaptureSpeciesMaxRank` ENABLE KEYS */;

-- *****************************************
-- DOWNLOAD FILES
-- *****************************************
/*!40000 ALTER TABLE `speciesDataGroup` DISABLE KEYS */;
alter table speciesDataGroup
modify speciesDataGroupId mediumint unsigned not null auto_increment primary key,
add unique(speciesDataGroupOrder);
/*!40000 ALTER TABLE `speciesDataGroup` ENABLE KEYS */;

/*!40000 ALTER TABLE `speciesToDataGroup` DISABLE KEYS */;
alter table speciesToDataGroup
add primary key (speciesId, speciesDataGroupId);
/*!40000 ALTER TABLE `speciesToDataGroup` ENABLE KEYS */;

/*!40000 ALTER TABLE `downloadFile` DISABLE KEYS */;
alter table downloadFile
modify downloadFileId mediumint unsigned not null auto_increment primary key;
/*!40000 ALTER TABLE `downloadFile` ENABLE KEYS */;
