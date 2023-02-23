-- this file contains the foreign key constraints.

-- ****************************************************
-- GENERAL
-- ****************************************************
/*!40000 ALTER TABLE `dataSourceToSpecies` DISABLE KEYS */;
alter table dataSourceToSpecies
add foreign key (dataSourceId) references dataSource(dataSourceId) on delete cascade,
add foreign key (speciesId) references species(speciesId) on delete cascade;
/*!40000 ALTER TABLE `dataSourceToSpecies` ENABLE KEYS */;

--  ****************************************************
--  TAXONOMY
--  ****************************************************

/*!40000 ALTER TABLE `species` DISABLE KEYS */;
alter table species
add foreign key (taxonId) references taxon(taxonId) on delete cascade,
add foreign key (dataSourceId) references dataSource(dataSourceId) on delete cascade;
/*!40000 ALTER TABLE `species` ENABLE KEYS */;

/*!40000 ALTER TABLE `speciesToSex` DISABLE KEYS */;
alter table speciesToSex
add foreign key (speciesId) references species(speciesId) on delete cascade;
/*!40000 ALTER TABLE `speciesToSex` ENABLE KEYS */;

/*!40000 ALTER TABLE `speciesToKeyword` DISABLE KEYS */;
alter table speciesToKeyword
add foreign key (speciesId) references species(speciesId) on delete cascade,
add foreign key (keywordId) references keyword(keywordId) on delete cascade;
/*!40000 ALTER TABLE `speciesToKeyword` ENABLE KEYS */;

--  ****************************************************
--  CONFIDENCE AND EVIDENCE ONTOLOGIES
--  ****************************************************

--  ****************************************************
--  ANATOMY AND DEVELOPMENT
--  ****************************************************
/*!40000 ALTER TABLE `stageTaxonConstraint` DISABLE KEYS */;
alter table stageTaxonConstraint
add foreign key (stageId) references stage(stageId) on delete cascade,
add foreign key (speciesId) references species(speciesId) on delete cascade;
/*!40000 ALTER TABLE `stageTaxonConstraint` ENABLE KEYS */;

/*!40000 ALTER TABLE `stageNameSynonym` DISABLE KEYS */;
alter table stageNameSynonym
add foreign key (stageId) references stage(stageId) on delete cascade;
/*!40000 ALTER TABLE `stageNameSynonym` ENABLE KEYS */;

/*!40000 ALTER TABLE `stageXRef` DISABLE KEYS */;
alter table stageXRef
add foreign key (stageId) references stage(stageId) on delete cascade;
/*!40000 ALTER TABLE `stageXRef` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntity` DISABLE KEYS */;
alter table anatEntity
add foreign key (startStageId) references stage(stageId),
add foreign key (endStageId) references stage(stageId);
/*!40000 ALTER TABLE `anatEntity` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntityTaxonConstraint` DISABLE KEYS */;
alter table anatEntityTaxonConstraint
add foreign key (anatEntityId) references anatEntity(anatEntityId) on delete cascade,
add foreign key (speciesId) references species(speciesId) on delete cascade;
/*!40000 ALTER TABLE `anatEntityTaxonConstraint` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntityNameSynonym` DISABLE KEYS */;
alter table anatEntityNameSynonym
add foreign key (anatEntityId) references anatEntity(anatEntityId) on delete cascade;
/*!40000 ALTER TABLE `anatEntityNameSynonym` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntityXRef` DISABLE KEYS */;
alter table anatEntityXRef
add foreign key (anatEntityId) references anatEntity(anatEntityId) on delete cascade;
/*!40000 ALTER TABLE `anatEntityXRef` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntityRelation` DISABLE KEYS */;
alter table anatEntityRelation
add foreign key (anatEntitySourceId) references anatEntity(anatEntityId) on delete cascade,
add foreign key (anatEntityTargetId) references anatEntity(anatEntityId) on delete cascade;
/*!40000 ALTER TABLE `anatEntityRelation` ENABLE KEYS */;

/*!40000 ALTER TABLE `anatEntityRelationTaxonConstraint` DISABLE KEYS */;
alter table anatEntityRelationTaxonConstraint
add foreign key (anatEntityRelationId) references anatEntityRelation(anatEntityRelationId) on delete cascade,
add foreign key (speciesId) references species(speciesId) on delete cascade;
/*!40000 ALTER TABLE `anatEntityRelationTaxonConstraint` ENABLE KEYS */;

--  ****************************************************
--  SIMILARITY ANNOTATIONS
--  ****************************************************
/*!40000 ALTER TABLE `summarySimilarityAnnotation` DISABLE KEYS */;
alter table summarySimilarityAnnotation
add foreign key (taxonId) references taxon(taxonId) on delete cascade,
add foreign key (CIOId) references CIOStatement(CIOId) on delete cascade;
/*!40000 ALTER TABLE `summarySimilarityAnnotation` ENABLE KEYS */;

/*!40000 ALTER TABLE `similarityAnnotationToAnatEntityId` DISABLE KEYS */;
alter table similarityAnnotationToAnatEntityId
add foreign key (summarySimilarityAnnotationId) references summarySimilarityAnnotation(summarySimilarityAnnotationId) on delete cascade,
add foreign key (anatEntityId) references anatEntity(anatEntityId) on delete cascade;
/*!40000 ALTER TABLE `similarityAnnotationToAnatEntityId` ENABLE KEYS */;

/*!40000 ALTER TABLE `rawSimilarityAnnotation` DISABLE KEYS */;
alter table rawSimilarityAnnotation
add foreign key (summarySimilarityAnnotationId) references summarySimilarityAnnotation(summarySimilarityAnnotationId) on delete cascade,
add foreign key (ECOId) references evidenceOntology(ECOId) on delete cascade,
add foreign key (CIOId) references CIOStatement(CIOId) on delete cascade;
/*!40000 ALTER TABLE `rawSimilarityAnnotation` ENABLE KEYS */;

-- ****************************************************
-- GENE AND TRANSCRIPT INFO
-- ****************************************************
/*!40000 ALTER TABLE `OMAHierarchicalGroup` DISABLE KEYS */;
alter table OMAHierarchicalGroup
add foreign key (taxonId) references taxon(taxonId) on delete set null;
/*!40000 ALTER TABLE `OMAHierarchicalGroup` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneOntologyTermAltId` DISABLE KEYS */;
alter table geneOntologyTermAltId
add foreign key (goId) references geneOntologyTerm(goId) on delete cascade;
/*!40000 ALTER TABLE `geneOntologyTermAltId` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneOntologyRelation` DISABLE KEYS */;
alter table geneOntologyRelation
add foreign key (goAllTargetId) references geneOntologyTerm(goId) on delete cascade,
add foreign key (goAllSourceId) references geneOntologyTerm(goId) on delete cascade;
/*!40000 ALTER TABLE `geneOntologyRelation` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneOrthologs` DISABLE KEYS */;
alter table geneOrthologs
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (targetGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (taxonId) references taxon(taxonId) on delete cascade;
/*!40000 ALTER TABLE `geneOrthologs` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneParalogs` DISABLE KEYS */;
alter table geneParalogs
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (targetGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (taxonId) references taxon(taxonId) on delete cascade;
/*!40000 ALTER TABLE `geneParalogs` ENABLE KEYS */;

/*!40000 ALTER TABLE `gene` DISABLE KEYS */;
alter table gene
add foreign key (speciesId) references species(speciesId) on delete cascade,
add foreign key (geneBioTypeId) references geneBioType(geneBioTypeId) on delete set null,
add foreign key (OMAParentNodeId) references OMAHierarchicalGroup(OMANodeId) on delete set null;
/*!40000 ALTER TABLE `gene` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneNameSynonym` DISABLE KEYS */;
alter table geneNameSynonym
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade;
/*!40000 ALTER TABLE `geneNameSynonym` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneXRef` DISABLE KEYS */;
alter table geneXRef
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (dataSourceId) references dataSource(dataSourceId) on delete cascade;
/*!40000 ALTER TABLE `geneXRef` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneToTerm` DISABLE KEYS */;
alter table geneToTerm
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade;
/*!40000 ALTER TABLE `geneToTerm` ENABLE KEYS */;

/*!40000 ALTER TABLE `geneToGeneOntologyTerm` DISABLE KEYS */;
alter table geneToGeneOntologyTerm
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (goId) references geneOntologyTerm(goId) on delete cascade;
/*!40000 ALTER TABLE `geneToGeneOntologyTerm` ENABLE KEYS */;

/*!40000 ALTER TABLE `transcript` DISABLE KEYS */;
alter table transcript
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade;
/*!40000 ALTER TABLE `transcript` ENABLE KEYS */;

-- ****************************************************
-- CONDITIONS
-- ****************************************************
/*!40000 ALTER TABLE `cond` DISABLE KEYS */;
alter table cond
add foreign key (exprMappedConditionId) references cond(conditionId) on delete cascade,
add foreign key (anatEntityId) references anatEntity(anatEntityId) on delete cascade,
add foreign key (cellTypeId) references anatEntity(anatEntityId) on delete cascade,
add foreign key (stageId) references stage(stageId) on delete cascade,
add foreign key (speciesId) references species(speciesId) on delete cascade;
/*!40000 ALTER TABLE `cond` ENABLE KEYS */;

/*!40000 ALTER TABLE `remapCond` DISABLE KEYS */;
alter table remapCond
add foreign key (remappedConditionId) references cond(conditionId) on delete cascade;
/*!40000 ALTER TABLE `remapCond` ENABLE KEYS */;

/*!40000 ALTER TABLE `remapExpression` DISABLE KEYS */;
alter table remapExpression
add foreign key (remappedExpressionId) references expression(expressionId) on delete cascade;
/*!40000 ALTER TABLE `remapExpression` ENABLE KEYS */;

/*!40000 ALTER TABLE `globalCond` DISABLE KEYS */;
alter table globalCond
add foreign key (anatEntityId) references anatEntity(anatEntityId) on delete cascade,
add foreign key (cellTypeId) references anatEntity(anatEntityId) on delete cascade,
add foreign key (stageId) references stage(stageId) on delete cascade,
add foreign key (speciesId) references species(speciesId) on delete cascade;
/*!40000 ALTER TABLE `globalCond` ENABLE KEYS */;

/*!40000 ALTER TABLE `globalCondToCond` DISABLE KEYS */;
alter table globalCondToCond
add foreign key (conditionId) references cond(conditionId) on delete cascade,
add foreign key (globalConditionId) references globalCond(globalConditionId) on delete cascade;
/*!40000 ALTER TABLE `globalCondToCond` ENABLE KEYS */;

-- ****************************************************
-- EXPRESSION DATA
-- ****************************************************
/*!40000 ALTER TABLE `expression` DISABLE KEYS */;
alter table expression
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (conditionId) references cond(conditionId) on delete cascade;
/*!40000 ALTER TABLE `expression` ENABLE KEYS */;

/*!40000 ALTER TABLE `globalExpression` DISABLE KEYS */;
alter table globalExpression
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (globalConditionId) references globalCond(globalConditionId) on delete cascade
-- disable this foreign key constraint, this creates way too many indexes that are not really needed
-- ,
-- add foreign key(gCIdPValBDAffy) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDEst) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDInSitu) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDRnaSeq) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyEst) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyInSitu) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyRnaSeq) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDEstInSitu) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDEstRnaSeq) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDEstScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDInSituRnaSeq) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDInSituScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDRnaSeqScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyEstInSitu) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyEstRnaSeq) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyEstScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyInSituRnaSeq) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyInSituScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyRnaSeqScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDEstInSituRnaSeq) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDEstInSituScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDEstRnaSeqScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDInSituRnaSeqScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyEstInSituRnaSeq) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyEstInSituScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyEstRnaSeqScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyInSituRnaSeqScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDEstInSituRnaSeqScRnaSeqFL) references globalCond(globalConditionId) on delete set null,
-- add foreign key(gCIdPValBDAffyEstInSituRnaSeqScRnaSeqFL) references globalCond(globalConditionId) on delete set null
;
/*!40000 ALTER TABLE `globalExpression` ENABLE KEYS */;

-- ****************************************************
-- RAW EST DATA
-- ****************************************************
/*!40000 ALTER TABLE `estLibrary` DISABLE KEYS */;
alter table estLibrary
add foreign key (conditionId) references cond(conditionId) on delete cascade,
add foreign key (dataSourceId) references dataSource(dataSourceId);
/*!40000 ALTER TABLE `estLibrary` ENABLE KEYS */;

/*!40000 ALTER TABLE `estLibraryToKeyword` DISABLE KEYS */;
alter table estLibraryToKeyword
add foreign key (estLibraryId) references estLibrary(estLibraryId) on delete cascade,
add foreign key (keywordId) references keyword(keywordId) on delete cascade;
/*!40000 ALTER TABLE `estLibraryToKeyword` ENABLE KEYS */;

/*!40000 ALTER TABLE `expressedSequenceTag` DISABLE KEYS */;
alter table expressedSequenceTag
add foreign key (estLibraryId) references estLibrary(estLibraryId) on delete cascade,
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (expressionId) references expression(expressionId) on delete set null;
/*!40000 ALTER TABLE `expressedSequenceTag` ENABLE KEYS */;

/*!40000 ALTER TABLE `estLibraryExpression` DISABLE KEYS */;
alter table estLibraryExpression
add foreign key (expressionId) references expression(expressionId) on delete cascade,
add foreign key (estLibraryId) references estLibrary(estLibraryId) on delete cascade;
/*!40000 ALTER TABLE `estLibraryExpression` ENABLE KEYS */;
--  ****************************************************
--  RAW AFFYMETRIX DATA
--  ****************************************************
/*!40000 ALTER TABLE `microarrayExperiment` DISABLE KEYS */;
alter table microarrayExperiment
add foreign key (dataSourceId) references dataSource(dataSourceId);
/*!40000 ALTER TABLE `microarrayExperiment` ENABLE KEYS */;

/*!40000 ALTER TABLE `microarrayExperimentToKeyword` DISABLE KEYS */;
alter table microarrayExperimentToKeyword
add foreign key (microarrayExperimentId) references microarrayExperiment(microarrayExperimentId) on delete cascade,
add foreign key (keywordId) references keyword(keywordId) on delete cascade;
/*!40000 ALTER TABLE `microarrayExperimentToKeyword` ENABLE KEYS */;

/*!40000 ALTER TABLE `affymetrixChip` DISABLE KEYS */;
alter table affymetrixChip
add foreign key (microarrayExperimentId) references microarrayExperiment(microarrayExperimentId) on delete cascade,
add foreign key (chipTypeId) references chipType(chipTypeId) on delete set null,
add foreign key (conditionId) references cond(conditionId) on delete cascade;
/*!40000 ALTER TABLE `affymetrixChip` ENABLE KEYS */;

/*!40000 ALTER TABLE `affymetrixProbeset` DISABLE KEYS */;
alter table affymetrixProbeset
add foreign key (bgeeAffymetrixChipId) references affymetrixChip(bgeeAffymetrixChipId) on delete cascade,
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (expressionId) references expression(expressionId) on delete set null;
/*!40000 ALTER TABLE `affymetrixProbeset` ENABLE KEYS */;

/*!40000 ALTER TABLE `microarrayExperimentExpression` DISABLE KEYS */;
alter table microarrayExperimentExpression
add foreign key (expressionId) references expression(expressionId) on delete cascade,
add foreign key (microarrayExperimentId) references microarrayExperiment(microarrayExperimentId) on delete cascade;
/*!40000 ALTER TABLE `microarrayExperimentExpression` ENABLE KEYS */;

--  ****************************************************
--  RAW IN SITU DATA
--  ****************************************************
/*!40000 ALTER TABLE `inSituExperiment` DISABLE KEYS */;
alter table inSituExperiment
add foreign key (dataSourceId) references dataSource(dataSourceId);
/*!40000 ALTER TABLE `inSituExperiment` ENABLE KEYS */;

/*!40000 ALTER TABLE `inSituExperimentToKeyword` DISABLE KEYS */;
alter table inSituExperimentToKeyword
add foreign key (inSituExperimentId) references inSituExperiment(inSituExperimentId) on delete cascade,
add foreign key (keywordId) references keyword(keywordId) on delete cascade;
/*!40000 ALTER TABLE `inSituExperimentToKeyword` ENABLE KEYS */;

/*!40000 ALTER TABLE `inSituEvidence` DISABLE KEYS */;
alter table inSituEvidence
add foreign key (inSituExperimentId) references inSituExperiment(inSituExperimentId) on delete cascade;
/*!40000 ALTER TABLE `inSituEvidence` ENABLE KEYS */;

/*!40000 ALTER TABLE `inSituSpot` DISABLE KEYS */;
alter table inSituSpot
add foreign key (inSituEvidenceId) references inSituEvidence(inSituEvidenceId) on delete cascade,
add foreign key (conditionId) references cond(conditionId) on delete cascade,
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (expressionId) references expression(expressionId) on delete set null;
/*!40000 ALTER TABLE `inSituSpot` ENABLE KEYS */;

/*!40000 ALTER TABLE `inSituExperimentExpression` DISABLE KEYS */;
alter table inSituExperimentExpression
add foreign key (expressionId) references expression(expressionId) on delete cascade,
add foreign key (inSituExperimentId) references inSituExperiment(inSituExperimentId) on delete cascade;
/*!40000 ALTER TABLE `inSituExperimentExpression` ENABLE KEYS */;

--  ****************************************************
--  RAW RNA-SEQ DATA
--  ****************************************************
/*!40000 ALTER TABLE `rnaSeqExperiment` DISABLE KEYS */;
alter table rnaSeqExperiment
add foreign key (dataSourceId) references dataSource(dataSourceId);
/*!40000 ALTER TABLE `rnaSeqExperiment` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqExperimentToKeyword` DISABLE KEYS */;
alter table rnaSeqExperimentToKeyword
add foreign key (rnaSeqExperimentId) references rnaSeqExperiment(rnaSeqExperimentId) on delete cascade,
add foreign key (keywordId) references keyword(keywordId) on delete cascade;
/*!40000 ALTER TABLE `rnaSeqExperimentToKeyword` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqLibrary` DISABLE KEYS */;
alter table rnaSeqLibrary
add foreign key (rnaSeqExperimentId) references rnaSeqExperiment(rnaSeqExperimentId) on delete cascade,
add foreign key (rnaSeqPopulationCaptureId) references rnaSeqPopulationCapture(rnaSeqPopulationCaptureId) on delete cascade;
/*!40000 ALTER TABLE `rnaSeqLibrary` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqRun` DISABLE KEYS */;
alter table rnaSeqRun
add foreign key (rnaSeqLibraryId) references rnaSeqLibrary(rnaSeqLibraryId) on delete cascade;
/*!40000 ALTER TABLE `rnaSeqRun` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqLibraryAnnotatedSample` DISABLE KEYS */;
alter table rnaSeqLibraryAnnotatedSample
add foreign key (rnaSeqLibraryId) references rnaSeqLibrary(rnaSeqLibraryId) on delete cascade,
add foreign key (conditionId) references cond(conditionId) on delete cascade;
/*!40000 ALTER TABLE `rnaSeqLibraryAnnotatedSample` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqLibraryAnnotatedSampleGeneResult` DISABLE KEYS */;
alter table rnaSeqLibraryAnnotatedSampleGeneResult
add foreign key (rnaSeqLibraryAnnotatedSampleId) references rnaSeqLibraryAnnotatedSample(rnaSeqLibraryAnnotatedSampleId) on delete cascade,
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade,
add foreign key (expressionId) references expression(expressionId) on delete set null;
/*!40000 ALTER TABLE `rnaSeqLibraryAnnotatedSampleGeneResult` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqLibraryIndividualSample` DISABLE KEYS */;
alter table rnaSeqLibraryIndividualSample
add foreign key (rnaSeqLibraryAnnotatedSampleId) references rnaSeqLibraryAnnotatedSample(rnaSeqLibraryAnnotatedSampleId) on delete cascade;
/*!40000 ALTER TABLE `rnaSeqLibraryIndividualSample` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqLibraryIndividualSampleGeneResult` DISABLE KEYS */;
alter table rnaSeqLibraryIndividualSampleGeneResult
add foreign key (rnaSeqLibraryIndividualSampleId) references rnaSeqLibraryIndividualSample(rnaSeqLibraryIndividualSampleId) on delete cascade,
add foreign key (bgeeGeneId) references gene(bgeeGeneId) on delete cascade;
/*!40000 ALTER TABLE `rnaSeqLibraryIndividualSampleGeneResult` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls` DISABLE KEYS */;
alter table rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls
add foreign key (rnaSeqPopulationCaptureId) references rnaSeqPopulationCapture(rnaSeqPopulationCaptureId) on delete cascade,
add foreign key (geneBioTypeId) references geneBioType(geneBioTypeId) on delete cascade;
/*!40000 ALTER TABLE `rnaSeqPopulationCaptureToBiotypeExcludedAbsentCalls` ENABLE KEYS */;

/*!40000 ALTER TABLE `rnaSeqPopulationCaptureSpeciesMaxRank` DISABLE KEYS */;
alter table rnaSeqPopulationCaptureSpeciesMaxRank
add foreign key (rnaSeqPopulationCaptureId) references rnaSeqPopulationCapture(rnaSeqPopulationCaptureId) on delete cascade,
add foreign key (speciesId) references species(speciesId) on delete cascade;
/*!40000 ALTER TABLE `rnaSeqPopulationCaptureSpeciesMaxRank` ENABLE KEYS */;

/*!40000 ALTER TABLE `downloadFile` DISABLE KEYS */;
alter table downloadFile
add foreign key (speciesDataGroupId) references speciesDataGroup(speciesDataGroupId) on delete cascade;
/*!40000 ALTER TABLE `downloadFile` ENABLE KEYS */;

/*!40000 ALTER TABLE `speciesToDataGroup` DISABLE KEYS */;
alter table speciesToDataGroup
add foreign key (speciesId) references species(speciesId) on delete cascade,
add foreign key (speciesDataGroupId) references speciesDataGroup(speciesDataGroupId) on delete cascade;
/*!40000 ALTER TABLE `speciesToDataGroup` ENABLE KEYS */;


