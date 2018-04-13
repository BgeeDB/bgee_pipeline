-- SQL file to create the Bgee database. Primary keys and other constraints
-- such as unique indexes are defined in bgeeConstraint.sql. Indexes defined solely
-- for performance issues are defined in bgeeIndex.sql. Foreign key constraints
-- are defined in bgeeForeignKey.sql.
--
-- To load a dump into the database, you should typically do:
-- mysql -u root -p -e "create database bgee_vXX"
-- mysql -u root -p bgee_vXX < bgeeSchema.sql
-- mysql -u root -p bgee_vXX < myDumpFile.sql
-- mysql -u root -p bgee_vXX < bgeeConstraint.sql
-- mysql -u root -p bgee_vXX < bgeeIndex.sql
-- mysql -u root -p bgee_vXX < bgeeForeignKey.sql
--
-- Altering a table after data insertion, to add indexes and foreign key constraints,
-- can fail if the table is very large, with the error 1206: "ERROR 1206 (HY000):
-- The total number of locks exceeds the lock table size". To solve this problem,
-- you have to increase the buffer pool size, or you have to insert the data AFTER
-- indexes and foreign key constraints generation.
-- The foreign key insertion should be done after the indexes creation to avoid
-- generating redundant indexes (as foreign key constraints require indexes and
-- create them if needed).

ALTER DATABASE CHARACTER SET utf8 COLLATE utf8_general_ci;

-- ****************************************************
-- GENERAL
-- ****************************************************
create table author (
    authorId smallInt unsigned not null,
    authorName varchar(255) not null COMMENT 'Bgee team author names'
) engine = innodb;

create table dataSource (
    dataSourceId smallInt unsigned not null,
    dataSourceName varchar(255) not null COMMENT 'Data source name',
    XRefUrl varchar(255) not null default '' COMMENT 'URL for cross-references to data sources',
-- path to experiment for expression data sources (ArrayExpress, GEO, NCBI, in situ databases, ...)
-- parameters such as experimentId are defined by the syntax [experimentId] for instance
    experimentUrl varchar(255) not null default '' COMMENT 'URL to experiment for expression data sources',
-- path to in situ evidence for in situ databases,
-- to Affymetrix chips for affymetrix data
-- parameters such as experimentId are defined by the syntax [experimentId] for instance
    evidenceUrl varchar(255) not null default '' COMMENT 'URL to evidence for expression data sources',
-- url to the home page of the ressource
    baseUrl varchar(255) not null default '' COMMENT 'URL to the home page of data sources',
    releaseDate timestamp null COMMENT 'Date of data source used',
-- e.g.: Ensembl 67, cvs version xxx
    releaseVersion varchar(255) not null default '' COMMENT 'Version of data source used',
    dataSourceDescription TEXT COMMENT 'Description of data source',
-- to define if this dataSource should be displayed on the page listing data sources
    toDisplay boolean not null default 0 COMMENT 'Display this data source in listing data source page?',
-- a cat to organize the display
    category enum('', 'Genomics database', 'Proteomics database',
        'In situ data source', 'Affymetrix data source', 'EST data source', 'RNA-Seq data source',
        'Ontology') COMMENT 'Data source category to organize the display',
-- to organize the display. Default value is the highest value, so that this field is the last to be displayed
    displayOrder tinyint unsigned not null default 255 COMMENT 'Data source display ordering'
) engine = innodb;

create table dataSourceToSpecies (
    dataSourceId smallInt unsigned not null COMMENT 'Data source id',
    speciesId mediumint unsigned not null COMMENT 'NCBI species taxon id',
    dataType enum('affymetrix', 'est', 'in situ', 'rna-seq') not null COMMENT 'Data type',
    infoType enum('data', 'annotation') not null COMMENT 'Information type'
) engine = innodb;

create table keyword (
    keywordId int unsigned not null,
    keyword varchar(255) not null COMMENT 'Aggregate keywords seen in all tables'
) engine = innodb;


-- ****************************************************
-- TAXONOMY
-- ****************************************************

-- The NCBI taxonomy, stored as a nested set model. This does not include species,
-- that are stored in a different table. This is because, while a species is a taxon,
-- we store some additional information for them.
--
-- Only taxa that are ancestors of a species included in Bgee are stored. The column
-- "bgeeSpeciesLCA" specifies if they are moreover a least common ancestor of at least
-- two species used in Bgee. For instance: if Bgee was using zebrafish, mouse and human,
-- "Euarchontoglires" would be the most common ancestor of human and mouse,
-- and "Euteleostomi" the most common ancestor of human, mouse and zebrafish.
-- This allows to provide a simplified display to the users, where only these relevant
-- branchings are used.
-- We neverthless also store all the ancestors of the species used in Bgee
-- (for instance, we would still store "Eutheria", "Theria", etc), as they are used
-- for the gene hierarchical groups, in case users want to have a finer control
-- on the paralogous/orthologous genes to retrieve and compare, and also for
-- the transitive evolutionary relations, in case users want to specify in which
-- common ancestor a structure should have existed.
create table taxon (
    taxonId mediumint unsigned not null COMMENT 'NCBI taxon id',
    taxonScientificName varchar(255) not null COMMENT 'NCBI taxon scientific name',
    taxonCommonName varchar(255) COMMENT 'NCBI taxon common name',
    taxonLeftBound int unsigned not null COMMENT 'Left bound taxon id',
    taxonRightBound int unsigned not null COMMENT 'Right bound taxon id',
    taxonLevel mediumint unsigned not null COMMENT 'How deep is the taxon node',
-- bgeeSpeciesLCA defines whether this taxon is the Least Common Ancestor of at least
-- two species used in Bgee. This allows to easily identify important branching.
    bgeeSpeciesLCA boolean not null COMMENT 'Is the Least Common Ancestor of at least two species used in Bgee?'
) engine = innodb;

create table species (
    speciesId mediumint unsigned not null COMMENT 'NCBI species taxon id',
-- example: homo
    genus varchar(70) not null COMMENT 'Genus name',
-- example: sapiens
    species varchar(70) not null COMMENT 'Species name',
-- exemple: human
    speciesCommonName varchar(70) not null COMMENT 'NCBI species common name',
-- integer allowing to sort the species in preferred display order
    speciesDisplayOrder smallint unsigned not null,
-- ID of the taxon which this species belongs to, present in the table `taxonomy`.
-- For instance, if this species is `human`, it belongs to the taxon `homo` (taxon ID 9605).
    taxonId mediumint unsigned not null COMMENT 'NCBI taxon id this species belongs to (most of the time genus taxon id)',
-- Path to retrieve the genome file we use for this species, from the GTF directory
-- of the Ensembl FTP, without the Ensembl version suffix, nor the file type suffixes.
-- For instance, for human, the GTF file in Ensembl 75 is stored at:
-- ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
-- This field would then contain: homo_sapiens/Homo_sapiens.GRCh37
-- This field is needed because we use for some species the genome of another species
-- (for instance, chimp genome for bonobo species).
    genomeFilePath varchar(100) not null COMMENT 'GTF annotation path used to map this species in Ensembl FTP',
    genomeVersion varchar(50) not null,
    dataSourceId smallInt unsigned not null COMMENT 'source for genome information',
-- ID of the species whose the genome was used for this species. This is used
-- when a genome is not in Ensembl. For instance, for bonobo (ID 9597), we use the chimp
-- genome (ID 9598), because bonobo is not in Ensembl.
-- We don't use a foreign key constraint here, because maybe the species whose the genome
-- was used does not have any data in Bgee, and thus is not in the taxon table.
-- If the correct genome of the species was used, the value of this field is 0.
    genomeSpeciesId mediumint unsigned not null default 0 COMMENT 'NCBI species taxon id used for mapping (0 if the same species)'
) engine = innodb;

-- which sex values are permitted for each species.
-- each species will usually have several entries in this table
create table speciesToSex (
    speciesId mediumint unsigned not null,
-- values correspond to some of the existing values in the `cond` table
-- XXX: maybe we actually need a `sex` table?
-- XXX: maybe we'll need an "asexual" value for some species? Or can we always assign kind of a "sex"?
    sex enum('hermaphrodite', 'female', 'male') not null
) engine = innodb;

-- represent mainly alternative common names (for instance, 'rhesus monkey', 'roundworm'),
-- or alternative taxon names related to a species.
create table speciesToKeyword (
    speciesId mediumint unsigned not null,
    keywordId int unsigned not null
) engine = innodb;

-- ****************************************************
-- CONFIDENCE AND EVIDENCE ONTOLOGIES
-- ****************************************************

-- Branch 'confidence information statement' of the CIO
-- (see https://github.com/BgeeDB/confidence-information-ontology).
-- We only use CIO statements for annotations, so we do not insert terms from the branch
-- 'confidence information element'. Moreover, we only insert terms associated to
-- a 'confidence level' term.
-- Also, we do not use relations between CIO terms yet, so they are not inserted for now;
-- if they were to be inserted, we coud use a nested set model (as for the tables stage,
-- taxon, etc), as there is a single is_a inheritance between CIO statements.
-- TODO All confidence levels used in Bgee should use these CIO statements,
-- rather than the enum fields 'low quality'/'high quality'.
-- In order to not use terms from the branch 'confidence information element',
-- this table has 3 columns capturing the 3 different types of CI element a statement
-- can be associatd to. This is highly dependent on the current state of the ontology,
-- this table should be changed if the ontology changed.
create table CIOStatement (
    CIOId varchar(20) not null COMMENT 'Confidence Information Ontology id',
    CIOName varchar(255) not null COMMENT 'Confidence Information Ontology name',
    CIODescription TEXT COMMENT 'Confidence Information Ontology description',
-- define whether this CIO term is used to capture a trusted evidence line (= 1), or whether
-- it indicates that the evidence should not be trusted (= 0).
    trusted tinyint unsigned not null default 0 COMMENT 'Trusted evidence (= 1) or not (default) (= 0)',
-- represent the level of confidence that can be put in a CI statement.
-- These enum fields correspond exactly to the labels of the relevant classes,
-- leaves of the branch 'confidence level'.
-- can be null when the evidenceConcordance is 'strongly conflicting'
    confidenceLevel enum('low confidence level', 'medium confidence level', 'high confidence level')
                                                                 COMMENT 'Confidence level: low/medium/high or null',
-- capture whether there are multiple evidence lines available related to an assertion,
-- and whether they are congruent or conflicting.
-- These enum fields correspond exactly to the labels of the relevant classes,
-- leaves of the branch 'evidence concordance'.
    evidenceConcordance enum('single evidence', 'congruent', 'weakly conflicting', 'strongly conflicting') not null
                                                                 COMMENT 'Evidence concordance: single/congruent/weakly conflicting/strongly conflicting',
-- capture, when there are several evidence lines available related to a same assertion,
-- whether there are of a same or different experimental or computational types.
-- These enum fields correspond exactly to the labels of the relevant classes,
-- leaves of the branch 'evidence type concordance'.
-- It is only applicable when a statement doesn't have an evienceConcordance = 'single evidence'
-- (so this field is null for, and only for, confidence from single evidence)
    evidenceTypeConcordance enum('same type', 'different type') COMMENT 'Evidence type concordance for not "single" evidence: same/different types'
) engine = innodb;

-- Evidence Ontology (see http://www.evidenceontology.org/).
-- Note that we do not insert pre-composed terms used to distinguish between
-- evidence based on manual or automatic assertion (such terms have a relation 'used in'
-- to either 'manual assertion' or 'automatic assertion', and are subclasses of either
-- 'evidence used in manual assertion' or 'evidence used in automatic assertion').
-- So, for instance, we will insert the term 'genetic similarity evidence', not the terms
-- 'genetic similarity evidence used in automatic assertion' and
-- 'genetic similarity evidence used in manual assertion'.
-- Also, we do not use relations between ECO terms yet, so they are not inserted for now;
-- if they were to be inserted, we coud use a nested set model (as for the tables stage,
-- taxon, etc), as there is a single is_a inheritance between these terms.
create table evidenceOntology (
    ECOId varchar(20) not null COMMENT 'Evidence Ontology id',
    ECOName varchar(255) not null COMMENT 'Evidence Ontology name',
    ECODescription TEXT COMMENT 'Evidence Ontology description'
) engine = innodb;


-- ****************************************************
-- ANATOMY AND DEVELOPMENT
-- ****************************************************

create table stage (
    stageId varchar(20) not null COMMENT 'Developmental stage id',
    stageName varchar(255) not null COMMENT 'Developmental stage name',
    stageDescription TEXT COMMENT 'Developmental stage description',
    stageLeftBound int unsigned not null COMMENT '???',
    stageRightBound int unsigned not null COMMENT '???',
    stageLevel int unsigned not null COMMENT 'How deep is the developmental stage',
    tooGranular tinyint unsigned not null default 0 COMMENT '??? Stage is too granular (= 1), or not (default) (= 0)',
    groupingStage tinyint unsigned not null default 0 COMMENT 'Stage to be grouped (= 1) or not (default) (= 0)'
) engine = innodb;

create table stageTaxonConstraint (
    stageId varchar(20) not null COMMENT 'Developmental stage id',
-- if speciesId is null, it means that the stage exists in all species.
-- The aim is to have an entry in this table for each stage,
-- to avoid looking for stages not present, and to avoid creating an entry for each species
-- when the stage exists in all species.
    speciesId mediumint unsigned COMMENT 'NCBI species taxon id on which stage is constrained, or null if stage exists for all species'
) engine = innodb;

create table stageNameSynonym (
    stageId varchar(20) not null COMMENT 'Developmental stage id',
    stageNameSynonym varchar(255) not null COMMENT 'Developmental stage name synonym'
) engine = innodb;

-- XRefs of developmental terms in the Uberon ontology
create table stageXRef (
    stageId varchar(20) not null COMMENT 'Developmental stage id',
    stageXRefId varchar(20) not null COMMENT 'Developmental stage cross-reference id'
) engine = innodb;

create table anatEntity (
    anatEntityId varchar(20) not null COMMENT 'Anatomical entity id',
    anatEntityName varchar(255) not null COMMENT 'Anatomical entity name',
    anatEntityDescription TEXT COMMENT 'Anatomical entity description',
    startStageId varchar(20) not null COMMENT 'Start to exist at this developmental stage',
    endStageId varchar(20) not null COMMENT 'Finish to exist at this developmental stage',
-- a boolean defining whether this anatomical entity is part of
-- a non-informative subset in Uberon, as, for instance,
-- 'upper_level "abstract upper-level terms not directly useful for analysis"'
    nonInformative boolean not null default 0 COMMENT 'Is non-informative (e.g. too broad) (= 1), or not (default) (= 0)'
) engine = innodb;

create table anatEntityTaxonConstraint (
    anatEntityId varchar(20) not null COMMENT 'Anatomical entity id',
-- if speciesId is null, it means that the anatEntity exists in all species.
-- The aim is to have an entry in this table for each anatEntity,
-- to avoid looking for anatEntities not present, and to avoid creating an entry for each species
-- when the anatEntity exists in all species.
    speciesId mediumint unsigned COMMENT 'NCBI species taxon id on which anatomical entity is constrained, or null if anatomical entity exists for all species'
) engine = innodb;

-- XRefs of anatomical terms in the Uberon ontology
create table anatEntityXRef (
    anatEntityId varchar(20) not null COMMENT 'Anatomical entity id',
    anatEntityXRefId varchar(20) not null COMMENT 'Anatomical entity cross-reference id'
) engine = innodb;

create table anatEntityNameSynonym (
    anatEntityId varchar(20) not null COMMENT 'Anatomical entity id',
    anatEntityNameSynonym varchar(255) not null COMMENT 'Anatomical entity name synonym'
) engine = innodb;

create table anatEntityRelation (
    anatEntityRelationId int unsigned not null,
    anatEntitySourceId varchar(20) not null COMMENT 'Anatomical entity source id',
    anatEntityTargetId varchar(20) not null COMMENT 'Anatomical entity target id',
-- there is no distinction made in Bgee between is_a and part_of
    relationType enum('is_a part_of', 'develops_from', 'transformation_of')
                                                COMMENT 'Relation type between anatEntitySourceId & its anatEntityTargetId',
-- relationStatus - direct: the relation is direct between anatEntityParentId and
-- anatEntityDescentId; indirect: this is an indirect relation between two
-- anatomical entities, that have been composed (e.g., part_of o is_a -> part_of);
-- reflexive: a special line added for each anatomical entity, where anatEntityTargetId
-- is equal to anatEntitySourceId. This is useful to get in one join all the descendants
-- of an antomical entity, plus itself (otherwise it requires a 'or' in the join clause,
-- which is non-optimal)
    relationStatus enum('direct', 'indirect', 'reflexive')
                                                COMMENT 'Relation status between anatEntitySourceId & its anatEntityTargetId'
) engine = innodb;

create table anatEntityRelationTaxonConstraint (
    anatEntityRelationId int unsigned not null COMMENT 'Anatomical entity relation id',
-- if speciesId is null, it means that the anatEntityRelation exists in all species.
-- The aim is to have an entry in this table for each anatEntityRelation,
-- to avoid looking for anatEntityRelations not present, and to avoid creating an entry for each species
-- when the anatEntityRelation exists in all species.
    speciesId mediumint unsigned COMMENT 'NCBI species taxon id on which anatomical entity relation is constrained, or null if anatomical entity relation exists for all species'
) engine = innodb;


-- ****************************************************
-- SIMILARITY ANNOTATIONS
-- (See https://github.com/BgeeDB/anatomical-similarity-annotations/)
-- ****************************************************

-- This table captures 'summary' similarity annotations: when several evidence lines
-- are available related to an assertion (same HOM ID, taxon ID, anatEntity IDs),
-- they are summarized into a single summary annotation, providing a global
-- confidence level, emerging from all evidence lines available.
-- See table 'rawSimilarityAnnotation' to retrieve associated single evidence
-- with single confidence level.
-- For convenience, all annotations are inserted in this table, even when only
-- a single evidence is available related to an annotation.
create table summarySimilarityAnnotation (
    summarySimilarityAnnotationId mediumint unsigned not null COMMENT 'Summary similarity annotation id',
-- for now, we only capture annotations of 'historical homology' (HOM:0000007),
-- so we do not use a field 'HOMId'. We should, if in the future we captured other types
-- of similarity annotations.
-- HOMId varchar(20) not null COMMENT 'Historical homology id',
-- the taxon targeted by the similarity annotation
-- (note that the similarity annotation file lets open the possibility of capturing
-- several taxon IDs, for instance to define in which taxa a structure is
-- functionally equivalent, as this type of relation would not originate from
-- a common ancestor; but, this is not yet done, so we use this field, and not a link table).
    taxonId mediumint unsigned not null COMMENT 'NCBI species taxon id targeted by the similarity annotation',
-- define whether this annotation is negated (using the NOT qualifier of the similarity
-- annotation file); this would mean that there existed only negative evidence lines
-- related to this annotation (when evidence lines are conflicting, the summary annotation
-- is considered positive, because we are primarly interested in positive annotations).
    negated boolean not null default 0 COMMENT 'Is this annotation negated (= 1), or not (default) (= 0)',
-- the ID of the confidence statement associated to this summary annotation;
-- allows to capture the global confidence level, whether evidence lines were congruent
-- or conflicting, etc.
-- If this summary annotation corresponds to an annotation supported by
-- a single evidence line, then this CIO term will be the same as the one used in the table
-- 'rawSimilarityAnnotation' for the related single annotation. Otherwise,
-- it will be a CIO statement from the 'multiple evidence lines' branch.
    CIOId varchar(20) not null COMMENT 'Confidence Information Ontology id'
) engine = innodb;

-- similarity annotations can target several anatomical entities (e.g., to capture
-- the homology between 'lung' and 'swim bladder'), and an antomical entity can be targeted
-- by several annotations (e.g., to capture multiple homology hypotheses);
-- this is why we need this link table.
create table similarityAnnotationToAnatEntityId (
    summarySimilarityAnnotationId mediumint unsigned not null COMMENT 'Summary similarity annotation id',
    anatEntityId varchar(20) not null COMMENT 'Anatomical entity id'
) engine = innodb;

-- Represent raw similarity annotations, capturing one single evidence line,
-- which corresponds to the GO guidelines to capture sources of annotations (see
-- http://geneontology.org/page/guide-go-evidence-codes). When several evidence lines are available
-- related to a same assertion, they are captured in an summary annotation, summarizing
-- all evidence lines available, see table 'summarySimilarityAnnotation';
-- for convenience, annotations are all present in the table 'summarySimilarityAnnotation'
-- anyway, even when they capture a single evidence.
-- So, this table provides information about single evidence related to an annotation:
-- the individual evidence and confidence codes, the reference ID, the supporting text,
-- the annotator who made the annotation, etc... Other "global" information are present
-- in the table 'summarySimilarityAnnotation' (e.g., the HOM ID, the taxon ID).
-- Targeted anatomical entities are stored in the table 'similarityAnnotationToAnatEntityId',
-- and can be retrieved through the table 'summarySimilarityAnnotation'.
create table rawSimilarityAnnotation (
-- the associated 'summary' annotation
    summarySimilarityAnnotationId mediumint unsigned not null COMMENT 'Summary similarity annotation id',
-- define whether this annotation is negated (using the NOT qualifier of the similarity
-- annotation file: used to capture an information rejecting a putative relation
-- between structures, that could otherwise seem plausible).
    negated boolean not null default 0 COMMENT 'Is this annotation negated (= 1), or not (default) (= 0)',
-- capture how the annotation is supported
    ECOId varchar(20) not null COMMENT 'Evidence Ontology id',
-- the ID of the confidence statement associated to this annotation;
-- it can only be a confidence statement from the branch
-- 'confidence statement from single evidence'
-- XXX: maybe we could have a trigger to check that it is a term from the correct branch;
-- it would require to store branch information, or relations between terms, in the table
-- 'CIOStatement'.
    CIOId varchar(20) not null COMMENT 'Confidence Information Ontology id',
-- Unique identifier of a single source, cited as an authority for asserting the relation.
-- Can be a DOI, a Pubmed ID, an ISBN, an URL.
-- XXX: should it be a TEXT field, to store long URL? It does not seem to be necessary for now.
    referenceId varchar(255) not null COMMENT 'Any unique id cited as an authority for asserting the relation',
-- information provided for convenience, manually captured.
    referenceTitle TEXT not null COMMENT 'Refence id description',
-- A quote from the reference, supporting the annotation. If possible, it should
-- also support the choice of the ECO and CIO IDs.
    supportingText TEXT not null COMMENT 'Annotation support description',
-- The database which made the annotation. Used for tracking the source of
-- an individual annotation. Currently, only Bgee is working on this file.
    assignedBy varchar(20) not null COMMENT 'Database source for the annotation',
-- A code allowing to identify the curator who made the annotation, from the database
-- defined above.
-- XXX: if we assume that these annotations are internal to Bgee, then this field
-- should be named 'authorId', and have a foreign key constraint to the table 'author';
-- if we assume these annotations are a community effort, then this simple varchar field
-- is fine. Let's stick to the community effort.
    curator varchar(20) not null COMMENT 'A code allowing to identify the curator who made the annotation',
-- Date when the annotation was made (AAAA-MM-JJ)
-- XXX: other tables store date in different format, which one is better really?
    annotationDate date COMMENT 'Date when the annotation was made'
) engine = innodb;


-- ****************************************************
-- GENE AND TRANSCRIPT INFO
-- ****************************************************
-- Hierarchical Orthologous Groups from OMA.

-- All the nodes of a particular group are stored in a nested set model.
-- A node in the tree could be a speciation node or a duplication node.
-- The OMANodeLeftBound and OMANodeRightBound correspond to the left and right bound IDs of the nested set model.
-- Note: to use the nested set model, we often need to join this table to itself,
-- using a range condition on left and right bounds for the join clause; sadly,
-- there is a performance issue for such queries in MySQL, see
-- http://www.percona.com/blog/2010/05/17/joining-on-range-wrong/
create table OMAHierarchicalGroup (
    -- A unique ID for each node inside an OMA Hierarchical Orthologous Group.
    -- Auto generated by us, unique over all groups (use as primary key)
    OMANodeId int unsigned not null COMMENT 'OMA Hierarchical Orthologous node id',
    -- The ID of Hierarchical Orthologous Group as provided by OMA.
    -- Only for Xref purpose.
    OMAGroupId varchar(255) not null COMMENT 'OMA Hierarchical Orthologous Group id',
    -- Bounds generated over all groups.
    OMANodeLeftBound int unsigned not null COMMENT 'OMA left bound id in the nested set model',
    OMANodeRightBound int unsigned not null COMMENT 'OMA right bound id in the nested set model',
    -- The ID corresponding to the level of taxonomy as in NCBI.
    -- Some nodes have no taxonomy ID because they correspond to a duplication node (paralogy group).
    taxonId mediumint unsigned COMMENT 'NCBI taxon id corresponding to the level of taxonomy'
) engine = innodb;

create table geneOntologyTerm (
    goId char(10) not null COMMENT 'Gene Ontology id',
    goTerm varchar(255) not null COMMENT 'Gene Ontology term',
    goDomain enum ('biological process', 'cellular component', 'molecular function')
                                    COMMENT 'Gene Ontology domain'
) engine = innodb;

-- link a GO ID to its alternative IDs
create table geneOntologyTermAltId (
    goId char(10) not null COMMENT 'Gene Ontology id',
    goAltId char(10) not null COMMENT 'Gene Ontology alternative id'
) engine = innodb;

-- list all is_a or part_of relations between GO terms, even indirect.
-- Relations other than is_a or part_of are not considered.
create table geneOntologyRelation (
    goAllTargetId char(10) not null COMMENT 'Gene Ontology target id for is_a or part_of relations from goAllSourceId',
    goAllSourceId char(10) not null COMMENT 'Gene Ontology source id'
) engine = innodb;

create table geneBioType (
    geneBioTypeId smallint unsigned not null COMMENT 'Gene BioType id (type of gene)',
    geneBioTypeName varchar(255) not null default '' COMMENT 'Gene BioType name'
) engine = innodb;

create table gene (
-- warning, maybe this bgeeGeneId will need to be changed to an 'int' when we reach around 200 species
    bgeeGeneId mediumint unsigned not null COMMENT 'Numeric internal gene ID used for improving performances',
    geneId varchar(20) not null COMMENT 'Real gene id',
    geneName varchar(255) not null default '' COMMENT 'Gene name',
    geneDescription TEXT COMMENT 'Gene description',
    speciesId mediumint unsigned not null COMMENT 'NCBI species taxon id this gene belongs to',
    geneBioTypeId smallint unsigned COMMENT 'Gene BioType id (type of gene)',
-- can be null if the gene does not belong to a hierarchical group
-- a gene can belong to one and only one group
-- OMA parent node ID instead of OMA node ID to avoid create group for all genes
    OMAParentNodeId int unsigned default null COMMENT 'OMA Hierarchical Orthologous parent node id',
-- defines whether the gene ID is present in Ensembl. For some species, they are not
-- (for instance, bonobo; we use chimp genome)
    ensemblGene boolean not null default 1 COMMENT 'Is the gene in Ensembl (default) (= 1), if not (= 0)',
    geneMappedToGeneIdCount tinyint unsigned not null default 1 COMMENT 'number of genes in the Bgee database with the same Ensembl gene ID. In Bgee, for some species with no genome available, we use the genome of a closely-related species, such as chimpanzee genome for analyzing bonobo data. For this reason, a same Ensembl gene ID can be mapped to several species in Bgee. The value returned here is equal to 1 when the Ensembl gene ID is uniquely used in the Bgee database.'
) engine = innodb;

create table geneToOma (
    bgeeGeneId mediumint unsigned not null,
    OMANodeId int unsigned not null COMMENT 'OMA Hierarchical Orthologous node id',
    taxonId mediumint unsigned not null
) engine = innodb;

create table geneNameSynonym (
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
    geneNameSynonym varchar(255) not null COMMENT 'Gene name synonym'
) engine = innodb;

create table geneXRef (
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
    XRefId varchar(20) not null COMMENT 'Cross-reference id',
    XRefName varchar(255) not null default '' COMMENT 'Cross-reference name',
    dataSourceId smallInt unsigned not null COMMENT 'Data Source id the cross-reference comes from'
) engine = innodb;

create table geneToTerm (
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
    term varchar(255) not null COMMENT '???Gene alias name'
) engine = innodb;

-- TODO: use IDs from the Evidence Ontology, rather than the Evidence Codes used by
-- the GO consortium. The field 'goEvidenceCode' would be replaced by a field 'ECOId',
-- with a foreign key constraint to the table 'evidenceOntology'. This would require
-- for the application inserting data in this table to retrieve the mapping between
-- ECO IDs and Evidence Codes from the Evidence Ontology.
create table geneToGeneOntologyTerm (
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
    goId char(10) not null COMMENT 'Gene Ontology id',
    goEvidenceCode varchar(20) not null default '' COMMENT 'Gene Ontology Evidence Code'
) engine = innodb;

create table transcript (
    bgeeTranscriptId int unsigned not null COMMENT 'Numeric internal transcript ID used for improving performances',
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID this transcript is mapped to',
    transcriptId varchar(40) not null COMMENT 'Real transcript ID',
    transcriptName varchar(255) not null default '',
    transcriptDescription TEXT,
    transcriptLength mediumint unsigned not null,
    effectiveTranscriptLength mediumint unsigned not null
) engine = innodb;


-- ****************************************************
-- CONDITIONS
-- ****************************************************
-- 'condition' is a reserved keyword in MySQL, we can't use it as table name
create table cond (
    conditionId           mediumint unsigned not null COMMENT 'Internal condition ID. Each condition is species-specific',
    exprMappedConditionId mediumint unsigned not null COMMENT 'the condition ID that should be used for insertion into the expression table: too-granular conditions (e.g., 43 yo human stage, or sexInferred=1) are mapped to less granular conditions for summary. Equal to conditionId if condition is not too granular.',
    anatEntityId          varchar(20)  not null       COMMENT 'Uberon anatomical entity ID',
    stageId               varchar(20)  not null       COMMENT 'Uberon stage ID',
    speciesId             mediumint unsigned not null COMMENT 'NCBI species taxon ID',
-- NA: not available from source information
-- not annotated: information not captured by Bgee
-- If an ENUM column is declared NOT NULL, its default value is the first element of the list
    sex enum('not annotated', 'hermaphrodite', 'female', 'male', 'mixed', 'NA') not null
    COMMENT 'Sex information. NA: not available from source information; not annotated: information not captured by Bgee. Note that all conditions used in the expression tables have "NA" replaced with "not annotated".',
    sexInferred boolean not null default 0
    COMMENT 'Whether sex information was retrieved from annotation (false), or inferred from information in Uberon (true). Note that all conditions used in the expression tables use a "0" value.',
-- For now, strains are captured as free-text format, only 4 term are "standardized":
-- 'NA', 'not annotated', 'wild-type', 'confidential_restricted_data'.
-- This should be improved in a further release (free-text is hardly satisfiable).
    strain varchar(100) not null default 'not annotated'
    COMMENT 'Strain information. NA: not available from source information; not annotated: information not captured by Bgee; confidential_restricted_data: information cannot be disclosed publicly. Note that all conditions used in the expression tables have "NA", "not annotated" and "confidential_restricted_data" replaced with "wild-type"'
) engine = innodb COMMENT 'This table stores the "raw" conditions used to annotate data and used in the "raw" expression table, where data are not propagated nor precomputed';

create table globalCond (
    globalConditionId mediumint unsigned not null,
    anatEntityId          varchar(20)  COMMENT 'Uberon anatomical entity ID. Can be null in this table if this condition aggregates data according to other condition parameters (e.g., grouping all data in a same stage whatever the organ is).',
    stageId               varchar(20)  COMMENT 'Uberon stage ID. Can be null in this table if this condition aggregates data according to other condition parameters (e.g., grouping all data in a same organ whatever the dev. stage is).',
    speciesId             mediumint unsigned not null COMMENT 'NCBI species taxon ID',
-- NA: not available from source information
-- not annotated: information not captured by Bgee
-- If an ENUM column is declared NOT NULL, its default value is the first element of the list
-- In this table, only 'not annotated' is used to replace 'NA', as for conditions
-- used in expression table
    sex enum('not annotated', 'hermaphrodite', 'female', 'male', 'mixed')
    COMMENT 'Sex information. NA: not available from source information; not annotated: not used in this table, since all conditions used in the expression tables have "NA" replaced with "not annotated". Can be null in this table if this condition aggregates data according to other condition parameters (e.g., grouping all data in a same organ whatever the sex is).',
-- For now, strains are captured as free-text format, only 4 term are "standardized":
-- 'NA', 'not annotated', 'wild-type', 'confidential_restricted_data'.
-- In this table, only 'wild-type' is used to replace 'NA', 'not annotated', and
-- 'confidential_restricted_data', as for conditions used in expression table.
    strain varchar(100)
    COMMENT 'Strain information. NA: not available from source information; not annotated: information not captured by Bgee; confidential_restricted_data: information cannot be disclosed publicly',

-- ** RANKS **
-- max ranks in each data type and condition, notably used to allow normalization
-- between data types and conditions. For EST and in situ data, they are also used for computation
-- of weighted mean between data types: for these data types, because we pool together all data
-- in a same condition, instead of computing a mean between samples, and because we use "dense ranking"
-- instead of fractional ranking (so that the max rank is equal to the number of distinct ranks),
-- it is irrelevant to consider a sum of the number of distinct ranks in each sample for weighting
-- the mean, as for Affymetrix and EST data.
-- Note: these values are the same for all genes in a condition-species, this is why they are stored in this table.
    affymetrixMaxRank decimal(9,2) unsigned,
    rnaSeqMaxRank decimal(9,2) unsigned,
    estMaxRank decimal(9,2) unsigned,
    inSituMaxRank decimal(9,2) unsigned,

    affymetrixGlobalMaxRank decimal(9,2) unsigned COMMENT 'This max rank is computed by taking into account all data in this condition, but also in all child conditions.',
    rnaSeqGlobalMaxRank decimal(9,2) unsigned COMMENT 'This max rank is computed by taking into account all data in this condition, but also in all child conditions.',
    estGlobalMaxRank decimal(9,2) unsigned COMMENT 'This max rank is computed by taking into account all data in this condition, but also in all child conditions.',
    inSituGlobalMaxRank decimal(9,2) unsigned COMMENT 'This max rank is computed by taking into account all data in this condition, but also in all child conditions.'
) engine = innodb COMMENT 'This table contains all condition used in the globalExpression table. It thus includes "real" conditions used in the raw expression table, but mostly conditions resulting from the propagation of expression calls in the globalExpression table. It results from the computation of propagated calls according to different condition parameters combination (e.g., grouping all data in a same anat. entity, or all data in a same anat. entity - stage, or data in anat. entity - sex). This is why the fields anatEntityId, stageId, sex, strain, can be null in this table (but not all of them at the same time).';

create table globalCondToCond (
    globalConditionId mediumint unsigned not null,
    conditionId mediumint unsigned not null,
    conditionRelationOrigin enum('self', 'descendant', 'parent') not null COMMENT 'Define whether the data from the raw conditions used for production of global calls in this global condition comes from raw conditions mapped to the globalCondition itself, a descendant global condition, or a parent global condition.'
) engine = innodb
comment = 'this table allows to link globalConditions to the raw conditions that were aggregated to produce global expression calls in the globalExpression table.';

-- ****************************************************
-- RAW EST DATA
-- ****************************************************
create table estLibrary (
    estLibraryId varchar(50) not null,
    estLibraryName varchar(255) not null,
    estLibraryDescription text,
    conditionId mediumint unsigned not null,
    dataSourceId smallInt unsigned not null
) engine = innodb;

create table estLibraryToKeyword (
    estLibraryId varchar(50) not null,
    keywordId int unsigned not null
) engine = innodb;

create table expressedSequenceTag (
    estId varchar(50) not null,
-- ESTs have two IDs in Unigene
    estId2 varchar(50) not null default '',
    estLibraryId varchar(50) not null,
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
    UniGeneClusterId varchar(70) not null default '',
    expressionId int unsigned,
-- Warning, qualities must be ordered, the index in the enum is used in many queries
    estData enum('no data', 'poor quality', 'high quality') default 'no data'
) engine = innodb;

create table estLibraryExpression (
    expressionId int unsigned not null,
    estLibraryId varchar(50) not null,
    estCount mediumint unsigned not null default 0
        comment 'number of ESTs in this library mapped to the gene associated to this expressionId',
-- no 'callDirection' column for ESTs, only 'present' calls are generated from ESTs
    estLibraryCallQuality enum('poor quality', 'high quality') not null
        comment 'Inferred quality for this call based on this library (based on the number of ESTs mapped to gene, see Audic and Claverie 1997). Value "poor quality" instead of "low quality" for historical reasons.'
) engine = innodb
comment = 'This table stores information about expression calls produced from EST libraries, that is then used in Bgee to compute global summary expression calls and qualities. Only "present" calls are generated from ESTs (no "absent" calls).';

-- ****************************************************
-- RAW AFFYMETRIX DATA
-- ****************************************************
create table microarrayExperiment (
    microarrayExperimentId varchar(70) not null,
    microarrayExperimentName varchar(255) not null default '',
    microarrayExperimentDescription text,
    dataSourceId smallInt unsigned not null
) engine = innodb;

create table microarrayExperimentToKeyword (
    microarrayExperimentId varchar(70) not null,
    keywordId int unsigned not null
) engine = innodb;

create table chipType (
    chipTypeId varchar(70) not null,
    chipTypeName varchar(255) not null,
    cdfName varchar(255) not null,
    isCompatible tinyint(1) not null default 1,
    qualityScoreThreshold decimal(10, 2) unsigned not null default 0,
-- percentage of present probesets
-- 100.00
    percentPresentThreshold decimal(5, 2) unsigned not null default 0,

-- this field is used for rank computations, and is set after all expression data insertion,
-- this is why null value is permitted.
    chipTypeMaxRank decimal(9,2) unsigned COMMENT 'The max fractional rank in this chip type (see `rank` field in affymetrixProbeset table)'
) engine = innodb;

-- this table represents mapping of affymetrix probesets in general,
-- not constrainted by the tables chipType and afymetrixProbeset
-- (that means for instance that you can insert in this table a mapping
-- for a probeset not present in the table affymetrixProbeset)
-- => so, NO foreign keys to the tables affymetrixProbeset and chipType.
-- moreover, the probeset mapping can be use for other tables
-- (deaAffymetrixProbesetGroups)
-- create table affymetrixProbesetMapping(
-- chipTypeId varchar(70) not null,
-- affymetrixProbesetId varchar(70) not null,
-- bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID'
-- ) engine = innodb;

create table affymetrixChip (
-- affymetrixChipId are not unique (couple affymetrixChipId - microarrayExperimentId is)
-- then we need an internal ID to link to affymetrixProbeset
-- warning, SMALLINT UNSIGNED only allows for 65535 chips to be inserted (we have 12,996 as of Bgee 14)
    bgeeAffymetrixChipId smallint unsigned not null,
    affymetrixChipId varchar(255) not null,
    microarrayExperimentId varchar(70) not null,
-- define only if CEL file available, normalization gcRMA, detection schuster
    chipTypeId varchar(70),
    scanDate varchar(70) not null default '',
-- An <code>enum</code> listing the different methods used ib Bgee
-- to normalize Affymetrix data:
-- * MAS5: normalization using the MAS5 software. Using
-- this naormalization usually means that only the processed MAS5 files
-- were available, otherwise another method would be used.
-- * RMA: normalization by RMA method.
-- * gcRMA: normalization by gcRMA method. This is the default
-- method in Bgee when raw data are available.
    normalizationType enum('MAS5', 'RMA', 'gcRMA') not null,
-- An <code>enum</code> listing the different methods to generate expression calls
-- on Affymetrix chips:
-- * MAS5: expression calls from the MAS5 software. Such calls
-- are usually taken from a processed MAS5 file, and imply that the data
-- were also normalizd using MAS5.
-- * Schuster: Wilcoxon test on the signal of probesets
-- against a subset of weakly expressed probesets, to generate expression calls
-- (see https://www.ncbi.nlm.nih.gov/pubmed/17594492). Such calls usually implies
-- that raw data were available, and were normalized using gcRMA.
    detectionType enum('MAS5', 'Schuster') not null,
    conditionId mediumint unsigned not null,
-- arIQR_score Marta score
-- can be set to 0 if it is a MAS5 file
-- 99999999.99
    qualityScore decimal(10, 2) unsigned not null default 0,
-- percentage of present probesets
-- 100.00
    percentPresent decimal(5, 2) unsigned not null,

-- the following fields are used for rank computations, and are set after all expression data insertion,
-- this is why null value is permitted.
    chipMaxRank decimal(9,2) unsigned COMMENT 'The max fractional rank in this chip (see `rank` field in affymetrixProbeset table)',
    chipDistinctRankCount mediumint unsigned COMMENT 'The count of distinct rank in this chip (see `rank` field in affymetrixProbeset table, used for weighted mean rank computations)'
) engine = innodb;

create table affymetrixProbeset (
    affymetrixProbesetId varchar(70) not null,
    bgeeAffymetrixChipId smallint unsigned not null,
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
    normalizedSignalIntensity decimal(13,5) unsigned not null default 0,
-- Warning, flags must be ordered, the index in the enum is used in many queries
    detectionFlag enum('undefined', 'absent', 'marginal', 'present') not null default 'undefined',
    expressionId int unsigned,
-- rank is not "not null" because we update this information afterwards.
-- note that this corresponds to the rank of the gene, not of the probeset
-- (so, all probesets mapped to a same gene have the same rank, based on its highest signal intensity)
    rank decimal(9, 2) unsigned,
-- Warning, qualities must be ordered, the index in the enum is used in many queries
    affymetrixData enum('no data', 'poor quality', 'high quality') not null default 'no data',
-- When expressionId is null, the result is not used for the summary of expression.
-- Reasons are:
-- * pre filtering: Probesets always seen as "absent" or "marginal" over the whole dataset are removed
-- * noExpression conflict: a "noExpression" result has been removed because of expression in a sub-condition.
-- Note: as of Bgee 14, we haven't remove this reason for exclusion, but we don't use it for now,
-- as we might want to take into account noExpression in parent conditions for generating
-- a global expression calls, where there is expression in a sub-condition.
-- Maybe we'll discard them again, but I don't think so, it'll allow to present absolutely
-- all data available about a call to users.
-- * undefined: only 'undefined' calls have been seen
--
-- Note that, as of Bgee 14, 2 reasons for exclusion were removed: 'bronze quality' and 'absent low quality'.
-- 'bronze quality' exclusion was removed, because now we always propagate expression evidence,
-- so a 'bronze quality' call can provide additional evidence to a parent structure.
-- 'bronze quality' used to be: for a gene/condition, no "present high" and mix of "present low" and "absent".
-- 'absent low quality' was removed, because we now use a same consistent mechanism for present/absent calls,
-- taking also into account 'absent low quality' evidence.
-- 'absent low quality' used to be: probesets always "absent" for this gene/condition,
-- but only seen by MAS5 (that we do not trust = "low quality" - "noExpression" should always be "high quality").
    reasonForExclusion enum('not excluded', 'pre-filtering',
        'noExpression conflict', 'undefined') not null default 'not excluded'
) engine = innodb;

create table microarrayExperimentExpression (
    expressionId int unsigned not null,
    microarrayExperimentId varchar(70) not null,
    presentHighMicroarrayChipCount smallint unsigned not null default 0
        comment 'number of chips in this experiment that produced this call as present high quality',
    presentLowMicroarrayChipCount  smallint unsigned not null default 0
        comment 'number of chips in this experiment that produced this call as present low quality',
    absentHighMicroarrayChipCount  smallint unsigned not null default 0
        comment 'number of chips in this experiment that produced this call as absent high quality',
    absentLowMicroarrayChipCount   smallint unsigned not null default 0
        comment 'number of chips in this experiment that produced this call as absent low quality',
    microarrayExperimentCallDirection enum('present', 'absent') not null
        comment 'Inferred direction for this call based on this experiment ("present" chips always win over "absent" chips)',
    microarrayExperimentCallQuality enum('poor quality', 'high quality') not null
        comment 'Inferred quality for this call based on this experiment (from all chips, "present high" > "present low" > "absent high" > "absent low"). Value "poor quality" instead of "low quality" for historical reasons.'
) engine = innodb
comment = 'This table stores information about expression calls produced from microarray experiments, that is then used in Bgee to compute global summary expression calls and qualities.';

-- ****************************************************
-- IN SITU HYBRIDIZATION DATA
-- ****************************************************
create table inSituExperiment (
    inSituExperimentId varchar(70) not null,
    inSituExperimentName varchar(255) not null default '',
    inSituExperimentDescription text,
    dataSourceId smallInt unsigned not null
) engine = innodb;

create table inSituExperimentToKeyword (
    inSituExperimentId varchar(70) not null,
    keywordId int unsigned not null
) engine = innodb;

-- evidence: picture, figure, paper, ...
create table inSituEvidence (
    inSituEvidenceId varchar(70) not null,
    inSituExperimentId varchar(70) not null,
-- some databases do not allow to distinguish different samples used in an experiment,
-- all results are merged into one "fake" sample. In that case, this boolean is false.
    evidenceDistinguishable boolean not null default 1,
-- an information used to generate URLs to this sample, taht can be used in the evidenceUrl
-- of the related DataSource. For instance, in MGI this represents
-- the ID of the image to link to (but as an image is not always available, we cannot
-- use it as the inSituEvidenceId)
    inSituEvidenceUrlPart varchar(255) not null default ''
) engine = innodb;

-- Absent spots can be associated to an expressionId, if there is other data
-- showing expression for the same gene/organ/stage
create table inSituSpot (
    inSituSpotId varchar(70) not null,
    inSituEvidenceId varchar(70) not null,
    -- for control purpose only (used in other databases)
    inSituExpressionPatternId varchar(70) not null,
    conditionId mediumint unsigned not null,
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
-- Warning, tags must be ordered, the index in the enum is used in many queries
    detectionFlag enum('undefined', 'absent', 'present') default 'undefined',
    expressionId int unsigned,
-- Warning, qualities must be ordered, the index in the enum is used in many queries
    inSituData enum('no data', 'poor quality', 'high quality') default 'no data',
-- When expressionId is null, the result is not used for the summary of expression.
-- Reasons are:
-- * pre filtering: Probesets always seen as "absent" or "marginal" over the whole dataset are removed
-- * noExpression conflict: a "noExpression" result has been removed because of expression in a sub-condition.
-- Note: as of Bgee 14, we haven't remove this reason for exclusion, but we don't use it for now,
-- as we might want to take into account noExpression in parent conditions for generating
-- a global expression calls, where there is expression in a sub-condition.
-- Maybe we'll discard them again, but I don't think so, it'll allow to present absolutely
-- all data available about a call to users.
-- * undefined: only 'undefined' calls have been seen
--
-- Note that, as of Bgee 14, 2 reasons for exclusion were removed: 'bronze quality' and 'absent low quality'.
-- 'bronze quality' exclusion was removed, because now we always propagate expression evidence,
-- so a 'bronze quality' call can provide additional evidence to a parent structure.
-- 'bronze quality' used to be: for a gene/condition, no "present high" and mix of "present low" and "absent".
-- 'absent low quality' was removed, because we now use a same consistent mechanism for present/absent calls,
-- taking also into account 'absent low quality' evidence.
-- 'absent low quality' used to be: probesets always "absent" for this gene/condition,
-- but only seen by MAS5 (that we do not trust = "low quality" - "noExpression" should always be "high quality").
    reasonForExclusion enum('not excluded', 'pre-filtering',
        'noExpression conflict', 'undefined') not null default 'not excluded'
) engine = innodb;

create table inSituExperimentExpression (
    expressionId int unsigned not null,
    inSituExperimentId varchar(70) not null,
    presentHighInSituSpotCount smallint unsigned not null default 0
        comment 'number of spots in this experiment that produced this call as present high quality',
    presentLowInSituSpotCount  smallint unsigned not null default 0
        comment 'number of spots in this experiment that produced this call as present low quality',
    absentHighInSituSpotCount  smallint unsigned not null default 0
        comment 'number of spots in this experiment that produced this call as absent high quality',
    absentLowInSituSpotCount   smallint unsigned not null default 0
        comment 'number of spots in this experiment that produced this call as absent low quality',
    inSituExperimentCallDirection enum('present', 'absent') not null
        comment 'Inferred direction for this call based on this experiment ("present" spots always win over "absent" spots)',
    inSituExperimentCallQuality enum('poor quality', 'high quality') not null
        comment 'Inferred quality for this call based on this experiment (from all spots, "present high" > "present low" > "absent high" > "absent low"). Value "poor quality" instead of "low quality" for historical reasons.'
) engine = innodb
comment = 'This table stores information about expression calls produced from in situ hybridization experiments, that is then used in Bgee to compute global summary expression calls and qualities.';

-- ****************************************************
-- RNA-Seq DATA
-- ****************************************************
create table rnaSeqExperiment (
-- primary exp ID, from GEO, patterns GSExxx
    rnaSeqExperimentId varchar(70) not null,
    rnaSeqExperimentName varchar(255) not null default '',
    rnaSeqExperimentDescription text,
    dataSourceId smallInt unsigned not null
) engine = innodb;

create table rnaSeqExperimentToKeyword (
    rnaSeqExperimentId varchar(70) not null,
    keywordId int unsigned not null
) engine = innodb;

create table rnaSeqPlatform (
    rnaSeqPlatformId varchar(255) not null,
    rnaSeqPlatformDescription text
) engine = innodb;

-- corresponds to one sample
-- uses to produce several runs
create table rnaSeqLibrary (
-- primary ID, from GEO, pattern GSMxxx
    rnaSeqLibraryId varchar(70) not null,
    rnaSeqExperimentId varchar(70) not null,
    rnaSeqPlatformId varchar(255) not null,
    conditionId mediumint unsigned not null,
-- TMM normalization factor
    tmmFactor decimal(8, 6) not null default 1.0,
-- FPKM threshold to consider a gene as expressed
    fpkmThreshold decimal(16, 6) not null,
-- TPM threshold to consider a gene as expressed
    tpmThreshold decimal(16, 6) not null,
    allGenesPercentPresent decimal(5, 2) unsigned not null default 0,
    proteinCodingGenesPercentPresent decimal(5, 2) unsigned not null default 0,
    intergenicRegionsPercentPresent decimal(5, 2) unsigned not null default 0,
    thresholdRatioIntergenicCodingPercent decimal(5, 2) unsigned not null default 0
            COMMENT 'Proportion intergenic/coding region used to define the threshold to consider a gene as expressed (should always be 5%, but some libraries do not allow to reach this value)',
-- total number of reads in library, including those not mapped.
-- In case of paired-end libraries, it's the number of pairs of reads;
-- In case of single read, it's the total number of reads
    allReadsCount int unsigned not null default 0,
-- total number of reads in library that were mapped to anything.
-- if it is not a paired-end library, this number is equal to leftMappedReadsCount
    mappedReadsCount int unsigned not null default 0,
-- a library is an assembly of different runs, and the runs can have different read lengths,
-- so we store the min and max read lengths
    minReadLength int unsigned not null default 0,
    maxReadLength int unsigned not null default 0,
-- Is the library built using paired end?
-- NA: info not used for pseudo-mapping. Default value in an enum is the first one.
    libraryType enum('NA', 'single', 'paired') not null,
    libraryOrientation enum('NA', 'forward', 'reverse', 'unstranded') not null,

-- the following fields are used for rank computations, and are set after all expression data insertion,
-- this is why null value is permitted.
    libraryMaxRank decimal(9,2) unsigned COMMENT 'The max fractional rank in this library (see `rank` field in rnaSeqResult table)',
    libraryDistinctRankCount mediumint unsigned COMMENT 'The count of distinct rank in this library (see `rank` field in rnaSeqResult table, used for weighted mean rank computations)'
) engine = innodb;

-- Store the information of runs used, pool together to generate the results
-- for a given library.
create table rnaSeqRun (
-- same ID in GEO and SRA, pattern SRR...
    rnaSeqRunId varchar(70) not null,
    rnaSeqLibraryId varchar(70) not null
) engine = innodb;

-- We sometimes discard some runs associated to a library, because of low mappability.
-- We keep track of these discarded runs in this table.
-- UPDATE Bgee 14: for pseudo-mapping using Kallisto, runs are pooled, so we can only exclude libraries,
-- not specific runs.
create table rnaSeqLibraryDiscarded (
    rnaSeqLibraryId varchar(70) not null
) engine = innodb;

-- This table contains TPM/RPKM/read count values for each gene for each library
-- and link them to an expressionId
create table rnaSeqResult (
    rnaSeqLibraryId varchar(70) not null,
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
    fpkm decimal(16, 6) not null,
    tpm decimal(16, 6) not null,
-- rank is not "not null" because we update this information afterwards
    rank decimal(9, 2) unsigned,
-- for information, measure not normalized for reads or genes lengths
    readsCount decimal(16, 6) unsigned not null,
    expressionId int unsigned,
    detectionFlag enum('undefined', 'absent', 'present') default 'undefined',
-- Warning, qualities must be ordered, the index in the enum is used in many queries.
-- We should only see genes with 'high quality' here
    rnaSeqData enum('no data', 'poor quality', 'high quality') default 'no data',
-- When expressionId is null, the result is not used for the summary of expression.
-- Reasons are:
-- * pre filtering: Probesets always seen as "absent" or "marginal" over the whole dataset are removed
-- * noExpression conflict: a "noExpression" result has been removed because of expression in a sub-condition.
-- Note: as of Bgee 14, we haven't remove this reason for exclusion, but we don't use it for now,
-- as we might want to take into account noExpression in parent conditions for generating
-- a global expression calls, where there is expression in a sub-condition.
-- Maybe we'll discard them again, but I don't think so, it'll allow to present absolutely
-- all data available about a call to users.
-- * undefined: only 'undefined' calls have been seen
--
-- Note that, as of Bgee 14, 2 reasons for exclusion were removed: 'bronze quality' and 'absent low quality'.
-- 'bronze quality' exclusion was removed, because now we always propagate expression evidence,
-- so a 'bronze quality' call can provide additional evidence to a parent structure.
-- 'bronze quality' used to be: for a gene/condition, no "present high" and mix of "present low" and "absent".
-- 'absent low quality' was removed, because we now use a same consistent mechanism for present/absent calls,
-- taking also into account 'absent low quality' evidence.
-- 'absent low quality' used to be: probesets always "absent" for this gene/condition,
-- but only seen by MAS5 (that we do not trust = "low quality" - "noExpression" should always be "high quality").
    reasonForExclusion enum('not excluded', 'pre-filtering',
        'noExpression conflict', 'undefined') not null default 'not excluded'
) engine = innodb;

-- This table contains TPM/RPKM/read count values for each transcript for each library
-- NOTE Bgee 14: as of Bgee 14 this table is not filled
create table rnaSeqTranscriptResult (
    rnaSeqLibraryId varchar(70) not null,
    bgeeTranscriptId int unsigned not null COMMENT 'Internal transcript ID',
    fpkm decimal(16, 6) not null,
    tpm decimal(16, 6) not null,
-- for information, measure not normalized for reads or genes lengths
    readsCount int unsigned not null
) engine = innodb;

create table rnaSeqExperimentExpression (
    expressionId int unsigned not null,
    rnaSeqExperimentId varchar(70) not null,
    presentHighRNASeqLibraryCount smallint unsigned not null default 0
        comment 'number of RNA-Seq libraries in this experiment that produced this call as present high quality',
    presentLowRNASeqLibraryCount  smallint unsigned not null default 0
        comment 'number of RNA-Seq libraries in this experiment that produced this call as present low quality',
    absentHighRNASeqLibraryCount  smallint unsigned not null default 0
        comment 'number of RNA-Seq libraries in this experiment that produced this call as absent high quality',
    absentLowRNASeqLibraryCount   smallint unsigned not null default 0
        comment 'number of RNA-Seq libraries in this experiment that produced this call as absent low quality',
    rnaSeqExperimentCallDirection enum('present', 'absent') not null
        comment 'Inferred direction for this call based on this experiment ("present" libraries always win over "absent" libraries)',
    rnaSeqExperimentCallQuality enum('poor quality', 'high quality') not null
        comment 'Inferred quality for this call based on this experiment (from all libraries, "present high" > "present low" > "absent high" > "absent low"). Value "poor quality" instead of "low quality" for historical reasons.'
) engine = innodb
comment = 'This table stores information about expression calls produced from RNA-Seq experiments, that is then used in Bgee to compute global summary expression calls and qualities.';

-- ****************************************************
-- RAW DIFFERENTIAL EXPRESSION ANALYSES
-- Note: dea = Differential Expression Analyses ;)
-- ****************************************************

-- several differential expression analyses can be performed
-- on the same experiment
create table differentialExpressionAnalysis (
    deaId smallint unsigned not null,
    detectionType enum('Limma - MCM'),
-- defines whether different organs at a same (broad) developmental stage
-- were compared ('anatomy'), or a same organ at different developmental stages
-- ('development')
    comparisonFactor enum('anatomy', 'development'),
-- microarrayExperimentId and rnaSeqExperimentId cannot be both null, ot both not null
-- at the same time. We use these fields rather than an association table,
-- because a DEA can belong to only one experiment, and because this would make
-- one join less needed.
    microarrayExperimentId varchar(70) default null,
    rnaSeqExperimentId varchar(70) default null
) engine = innodb;

-- a DEA can only be performed by comparing different conditions
-- (a condition being an organ at a developmental stage), with each condition
-- represented by several replicates. Such a group of replicates of a same condition
-- in a same DEA is a 'deaSampleGroup'.
-- While it would be possible to determine the condition (anatEntityId + stageId)
-- of a deaSampleGroup by looking at the individual samples (for instance,
-- looking at the condition of an affymetrixChip member of a deaSampleGroup),
-- this information is also present in this table (see anatEntityId and stageId fields).
-- This is because, for the sake of performing the analyses, too granular
-- developmental stages can be mapped to a broader parent stage (for instance,
-- mapping '24 yo human' to 'young adult'), otherwise the analyses could be
-- meaningless (e.g., performing a DEA on '24 yo human' vs. '25 yo human').
-- So the anatEntityId and stageId in this table can actually be different than
-- the annotated anatDevId and stageId of the samples (meaning, different than
-- in the table affymetrixChip or rnaSeqLibrary).
-- As of Bgee 13, a deaSampleGroup can either be a group of affymetrixChips,
-- or a group of rnaSeqLibraries. Their related samples will then be find
-- respectively in deaSampleGroupToAffymetrixChip, or deaSampleGroupToRnaSeqLibrary.
-- this can be determined by checking in the table differentialExpressionAnalysis
-- the fields microarrayExperimentId and rnaSeqExperimentId, to determine whether
-- the DEA was using Affymetrix, or RNA-Seq.
create table deaSampleGroup (
    deaSampleGroupId mediumint unsigned not null,
    deaId smallint unsigned not null,
    conditionId mediumint unsigned not null
) engine = innodb;

-- An association table to link an affymetrixChip to the deaSampleGroup it belongs to.
-- A same chip can be part of several groups, for instance if it was use for DEAs
-- with different comparisonFactors. But all the affymetrixChips inside a deaSampleGroup
-- are unique
create table deaSampleGroupToAffymetrixChip (
    deaSampleGroupId mediumint unsigned not null,
    bgeeAffymetrixChipId smallint unsigned not null
) engine = innodb;

-- An association table to link a rnaSeqLibrary to the deaSampleGroup it belongs to.
-- A same library can be part of several groups, for instance if it was use for DEAs
-- with different comparisonFactors. But all the rnaSeqLibraries inside a deaSampleGroup
-- are unique
create table deaSampleGroupToRnaSeqLibrary (
    deaSampleGroupId mediumint unsigned not null,
    rnaSeqLibraryId varchar(70) not null
) engine = innodb;

-- differentialExpressionAnalysisProbesetsSummary
-- a line in this table is a summary of a set of probesets, used for the
-- differential expression analysis, belonging to different
-- affymetrix chips, corresponding to one group of chips
create table deaAffymetrixProbesetSummary (
-- deaAffymetrixProbesetSummaryId corresponds to the IDs of the probesets used for this summary
-- (all of them have the same of course). These probesets belong to the affymetrix chips, retrieved using the field `deaChipsGroupId`
-- and the table `deaChipsGroupToAffymetrixChip`
    deaAffymetrixProbesetSummaryId varchar(70) not null,
    deaSampleGroupId mediumint unsigned not null,
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
    foldChange decimal(7,2) not null default 0,
    differentialExpressionId int unsigned,
-- Warning, qualities must be ordered, the index in the enum is used in many queries
-- 'not expressed' = gene never seen as 'expressed' in the conditions studied ('marginal' is not considered)
-- 'no diff expression' = gene has expression, but no significant fold change observed
    differentialExpressionAffymetrixData enum('no data', 'not expressed', 'no diff expression', 'poor quality', 'high quality') default 'no data',
-- p-value adjusted by Benjamini-Hochberg procedure
-- "number of digits to the right of the decimal point (the scale). It has a range of 0 to 30"
    deaRawPValue decimal(31, 30) unsigned not null default 1,
-- excluded if not expressed in ALL samples in a given analysis
-- (it is not excluded if expressed in at least one condition)
    reasonForExclusion enum('not excluded', 'not expressed') not null default 'not excluded'
) engine = innodb;

-- deaRNASeqSummary
-- a line in this table is a summary of a set of RNA-Seq results, used for the
-- differential expression analysis, belonging to different runs, corresponding to one group of runs
create table deaRNASeqSummary (
    geneSummaryId mediumint unsigned not null,
    deaSampleGroupId mediumint unsigned not null,
    foldChange decimal(7,2) not null default 0,
    differentialExpressionId int unsigned,
-- Warning, qualities must be ordered, the index in the enum is used in many queries
-- 'not expressed' = gene never seen as 'expressed' in the conditions studied ('marginal' is not considered)
-- 'no diff expression' = gene has expression, but no significant fold change observed
    differentialExpressionRNASeqData enum('no data', 'not expressed', 'no diff expression', 'poor quality', 'high quality') default 'no data',
-- p-value adjusted by Benjamini-Hochberg procedure
-- "number of digits to the right of the decimal point (the scale). It has a range of 0 to 30"
    deaRawPValue decimal(31, 30) unsigned not null default 1,
-- excluded if not expressed in ALL samples in a given analysis
-- (it is not excluded if expressed in at least one condition)
    reasonForExclusion enum('not excluded', 'not expressed') not null default 'not excluded'
) engine = innodb;

-- ****************************************************
-- SUMMARY EXPRESSION CALLS
-- ****************************************************

-- This table is a summary of expression calls for a given gene-condition
-- gene - anatomical entity - developmental stage - sex- strain, over all the experiments
-- for all data types with no propagation nor experiment expression summary.
create table expression (
    expressionId int unsigned not null COMMENT 'Internal expression ID, not stable between releases.',
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID, not stable between releases.',
    conditionId mediumint unsigned not null COMMENT 'ID of condition in the related condition table ("cond"), not stable between releases.'
) engine = innodb
comment = 'This table is a summary of expression calls for a given gene-condition (anatomical entity - developmental stage - sex- strain), over all the experiments and data types, with no propagation nor experiment expression summary.';

-- This table is a summary of expression calls for a given gene-condition
-- gene - anatomical entity - developmental stage - sex- strain, over all the experiments
-- for all data types, with all data propagated and reconciled, with experiment expression summaries computed.
-- DESIGN note: this table uses an ugly design with enumerated columns. For a discussion about this decision,
-- see http://stackoverflow.com/q/42781299/1768736
create table globalExpression (
    globalExpressionId int unsigned not null COMMENT 'Internal expression ID, not stable between releases.',
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID, not stable between releases.',
    globalConditionId mediumint unsigned not null COMMENT 'ID of condition in the related condition table ("cond"), not stable between releases.',

-- ** OBSERVED DATA STATES ** --
-- It is not enough to only know that there are some data in the condition itself
-- (see 'SELF' expression summaries below), because we need to distinguish between
-- data propagated along a condition parameter (e.g., a developmental stage),
-- but observed in another (e.g., an anat. entity).
-- Note that these enum values must stay in sync with org.bgee.model.dao.api.expressiondata.DAOPropagationState.
-- And these enum are null when the expression calls were propagated without taking into account
-- the related condition parameter.
-- Note that for EST data, there is no propagation of ABSENT calls from parent conditions,
-- so this simplifies its enum.
    estAnatEntityPropagationState ENUM('self', 'descendant', 'self and descendant') COMMENT 'The origin of the propagated EST data related to the anatomical entity of the related condition. If null, it means that anatomical entities were not considered in this propagated call, or that the call is not supported by any EST data.',
    estStagePropagationState ENUM('self', 'descendant', 'self and descendant') COMMENT 'The origin of the propagated EST data related to the developmental stage of the related condition. If null, it means that developmental stages were not considered in this propagated call, or that the call is not supported by any EST data.',
    estConditionObservedData BOOLEAN COMMENT 'Whether some EST data were observed in this condition itself. If null, it means that the call is not supported by any EST data. This field is redundant as compared to the "self" experiment counts below, but is more practical to use.',

    affymetrixAnatEntityPropagationState ENUM('all', 'self', 'ancestor', 'descendant', 'self and ancestor',
    'self and descendant', 'ancestor and descendant') COMMENT 'The origin of the propagated Affymetrix data related to the anatomical entity of the related condition. If null, it means that anatomical entities were not considered in this propagated call, or that the call is not supported by any Affymetrix data.',
    affymetrixStagePropagationState ENUM('all', 'self', 'ancestor', 'descendant', 'self and ancestor',
    'self and descendant', 'ancestor and descendant') COMMENT 'The origin of the propagated Affymetrix data related to the developmental stage of the related condition. If null, it means that developmental stages were not considered in this propagated call, or that the call is not supported by any Affymetrix data.',
    affymetrixConditionObservedData BOOLEAN COMMENT 'Whether some Affymetrix data were observed in this condition itself. If null, it means that the call is not supported by any Affymetrix data. This field is redundant as compared to the "self" experiment counts below, but is more practical to use.',

    inSituAnatEntityPropagationState ENUM('all', 'self', 'ancestor', 'descendant', 'self and ancestor',
    'self and descendant', 'ancestor and descendant') COMMENT 'The origin of the propagated in situ hybridization data related to the anatomical entity of the related condition. If null, it means that anatomical entities were not considered in this propagated call, or that the call is not supported by any in situ hybridization data.',
    inSituStagePropagationState ENUM('all', 'self', 'ancestor', 'descendant', 'self and ancestor',
    'self and descendant', 'ancestor and descendant') COMMENT 'The origin of the propagated in situ hybridization data related to the developmental stage of the related condition. If null, it means that developmental stages were not considered in this propagated call, or that the call is not supported by any in situ hybridization data.',
    inSituConditionObservedData BOOLEAN COMMENT 'Whether some in situ hybridization data were observed in this condition itself. If null, it means that the call is not supported by any in situ hybridization data. This field is redundant as compared to the "self" experiment counts below, but is more practical to use.',

    rnaSeqAnatEntityPropagationState ENUM('all', 'self', 'ancestor', 'descendant', 'self and ancestor',
    'self and descendant', 'ancestor and descendant') COMMENT 'The origin of the propagated RNA-Seq data related to the anatomical entity of the related condition. If null, it means that anatomical entities were not considered in this propagated call, or that the call is not supported by any RNA-Seq data.',
    rnaSeqStagePropagationState ENUM('all', 'self', 'ancestor', 'descendant', 'self and ancestor',
    'self and descendant', 'ancestor and descendant') COMMENT 'The origin of the propagated RNA-Seq data related to the developmental stage of the related condition. If null, it means that developmental stages were not considered in this propagated call, or that the call is not supported by any RNA-Seq data.',
    rnaSeqConditionObservedData BOOLEAN COMMENT 'Whether some RNA-Seq data were observed in this condition itself. If null, it means that the call is not supported by any RNA-Seq data. This field is redundant as compared to the "self" experiment counts below, but is more practical to use.',

-- ** EXPRESSION SUMMARIES **
-- Note: EST data are not used to produce no-expression calls
    estLibPresentHighSelfCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of EST libraries showing expression of this gene in this condition (not taking into account sub-conditions) with a high quality.',
    estLibPresentLowSelfCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of EST libraries showing expression of this gene in this condition (not taking into account sub-conditions) with a low quality.',
    estLibPresentHighDescendantCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of EST libraries showing expression of this gene, solely in the sub-conditions of this condition, with a high quality.',
    estLibPresentLowDescendantCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of EST libraries showing expression of this gene, solely in the sub-conditions of this condition, with a low quality.',
    estLibPresentHighTotalCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of EST libraries showing expression of this gene in this condition or in sub-conditions with a high quality.',
    estLibPresentLowTotalCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of EST libraries showing expression of this gene in this condition or in sub-conditions with a low quality.',
    estLibPropagatedCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of EST libraries used to show presence of expression (low or high) in sub-conditions of this condition.',

    affymetrixExpPresentHighSelfCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing expression of this gene in this condition (not taking into account sub-conditions) with a high quality.',
    affymetrixExpPresentLowSelfCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing expression of this gene in this condition (not taking into account sub-conditions) with a low quality.',
    affymetrixExpAbsentHighSelfCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing absence of expression of this gene in this condition (not taking into account parent conditions) with a high quality.',
    affymetrixExpAbsentLowSelfCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing absence of expression of this gene in this condition (not taking into account parent conditions) with a low quality.',
    affymetrixExpPresentHighDescendantCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing expression of this gene, solely in the sub-conditions of this condition, with a high quality.',
    affymetrixExpPresentLowDescendantCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing expression of this gene, solely in the sub-conditions of this condition, with a low quality.',
    affymetrixExpAbsentHighParentCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing absence of expression of this gene, solely in the valid parent conditions of this condition, with a high quality.',
    affymetrixExpAbsentLowParentCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing absence of expression of this gene, solely in the valid parent conditions of this condition, with a low quality.',
    affymetrixExpPresentHighTotalCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing expression of this gene in this condition or in sub-conditions with a high quality.',
    affymetrixExpPresentLowTotalCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing expression of this gene in this condition or in sub-conditions with a low quality.',
    affymetrixExpAbsentHighTotalCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing absence of expression of this gene in this condition or valid parent conditions with a high quality.',
    affymetrixExpAbsentLowTotalCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments showing absence of expression of this gene in this condition or valid parent conditions with a low quality.',
    affymetrixExpPropagatedCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of Affymetrix experiments used either to show presence of expression (low or high) in sub-conditions of this condition, or absence of expression (low or high) in parent conditions of this condition.',

    inSituExpPresentHighSelfCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing expression of this gene in this condition (not taking into account sub-conditions) with a high quality.',
    inSituExpPresentLowSelfCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing expression of this gene in this condition (not taking into account sub-conditions) with a low quality.',
    inSituExpAbsentHighSelfCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing absence of expression of this gene in this condition (not taking into account parent conditions) with a high quality.',
    inSituExpAbsentLowSelfCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing absence of expression of this gene in this condition (not taking into account parent conditions) with a low quality.',
    inSituExpPresentHighDescendantCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing expression of this gene, solely in the sub-conditions of this condition, with a high quality.',
    inSituExpPresentLowDescendantCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing expression of this gene, solely in the sub-conditions of this condition, with a low quality.',
    inSituExpAbsentHighParentCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing absence of expression of this gene, solely in the valid parent conditions of this condition, with a high quality.',
    inSituExpAbsentLowParentCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing absence of expression of this gene, solely in the valid parent conditions of this condition, with a low quality.',
    inSituExpPresentHighTotalCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing expression of this gene in this condition or in sub-conditions with a high quality.',
    inSituExpPresentLowTotalCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing expression of this gene in this condition or in sub-conditions with a low quality.',
    inSituExpAbsentHighTotalCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing absence of expression of this gene in this condition or valid parent conditions with a high quality.',
    inSituExpAbsentLowTotalCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments showing absence of expression of this gene in this condition or valid parent conditions with a low quality.',
    inSituExpPropagatedCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of in situ hybridization experiments used either to show presence of expression (low or high) in sub-conditions of this condition, or absence of expression (low or high) in parent conditions of this condition.',

    rnaSeqExpPresentHighSelfCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing expression of this gene in this condition (not taking into account sub-conditions) with a high quality.',
    rnaSeqExpPresentLowSelfCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing expression of this gene in this condition (not taking into account sub-conditions) with a low quality.',
    rnaSeqExpAbsentHighSelfCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing absence of expression of this gene in this condition (not taking into account parent conditions) with a high quality.',
    rnaSeqExpAbsentLowSelfCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing absence of expression of this gene in this condition (not taking into account parent conditions) with a low quality.',
    rnaSeqExpPresentHighDescendantCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing expression of this gene, solely in the sub-conditions of this condition, with a high quality.',
    rnaSeqExpPresentLowDescendantCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing expression of this gene, solely in the sub-conditions of this condition, with a low quality.',
    rnaSeqExpAbsentHighParentCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing absence of expression of this gene, solely in the valid parent conditions of this condition, with a high quality.',
    rnaSeqExpAbsentLowParentCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing absence of expression of this gene, solely in the valid parent conditions of this condition, with a low quality.',
    rnaSeqExpPresentHighTotalCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing expression of this gene in this condition or in sub-conditions with a high quality.',
    rnaSeqExpPresentLowTotalCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing expression of this gene in this condition or in sub-conditions with a low quality.',
    rnaSeqExpAbsentHighTotalCount SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing absence of expression of this gene in this condition or valid parent conditions with a high quality.',
    rnaSeqExpAbsentLowTotalCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments showing absence of expression of this gene in this condition or valid parent conditions with a low quality.',
    rnaSeqExpPropagatedCount  SMALLINT UNSIGNED NOT NULL DEFAULT 0 COMMENT 'Count of RNA-Seq experiments used either to show presence of expression (low or high) in sub-conditions of this condition, or absence of expression (low or high) in parent conditions of this condition.',

-- ** RANKS **
-- For RNA-Seq data: mean ranks before normalization between data types and conditions.
-- Used for convenience during rank computations. It corresponds to the following:
-- gene ranks are computed for each sample, then a mean is computed for each gene and condition
-- of the expression table, weighted by the number of distinct ranks in each sample.
    rnaSeqMeanRank decimal(9, 2) unsigned COMMENT 'RNA-Seq mean rank for this gene-condition before normalization over all data types, conditions and species.',
-- For Affymetrix:
-- mean ranks *after within-datatype normalization*, before normalization between data types and conditions.
-- Used for convenience during rank computations. It corresponds to the following:
-- ranks are computed for each sample, then "normalized" between samples in a same condition
-- of the expression table ("within-datatype normalization", based on the genomic coverage of each chip type);
-- then a mean is computed for each gene and condition, weighted by the number of distinct ranks
-- in each sample.
    affymetrixMeanRank decimal(9, 2) unsigned COMMENT 'Affymetrix mean rank for this gene-condition, after normalization between different chip types, but before normalization over all data types, conditions and species.',
-- For EST and in situ data: ranks before normalization between data types and conditions.
-- Used for convenience during rank computations. It corresponds to the following:
-- For each condition of the expression table, all data are pooled together; they are not first
-- analyzed independently per libraries or experiments, as for Affymetrix and RNA-Seq data.
-- This is because the genomic coverage of EST or in situ experiments is usually very low,
-- and highly variable. Genes are ranked based on number of ESTs or of in situ evidence in each condition.
-- They are ranked using "dense ranking" instead of fractional ranking.
    estRank decimal(9, 2) unsigned COMMENT 'EST rank for this gene-condition before normalization over all data types, conditions and species. All EST libraries in a same condition are pulled together, so there is no concept of "mean", only a single rank is computed from EST data for each gene-condition.',
    inSituRank decimal(9, 2) unsigned COMMENT 'In situ hybridization rank for this gene-condition before normalization over all data types, conditions and species. All in situ evidence in a same condition are pulled together, so there is no concept of "mean", only a single rank is computed from in situ data for each gene-condition.',

-- All ranks are normalized between all data types and conditons, this is what we use to compute
-- the global mean rank of a gene in a condition. Basically, the max rank over all data types
-- and all conditions is retrieved, and used to normalize all ranks.
-- normalized rank = rank * (max of max rank over all conditions and data types) / (max rank for this condition and data type)
    rnaSeqMeanRankNorm decimal(9, 2) unsigned COMMENT 'RNA-Seq normalized mean rank for this gene-condition after normalization over all data types, conditions and species, computed from the field rnaSeqMeanRank, and rnaSeqMaxRank in the related condition table.',
    affymetrixMeanRankNorm decimal(9, 2) unsigned COMMENT 'Affymetrix normalized mean rank for this gene-condition after normalization over all data types, conditions and species, computed from the field affymetrixMeanRank, and affymetrixMaxRank in the related condition table.',
-- For EST and in situ, the rank is not a mean
    estRankNorm decimal(9, 2) unsigned COMMENT 'EST normalized rank for this gene-condition after normalization over all data types, conditions and species, computed from the field estRank, and estMaxRank in the related condition table.',
    inSituRankNorm decimal(9, 2) unsigned COMMENT 'In situ hybridization normalized rank for this gene-condition after normalization over all data types, conditions and species, computed from the field inSituRank, and inSituMaxRank in the related condition table.',

-- For Affymetrix and RNA-Seq data: sum of the number of distinct ranks in each sample
-- where this gene is considered, in this condition and data type (for RNA-Seq: the same set of genes
-- is considered in all conditions, so these values are all the same for all genes in a same condition-species;
-- for Affymetrix, it depends on the chip types, so it can vary between genes of same condition-species ).
-- Distinct ranks in samples are used to weight the mean rank of genes for each data type and condition.
-- By storing the sum of the distinct rank count, we will be able to compute the weighted mean
-- over all data types in a condition.
-- XXX: shoud we store this information in the condition table for RNA-Seq?
-- Or maybe we shouldn't constrain to have the same genomic coverage in all libraries of a condition?
--
-- For EST and in situ data, this is irrelevant as we pool all data for a same condition together,
-- and use dense ranking instead of fractional ranking. As a result, the max rank in each condition
-- is used for weighted mean computation between data types.
    rnaSeqDistinctRankSum int unsigned COMMENT 'Factor used to weight the RNA-Seq normalized mean rank (rnaSeqMeanRankNorm), to compute a global weighted mean rank between all data types. Corresponds to the sum of distinct ranks in each library mapped to this condition. Note that for EST and in situ data, the max rank found in the related condition table is instead used to compute the weighted mean between data types.',
    affymetrixDistinctRankSum int unsigned COMMENT 'Factor used to weight the Affymetrix normalized mean rank (affymetrixMeanRankNorm), to compute a global weighted mean rank between all data types. Corresponds to the sum of distinct ranks in each chip mapped to this condition. Note that for EST and in situ data, the max rank found in the related condition table is instead used to compute the weighted mean between data types.',

-- Same fields, but dedicated to "global" ranks, computed by taking into account
-- all data in a condition, but also all data in its descendant conditions.
    rnaSeqGlobalMeanRank decimal(9, 2) unsigned COMMENT 'RNA-Seq global mean rank for this gene in this condition and all its descendant conditions, before normalization over all data types, conditions and species.',
    affymetrixGlobalMeanRank decimal(9, 2) unsigned COMMENT 'Affymetrix global mean rank for this gene in this condition and all its descendant conditions, after normalization between different chip types, but before normalization over all data types, conditions and species.',
    estGlobalRank decimal(9, 2) unsigned COMMENT 'EST global rank for this gene in this condition and all its descendant conditions, before normalization over all data types, conditions and species. All EST libraries in a same condition in this condition and its descendant conditions are pulled together, so there is no concept of "mean", only a single rank is computed from EST data for each gene-condition.',
    inSituGlobalRank decimal(9, 2) unsigned COMMENT 'In situ hybridization global rank for this gene in this condition and all its descendant conditions, before normalization over all data types, conditions and species. All in situ evidence in a same condition and its descendant conditions are pulled together, so there is no concept of "mean", only a single rank is computed from in situ data for each gene-condition.',

    rnaSeqGlobalMeanRankNorm decimal(9, 2) unsigned COMMENT 'RNA-Seq normalized mean rank for this gene in this condition and all its descendant conditions, after normalization over all data types, conditions and species, computed from the field rnaSeqMeanRank, and rnaSeqMaxRank in the related condition table.',
    affymetrixGlobalMeanRankNorm decimal(9, 2) unsigned COMMENT 'Affymetrix normalized mean rank for this gene in this condition and all its descendant conditions, after normalization over all data types, conditions and species, computed from the field affymetrixMeanRank, and affymetrixMaxRank in the related condition table.',
    estGlobalRankNorm decimal(9, 2) unsigned COMMENT 'EST normalized rank for this gene in this condition and all its descendant conditions, after normalization over all data types, conditions and species, computed from the field estRank, and estMaxRank in the related condition table.',
    inSituGlobalRankNorm decimal(9, 2) unsigned COMMENT 'In situ hybridization normalized rank for this gene in this condition and all its descendant conditions, after normalization over all data types, conditions and species, computed from the field inSituRank, and inSituMaxRank in the related condition table.',

    rnaSeqGlobalDistinctRankSum int unsigned COMMENT 'Factor used to weight the RNA-Seq normalized global mean rank (rnaSeqGlobalMeanRankNorm), to compute a global weighted mean rank between all data types. Corresponds to the sum of distinct ranks in each library mapped to this condition and all its descendant conditions. Note that for EST and in situ data, the global max rank found in the related condition table is instead used to compute the weighted mean between data types.',
    affymetrixGlobalDistinctRankSum int unsigned COMMENT 'Factor used to weight the Affymetrix normalized global mean rank (affymetrixGlobalMeanRankNorm), to compute a global weighted mean rank between all data types. Corresponds to the sum of distinct ranks in each chip mapped to this condition and all its descendant conditions. Note that for EST and in situ data, the global max rank found in the related condition table is instead used to compute the weighted mean between data types.'
) engine = innodb
comment = 'This table is a summary of expression calls for a given gene-condition (anatomical entity - developmental stage - sex- strain), over all the experiments and data types, with all data propagated and reconciled, and with experiment expression summaries computed.';

-- ****************************************************
-- SUMMARY DIFF EXPRESSION CALLS
-- ****************************************************

create table differentialExpression (
    differentialExpressionId int unsigned not null,
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
    conditionId mediumint unsigned not null,
-- defines whether different organs at a same (broad) developmental stage
-- were compared ('anatomy'), or a same organ at different developmental stages
-- ('development')
    comparisonFactor enum('anatomy', 'development'),
-- *** Affymetrix ***
-- the diff expression call generated by Affymetrix
-- 'not expressed' = gene never seen as 'expressed' in the conditions studied ('marginal' is not considered)
-- 'no diff expression' = gene has expression, but no significant fold change observed
    diffExprCallAffymetrix enum('no data', 'not expressed', 'no diff expression', 'under-expression', 'over-expression') not null default 'no data',
-- confidence in the call generated by Affymetrix data
-- 'no data' is redundant but it is kept to keep the same indexes for all data states (for instance, rnaSeqData in expression table)
    diffExprAffymetrixData enum('no data', 'poor quality', 'high quality') default 'no data',
-- among all the analyses using Affymetrix comparing this condition, best p-value associated to this call
-- "number of digits to the right of the decimal point (the scale). It has a range of 0 to 30"
    bestPValueAffymetrix decimal(31, 30) unsigned not null default 1,
-- number of analyses using Affymetrix data where the same call is found
    consistentDEACountAffymetrix smallint unsigned not null default 0,
-- number of analyses using Affymetrix data where a different call is found
    inconsistentDEACountAffymetrix smallint unsigned not null default 0,
-- *** RNA-Seq ***
-- the diff expression call generated by RNA-Seq
-- 'not expressed' = gene never seen as 'expressed' in the conditions studied ('marginal' is not considered)
-- 'no diff expression' = gene has expression, but no significant fold change observed
    diffExprCallRNASeq enum('no data','not expressed', 'no diff expression', 'under-expression', 'over-expression') not null default 'no data',
-- confidence in the call generated by RNA-Seq data
-- 'no data' is redundant but it is kept to keep the same indexes for all data states (for instance, rnaSeqData in expression table)
    diffExprRNASeqData enum('no data', 'poor quality', 'high quality') default 'no data',
-- among all the analyses using RNA-Seq comparing this condition, best p-value associated to this call
-- "number of digits to the right of the decimal point (the scale). It has a range of 0 to 30"
    bestPValueRNASeq decimal(31, 30) unsigned not null default 1,
-- number of analyses using RNA-Seq data where the same call is found
    consistentDEACountRNASeq smallint unsigned not null default 0,
-- number of analyses using RNA-Seq data where a different call is found
    inconsistentDEACountRNASeq smallint unsigned not null default 0
) engine = innodb;

-- this version of the diff expression table is not considered as of Bgee 13
/*create table differentialExpression (
    differentialExpressionId int unsigned not null,
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
    conditionId mediumint unsigned not null,
-- defines whether different organs at a same (broad) developmental stage
-- were compared ('anatomy'), or a same organ at different developmental stages
-- ('development')
    comparisonFactor enum('anatomy', 'development'),
-- Warning, differentialExpressionCall must be ordered this way, the index in the enum
-- is used in many queries
    differentialExpressionCall enum('no diff expression', 'under-expression', 'over-expression'),
-- the maximum number of conditions compared for which this differential expression call
-- is valid. For instance, if a differential expression analysis comparing 3 conditions
-- generated a call for a given gene-organ-stage, and another analysis comparing
-- 6 conditions generated another call for the same gene-organ-stage (with different
-- direction and/or qualities), then maxNumberOfConditions will be 3 for the call
-- generated by the first analysis, and 6 for the other call.
--
-- But if the two analyses were generating the same call, then they would be only
-- one call in this table for the given gene-organ-stade, with a maxNumberOfConditions
-- equals to 6
--
-- default is 3, because as of Bgee 13, this is the minimum number of conditions
-- to perform a diff expression analysis
--
-- Examples of queries:
--
-- query to retrieve diff expression calls for a given Gene with no minimum number
-- of conditions compared requested:
-- select * from differentialExpression as t1 inner join
-- (
-- select geneId, organId, stageId, comparisonFactor, min(maxNumberOfConditions) as min
-- from differentialExpression where geneId = ? group by geneId, organId, stageId,
-- comparisonFactor)
-- ) as t2 on t1.geneId = t2.geneId and t1.organId = t2.organId and t1.stageId = t2.stageId
-- and t1.comparisonFactor = t2.comparisonFactor and t1.maxNumberOfConditions = t2.min
-- where t1.geneId = ?;
--
-- Alternatively: TO TEST, IT SEEMS WRONG
--
-- select * from differentialExpression as t1
-- where t1.geneId = ? and t1.maxNumberOfConditions =
-- (select min(maxNumberOfConditions) from differentialExpression as t2 where
-- t2.geneId = t1.geneId and t2.organId = t1.organId and t2.stageId = t1.geneId and
-- t2.comparisonFactor = t1.comparisonFactor);
--
-- Example of query to select the calls with the maximum number of conditions compared
-- for a given gene-organ-stage, with no minimum defined (select only the "best" calls):
--
-- select * from differentialExpression as t1 inner join
-- (
-- select geneId, organId, stageId, comparisonFactor, max(maxNumberOfConditions) as max
-- from differentialExpression where geneId = ? group by geneId, organId, stageId,
-- comparisonFactor)
-- ) as t2 on t1.geneId = t2.geneId and t1.organId = t2.organId and t1.stageId = t2.stageId
-- and t1.comparisonFactor = t2.comparisonFactor and t1.maxNumberOfConditions = t2.max
-- where t1.geneId = ?;
   maxNumberOfConditions smallint unsigned not null default 3,
-- Warning, qualities must be ordered this way, the index in the enum is used in many queries
    differentialExpressionAffymetrixData enum('no data', 'poor quality', 'high quality') default 'no data',
    differentialExpressionRnaSeqData enum('no data', 'poor quality', 'high quality') default 'no data'
) engine = innodb;*/


-- select((select count(1) from rnaSeqExperiment) + (select count(1) from rnaSeqLibrary) + (select count(1) from rnaSeqResults) + (select count(1) from rnaSeqExperimentToKeyword) + (select count(1) from affymetrixChip) + (select count(1) from affymetrixProbeset) + (select count(1) from author) + (select count(1) from chipType) + (select count(1) from dataSource) + (select count(1) from dataType) + (select count(1) from deaAffymetrixProbesetSummary) + (select count(1) from deaChipsGroup) + (select count(1) from deaChipsGroupToAffymetrixChip) + (select count(1) from detectionType) + (select count(1) from differentialExpression) + (select count(1) from differentialExpressionAnalysis) + (select count(1) from differentialExpressionAnalysisType) + (select count(1) from estLibrary) + (select count(1) from estLibraryToKeyword) + (select count(1) from expressedSequenceTag) + (select count(1) from expression) + (select count(1) from gene) + (select count(1) from geneBioType) + (select count(1) from geneFamily) + (select count(1) from geneFamilyPredictionMethod) + (select count(1) from geneNameSynonym) + (select count(1) from geneOntologyDescendants) + (select count(1) from geneOntologyTerm) + (select count(1) from geneToTerm) + (select count(1) from geneXRef) + (select count(1) from globalExpression) + (select count(1) from globalExpressionToExpression) + (select count(1) from hogDescendants) + (select count(1) from hogExpression) + (select count(1) from hogExpressionSummary) + (select count(1) from hogExpressionToExpression) + (select count(1) from hogNameSynonym) + (select count(1) from hogRelationship) + (select count(1) from hogXRef) + (select count(1) from homologousOrgansGroup) + (select count(1) from inSituEvidence) + (select count(1) from inSituExperiment) + (select count(1) from inSituExperimentToKeyword) + (select count(1) from inSituSpot) + (select count(1) from keyword) + (select count(1) from metaStage) + (select count(1) from metaStageNameSynonym) + (select count(1) from microarrayExperiment) + (select count(1) from microarrayExperimentToKeyword) + (select count(1) from normalizationType) + (select count(1) from organ) + (select count(1) from organDescendants) + (select count(1) from organNameSynonym) + (select count(1) from organRelationship) + (select count(1) from species) + (select count(1) from stage) + (select count(1) from stageNameSynonym) + (select count(1) from stageXRef));

-- ******************************************
-- AVAILABLE FILES FOR DOWNLOAD
-- ******************************************
-- see (https://gitlab.sib.swiss/Bgee/bgee_apps/issues/31)
create table downloadFile (
  downloadFileId mediumint unsigned not null,
-- path relative to the root of the download file directory, including file name
  downloadFileRelativePath varchar(255) not null,
-- currently, just the name of the file
  downloadFileName varchar(255) not null,
  downloadFileDescription text,
  downloadFileCategory enum("expr_simple", "expr_complete", "diff_expr_anatomy_complete", "diff_expr_anatomy_simple"
   , "diff_expr_dev_complete", "diff_expr_dev_simple", "ortholog",
   "affy_annot","rnaseq_annot","affy_data","rnaseq_data"),
  speciesDataGroupId mediumint unsigned not null,
  downloadFileSize int unsigned not null,
  downloadFileConditionParameters set('anatomicalEntity', 'developmentalStage')
) engine = innodb;

-- *****************************************
-- SPECIES CONFIG
-- *****************************************
-- a set of species (containing at least one element) for which a file was generated

create table speciesDataGroup(
  speciesDataGroupId mediumint unsigned not null,
  speciesDataGroupName varchar(255) not null,
  speciesDataGroupDescription text,
-- preferred order to display speciesDataGroups
  speciesDataGroupOrder tinyint unsigned not null default 255
) engine = innodb;

create table speciesToDataGroup(
  speciesDataGroupId mediumint unsigned not null,
  speciesId mediumint unsigned not null
) engine = innodb;

