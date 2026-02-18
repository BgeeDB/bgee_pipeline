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
    authorId   smallInt unsigned not null,
    authorName varchar(255)      not null COMMENT 'Bgee team author names'
) engine = innodb;

create table dataSource (
    dataSourceId          smallInt unsigned not null,
    dataSourceName        varchar(55)       not null             COMMENT 'Data source name',
    XRefUrl               varchar(255)      not null default ''  COMMENT 'URL for cross-references to data sources',
-- path to experiment for expression data sources (ArrayExpress, GEO, NCBI, in situ databases, ...)
-- parameters such as experimentId are defined by the syntax [experimentId] for instance
    experimentUrl         varchar(100)      not null default ''  COMMENT 'URL to experiment for expression data sources',
-- path to in situ evidence for in situ databases
-- parameters such as experimentId are defined by the syntax [experimentId] for instance
    evidenceUrl           varchar(100)      not null default ''  COMMENT 'URL to evidence for expression data sources',
-- url to the home page of the ressource
    baseUrl               varchar(100)      not null default ''  COMMENT 'URL to the home page of data sources',
    releaseDate           date                  null             COMMENT 'Date of data source used',
-- e.g.: Ensembl 87, git version xxx
    releaseVersion        varchar(35)       not null default ''  COMMENT 'Version of data source used',
    dataSourceDescription varchar(200)                           COMMENT 'Description of data source',
-- to define if this dataSource should be displayed on the page listing data sources
    toDisplay             boolean           not null default 0   COMMENT 'Display this data source in listing data source page?',
-- a cat to organize the display
    category              enum('', 'Genomics database', 'Proteomics database',
                               'In situ data source', 'RNA-Seq data source',
                               'Single-cell RNA-Seq data source', 'Ontology') COMMENT 'Data source category to organize the display',
-- to organize the display. Default value is the highest value, so that this field is the last to be displayed
    displayOrder          tinyint unsigned  not null default 255 COMMENT 'Data source display ordering'
) engine = innodb;

create table dataSourceToSpecies (
    dataSourceId smallInt  unsigned  not null COMMENT 'Data source id',
    speciesId    mediumint unsigned  not null COMMENT 'NCBI species taxon id',
    dataType     enum('in situ', 'rna-seq', 'single-cell RNA-Seq') not null COMMENT 'Data type',
    infoType     enum('data', 'annotation') not null COMMENT 'Information type'
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
    genomeAssemblyXRef varchar(250) not null default '' COMMENT 'XRef to the genome assembly',
    dataSourceId smallInt unsigned not null COMMENT 'source for genome information',
-- ID of the species whose the genome was used for this species. This is used
-- when a genome is not in Ensembl. For instance, for bonobo (ID 9597), we use the chimp
-- genome (ID 9598), because bonobo is not in Ensembl.
-- We don't use a foreign key constraint here, because maybe the species whose the genome
-- was used does not have any data in Bgee, and thus is not in the taxon table.
-- If the correct genome of the species was used, the value of this field is 0.
    genomeSpeciesId mediumint unsigned not null default 0 COMMENT 'NCBI species taxon id used for mapping (0 if the same species)',
    devOntologyXRef varchar(250) not null default '' COMMENT 'XRef to the developmental stage ontology of that species'
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

-- Note:
-- * query to obtain list of tissues:
-- we select the list of terms that are descendants of 'UBERON:0001062 anatomical entity',
-- and that are not part of the cell type graph. Indeed, a term such as
-- 'CL:0002252 epithelial cell of esophagus' is both a cell type,
-- and part of esophagus. We don't want to retrieve those. The query is:
-- ```
-- select distinct t1.anatEntityId from anatEntity as t1
-- inner join anatEntityRelation as t2 on t1.anatEntityId = t2.anatEntitySourceId
--     and t2.anatEntityTargetId = 'UBERON:0001062' and t2.relationType = 'is_a part_of'
-- left outer join anatEntityRelation as t3 on t1.anatEntityId = t3.anatEntitySourceId
--     and t3.anatEntityTargetId = 'GO:0005575' and t3.relationType = 'is_a part_of'
-- where t3.anatEntitySourceId is null;
-- ```
-- * query to obtain list of cell types:
-- simply retrieve the descendants of 'GO:0005575 cellular component'. The query is:
-- ```
-- select distinct t1.anatEntityId from anatEntity as t1
-- inner join anatEntityRelation as t2 on t1.anatEntityId = t2.anatEntitySourceId
-- where t2.anatEntityTargetId = 'GO:0005575' and relationType = 'is_a part_of';
-- ```
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

create table geneOrthologs (
    bgeeGeneId mediumint unsigned not null COMMENT 'Numeric internal gene ID used for improving performances',
    targetGeneId mediumint unsigned not null COMMENT 'Numeric internal gene ID of the orthologous gene',
    taxonId mediumint unsigned not null COMMENT 'NCBI taxon id at which orthology relation had been identified'
) engine = innodb;

create table geneParalogs (
    bgeeGeneId mediumint unsigned not null COMMENT 'Numeric internal gene ID used for improving performances',
    targetGeneId mediumint unsigned not null COMMENT 'Numeric internal gene ID of the orthologous gene',
    taxonId mediumint unsigned not null COMMENT 'NCBI taxon id of the closest parent speciation of this duplication'
) engine = innodb;

create table gene (
-- warning, maybe this bgeeGeneId will need to be changed to an 'int' when we reach around 200 species
    bgeeGeneId mediumint unsigned not null COMMENT 'Numeric internal gene ID used for improving performances',
    geneId varchar(64) not null COMMENT 'Real gene id',
    geneName varchar(255) not null default '' COMMENT 'Gene name',
    geneDescription TEXT COMMENT 'Gene description',
    speciesId mediumint unsigned not null COMMENT 'NCBI species taxon id this gene belongs to',
-- TODO: check if we should add 'not null' to geneBioTypeId.
-- This depends on pipeline. If we update biotype after insertion of gene, it's not possible to set 'not null'.
    geneBioTypeId smallint unsigned COMMENT 'Gene BioType id (type of gene)',
-- defines whether the gene ID is present in Ensembl. For some species, they are not
-- (for instance, bonobo; we use chimp genome)
    ensemblGene boolean not null default 1 COMMENT 'Is the gene in Ensembl (default) (= 1), if not (= 0)',
    seqRegionName varchar(255) not null default '' COMMENT 'Chromosomal or assembly name where this gene is',
    geneMappedToGeneIdCount tinyint unsigned not null default 1 COMMENT 'number of genes in the Bgee database with the same Ensembl gene ID. In Bgee, for some species with no genome available, we use the genome of a closely-related species, such as chimpanzee genome for analyzing bonobo data. For this reason, a same Ensembl gene ID can be mapped to several species in Bgee. The value returned here is equal to 1 when the Ensembl gene ID is uniquely used in the Bgee database.',
    expressionSummary varchar(255) not null default '' COMMENT 'Sentence generated from propagated expression calls that summarizes the expression for the anatomical entity and the celltype'
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
    conditionId           mediumint unsigned not null      COMMENT 'Internal condition ID. Each condition is species-specific',
    exprMappedConditionId mediumint unsigned not null      COMMENT 'the condition ID that should be used for insertion into the expression table: too-granular conditions (e.g., 43 yo human stage, or sexInferred=1) are mapped to less granular conditions for summary. Equal to conditionId if condition is not too granular.',
    anatEntityId          varchar(20)        not null      COMMENT 'Uberon anatomical entity ID',
    cellTypeId            varchar(20)        default null  COMMENT 'A second uberon anatomical entity ID used to manage composition of anatomical entities. Used only for single cell data for postcomposition of anatomical entity ID and cell type ID',
    stageId               varchar(20)        not null      COMMENT 'Uberon stage ID',
    speciesId             mediumint unsigned not null      COMMENT 'NCBI species taxon ID',
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

create table remapCond (
    incorrectConditionId mediumint unsigned not null,
    remappedConditionId mediumint unsigned not null
) engine = innodb COMMENT 'This table is used as an intermediary step for condition remapping, see remap_conditions.pl';

create table remapExpression (
    incorrectExpressionId int unsigned not null,
    remappedExpressionId int unsigned not null
) engine = innodb COMMENT 'This table is used as an intermediary step for condition remapping, see remap_conditions.sql';

create table globalCond (
    globalConditionId mediumint unsigned not null,
    anatEntityId          varchar(20)  COMMENT 'Uberon anatomical entity ID. Can be null in this table if this condition aggregates data according to other condition parameters (e.g., grouping all data in a same stage whatever the organ is).',
    cellTypeId            varchar(20)  default null COMMENT 'A second uberon anatomical entity ID used to manage composition of anatomical entities. Used only for single cell data for postcomposition of anatomical entity ID and cell type ID',
    stageId               varchar(20)  COMMENT 'Uberon stage ID. Can be null in this table if this condition aggregates data according to other condition parameters (e.g., grouping all data in a same organ whatever the dev. stage is).',
    speciesId             mediumint unsigned not null COMMENT 'NCBI species taxon ID',
-- NA: not available from source information
-- not annotated: information not captured by Bgee
-- If an ENUM column is declared NOT NULL, its default value is the first element of the list
-- In this table, only 'any' is used to replace 'not annotated', 'NA', 'mixed'
-- and also represents the propagation of calls along the sex 'ontology'.
    sex enum('any', 'hermaphrodite', 'female', 'male'),
-- For now, strains are captured as free-text format, only 4 term are "standardized":
-- 'NA', 'not annotated', 'wild-type', 'confidential_restricted_data'.
-- In this table, only 'wild-type' is used to replace 'NA', 'not annotated', and
-- 'confidential_restricted_data', as for conditions used in expression table.
    strain varchar(100)
    COMMENT 'Strain information. NA: not available from source information; not annotated: information not captured by Bgee; confidential_restricted_data: information cannot be disclosed publicly',

-- ** RANKS **
-- max ranks in each data type and condition, notably used to allow normalization
-- between data types and conditions. For in situ data, they are also used for computation
-- of weighted mean between data types: for these data type, because we pool together all data
-- in a same condition, instead of computing a mean between samples, and because we use "dense ranking"
-- instead of fractional ranking (so that the max rank is equal to the number of distinct ranks),
-- it is irrelevant to consider a sum of the number of distinct ranks in each sample for weighting
-- the mean.
-- Note: these values are the same for all genes in a condition-species, this is why they are stored in this table.
    bulkMaxRank decimal(9,2) unsigned,
    singleCellMaxRank decimal(9,2) unsigned,
    inSituMaxRank decimal(9,2) unsigned
) engine = innodb COMMENT 'This table includes "real" conditions used in the raw expression table, but mostly conditions resulting from the propagation of expression calls. It results from the computation of propagated calls according to different condition parameters combination (e.g., grouping all data in a same anat. entity, or all data in a same anat. entity - stage, or data in anat. entity - sex). This is why the fields anatEntityId, stageId, sex, strain, can be null in this table (but not all of them at the same time).';

CREATE TABLE globalCondRelation (
    sourceGlobalConditionId mediumint unsigned NOT NULL,
    targetGlobalConditionId mediumint unsigned NOT NULL
) engine = innodb COMMENT 'This table stores the relations between global conditions, allowing to reconstruct the global condition graph used for call propagation. A relation exists between a source global condition and a target global condition, when the target is a parent of the source in the global condition graph (e.g., when the anatEntityId of the target is a parent of the anatEntityId of the source in the anatomical entity ontology).';

create table condToSelfGlobalCond (
    conditionId mediumint unsigned not null,
    globalConditionId mediumint unsigned not null,
    -- subsetMask indicates which condition parameters were used to map an observed condition to the corresponding globalCondition. For instance, a subset mask of 3 (binary 11000) indicates that only anatEntityId and stageId were used to define the globalCondition for this mapping. A subset Mask of 7 (binary 11100) indicates that anatEntityId, stageId and cellTypeId were used, and so on. It is really useful to retrieve all observed global condition for a given subset of condition parameters and then be able to subset the condition graph only for these global conditions and their parents.
    subsetMask tinyint unsigned NOT NULL COMMENT '5-bit mask, values 1..31. bit 1: anatEntityId, bit 2: stageId, bit 3: celltypeId, bit 4: sex, bit 5: strain'
) engine = innodb
comment = 'this table allows to link cond to their self globalCondition depending on the subset of condition parameters used to define the globalCondition';

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
    pValue decimal(31, 30) unsigned default null,
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
    reasonForExclusion enum('not excluded', 'pre-filtering', 'undefined') not null default 'not excluded'
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

-- this table contains counts and abundance level for each gene at the level of an annotated
-- sample. Each pair of bgeeGeneId and rnaSeqLibraryAnnotatedSampleId is unique.
-- * for bulk RNA-Seq one result corresponds to one gene at one organ level.
-- * for BRB-Seq one result corresponds to one gene at one organ level (after demultiplexing of pooled libraries)
-- * for full length single cell RNA-Seq one result corresponds to one gene at one cell level
-- * for droplet base single cell RNA-Seq one result corresponds to one gene at one cell-type population level (combine all counts of same cell-type per library)
create table rnaSeqLibraryAnnotatedSampleGeneResult (
    rnaSeqLibraryAnnotatedSampleId mediumint unsigned not null COMMENT 'Internal ID used to define one library at one annotated condition',
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
--  abundance values inserted here are NOT TMM normalized,
--  these are the raw data before any normalization
    abundanceUnit enum ('tpm','cpm'),
    abundance decimal(16, 6) not null COMMENT 'abundance values, NOT log transformed',
--  rawRank is not "not null" because we update this information afterwards
    rawRank decimal(9, 2) unsigned,
--  for information, measure not normalized for reads or genes lengths
    readsCount decimal(16, 6) unsigned not null COMMENT 'As of Bgee 14, read counts are "estimated counts" produced using the Kallisto software. They are not normalized for read or gene lengths.',
    UMIsCount decimal(16, 6) unsigned not null COMMENT 'As of Bgee 15, UMI counts are "estimated counts" produced using the Kallisto software. They are not normalized for read or gene lengths.',
-- zScore can be negative
    zScore decimal(35, 30),
    pValue decimal(31, 30) unsigned default null COMMENT 'present calls are based on the pValue',
    expressionId int unsigned,
-- TODO: to remove as not used anymore since Bgee 15. pValues are now used to consider a call as present/absent.
    detectionFlag enum('undefined', 'absent', 'present') default 'undefined',
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
    rnaSeqData enum('no data','poor quality','high quality') default 'no data',
    reasonForExclusion enum('not excluded', 'pre-filtering', 'biotype not targeted',
    'undefined') not null default 'not excluded'
) engine = innodb;

create table rnaSeqExperiment (
-- primary exp ID, from GEO, patterns GSExxx
    rnaSeqExperimentId varchar(70) not null,
    rnaSeqExperimentName varchar(255) not null default '',
    rnaSeqExperimentDescription text,
    dataSourceId smallInt unsigned not null,
    numberOfAnnotatedCells int unsigned not null default 0,
    DOI varchar(255)
) engine = innodb;

create table rnaSeqExperimentToKeyword (
    rnaSeqExperimentId varchar(70) not null,
    keywordId int unsigned not null
) engine = innodb;

-- Corresponds to one library in the sense of one sequencing library. It can contains several
-- sample libraries each one of them potentially having different condition in case of library (e.g BRB-Seq)
-- or sample (e.g 10x) multiplexing
-- uses to produce several runs
create table rnaSeqLibrary (
-- primary ID, from GEO, pattern GSMxxx
    rnaSeqLibraryId varchar(70) not null,
    rnaSeqExperimentId varchar(70) not null,
    rnaSeqSequencerName varchar(255) not null,
    rnaSeqTechnologyName varchar(255) not null,
    rnaSeqTechnologyIsSingleCell tinyint(1) not null,
    sampleMultiplexing boolean not null default 0,
    libraryMultiplexing boolean not null default 0,
--  **** Columns related to the sampling protocol ***
--  TODO: check validity of enum
    strandSelection enum ('NA', 'forward', 'revert', 'unstranded'),
    cellCompartment enum('NA', 'nucleus', 'cell'),
    sequencedTranscriptPart enum ('NA', '3prime', '5prime', 'full length'),
    fragmentation smallint unsigned not null default 0 COMMENT 'corresponds to fragmentation of cDNA. 0 for long reads',
    rnaSeqPopulationCaptureId varchar(255) not null,
    genotype varchar(70),
    -- In case of single read, it's the total number of reads
    allReadsCount bigint unsigned not null default 0,
-- total number of reads in library that were mapped to anything.
-- if it is not a paired-end library, this number is equal to leftMappedReadsCount
    mappedReadsCount bigint unsigned not null default 0,
-- a library is an assembly of different runs, and the runs can have different read lengths,
-- so we store the min and max read lengths
    minReadLength int unsigned not null default 0,
    maxReadLength int unsigned not null default 0,
    -- Is the library built using paired end?
-- NA: info not used for pseudo-mapping. Default value in an enum is the first one.
    libraryType enum('NA', 'single', 'paired') not null,
    usedInPropagatedCalls tinyint(1) not null default 0
) engine = innodb;

-- XXX should we keep discarded info at rnaSeqLibrary level, at rnaSeqLibraryAnnotatedSample level,
-- or at both levels? IfrnaSeqLibraryAnnotatedSample level or both do we want to provide condition ?
--
-- We sometimes discard some runs associated to a library, because of low mappability.
-- We keep track of these discarded runs in this table.
-- UPDATE Bgee 14: for pseudo-mapping using Kallisto, runs are pooled, so we can only exclude libraries,
-- not specific runs.
create table rnaSeqLibraryDiscarded (
    rnaSeqLibraryId varchar(70) not null,
    rnaSeqLibraryDiscardReason varchar(255) not null default ''
) engine = innodb;

-- Store the information of runs used, pool together to generate the results
-- for a given library.
create table rnaSeqRun (
-- same ID in GEO and SRA, pattern SRR...
    rnaSeqRunId varchar(70) not null,
    rnaSeqLibraryId varchar(70) not null
) engine = innodb;

-- corresponds to one library as it was annotated
-- * for bulk RNA-Seq one library corresponds to one sample. For such data there is a
--   1-to-1 relation between rnaSeqLibrary and rnaSeqLibraryAnnotatedSample.
-- * for BRB-Seq one library contains several libraries pooled together. Each pooled library has its
--   own annotation. For such data there will be 1-to-many relation between rnaSeqLibrary and
--   rnaSeqLibraryAnnotatedSample.
-- * for full length single cell RNA-Seq one library corresponds to one sample. For such data
--   there is a 1-to-1 relation between rnaSeqLibrary and rnaSeqLibraryAnnotatedSample.
-- * for droplet base single cell RNA-Seq one library corresponds to a cell population. Each cell has
--   been annotated with a different barcode and each barcode has its own annotation. For such data
--   there will be 1-to-many relation between rnaSeqLibrary and rnaSeqLibraryAnnotatedSample.
create table rnaSeqLibraryAnnotatedSample (
    rnaSeqLibraryAnnotatedSampleId mediumint unsigned not null,
    rnaSeqLibraryId varchar(70) not null,
    conditionId mediumint unsigned not null,
--  all *AuthorAnnotation columns correspond to free text retrieved from paper by Bgee curators.
--  anatEntityAuthorAnnotation and stageAuthorAnnotation are at the library level they are unique
--  for a given combination of rnaSeqLibraryId and conditionId. However, it is possible to have different
--  cellTypeAuthorAnnotation for a given combination of rnaSeqLibraryId and conditionId as different
--  cellTypeAuthorAnnotation can be annotated with the same cellTypeId, especially when the cell ontology does
--  not contain terms precise enough.
    cellTypeAuthorAnnotation varchar(255) not null,
    anatEntityAuthorAnnotation varchar(255) not null,
    stageAuthorAnnotation varchar(255) not null,
    abundanceUnit enum('tpm', 'cpm'),
    meanAbundanceReferenceIntergenicDistribution decimal(16, 6) not null default -1 COMMENT 'mean TPM of the distribution of the reference intergenics regions in this library, NOT log transformed',
    sdAbundanceReferenceIntergenicDistribution decimal(16, 6) not null default -1 COMMENT 'standard deviation in TPM of the distribution of the reference intergenics regions in this library, NOT log transformed',
-- TMM normalization factor
    tmmFactor decimal(8, 6) not null default 1.0,
--  abundance threshold to consider a gene as expressed
    abundanceThreshold decimal(16, 6) not null default -1,
    allGenesPercentPresent decimal(5, 2) unsigned not null default 0,
    proteinCodingGenesPercentPresent decimal(5, 2) unsigned not null default 0,
    intergenicRegionsPercentPresent decimal(5, 2) unsigned not null default 0,
    pValueThreshold decimal(5, 4) unsigned not null default 0 COMMENT 'pValue threshold used to consider genes present/absent. (for Bgee15 this threshold should always be 0.05)',
-- total number of reads in library, including those not mapped.
-- In case of paired-end libraries, it's the number of pairs of reads;
-- total number of UMIs in library, including those not mapped.
    allUMIsCount int unsigned not null default 0,
-- total number of UMIs in library that were mapped to anything.
    mappedUMIsCount int unsigned not null default 0,
-- the following fields are used for rank computations, and are set after all expression data insertion,
-- this is why null value is permitted.
    rnaSeqLibraryAnnotatedSampleMaxRank decimal(9,2) unsigned COMMENT 'The max fractional rank in this sample (see `rank` field in rnaSeqLibraryAnnotatedSampleGeneResult table)',
    rnaSeqLibraryAnnotatedSampleDistinctRankCount mediumint unsigned COMMENT 'The count of distinct rank in this sample (see `rank` field in rnaSeqLibraryAnnotatedSampleGeneResult table, used for weighted mean rank computations)',
    multipleLibraryIndividualSample boolean not null default 0 COMMENT 'boolean true if the annotated sample contains several individual samples. e.g true for 10x as one annotated sample corresponds to one cell population and individual sample will correspond to each cell of this cell population',
    --  can be null as it is applicable only to pooled bulk samples like BRB-Seq
    barcode varchar(70) COMMENT 'barcode used to pool several samples in the same library',
-- these 3 columns have been added to be able to insert precise Salmon condition information
    time decimal(5, 2) unsigned default null,
    timeUnit varchar(35) default null,
    physiologicalStatus varchar(255) default null
) engine = innodb;

--  TO CLARIFY:
--  * comment from Fred :comes from sample and library demultiplexing. In scRNA-Seq, 1 sample = 1 cell. In bulk, 1 sample = 1 organ for instance)
--  * my feeling : comes only from sample demultiplexing with barcodes describing each cell. For library demultiplexing like BRB-Seq all librariesq
--    are already described in the table `rnaSeqLibraryAnnotatedSample` So for me for BRB-Seq
--    rnaSeqLibraryAnnotatedSample.multipleLibraryIndividualSample == 0.
create table rnaSeqLibraryIndividualSample (
    rnaSeqLibraryIndividualSampleId int unsigned not null,
    rnaSeqLibraryAnnotatedSampleId mediumint unsigned not null,
    barcode varchar(70) COMMENT 'barcode used to pool several samples in the same library',
    sampleName varchar(70),
    -- total number of UMIs in library that were mapped to this individual sample
    mappedUMIsCount int unsigned not null default 0
) engine = innodb;

-- gene result at individual sample level (e.g for each cell for 10x)
create table rnaSeqLibraryIndividualSampleGeneResult (
    rnaSeqLibraryIndividualSampleId int unsigned not null COMMENT 'Internal ID used to define one individual sample',
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID',
    abundanceUnit enum ('tpm','cpm'),
    abundance decimal(16, 6) not null COMMENT 'abundance values, NOT log transformed',
    readsCount decimal(16, 6) unsigned not null COMMENT 'As of Bgee 14, read counts are "estimated counts" produced using the Kallisto software. They are not normalized for read or gene lengths.',
    UMIsCount decimal(16, 6) unsigned not null ,
    rnaSeqData enum('no data','poor quality','high quality') default 'no data',
    reasonForExclusion enum('not excluded', 'pre-filtering', 'biotype not targeted',
    'undefined') not null default 'not excluded'
) engine = innodb;

-- called protocol until Bgee 15 but updated the name as protocol now regroup
-- a lot of different parameters (e.g population captured, strand, fragmentation size, ...)
create table rnaSeqPopulationCapture (
    rnaSeqPopulationCaptureId varchar(255) not null
) engine = innodb;

-- this table stores the biotypes for which calls will be generated in a given population captured (e.g. polyA, lncRNA, etc). For instance, for polyA, we might want to generate calls only for biotype having a polyA tail. This table allows to manage this information.
create table rnaSeqPopulationCaptureToBiotype (
    rnaSeqPopulationCaptureId varchar(255) not null COMMENT 'protocol ID for which a biotype will not be used to generate absent calls',
    geneBioTypeId smallint unsigned not null COMMENT 'biotype ID for which absent calls will not be generated.'
) engine = innodb;

-- this table stores the max rank for each population captured , presence or not of multiplexing, is single cell, and species. This max rank is used for normalization.
create table rnaSeqPopulationCaptureSpeciesMaxRank (
    rnaSeqPopulationCaptureId varchar(255) not null,
    rnaSeqTechnologyIsSingleCell tinyint unsigned not null,
    sampleMultiplexing tinyint unsigned not null,
    speciesId mediumint unsigned not null,
    maxRank decimal(9,2) unsigned not null COMMENT 'The max fractional rank in this protocol and species (see `rank` field in rnaSeqLibraryAnnotatedSampleGeneResult table)'
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
    conditionId mediumint unsigned not null COMMENT 'ID of condition in the related condition table ("cond"), not stable between releases.',
    bulkScore decimal(9,2) unsigned not null COMMENT 'The score of the expression call in the bulk RNA-Seq protocol',
    bulkPValue decimal(31,30) unsigned not null COMMENT 'The p-value of the expression call in the bulk RNA-Seq protocol',
    bulkWeight bigint unsigned not null COMMENT 'The weight of the expression call in the bulk RNA-Seq protocol',
    bulkNumberObs smallint unsigned not null COMMENT 'The number of observations for the expression call in the bulk RNA-Seq protocol',
    fullLengthScore decimal(9,2) unsigned not null COMMENT 'The score of the expression call in the full-length single-cell RNA-Seq protocol',
    fullLengthPValue decimal(31,30) unsigned not null COMMENT 'The p-value of the expression call in the full-length single-cell RNA-Seq protocol',
    fullLengthWeight bigint unsigned not null COMMENT 'The weight of the expression call in the full-length single-cell RNA-Seq protocol',
    fullLengthNumberObs smallint unsigned not null COMMENT 'The number of observations for the expression call in the full-length single-cell RNA-Seq protocol',
    dropletScore decimal(9,2) unsigned not null COMMENT 'The score of the expression call in the droplet-based single-cell RNA-Seq protocol',
    dropletPValue decimal(31,30) unsigned not null COMMENT 'The p-value of the expression call in the droplet-based single-cell RNA-Seq protocol',
    dropletWeight bigint unsigned not null COMMENT 'The weight of the expression call in the droplet-based single-cell RNA-Seq protocol',
    dropletNumberObs smallint unsigned not null COMMENT 'The number of observations for the expression call in the droplet-based single-cell RNA-Seq protocol',
    inSituScore decimal(9,2) unsigned not null COMMENT 'The score of the expression call in the in situ hybridization protocol',
    inSituPValue decimal(31,30) unsigned not null COMMENT 'The p-value of the expression call in the in situ hybridization protocol',
    inSituWeight bigint unsigned not null COMMENT 'The weight of the expression call in the in situ hybridization protocol',
    inSituNumberObs smallint unsigned not null COMMENT 'The number of observations for the expression call in the in situ hybridization protocol'
) engine = innodb
comment = 'This table is a summary of expression calls for a given gene-condition (anatomical entity - developmental stage - sex- strain), over all the experiments and data types, with no propagation nor experiment expression summary.';


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
  downloadFileCategory enum('expr_simple', 'expr_complete', 'diff_expr_anatomy_complete', 'diff_expr_anatomy_simple'
   , 'diff_expr_dev_complete', 'diff_expr_dev_simple', 'ortholog',
   'affy_annot','rnaseq_annot','affy_data','rnaseq_data', 'full_length_annot', 'full_length_data', 'droplet_based_annot',
   'droplet_based_data', 'droplet_based_h5ad', 'full_length_h5ad'),
  speciesDataGroupId mediumint unsigned not null,
  downloadFileSize int unsigned not null,
  downloadFileConditionParameters set('anatomicalEntity', 'developmentalStage', 'sex', 'strain')
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

