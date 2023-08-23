ALTER DATABASE CHARACTER SET utf8 COLLATE utf8_general_ci;

create table anatEntity (
    anatEntityId varchar(20) not null COMMENT 'Anatomical entity id',
    anatEntityName varchar(255) not null COMMENT 'Anatomical entity name',
    anatEntityDescription TEXT COMMENT 'Anatomical entity description',
    PRIMARY KEY(anatEntityId)
) engine = innodb;

create table species (
    speciesId mediumint unsigned not null COMMENT 'NCBI species taxon id',
-- example: homo
    genus varchar(70) not null COMMENT 'Genus name',
-- example: sapiens
    species varchar(70) not null COMMENT 'Species name',
-- exemple: human
-- warning, this column in the bgee schema is defined as not null but contains empty values.
-- In bgeelite we added a default value corresponding to an empty string ''
    speciesCommonName varchar(70) default '' COMMENT 'NCBI species common name',
    genomeVersion varchar(50) not null,
-- ID of the species whose the genome was used for this species. This is used
-- when a genome is not in Ensembl. For instance, for bonobo (ID 9597), we use the chimp
-- genome (ID 9598), because bonobo is not in Ensembl.
-- We don't use a foreign key constraint here, because maybe the species whose the genome
-- was used does not have any data in Bgee, and thus is not in the taxon table.
-- If the correct genome of the species was used, the value of this field is 0.
    genomeSpeciesId mediumint unsigned not null default 0 COMMENT 'NCBI species taxon id used for mapping (0 if the same species)',
    PRIMARY KEY(speciesId),
	UNIQUE(species, genus)
) engine = innodb;

create table dataSource (
    dataSourceId          smallInt unsigned not null,
    dataSourceName        varchar(255)      not null            COMMENT 'Data source name',
    XRefUrl               varchar(255)      not null default '' COMMENT 'URL for cross-references to data sources',
    baseUrl               varchar(255)      not null default '' COMMENT 'URL to the home page of data sources',
    releaseVersion        varchar(255)      not null default '' COMMENT 'Version of data source used',
    dataSourceDescription TEXT                                  COMMENT 'Description of data source',
    PRIMARY KEY(dataSourceId)
) engine = innodb;

create table gene (
-- warning, maybe this bgeeGeneId will need to be changed to an 'int' when we reach around 200 species
    bgeeGeneId mediumint unsigned not null auto_increment COMMENT 'Numeric internal gene ID used for improving performances',
    geneId varchar(20) not null COMMENT 'Real gene id',
    geneName varchar(255) not null default '' COMMENT 'Gene name',
    geneDescription TEXT COMMENT 'Gene description',
    dataSourceId smallInt unsigned not null,
    speciesId mediumint unsigned not null COMMENT 'NCBI species taxon id this gene belongs to',
    PRIMARY KEY (bgeeGeneId),
    UNIQUE(geneId, speciesId),
    FOREIGN KEY(speciesId) REFERENCES species(speciesId) ON DELETE CASCADE,
    FOREIGN KEY(dataSourceId) REFERENCES dataSource(dataSourceId) ON DELETE CASCADE
) engine = innodb;

create table stage (
    stageId varchar(20) not null COMMENT 'Developmental stage id',
    stageName varchar(255) not null COMMENT 'Developmental stage name',
    stageDescription TEXT COMMENT 'Developmental stage description',
    PRIMARY KEY(stageId)
) engine = innodb;

create table globalCond (
    globalConditionId mediumint unsigned not null auto_increment,
    anatEntityId varchar(20)  COMMENT 'Uberon anatomical entity ID. Can be null in this table if this condition aggregates data according to other condition parameters (e.g., grouping all data in a same stage whatever the organ is).',
    stageId varchar(20)  COMMENT 'Uberon stage ID. Can be null in this table if this condition aggregates data according to other condition parameters (e.g., grouping all data in a same organ whatever the dev. stage is).',
    cellTypeId            varchar(20)  default null COMMENT 'A second uberon anatomical entity ID used to manage composition of anatomical entities. Used only for single cell data for postcomposition of anatomical entity ID and cell type ID',
    sex enum('any', 'hermaphrodite', 'female', 'male'),
	strain varchar(100),
	speciesId mediumint unsigned not null COMMENT 'NCBI species taxon ID',
    PRIMARY KEY(globalConditionId),
-- not a primary key as for table cond, because some field can be null
	UNIQUE(anatEntityId, cellTypeId, stageId, speciesId, sex, strain),
	FOREIGN KEY(anatEntityId) REFERENCES anatEntity(anatEntityId) ON DELETE CASCADE,
	FOREIGN KEY (cellTypeId) REFERENCES anatEntity(anatEntityId) ON DELETE CASCADE,
	FOREIGN KEY(stageId) REFERENCES stage(stageId) ON DELETE CASCADE,
	FOREIGN KEY(speciesId) REFERENCES species(speciesId) ON DELETE CASCADE
) engine = innodb;
-- COMMENT 'This table contains all condition used in the globalExpression table. It thus includes "real" conditions, but mostly conditions resulting from the propagation of expression calls in the globalExpression table. It results from the computation of propagated calls according to different condition parameters combination (e.g., grouping all data in a same anat. entity, or all data in a same anat. entity - stage). This is why the fields anatEntityId or stageId can be null in this table (but not all of them at the same time).';


create table globalExpression (
    globalExpressionId int unsigned not null auto_increment COMMENT 'Internal expression ID, not stable between releases.',
    bgeeGeneId mediumint unsigned not null COMMENT 'Internal gene ID, not stable between releases.',
    globalConditionId mediumint unsigned not null COMMENT 'ID of condition in the related condition table ("globalCond"), not stable between releases.',
    summaryQuality varchar(10) not null,
    globalRank decimal(9, 2) unsigned not null COMMENT 'Normalized rank for this gene-condition after normalization over all data types, conditions and species',
    score decimal(9, 5) unsigned not null COMMENT 'Use the minimum and maximum rank of the species to normalize the expression to a value between 0 and 100',
    pValue decimal(31, 30) unsigned default null,
	propagationOrigin varchar(20) not null COMMENT 'The origin of the propagated expression calls : self, self and descendant, self and ancestor, or all',
    callType varchar(20) not null COMMENT 'Type of the call. Can be EXPRESSED or NOT_EXPRESSED',
	PRIMARY KEY(bgeeGeneId, globalConditionId),
    UNIQUE(globalExpressionId),
    FOREIGN KEY(bgeeGeneId) REFERENCES gene(bgeeGeneId) ON DELETE CASCADE,
    FOREIGN KEY(globalConditionId) REFERENCES globalCond(globalConditionId) ON DELETE CASCADE
) engine = innodb;
-- COMMENT = 'This table is a summary of expression calls for a given gene-condition (anatomical entity - developmental stage), over all the experiments and data types, with all data propagated and reconciled.';
