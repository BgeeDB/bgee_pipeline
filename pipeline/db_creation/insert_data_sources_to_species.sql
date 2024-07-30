-- see insert_data_sources.sql to see list of data source IDs
-- species IDs are NCBI tax IDs

-- *** query to retrieve species and data sources for RNA-Seq data ***
-- select t2.dataSourceId,count(DISTINCT t1.rnaSeqExperimentId),(SELECT speciesId from gene where bgeeGeneId = (select bgeeGeneId from rnaSeqResult WHERE rnaSeqResult.rnaSeqLibraryId = t1.rnaSeqLibraryId LIMIT 1)) as speciesId FROM rnaSeqLibrary AS t1 INNER JOIN rnaSeqExperiment AS t2 ON t1.rnaSeqExperimentId = t2.rnaSeqExperimentId GROUP BY speciesId, dataSourceId;
-- *** query to retrieve species and data sources for Affymetrix data ***
-- select t2.dataSourceId, count(DISTINCT t1.microarrayExperimentId), (SELECT speciesId from gene where bgeeGeneId = (select bgeeGeneId from affymetrixProbeset WHERE affymetrixProbeset.bgeeAffymetrixChipId = t1.bgeeAffymetrixChipId LIMIT 1)) as speciesId FROM affymetrixChip AS t1 INNER JOIN microarrayExperiment AS t2 ON t1.microarrayExperimentId = t2.microarrayExperimentId GROUP BY speciesId, dataSourceId;
-- *** query to retrieve species and data sources for EST data ***
-- select t1.dataSourceId, count(DISTINCT estLibraryId), (SELECT speciesId from gene where bgeeGeneId = (select bgeeGeneId from expressedSequenceTag WHERE expressedSequenceTag.estLibraryId = t1.estLibraryId LIMIT 1)) as speciesId from estLibrary as t1 GROUP BY speciesId, dataSourceId;
-- ** in situ **
-- select t2.dataSourceId, count(DISTINCT t1.inSituExperimentId), (SELECT speciesId from gene where bgeeGeneId = (select bgeeGeneId from inSituSpot WHERE inSituSpot.inSituEvidenceId = t1.inSituEvidenceId LIMIT 1)) as speciesId FROM inSituEvidence AS t1 INNER JOIN inSituExperiment AS t2 ON t1.inSituExperimentId = t2.inSituExperimentId GROUP BY speciesId, dataSourceId;

DELETE FROM dataSourceToSpecies;

--  IN SITU ----------------------------------

-- ZFIN in situ in zebrafish, data and annot
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(8, 7955, 'in situ', 'data'),
(8, 7955, 'in situ', 'annotation');
-- MGI in situ mouse, data and annot
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(9, 10090, 'in situ', 'data'),
(9, 10090, 'in situ', 'annotation');
-- Flybase in situ droso, data and annot
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(10, 7227, 'in situ', 'data'),
(10, 7227, 'in situ', 'annotation');
-- BDGP source for in situ data and annot in drosophila
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(14, 7227, 'in situ', 'data'),
(14, 7227, 'in situ', 'annotation');
-- Xenbase source for in situ data and annot in xenopus
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(15, 8364, 'in situ', 'data'),
(15, 8364, 'in situ', 'annotation');


-- AFFYMETRIX --------------------------------

-- ArrayExpress source for Affymetrix data in drosophila, human, zebrafish, mouse
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(11, 7227,  'affymetrix', 'data'),
(11, 9606,  'affymetrix', 'data'),
(11, 7955,  'affymetrix', 'data'),
(11, 10090, 'affymetrix', 'data');
-- GEO, source for Affymetrix data in c. elegans, drosophila, human, zebrafish, mouse
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(17, 6239,  'affymetrix', 'data'),
(17, 7227,  'affymetrix', 'data'),
(17, 7955,  'affymetrix', 'data'),
(17, 9544,  'affymetrix', 'data'),
(17, 9606,  'affymetrix', 'data'),
(17, 10090, 'affymetrix', 'data'),
(17, 10116, 'affymetrix', 'data');


-- EST ---------------------------------------

-- Unigene source for EST data in droso, zebrafish, xenopus, human, mouse
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(12, 7227,  'est', 'data'),
(12, 7955,  'est', 'data'),
(12, 8364,  'est', 'data'),
(12, 9606,  'est', 'data'),
(12, 10090, 'est', 'data');
-- smiRNAdb source for EST data in zebrafish human mouse
-- we hide it because it's outdated
-- INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
-- (13, 7955,  'est', 'data'),
-- (13, 9606,  'est', 'data'),
-- (13, 10090, 'est', 'data');


-- RNA-SEQ -----------------------------------

-- SRA source of RNA-Seq data in lots of species :p
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(19, 6239,   'rna-seq', 'data'),
(19, 7227,   'rna-seq', 'data'),
(19, 7237,   'rna-seq', 'data'),
(19, 7240,   'rna-seq', 'data'),
(19, 7740,   'rna-seq', 'data'),
(19, 7897,   'rna-seq', 'data'),
(19, 7918,   'rna-seq', 'data'),
(19, 7936,   'rna-seq', 'data'),
(19, 7955,   'rna-seq', 'data'),
(19, 7994,   'rna-seq', 'data'),
(19, 8010,   'rna-seq', 'data'),
(19, 8030,   'rna-seq', 'data'),
(19, 8049,   'rna-seq', 'data'),
(19, 8081,   'rna-seq', 'data'),
(19, 8090,   'rna-seq', 'data'),
(19, 8154,   'rna-seq', 'data'),
(19, 8355,   'rna-seq', 'data'),
(19, 8364,   'rna-seq', 'data'),
(19, 9031,   'rna-seq', 'data'),
(19, 9103,   'rna-seq', 'data'),
(19, 9258,   'rna-seq', 'data'),
(19, 9483,   'rna-seq', 'data'),
(19, 9531,   'rna-seq', 'data'),
(19, 9541,   'rna-seq', 'data'),
(19, 9544,   'rna-seq', 'data'),
(19, 9545,   'rna-seq', 'data'),
(19, 9555,   'rna-seq', 'data'),
(19, 9593,   'rna-seq', 'data'),
(19, 9597,   'rna-seq', 'data'),
(19, 9598,   'rna-seq', 'data'),
(19, 9606,   'rna-seq', 'data'),
(19, 9615,   'rna-seq', 'data'),
(19, 9685,   'rna-seq', 'data'),
(19, 9796,   'rna-seq', 'data'),
(19, 9823,   'rna-seq', 'data'),
(19, 9913,   'rna-seq', 'data'),
(19, 9925,   'rna-seq', 'data'),
(19, 9940,   'rna-seq', 'data'),
(19, 9974,   'rna-seq', 'data'),
(19, 9986,   'rna-seq', 'data'),
(19, 10090,  'rna-seq', 'data'),
(19, 10116,  'rna-seq', 'data'),
(19, 10141,  'rna-seq', 'data'),
(19, 10181,  'rna-seq', 'data'),
(19, 13616,  'rna-seq', 'data'),
(19, 28377,  'rna-seq', 'data'),
(19, 30608,  'rna-seq', 'data'),
(19, 32507,  'rna-seq', 'data'),
(19, 52904,  'rna-seq', 'data'),
(19, 60711,  'rna-seq', 'data'),
(19, 69293,  'rna-seq', 'data'),
(19, 105023, 'rna-seq', 'data');

-- GEO source of RNA-Seq data in lots of species :p
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(100, 6239,   'rna-seq', 'data'),
(100, 7227,   'rna-seq', 'data'),
(100, 7237,   'rna-seq', 'data'),
(100, 7240,   'rna-seq', 'data'),
(100, 7740,   'rna-seq', 'data'),
(100, 7897,   'rna-seq', 'data'),
(100, 7918,   'rna-seq', 'data'),
(100, 7936,   'rna-seq', 'data'),
(100, 7955,   'rna-seq', 'data'),
(100, 7994,   'rna-seq', 'data'),
(100, 8010,   'rna-seq', 'data'),
(100, 8030,   'rna-seq', 'data'),
(100, 8049,   'rna-seq', 'data'),
(100, 8081,   'rna-seq', 'data'),
(100, 8090,   'rna-seq', 'data'),
(100, 8154,   'rna-seq', 'data'),
(100, 8355,   'rna-seq', 'data'),
(100, 8364,   'rna-seq', 'data'),
(100, 9031,   'rna-seq', 'data'),
(100, 9103,   'rna-seq', 'data'),
(100, 9258,   'rna-seq', 'data'),
(100, 9483,   'rna-seq', 'data'),
(100, 9531,   'rna-seq', 'data'),
(100, 9541,   'rna-seq', 'data'),
(100, 9544,   'rna-seq', 'data'),
(100, 9545,   'rna-seq', 'data'),
(100, 9555,   'rna-seq', 'data'),
(100, 9593,   'rna-seq', 'data'),
(100, 9597,   'rna-seq', 'data'),
(100, 9598,   'rna-seq', 'data'),
(100, 9606,   'rna-seq', 'data'),
(100, 9615,   'rna-seq', 'data'),
(100, 9685,   'rna-seq', 'data'),
(100, 9796,   'rna-seq', 'data'),
(100, 9823,   'rna-seq', 'data'),
(100, 9913,   'rna-seq', 'data'),
(100, 9925,   'rna-seq', 'data'),
(100, 9940,   'rna-seq', 'data'),
(100, 9974,   'rna-seq', 'data'),
(100, 9986,   'rna-seq', 'data'),
(100, 10090,  'rna-seq', 'data'),
(100, 10116,  'rna-seq', 'data'),
(100, 10141,  'rna-seq', 'data'),
(100, 10181,  'rna-seq', 'data'),
(100, 13616,  'rna-seq', 'data'),
(100, 28377,  'rna-seq', 'data'),
(100, 30608,  'rna-seq', 'data'),
(100, 32507,  'rna-seq', 'data'),
(100, 52904,  'rna-seq', 'data'),
(100, 60711,  'rna-seq', 'data'),
(100, 69293,  'rna-seq', 'data'),
(100, 105023, 'rna-seq', 'data');

-- wormbase source for in situ data and annot in c. elegans
-- also source of annotations for affymetrix and RNA-Seq
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(23, 6239, 'in situ',    'data'),
(23, 6239, 'in situ',    'annotation'),
(23, 6239, 'affymetrix', 'annotation'),
(23, 6239, 'rna-seq',    'annotation');


-- Bgee source of annotation for EST, Affymetrix and RNA-Seq data in all species
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(25, 6239,   'affymetrix', 'annotation'),
(25, 6239,   'est',        'annotation'),
(25, 6239,   'rna-seq',    'annotation'),
(25, 7227,   'affymetrix', 'annotation'),
(25, 7227,   'est',        'annotation'),
(25, 7227,   'rna-seq',    'annotation'),
(25, 7955,   'affymetrix', 'annotation'),
(25, 7955,   'est',        'annotation'),
(25, 7955,   'rna-seq',    'annotation'),
(25, 8364,   'affymetrix', 'annotation'),
(25, 8364,   'est',        'annotation'),
(25, 8364,   'rna-seq',    'annotation'),
(25, 9031,   'affymetrix', 'annotation'),
(25, 9031,   'est',        'annotation'),
(25, 9031,   'rna-seq',    'annotation'),
(25, 9258,   'affymetrix', 'annotation'),
(25, 9258,   'est',        'annotation'),
(25, 9258,   'rna-seq',    'annotation'),
(25, 9544,   'affymetrix', 'annotation'),
(25, 9544,   'est',        'annotation'),
(25, 9544,   'rna-seq',    'annotation'),
(25, 9593,   'affymetrix', 'annotation'),
(25, 9593,   'est',        'annotation'),
(25, 9593,   'rna-seq',    'annotation'),
(25, 9597,   'affymetrix', 'annotation'),
(25, 9597,   'est',        'annotation'),
(25, 9597,   'rna-seq',    'annotation'),
(25, 9598,   'affymetrix', 'annotation'),
(25, 9598,   'est',        'annotation'),
(25, 9598,   'rna-seq',    'annotation'),
(25, 9606,   'affymetrix', 'annotation'),
(25, 9606,   'est',        'annotation'),
(25, 9606,   'rna-seq',    'annotation'),
(25, 9823,   'affymetrix', 'annotation'),
(25, 9823,   'est',        'annotation'),
(25, 9823,   'rna-seq',    'annotation'),
(25, 9913,   'affymetrix', 'annotation'),
(25, 9913,   'est',        'annotation'),
(25, 9913,   'rna-seq',    'annotation'),
(25, 10090,  'affymetrix', 'annotation'),
(25, 10090,  'est',        'annotation'),
(25, 10090,  'rna-seq',    'annotation'),
(25, 10116,  'affymetrix', 'annotation'),
(25, 10116,  'est',        'annotation'),
(25, 10116,  'rna-seq',    'annotation'),
(25, 13616,  'affymetrix', 'annotation'),
(25, 13616,  'est',        'annotation'),
(25, 13616,  'rna-seq',    'annotation'),
(25, 28377,  'affymetrix', 'annotation'),
(25, 28377,  'est',        'annotation'),
(25, 28377,  'rna-seq',    'annotation'),
(25, 7237,   'rna-seq',    'annotation'),
(25, 7240,   'rna-seq',    'annotation'),
(25, 7740,   'rna-seq',    'annotation'),
(25, 7897,   'rna-seq',    'annotation'),
(25, 7918,   'rna-seq',    'annotation'),
(25, 7936,   'rna-seq',    'annotation'),
(25, 7994,   'rna-seq',    'annotation'),
(25, 8010,   'rna-seq',    'annotation'),
(25, 8030,   'rna-seq',    'annotation'),
(25, 8049,   'rna-seq',    'annotation'),
(25, 8081,   'rna-seq',    'annotation'),
(25, 8090,   'rna-seq',    'annotation'),
(25, 8154,   'rna-seq',    'annotation'),
(25, 8355,   'rna-seq',    'annotation'),
(25, 9103,   'rna-seq',    'annotation'),
(25, 9483,   'rna-seq',    'annotation'),
(25, 9531,   'rna-seq',    'annotation'),
(25, 9541,   'rna-seq',    'annotation'),
(25, 9545,   'rna-seq',    'annotation'),
(25, 9555,   'rna-seq',    'annotation'),
(25, 9615,   'rna-seq',    'annotation'),
(25, 9685,   'rna-seq',    'annotation'),
(25, 9796,   'rna-seq',    'annotation'),
(25, 9925,   'rna-seq',    'annotation'),
(25, 9940,   'rna-seq',    'annotation'),
(25, 9974,   'rna-seq',    'annotation'),
(25, 9986,   'rna-seq',    'annotation'),
(25, 10141,  'rna-seq',    'annotation'),
(25, 10181,  'rna-seq',    'annotation'),
(25, 30608,  'rna-seq',    'annotation'),
(25, 32507,  'rna-seq',    'annotation'),
(25, 52904,  'rna-seq',    'annotation'),
(25, 60711,  'rna-seq',    'annotation'),
(25, 69293,  'rna-seq',    'annotation'),
(25, 105023, 'rna-seq',    'annotation');

-- Bgee source of annotation for single cell data
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(25, 9606,   'single-cell RNA-Seq', 'annotation'),
(25, 10090,  'single-cell RNA-Seq', 'annotation'),
(25, 7227,  'single-cell RNA-Seq', 'annotation');

-- source of data for single cell data
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(102, 9606,   'single-cell RNA-Seq', 'data'),
(102, 10090,  'single-cell RNA-Seq', 'data'),
(102, 7227,  'single-cell RNA-Seq', 'data');

-- GTEx data source
INSERT INTO dataSourceToSpecies (dataSourceId, speciesId, dataType, infoType) VALUES
(31, 9606, 'rna-seq', 'data');

