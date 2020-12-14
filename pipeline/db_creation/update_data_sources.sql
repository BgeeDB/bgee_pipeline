-- see insert_data_sources.sql to see list of data source IDs

-- NCBI Taxonomy
update dataSource set releaseDate = '2020-12-14', releaseVersion = '',                             displayOrder = 1  where dataSourceId = 1;
-- Ensembl
update dataSource set releaseDate = '2020-12-01', releaseVersion = '102',                          displayOrder = 1  where dataSourceId = 2;
-- EMBL (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = ''                                                where dataSourceId = 3;
-- Uniprot/trEMBL (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = ''                                                where dataSourceId = 4;
-- Uniprot/SwissProt (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = ''                                                where dataSourceId = 5;
-- mirbase
update dataSource set releaseDate = '2014-06-23', releaseVersion = '21',                           displayOrder = 9  where dataSourceId = 6;
-- 4DXpress
update dataSource set releaseDate = null,         releaseVersion = ''                                                where dataSourceId = 7;
-- ZFIN
update dataSource set releaseDate = '2017-03-27', releaseVersion = '',                             displayOrder = 5  where dataSourceId = 8;
-- MGI
update dataSource set releaseDate = '2016-12-05', releaseVersion = '',                             displayOrder = 6  where dataSourceId = 9;
-- FlyBase
update dataSource set releaseDate = '2016-12-06', releaseVersion = ''                                                where dataSourceId = 10;
-- ArrayExpress
update dataSource set releaseDate = '2016-08-22', releaseVersion = '',                             displayOrder = 4  where dataSourceId = 11;
-- Unigene
update dataSource set releaseDate = '2016-12-05', releaseVersion = '',                             displayOrder = 9  where dataSourceId = 12;
-- smiRNAdb
update dataSource set releaseDate = '2009-04-23', releaseVersion = '2',                            displayOrder = 11 where dataSourceId = 13;
-- BDGP
update dataSource set releaseDate = '2016-12-04', releaseVersion = '',                             displayOrder = 8  where dataSourceId = 14;
-- Xenbase
update dataSource set releaseDate = '2016-11-15', releaseVersion = '',                             displayOrder = 7  where dataSourceId = 15;
-- neXtProt
update dataSource set releaseDate = null,         releaseVersion = ''                                                where dataSourceId = 16;
-- GEO (Affymetrix)
update dataSource set releaseDate = '2016-08-22', releaseVersion = '',                             displayOrder = 3  where dataSourceId = 17;
-- GO
update dataSource set releaseDate = '2018-09-19', releaseVersion = ''                                                where dataSourceId = 18;
-- SRA
update dataSource set releaseDate = '2020-08-07', releaseVersion = '',                             displayOrder = 2  where dataSourceId = 19;
-- HGNC (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = '',                             displayOrder = 1  where dataSourceId = 20;
-- CCDS (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = '',                             displayOrder = 1  where dataSourceId = 21;
-- RGD (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = '',                             displayOrder = 1  where dataSourceId = 22;
-- WormBase
update dataSource set releaseDate = '2016-08-29', releaseVersion = '',                             displayOrder = 1  where dataSourceId = 23;
-- EnsemblMetazoa
update dataSource set releaseDate = '2020-12-01', releaseVersion = '49',                           displayOrder = 1  where dataSourceId = 24;
-- Bgee
update dataSource set releaseDate = '2020-03-26', releaseVersion = '15.0',                         displayOrder = 1  where dataSourceId = 25;
-- Uberon
update dataSource set releaseDate = '2016-07-14', releaseVersion = '',                             displayOrder = 1  where dataSourceId = 26;
-- Developmental stage ontologies
update dataSource set releaseDate = '2016-12-08', releaseVersion = '',                             displayOrder = 1  where dataSourceId = 27;
-- OMA
update dataSource set releaseDate = '2017-08-21', releaseVersion = '',                             displayOrder = 1  where dataSourceId = 28;
-- Anatomical similarity annotations
update dataSource set releaseDate = '2016-07-14', releaseVersion = '',                             displayOrder = 1  where dataSourceId = 29;
-- CIO
update dataSource set releaseDate = '2015-03-10', releaseVersion = '',                             displayOrder = 1  where dataSourceId = 30;
-- GTEx-dbGAP
update dataSource set releaseDate = '2015-10-03', releaseVersion = '6.p1',                         displayOrder = 3  where dataSourceId = 31;
-- RefSeq
update dataSource set releaseDate = '2020-12-01', releaseVersion = '',                             displayOrder = 1  where dataSourceId = 37;
-- GenBank
update dataSource set releaseDate = '2020-12-01', releaseVersion = '',                             displayOrder = 1  where dataSourceId = 38;
-- GEO (RNASeq)
update dataSource set releaseDate = '2018-08-07', releaseVersion = '',                             displayOrder = 2  where dataSourceId = 100;

