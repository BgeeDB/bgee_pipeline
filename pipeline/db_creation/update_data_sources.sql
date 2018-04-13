-- see insert_data_sources.sql to see list of data source IDs

-- NCBI Taxonomy
update dataSource set releaseDate = '2016-07-20', releaseVersion = '', displayOrder = 1  where dataSourceId = 1;
-- Ensembl
update dataSource set releaseDate = '2016-03-18', releaseVersion = 'Release 84', displayOrder = 1  where dataSourceId = 2;
-- EnsemblMetazoa
update dataSource set releaseDate = '2016-03-18', releaseVersion = 'Release 30', displayOrder = 1  where dataSourceId = 24;
-- Uberon
update dataSource set releaseDate = '2016-07-14', releaseVersion = '', displayOrder = 1  where dataSourceId = 26;
-- SRA
update dataSource set releaseDate = null,         releaseVersion = '',                                             displayOrder = 2  where dataSourceId = 19;
-- EMBL
update dataSource set releaseDate = null,         releaseVersion = ''                                                                where dataSourceId = 3;
-- Uniprot/trEMBL
update dataSource set releaseDate = null,         releaseVersion = ''                                                                where dataSourceId = 4;
-- Uniprot/SwissProt
update dataSource set releaseDate = null,         releaseVersion = ''                                                                where dataSourceId = 5;
-- mirbase
update dataSource set releaseDate = '2012-08',    releaseVersion = 'Release 19',                                   displayOrder = 9  where dataSourceId = 6;
-- 4DXpress
update dataSource set releaseDate = null,         releaseVersion = ''                                                                where dataSourceId = 7;
-- ZFIN
update dataSource set releaseDate = '2012-11-07', releaseVersion = '',                                             displayOrder = 5  where dataSourceId = 8;
-- MGI
update dataSource set releaseDate = '2012-11-14', releaseVersion = 'MGI 5.10.03',                                  displayOrder = 6  where dataSourceId = 9;
-- FlyBase
update dataSource set releaseDate = null,         releaseVersion = ''                                                                where dataSourceId = 10;
-- ArrayExpress
update dataSource set releaseDate = null,         releaseVersion = '',                                             displayOrder = 4  where dataSourceId = 11;
-- Unigene
update dataSource set releaseDate = null,         releaseVersion = '',                                             displayOrder = 9  where dataSourceId = 12;
-- smiRNAdb
update dataSource set releaseDate = '2009-02',    releaseVersion = 'v2',                                           displayOrder = 11 where dataSourceId = 13;
-- BDGP
update dataSource set releaseDate = '2012-11-02', releaseVersion = 'Partial integration of "Release 3 + updates"', displayOrder = 8  where dataSourceId = 14;
-- Xenbase
update dataSource set releaseDate = '2012-11-17', releaseVersion = '2.8.1',                                        displayOrder = 7  where dataSourceId = 15;
-- neXtProt
update dataSource set releaseDate = null,         releaseVersion = ''                                                                where dataSourceId = 16;
-- GEO
update dataSource set releaseDate = null,         releaseVersion = '',                                             displayOrder = 3  where dataSourceId = 17;
-- GO
update dataSource set releaseDate = '2012-11-06', releaseVersion = 'data-version: 2012-11-07'                                        where dataSourceId = 18;
-- SRA
update dataSource set releaseDate = '2016-03-18', releaseVersion = '',                                             displayOrder = 2  where dataSourceId = 19;
-- HGNC
update dataSource set releaseDate = '2012-11-07', releaseVersion = '',                                             displayOrder = 1  where dataSourceId = 20;
-- CCDS
update dataSource set releaseDate = '2013-11-29', releaseVersion = '15',                                           displayOrder = 1  where dataSourceId = 21;
-- RGD
update dataSource set releaseDate = '2013-12-26', releaseVersion = '5.0',                                          displayOrder = 1  where dataSourceId = 22;
-- WormBase
update dataSource set releaseDate = '2013-10-11', releaseVersion = 'WS240',                                        displayOrder = 1  where dataSourceId = 23;
-- GTEx-dbGAP
update dataSource set releaseDate = '2015-10-03', releaseVersion = 'v6.p1',                                        displayOrder = 3  where dataSourceId = 31;
