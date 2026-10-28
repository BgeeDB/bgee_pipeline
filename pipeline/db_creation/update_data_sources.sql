-- see insert_data_sources.sql to see list of data source IDs

-- NCBI Taxonomy
update dataSource set releaseDate = '2020-12-14', releaseVersion = '',                             displayOrder = 1   where dataSourceId = 1;
-- Ensembl
update dataSource set releaseDate = '2024-10-17', releaseVersion = '113',                          displayOrder = 1   where dataSourceId = 2;
-- EMBL (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = ''                                                 where dataSourceId = 3;
-- Uniprot/trEMBL (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = ''                                                 where dataSourceId = 4;
-- Uniprot/SwissProt (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = ''                                                 where dataSourceId = 5;
-- mirbase (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = '',                             displayOrder = 9   where dataSourceId = 6;
-- ZFIN
update dataSource set releaseDate = '2025-08-05', releaseVersion = '',                             displayOrder = 5   where dataSourceId = 8;
-- MGI
update dataSource set releaseDate = '2025-07-21', releaseVersion = '',                             displayOrder = 6   where dataSourceId = 9;
-- FlyBase
update dataSource set releaseDate = '2025-07-23', releaseVersion = ''                                                 where dataSourceId = 10;
-- ArrayExpress RNA-Seq
update dataSource set releaseDate = '2023-04-04', releaseVersion = '',                             displayOrder = 255 where dataSourceId = 101;
-- BDGP
update dataSource set releaseDate = '2021-08-31', releaseVersion = '',                             displayOrder = 8   where dataSourceId = 14;
-- Xenbase
update dataSource set releaseDate = '2025-09-03', releaseVersion = '',                             displayOrder = 7   where dataSourceId = 15;
-- GEO (RNASeq)
update dataSource set releaseDate = '2023-04-04', releaseVersion = '',                             displayOrder = 2   where dataSourceId = 100;
-- GO
update dataSource set releaseDate = '2021-02-08', releaseVersion = ''                                                 where dataSourceId = 18;
-- SRA (RNA-Seq)
update dataSource set releaseDate = '2023-04-04', releaseVersion = '',                             displayOrder = 2   where dataSourceId = 19;
-- SRA (scRNA-Seq)
update dataSource set releaseDate = '2023-04-04', releaseVersion = '',                             displayOrder = 1   where dataSourceId = 102;
-- HGNC (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = '',                             displayOrder = 1   where dataSourceId = 20;
-- CCDS (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = '',                             displayOrder = 1   where dataSourceId = 21;
-- RGD (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = '',                             displayOrder = 1   where dataSourceId = 22;
-- WormBase
update dataSource set releaseDate = '2025-03-31', releaseVersion = '',                             displayOrder = 1   where dataSourceId = 23;
-- EnsemblMetazoa
update dataSource set releaseDate = '2024-10-17', releaseVersion = '60',                           displayOrder = 1   where dataSourceId = 24;
-- Bgee
update dataSource set releaseDate = '2026-11-21', releaseVersion = '16.0',                         displayOrder = 1   where dataSourceId = 25;
-- Uberon
update dataSource set releaseDate = '2020-09-16', releaseVersion = '',                             displayOrder = 1   where dataSourceId = 26;
-- Developmental stage ontologies
update dataSource set releaseDate = '2021-03-01', releaseVersion = '',                             displayOrder = 1   where dataSourceId = 27;
-- OMA
update dataSource set releaseDate = '2026-05-01', releaseVersion = '',                             displayOrder = 1   where dataSourceId = 28;
-- Anatomical similarity annotations
update dataSource set releaseDate = '2016-07-14', releaseVersion = '',                             displayOrder = 1   where dataSourceId = 29;
-- CIO
update dataSource set releaseDate = '2015-03-10', releaseVersion = '',                             displayOrder = 1   where dataSourceId = 30;
-- GTEx-dbGAP
update dataSource set releaseDate = '2015-10-03', releaseVersion = '6.p1',                         displayOrder = 3   where dataSourceId = 31;
-- RefSeq genomic
update dataSource set releaseDate = '2025-06-19', releaseVersion = '',                             displayOrder = 1   where dataSourceId = 37;
-- GenBank (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = '',                             displayOrder = 1   where dataSourceId = 38;
-- RefSeq nucleotide (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = '',                             displayOrder = 1   where dataSourceId = 41;
-- RefSeq protein (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = '',                             displayOrder = 1   where dataSourceId = 42;
-- NCBI Gene ID (taken from Ensembl xrefs)
update dataSource set releaseDate = null,         releaseVersion = '',                             displayOrder = 1   where dataSourceId = 43;
-- Alliance of Genome Resources
update dataSource set releaseDate = '2025-08-28', releaseVersion = '8.2.0',                        displayOrder = 1   where dataSourceId = 44;

