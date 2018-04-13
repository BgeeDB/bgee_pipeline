# Bgee pipeline

**General information:**

1. [Introduction](#introduction)
2. [List of pipeline steps](#list-of-pipeline-steps)
   * [Pipeline and database initialization](#pipeline-and-database-initialization)
   * [Taxa, genomes, ontologies (e.g., Uberon)](#taxa-genomes-ontologies)
   * [Raw expression data analyses and insertion](raw-data analyses-and-insertion): RNA-Seq, Affymetrix, *in situ* hybridization, EST analyses.
   * [Bgee post-processing steps](#bgee-post-processing-steps)

Shortcut note: for the RNA-Seq analysis pipeline, see [RNA_Seq/](RNA_Seq/).


**Developer guidelines**

1. [Developer guidelines](#developer-guidelines)

# General information

## Introduction

Through all the documentation, `RELEASE` will denote the current Bgee version
(e.g., if the current release number is `14`, `bgee_vRELEASE` means `bgee_v14`).

Each step in the Bgee pipeline is represented by a specific folder,
containing a Makefile, and related scripts. Variables common
to several steps are defined in the file `Makefile.common`.
Sensitive variables are stored in the file `Makefile.Config`.

Each Makefile ultimately generates an output file, called step_verification_RELEASE.txt, in the corresponding output folder.
This file is generated for the Makefile to determine whether a step should be re-run, and for developers to control that the step was correctly executed.
These files are committed to git, so that results can be compared between releases.
They are not meant to be the output of the Makefiles, but, rather, small files to be added to git, and to served as control of the procedures.


## List of pipeline steps

### Pipeline and database initialization

1. Pipeline initialization: see [init/](init/).
2. Database creation: see [db_creation/](db_creation/).

### Taxa, genomes, ontologies

1. Species and taxon information: see [species/](species/).
2. Genomes and gene-related information: see [genes/](genes/).
3. Anatomical ontology (Uberon) and developmental stage ontologies: see [uberon/](uberon/)

### Raw data analyses and insertion

1. RNA-Seq data analyses: see [RNA_Seq/](RNA_Seq/).
2. Affymetrix data analyses: see [Affymetrix](Affymetrix/).
3. *In situ* hybridization data analyses: see [In_situ/](In_situ/).
4. EST data analyses: see [ESTs/](ESTs/).
5. Differential expression analyses: see [Differential_expression/](Differential_expression/)

### Bgee post-processing steps

1. Annotation sanity checks: see [post_processing/](post_processing/).
2. Propagation/reconciliation of present/absent expression calls: see [post_processing/](post_processing/).
3. Computations of expression rank scores: see [post_processing/](post_processing/).
4. Generation of files containing data available for download: see [download_files/](download_files/).
5. Generation of XRefs to Uniprot: see [download_files/](download_files/).
6. insertion of information about versions of the data sources used: see [db_creation Makefile](db_creation/Makefile), target `update_data_sources.sql`.

# Developer guidelines

1. [Keeping track of data source versions](#keeping-track-of-data-source-versions)
2. [Pipeline configuration](#configuration)
3. [Running and re-running pipeline steps](#running-and-re-running-pipeline-steps)

## Keeping track of data source versions

* At each step of the pipeline, you will need to update the file `db_creation/update_data_sources.sql`, that keeps track of the version of the data sources used for the current release.
  This file will be used at the end of the pipeline run, to insert this information into the database.
  The reason why this information is not managed by the Makefiles, is that the ways to obtain this information are too disparate between data sources
  (sometimes you have to look at the home page of the website, sometimes to look at a specific file, sometimes you cannot use the modification date of the file, but need to look for a release date inside the file, etc).

## Configuration

Before running the pipeline on a specific machine, you need to perform some configurations:

* in `Makefile.Config`: edit this file with correct values of logins and passwords. The correct values
should not be versioned! (easier than to encrypt the file)

* in `Makefile.common`, edit the following variables as needed:
  * `RELEASE`: version of Bgee for which the pipeline is being run
  * `ENSRELEASE`: version used of Ensembl
  * `TMP DIR`: where to store (potentially large) TMP files
  * Servers and ports configuration: `DBHOST` and `DBPORT` for MySQL database; `ANNOTATORHOST`
  denoting the server storing Affymetrix raw data, and Ensembl local version; `DATAHOST`,
  an additional backup machine; `PIPEHOST`, name of the machine on which the pipeline
  is run;
  
## Running and re-running pipeline steps

* To re-run the last operation performed by a pipeline step, remove its step_verification_RELEASE.txt file.
  To re-run the step all from scratch, use the command make clean.
  In that case, data inserted in the database are not cleaned automatically, for safety, you would need to remove inserted data yourself.
  This documentation often explains how to do it.
  The command clean only takes care of the generated files.


