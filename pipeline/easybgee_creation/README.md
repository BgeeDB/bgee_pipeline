Goal: Create the Easy Bgee database.

## Information about the Easy Bgee database

### Why creating a Easy Bgee database ?

The aims of the Easy Bgee database are :
* used to generate Bgee triple store
* creating a dump lighter than the bgee database one
* allowing users to use a local database containing all expressed calls
* allowing users to use a database with an easy to understand schema
* usable in the BioSODA project

### Tables present in the schema of Easy Bgee

* species : all species present in Bgee database
* gene : all genes (the primary key is an internal bgeeGeneId that is not comparable in different releases of the DB)
* anatEntity : anatomical entities as described in the Uberon ontology (https://uberon.github.io/)
* stage : stages as described in the dev. stage ontology (https://github.com/obophenotype/developmental*stage*ontologies)
* globalExpression : calls as shown in the gene page of the Bgee 14 webapp. 
* globalCond : condition describing under which globalExpression were studied. One condition correspond to one dev. stage, one anat. entity and one species.  

### Information available only in Easy Bgee

* gene expression calls corresponding to presence/absence of expression with their quality (SILVER or GOLD), rank and expression score
* condition (developmental stage, anatomical entity, and species) under which the gene expression calls have been detected

### Information NOT available in Easy Bgee

* Evidence (raw data) used to generate the calls.
* transcript
* Xrefs, synonymes
* taxonomy, OMA HOGs
* structure of ontologies (relations, taxon constraints, ...)

## Description

This step will create a database containing a subset of information present in the Bgee database following the database schema defined in the file easyBgeeSchema.sql
The new database is named `easybgee_vRELEASE` (e.g., `easybgee_v14`).

## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`
  
## Detail all rules

* `make create_schema`: create the schema of the Easy Bgee database
* `make extract_data_from_bgee`: create a tsv file containing all data extracted from Bgee database that should be persisted in the Easy Bgee database
* `make import_to_easybgee`: use the previously generated tsv file and import all data it contains in the Easy Bgee database
* `make drop_easybgee`: drop the Easy Bgee database. 

## Data verification

* Check all tsv files in `generated_files/easybgee_creation/`. One file is created for each table of the Bgee lite database


