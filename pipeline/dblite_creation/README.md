Goal: Create the Bgee lite database.

## Information about the Bgee lite database

### Why creating a Bgee lite database ?

The aims of the Bgee lite database are :
* usable in the BioSODA project
* creating a dump lighter than the bgee database one
* allowing users to use a local database containing all expressed calls
* allowing users to use a database with an easy to understand schema

### Tables present in the schema of Bgee lite

* species : all species present in Bgee database
* gene : all genes (the primary key is an internal bgeeGeneId that is not comparable in different releases of the DB)
* anatEntity : anatomical entities as described in the Uberon ontology (https://uberon.github.io/)
* stage : stages as described in the dev. stage ontology (https://github.com/obophenotype/developmental*stage*ontologies)
* globalExpression : calls as shown in the gene page of the Bgee 14 webapp. 
* globalCond : condition describing under which globalExpression were studied. One condition correspond to one dev. stage, one anat. entity and one species.  

### Information available only in Bgee lite

* gene expression calls corresponding to presence of expression with their quality (SILVER or GOLD)
* condition (developmental stage, anatomical entity, and species) under which the gene expression calls have been detected

### Information NOT available in Bgee lite

* Bgee normalized expression rank
* calls corresponding to absence of expression
* Evidence (raw data) used to generate the calls.
* transcript
* Xrefs, synonymes
* taxonomy, OMA HOGs
* structure of ontologies (relations, taxon constraints, ...)

## Description

This step will create a lite version of Bgee following the database schema defined in the file bgeeLiteSchema.sql
The new database is named `bgeelite_vRELEASE` (e.g., `bgeelite_v14`).

## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`
  
## Detail all rules

* `make create_schema`: create the schema of the Bgee lite database
* `make extract_data_from_bgee`: create a tsv file containing all data extracted from Bgee database that should be persisted in the Bgee lite database
* `make import_to_bgeelite`: use the previously generated tsv file and import all data it contains in the Bgee lite database
* `make drop_bgee_lite`: drop the Bgee lite database. 

## Data verification

* Check all tsv files in `generated_files/dblite_creation/`. One file is created for each table of the Bgee lite database


