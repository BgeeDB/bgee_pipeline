Goal: Create the Bgee lite database.

## Description

This step will create a lite version of bgee following the DB schema defined in the file bgeeLiteSchema.sql
The new database is named `bgeelite_vRELEASE` (e.g., `bgeelite_v14`).

## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`
  
## Detail all rules

* `make create_schema`: create the schema of the bgeelite database
* `make extract_data_from_bgee`: create a tsv file containing all data extracted from bgee database that should be persisted in the bgeelite database
* `make import_to_bgeelite`: use the previously generated tsv file and import all data it contains in the bgeelite database
* `make drop_bgee_lite`: drop the bgeelite database. 

## Data verification

* Check all tsv files in `generated_files/dblite_creation/`. One file is created for each table of the bgeelite database

## Information about the files generated for OncoMX

### Specifications for the files containing calls of presence/absence of expression:

## Information about the bgeelite database

### Why creating a bgeelite database ?

The aims of the bgeelite database are :
* usable in the BioSoda project
* creating a dump lighter than the bgee database one
* allowing users to create a local database containing all expressed calls
* allowing users to create a database with an easy to understand schema

### Tables present in the schema of bgeelite

* species : all species present in bgee database
* gene : all genes (the primary key is an internal bgeeGeneId that is not comparable in different releases of the DB)
* anatEntity : anatomical entities as described in the UBERON ontology (https://uberon.github.io/)
* stage : stages as described in the dev. stage ontology (https://github.com/obophenotype/developmental*stage*ontologies)
* globalExpression : calls as shown in the gene page of the Bgee 14 webapp. 
* globalCond : condition describing under which globalExpression were studied. One condition correspond to one dev. stage, one anat. entity and one species.  

### Information available in bgeelite

* gene expression calls corresponding to presence of expression with their quality (SILVER or GOLD)
* condition (developmental stage, anatomical entity, and species) under which the gene expression calls have been detected

### Information NOT available in bgeelite

* Bgee normalized expression rank
* calls corresponding to absence of expression
* Evidence (raw data) used to generate the calls.
* transcript
* Xrefs, synonymes
* taxonomy, OMA HOGs
* structure of ontologies (relations, taxon constraints, ...)
