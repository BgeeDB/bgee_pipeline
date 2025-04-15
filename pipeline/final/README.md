**Requirements**: having successfully run all other insertion steps for all data types, genes, organs, stages, gene families

**Goals**:
* generate useful files to help web indexing and build Sphinx indexes

## Details

- Generate a **sitemap.xml** index file with all sub-sitemap files currently for main pages and gene pages.
- Generate or update Sphinx indexes
- Generate a dump of our schema.org metadata

## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`

## Data verification

## Error handling

## Other notable Makefile targets

