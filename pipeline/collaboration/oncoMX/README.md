Requirements: having successfully run all previous steps, this is a post-processing steps.

Goal: generate expression data for the oncoMX project (Analyze of gene expression of cancers for 28 different organs)

## Details

* `make oncoMX`: generates the processed expression files for oncoMX

## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`

## Data verification

* Check the files generated_files/colaboration_output_files/oncoMX/oncoMX_expression_'speciesId'.tsv

## Documentation to send to oncoMX 

* Last edition : 7th September 2017

* Content : This file contains expression data represented in 6 columns:

geneId : 		corresponds to UniProt/SwissProt ID of human genes (18807 genes)
anatEntityId : 	The developmental stage for which expression has been observed.
				Only contains genes expressed during the developmental stage called post-juvenile adult stage(UBERON:0000113). 
				It means that all expression data coming from following developmental stages have been removed :
					- Sexually immature stage (UBERON:0000112) => between 1 month and 12 years
					- All embryonic dev. stages
Expression data correspond to expression of a gene in some conidtion (anatomical entity, dev. stage and sex)

The file 
The score value is between 0 and 10. This score combines 2 information :
	- expression level : between 0 and 5 (0 => no expression 5 => best expression 



