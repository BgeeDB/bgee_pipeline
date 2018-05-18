# Download files with rank score information

The files present in these directories provide information about rank scores computed for the Bgee database (https://bgee.org/).

## Ranks and max ranks

The putative max rank for this release of Bgee, over all species, genes and conditions is: `47947.00`.

A lower rank score for a gene means that the gene is more highly expressed in a condition. We "normalize" all ranks so that they range to the same max value. However, the putative min rank varies between species.

## Directory structure

* `anat_entity/`: provide rank scores grouped by genes and anatomical entities. It means that only one best score is displayed for a gene expressed in an anatomical entity (e.g., the different developmental stages during which the gene is expressed in this anatomical structure are not shown)

* `condition/`: provide rank scores for all conditions where genes are expressed (e.g., developmental stage information is provided).

## Files

For each species and data grouping, 5 files are generated: one considering all data available, one considering Affymetrix data only, one EST data only, one _in situ_ hybridization data only, one RNA-Seq data only.

### File name pattern

File names follow the following pattern: `<NCBI_tax_id>_<grouping>_<data_type>_<species_name>.tsv.zip`. For instance: `9606_condition_all_data_Homo_sapiens.tsv.zip`.

* `<NCBI_tax_id>`: the NCBI taxonomy ID of the species (e.g., `9606` for human)
* `<grouping>`: the criteria which the grouping of the scores was based on, see [Directory structure](#directory-structure) above. Permitted values: 
  * `anat_entity`: rank scores grouped by genes and anatomical entities.
  * `condition`: rank scores provided for all conditions where genes are expressed.
* `<data_type>`: the data types that were considered to produce the rank scores. Permitted values: 
  * `all_data`: all data types were considered to produce a global mean rank score. This is the score we recommend, used on the Bgee website. 
  * `affymetrix`: scores based exclusively on Affymetrix data.
  * `est`: scores based exclusively on EST data.
  * `rna_seq`: scores based exclusively on RNA-Seq data.
  * `in_situ`: scores based exclusively on _in situ_ hybridization data.
* `<species_name>`: the scientific name of the species, with spaces replaced with `_`, e.g. `Homo_sapiens`.

### File header

* `Ensembl gene ID`
* `gene name`
* `anatomical entity ID`
* `anatomical entity name`
* If data grouped by all conditions: `developmental stage ID`
* If data grouped by all conditions: `developmental stage name`
* `rank score`: the information you are looking for ;)
* `XRefs to BTO`: cross-references to the BTO ontology, separated by `|` when several are available for a same anatomical entity. Please contact us <bgee@sib.swiss> if you need cross-references to other ontologies or controlled vocabularies. 

## Notes
The scores are exactly the same as those used on the Bgee website, although some slight differences might appear: 

* differences in ordering: on the Bgee website, conditions with identical rank scores are ordered so that the more precise conditions appear first. For instance, if a gene is expressed in a same organ at different developmental stages, and these calls of expression are associated to the same rank scores, the more precise stages will be displayed first, e.g.: `human late adulthood stage` displayed before `human adult stage`. This ordering is more computer-intensive, and is not performed for generating these files (please contact us <bgee@sib.swiss> if you need this information). Conditions with equal ranks are sorted by alphabetical order of their IDs. 

* differences in anatomical structures displayed: on the Bgee website, we detect when for a call of expression there exists a more precise call at a better rank. We consider these calls to be "redundant". For instance, if there exists a call of expression of a gene in `cerebellum` at `adult stage` with rank score of `1,000`, and another call for the same gene in `brain` at `adult stage` with rank score of `2,000`, then the expression call in `brain` is considered as redundant for web interface purposes. 

In Bgee, if an anatomical entity is supported only by "redundant calls", we hide it on the webinterface. This explains why more information is present in these download files. Detection of redundant calls is more computer-intensive, and is not performed for generating these files. Please contact us <bgee@sib.swiss> if you need this information. 
