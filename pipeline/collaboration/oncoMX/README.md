Requirements: having successfully run all previous steps, this is a post-processing steps.

Goal: generate expression data for the oncoMX project (Analyze of gene expression of cancers for 28 different organs)

## Details

* `make generateFiles`: generates the processed expression files for oncoMX

## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`

## Data verification

* Check the files in `generated_files/collaboration/oncoMX/`, notably the verification file `generateFiles`.

## Information about the files generated for OncoMX

### Specifications for the files containing calls of presence/absence of expression:

* generates one file for human, one file for mouse
* based only on bulk RNA-Seq data
* includes only data corresponding to OncoMX cancer types (see file
[human_doid_slim_uberon_mapping.csv](../../../source_files/collaboration/oncoMX/human_doid_slim_uberon_mapping.csv).
We retrieve data for the Uberon terms in this file, and all their children.
* includes only data for adult developmental stages (and their children)
* Compute a gene-centric qualitative expression level and an organ-centric qualitative expression level
with the categories "high/medium/low/absent" (see [Computation of the expression levels](#computation-of-the-expression-levels))

### Selection of the anatomical entities

We retrieve the anatomical entities selected by OncoMX (stored in the file
[human_doid_slim_uberon_mapping.csv](../../../source_files/collaboration/oncoMX/human_doid_slim_uberon_mapping.csv)),
column `UBERON_doid.owl`), and all the children of these terms in the Uberon ontology
by `part_of` and `is_a` relationships.

### Selection of the developmental stages

In order to be consistent between human and mouse data, we have selected the developmental stage
`UBERON:0000113 post-juvenile adult stage`. We retrieve data for this stage, and all its child stages
by `part_of` relations.

* In human, this covers the stages from adolescent (13 yo) to death
(see [report for human dev. stages](https://github.com/obophenotype/developmental-stage-ontologies/blob/master/external/bgee/report.md#homo-sapiens))
* In mouse, this covers the stages from early adult stage (6 weeks) to death
(see [report for mouse dev. stages](https://github.com/obophenotype/developmental-stage-ontologies/blob/master/external/bgee/report.md#mus-musculus))

Note: since OncoMX wanted only data in adult stages, maybe a better dev. stage for human would be `HsapDv:0000087 human adult stage (human)`
(from 19yo to death).

### Computation of the expression levels

We produce for OncoMX expression level categories. For `absent` expression calls, the category is alway `ABSENT`.
For `present` expression calls, we produce one expression level category relative to the expression levels
of the gene, and one expression level category relative to the expression levels in the anatomical entity.
The value for these categories is one of `HIGH`, `MEDIUM`, and `LOW`.


For `present` expression calls, the category is computed by comparing the rank score of the expression call
to the min. and max rank scores of an entity. So, for a same expression call, we compute two expression level categories:
* one relative to the gene (by comparing the expression rank score of the call to the min. and max ranks
of the gene, in any anatomical structure where it is expressed).
* one relative to the anatomical entity (by comparing to the min. and max ranks in the anatomical entity,
considering any gene expressed in it)


As an example, if a gene has the following expression calls:

```
Gene ID  Anat. entity ID Expression call Expression rank
Gene1    AnatEntity1     EXPRESSED       2500
Gene1    AnatEntity2     EXPRESSED       10000
Gene1    AnatEntity3     EXPRESSED       22500
```

And the following genes are expressed in `AnatEntity1`:

```
Gene ID  Anat. entity ID Expression call Expression rank
Gene2    AnatEntity1     EXPRESSED       1000
Gene1    AnatEntity1     EXPRESSED       2500
Gene3    AnatEntity1     EXPRESSED       4000
```

As a result:
* The min. and max ranks for `Gene1` are respectively `2500` and `22500`.
The expression level category for `Gene1` in `AnatEntity1`,
relative to the expression levels of `Gene1`, is `HIGH`
(threshold for HIGH: 9167; threshold for MEDIUM: 15833).
* The min. and max ranks in `AnatEntity1` are respectively `1000` and `4000`.
The expression level category for `Gene1` in `AnatEntity1`,
relative to the expression levels in `AnatEntity1`, is `MEDIUM`
(threshold for HIGH: 2000; threshold for MEDIUM: 3000).

As shown above, the category is computed by dividing the range between the min. and the max ranks
into three categories. Only `present` expression calls are considered for retrieving
the min. and max ranks.

**Important:** if the difference between the max and the min. ranks is less than `100`,
the category of the expression call is always `HIGH`.

### Header of the generated files

* Ensembl gene ID
* Gene name
* Anatomical entity ID
* Anatomical entity name
* Developmental stage ID
* Developmental stage name
* Expression level relative to gene
* Expression level relative to anatomical entity
* Call quality
* Expression rank score
* Expression score

## Former version of the file generated for OncoMX

**OUTDATED, NOT USED ANYMORE, AS OF FEB. 2019**

This version of the file was generated by the perl script `generate_oncoMX_expression_files.pl`

* This file contains expression data represented in 6 columns:

geneId :  corresponds to UniProt/SwissProt ID of human genes (18807 protein-coding genes).

anatEntityId : the anatomical entity where the expression has been observed. These anatomical entities come from the Uberon ontology (http://uberon.github.io/).

score : A value between 0 and 10. This score combines 2 levels of information :
	  - expression level : between 0 and 5 (0 => no expression 5 => highest expression),
	  - tissue specificity : +5 if of this gene expression is classified as tissue specific, only for the tissue of highest expression of the gene.

Expression Level : Expression is reported here based only on TPMs from RNA-seq of curated healthy samples. In the Bgee database we can have more than one expression experiment per gene in the same conditions. The expression level in each line of this file is calculated using the average of all Bgee expression for the same gene and condition, over experiments and samples. We only take into account 'present' expression values (see Pres/Abs column) when calculating the average. We used this approach because what we call condition doesn't take into account the full range of possible condition variation (e.g circadian cycle, preprandial/postprandial, tiredness, temperature). 
In Bgee average values ofTPM are between 0 and 460000 for each pair of gene and condition in human. In order to create an expression score between 0 and 5 we normalized the log2 of the average TPM value. 
				
Tissue Specificity : In this file we calculated the tissue specificity using Tau, which is the best ranked method to measure expression specificity (https://doi.org/10.1093/bib/bbw008). Calculation of the tissue specificity of each gene takes as input one TPM value for each anatomical entity where this gene is expressed. In Bgee we have more than one expression value for each pair gene â€“ anatomical entity, as we also take into account the sex and dev. stages. That is why, in order to calculate the tissue specificity, we took the maximum value of all average expression (see "Expression Level" part) for each pair gene - anatomical entity
				
	Score examples :
	0.9 -> low expression of a non tissue specific gene
	4.6 -> high expression of a non tissue specific gene
	5.9 -> low expression of a tissue specific gene
	9.6 -> high expression of a tissue specific gene

devStageId : The developmental stage for which expression has been observed. We only contains report expression during from the post-juvenile adult stage (UBERON:0000113) 
and its children (more specific stages included within it). 
This means that all expression data coming from the following developmental stages has been removed:
	- Sexually immature stage (UBERON:0000112) => between 1 month and 12 years
	- All embryonic developmental stages.
				
sex : The sex corresponding to samples where expression data come from, when known.

Pres/Abs : Define if we detected presence or absence of expression. Based on mapping RNA-seq reads both to genes and to intergenic regions; genes are called "present" if their TPM in one sample is significantly above the distribution of TPMs for intergenic regions in the same sample. An expression defined as 'present' means that expression is found present for this gene in this condition in at least one experiment, but it does not mean that all experiments detect the presence of expression. An expression defined as "absent" means that all experiments report absence of expression for this gene in this condition. 



FYI : In this document we call condition a unique combination of :
	- anatEntityId
	- sex
	- devStageId



