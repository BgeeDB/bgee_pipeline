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

## Description of the expression data file

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



