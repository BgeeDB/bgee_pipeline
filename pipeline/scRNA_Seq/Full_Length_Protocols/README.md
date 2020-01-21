# Bgee: Single Cell RNA-Seq data analysis pipeline


**General information:**

1. [Introduction](#introduction)
2. [Step 1: Data annotation](#step-1-data-annotation)
3. [Step 2: Verification, data download and preparation of info file](#step-2-verification-data-download-and-preparation-of-info-file)
   1. [Verification: metadata from source and quantity of cells](#verification-metadata-from-source-and-quantity-of-cells)
   2. [Data download](#data-download)
   3. [Preparation of informative file](#preparation-of-informative-file)
   
4. [Step 3: scRNA-Seq library analyses](#step-3-scrna-seq-library-analyses)
   1. [Data preparation](#data-preparation)
   2. [Pseudo-alignment](#pseudo-alignment)
   3. [Quality control for cell population](#quality-control-for-cell-population)  
   4. [Result processing at individual cell](#result-processing-at-individual-cell)
   5. [Result processing at cell population](#result-processing-at-cell-population)
   6. [Post-processing: expression calls and rank computation](#post-processing-expression-calls-and-rank-computation)

**Detailed guidelines:**

1. [Preparation steps](#preparation-steps)
2. [Mapping the libraries](#mapping-the-libraries)
3. [Validate cell-type and experiment](#validate-cell-type-and-experiment)
4. [Export global information per cell population ](#export-global-information-per-cell-population )
5. [Presence calls](#presence-calls)
  

## General information:

### Introduction

scRNA-Seq data are used to produce:

* baseline calls of presence of expression
* ranking of these baseline calls to identify the most important conditions with expression, for each gene
* calls of differential over-/under-expression

These results are then integrated in a consistent manner with all other results from other data types, to produce global calls of expression, and to compute gene expression ranks for calls of expression.

### Step 1: Data annotation

For scRNA-Seq we manualy annotated healthy WT data using information from GEO or papers, or provided by the Model Organism Database WormBase.
All the data treated are present in SRA repository. 
The protocols selected are full-length protocols, mainly `SMART-Seq`, `SMART-Seq2` and `SMARTer Ultra Low`.

The annotation of each cell-type is done by using the scientific information provided by the previous sources.


### Step 2: Verification, data download and preparation of info file

#### Verification: metadata from source and quantity of cells

After the annotation process, were each library correspond to an individual cell, we verify if the Bgee annotation are in concordance with repository metadata from where we download the data. 
In the next step we validate experiments based on a minimum quantity of cells per cell-type in each experiment and species.

We retrieve experiments to the next step of the pipeline if they have at least 100 cells per cell-type.

#### Data download

The download of the data is done after cleaning the experiments that not match the requirement in the previous step.
This allow us to save computation time in the dowloading process.
The cleaning data annotated are downloaded from SRA using wget function in R. All files extracted are FASTQ files.

GTF annotation files and genome sequence fasta files are retrieved from Ensembl and Ensembl metazoa for all species included in Bgee (see RNA-Seq pipeline). 
This information is used to identify sequences of genic regions, exonic regions, and intergenic regions, as described in the `RNA-Seq pipeline`. It is also used to generate indexed transcriptome files for all species in Bgee, using the `TopHat` and `Kallisto` software.

#### Preparation of informative file

An informative file is created in the following step of the pipeline by collecting information from the manual annotation and from the first quality control step using `FASTP` software for each FASTQ file that correspond to each individual cell. 

### Step 3: scRNA-Seq library analyses

For each independent library, we do:

#### Data preparation

* For each cell/library, check for presence of single-end FASTQ read file, or of the two FASTQ files for paired-end runs.
* For each cell/library, estimation of read length, by using the mean of all reads of FASTQ file determined by `FASTP`, in order to check which K-mer length should be applied.
* If the reads are too short, less 31 bp, for Kallisto indexing with default k-mer length, the k-mer length is set to 15 nucleotides. 
* A FASTP file is generated for each FASTQ file, this means each individual cell, to check for potential problems.
* FASTP file also provide stats information about possible trimmed sample.

#### Pseudo-alignment

The following parameters are used:

* No bootstrapping
* K-mer length for indexing: see `RNA-Seq pipeline`
* For single-end libraries, it is needed to provide the fragment length, we use by default: 180 bp and a value of 20bp for sd.
* For paired-end libraries, Kallisto can estimate the fragment length and sd, so it is not needed to provide this information.

#### Quality control for cell population

* For each cell-type population we compute the ratio of how many times a gene is detected across the number of cells (TPM > 0).
* Verify if the cell-type population follow a bimodal distribution.

#### Result processing at individual cell

From the Kallisto output, for each genomic feature, counts of pseudo-aligned reads are retrieved.
The pseudo-aligned read counts, and the genomic feature effective lengths, are used to compute TPM and FPKM values. 
We sum at the gene level the counts of pseudo-aligned reads computed at the transcript level by Kallisto.
The values provided in Bgee are computed on the basis of genic regions only, but we also compute this information considering other genomic features,this means for calling genes as present.


#### Result processing at cell population

* Sum calls

Computation of the ratio for each gene that belongs to the same cell-type population per experiment and species.

![Boxplot](img/CodeCogsEqn.png)

Bgee also provide information about the cut-off values applied to call genes as present for the cell-type population.

* Sum raw counts

From the Kallisto output, for each cell-type population that belongs to same experiment and species the counts of pseudo-aligned reads are retrieved and summed, the effective length are recalculated based on the weigth mean.
The pseudo-aligned summed read counts, and the genomic feature weigthed mean effective lengths, are used to compute TPM and FPKM values. 
We sum at the gene level the counts of pseudo-aligned reads computed at the transcript level by Kallisto.

#### Post-processing: expression calls and rank computation

##### Expression calls

* Per individual cell

To define the call of expression of a gene in a library as "present", we check whether its level of expression is over the background transcriptional noise in this library. To estimate the background transcriptional noise in each library, we use the level of expression of a set of intergenic regions (described in the RNA-Seq pipeline). 
How we define this set of intergenic regions is described in the developer documentation section of RNA-Seq pipeline.

* Per cell population

To define the call of expression of a gene in a cell population as "present", we compute the density ratio.
Based on the densoty ratio the gene can be classified as present with high confidence or low confidence of expression.
The classification is based on 2 main cut-offs: proportion (proportion intergenic/proportion coding = 0.05) and density (density ratio where the density of intergenic region is low then the density of protein coding).

##### Rank computation

To calculate the ranks we used the TPM values at cell population level (output file from the Sum Raw Counts).

## Detailed guidelines:

### Preparation steps

The preparation step is done by executing 4 main R scripts available in the folder [0Preparation/](0Preparation/)

* [pre_process_control_annotation.R](#pre-process-control-annotation-R)
   * Clean experiments that not pass requirement of minimum number of cells.
   
* [retrieve_metadata.R](#retrive-metadata-R)
   * Retrive metadata of the libraries that need to be downloaded.

* [download_cleaning_data.R](#download-cleaning-data-R)
   * Download all data annotated with minimum requirements.
   
* [prepare_scrna_seq_sample_info.R](#prepare-scrna-seq-sample-info-R)
   * Create an info file about all libraries by collecting information from manual annotation and FASTP software


### Mapping and analysis of the libraries

To run the mapping for each library a R script from the folder [1Run/](1Run/) should be executed:

* [kallisto_analysis.R](#kallisto-analysis-R)
   * The first step of this script allow to run the pseudo-alignment for each individual library.
   * The second step allow to perform the analysis where the transcript are summed to gene level and the TPM and FPKM recalculated, more details can be founded in the RNA-Seq pipeline for this second step.

In this step we have in consideration that exist a folder that was generated from RNA-Seq pipeline for all species with the following information:
   * transcriptome index (15 k-mer and 31 k-mer)
   * gene2transcript file
   * gene2biotype file

### Validate cell-type and experiment

In order to validate if an experiment should be integrated in Bgee, a validation createria is applied by using a quality control script in [1Run/QC_cellPopulation.R](1Run/QC_cellPopulation.R).
The criteria is to verify if a cell-type from an experiment and species follow a bimodal distribution.
For this we quantify how many times the TPM value is higher then zero for each gene that belongs to biotype protein coding across the number of the cells and then calculate the ratio. After that a R `package LaplacesDemon` is applied to determine if the ratio distribution is indeed bimodal. 

![Boxplot](img/validation_experiment.png)


### Export global information per cell population 

* Sum of raw counts

In order to provide to the user a global information per cell-type, experiment and species we provide a abundance file with sum of est_counts from kallisto for each transcript, followed by the weight mean of effective length across all cells and then the recalculation of TPM and FPKM. The information in the end is reported at gene level.
This is executed using the script: [1Run/Sum_RawCounts.R](1Run/Sum_RawCounts.R).

### Presence calls

In this step of the pipeline to call present genes (note that in single cell RNA-Seq data we don't call absent genes) we use the reference intergenic regions from the RNA-Seq pipeline. This comes from the fact that the density of deconvolute intergenic regions are tendencialy less noisy, as is showed in the graphics below, where the summed by species using just single cell data or by pulling all libraries together, this means RNA-Seq and scRNA-Seq, provide a more high overlap between intergenic and coding regions.

![Boxplot](img/sum_by_species_2.png)

* Per cell

In order to call present genes using the intergenic regions the script [1Run/scRNA-Seq_Callpresent.R](1Run/scRNA-Seq_Callpresent.R) should be executed.
After call expressed genes for each individual cell, the output is exported with information about present genes in each library, as well as, is collected a global information about each library in the file `All_samples.tsv`. In order to summarize the information in a visual way a global plot with the distribution of proportion of coding present is exported, as represented below for each species.
We should have in consideration, that because of high propotion of zeros in single cell RNA-Seq data, exist a proportional effect in the proportion of protein coding genes called present.

![Boxplot](img/ProportionCodingIndividualCell.png)

* Per cell-type population

In order to provide a call of expressed genes at cell-type population level, we introduze an approach based on the aggregation of the calls, by executing the following code [1Run/Sum_Calls_cellPopulation.R](1Run/Sum_Calls_cellPopulation.R).
The confidente approach at population level allow us to call expressed genes based in two criteria: threshold ratio cut-off based on the proportion of intergenic and coding regions and the second criteria is the density, this means the cut-off is applied when the density of intergenic regions is lower than the coding regions, as demonstrated in the graphic below.

![Boxplot](img/cellpopulation.png)

Based on this approach we are able to call "high confident" expressed genes at population level (all genes at the right of the threshold ratio cut-off) and "low confident" expressed genes (all genes between the density cut-off and threshold ratio cut-off).

At the population level a `Sum_Calls file` is exported for each cell-type population, as well as, the plot with cut-offs.
In a global way a `Stats_SumFile.tsv` is exported as well as the plot about all cells and species present in the annotation file, as represented below.


![Boxplot](img/SumCalls_allinfo.png)



