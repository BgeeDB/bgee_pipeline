# Bgee: Single Cell RNA-Seq data analysis pipeline


**General information:**

1. [Introduction](#introduction)
2. [Step 1: Data annotation](#step-1-data-annotation)
3. [Step 2: Metadata and data download](#step-2-metadata-and-data-download)
   1. [Metadata](#metadata)
   2. [Data download](#data-download)

4. [Step 3: scRNA-Seq library analyses](#step-3-scrna-seq-library-analyses)
   1. [Data preparation](#data-preparation)
   2. [Pseudo-alignment](#pseudo-alignment)
   3. [Processing the files](#processing-the-files)
   4. [Loading the libraries for analysis](#loading-the-libraries-for-analysis)
   5. [Result processing at individual cell](#result-processing-at-individual-cell)
   6. [Result processing at cell population](#result-processing-at-cell-population)
   7. [Post-processing: expression calls and rank computation](#post-processing-expression-calls-and-rank-computation)
   
**Detailed guidelines (Developer):**

## General information:

### Introduction

scRNA-Seq data are used in Bgee to produce:

* baseline calls of presence of expression
* ranking of these baseline calls to identify the most important conditions with expression, for each gene
* calls of differential over-/under-expression

These results are then integrated in a consistent manner with all other results from other data types, to produce global calls of expression, and to compute gene expression ranks for calls of expression.

### Step 1: Data annotation

For scRNA-Seq we manualy annotated healthy WT data using information from GEO or from papers. The annotation of each library is done by using the scientific information provided by these sources.
All the data treated are present in the SRA, EBI ArrayExpress or Human Cell Atlas repositories.
The protocols selected at present for target based are only from `10X platform`.

### Step 2: Metadata and data download

#### Metadata

After the annotation process, where each library corresponds to multiple number of cells, we retrieve information if possible from the repository metadata in order to verify if the Bgee annotation and the repository are in concordance (also to complete information, as for example SRR Id's), if yes a file is write to processed in next steps of the pipeline `metadata_info_10X.txt`. If they are not, we save this information in `metadata_notMatch_10X.txt` and do not further use these libraries.

#### Data download

The download of the data is done only for the experiments that match the requirements in the previous steps.
These data are downloaded from three main sources:

- SRA: The data is downloaded using wget function in R. For each library three fastq files are downloaded, R1, R2 and I, respectively. 

- EBI ArrayExpress: The data is downloaded using wget function in R. All downloaded files for each library are in BAM format. Then using the tool `bamTofastq` from 10X the files are converted to original FASTQ format files in order to make possible all analysis from the scratch. 

- HCA: The data is downloaded using the HCA Command Line Interface (HCA-CLI) (https://data.humancellatlas.org/guides/consumer-vignettes/intro-to-downloading-and-analyzing#step-two-downloading-files-from-a-file-manifest-with-the-hca-cli). In order to do the download of the target libraries a manisfest file need to be downloaded first from the Human Cell Atlas Data Portal. This file is retrieved by experiment. 
The downloaded of each library is done by retrieving directly the FASTQ file. 

GTF annotation files and genome sequence fasta files are retrieved from Ensembl and Ensembl metazoa for all species included in Bgee (see `RNA-Seq pipeline`).
This information is used to identify sequences of genic regions, exonic regions, and intergenic regions, as described in the `RNA-Seq pipeline`. It is also used to generate indexed transcriptome files for all species in Bgee, using the `TopHat` and `Kallisto` software.


### Step 3: scRNA-Seq library analyses

For each independent library, we do:

#### Data preparation

* Check for presence 3 fastq.gz files: R1, I1 and R2 which contains the sequence of the cell barcode + UMI, the sample index, and the cDNA sequence, respectively.
* Estimation of read length, by using the mean of all reads of FASTQ file determined by `FASTP`.
* A FASTP file is generated for each FASTQ file to check for potential problems; it also provide information about possible trimmed samples.

#### Pseudo-alignment

The following parameters are used:

* K-mer length for indexing use the default K-mer size from Kallisto (31 K-mer).
* Pseudoalign the reads with `kallisto bus`.
* The technology argument referent to the single-cell technology is provided depending on the annotation information to the target library, this is passed to the argument as: 10xV2 or 10xV3.
  
#### Processing the files

The output files from Kallisto bus (matrix.ec, output.bus, run_info.json and transcripts.txt), are processed using the `bustools` software to go from the BUS file to a gene-UMI count matrix. The procedue is done in the following way:

* Correcting the barcodes that are within one hamming distance of the barcodes using the whitelist from 10X platform.
* Sort the busfile by organizing the file by: barcode, UMI, set, and multiplicity.
* Generate the UMI count matrix.

#### Loading the libraries for analysis

For each library the gene-UMI count matrix is loaded directly into R for analysis by initialy creating a sparseMatrix. Then the following steps are subsequently performed:

* Using `SEURAT` R package to perform quality control metrics as: filtering cells based on the knee plot.
* Target cell-types in each library by using information of the barcodes list or gene markers provided by the authors to identify specific cells.

#### Result processing at individual cell

For unique cell-type in each individual library a table is exported containing all barcodes (cells per column) and the genomic feature (gene Id per row) with correspondent UMI counts per individual barcode/feature.
The UMI counts are used to compute the normalized values per gene / per cell, this means CPM values.
For calling genes present, we compute genic regions and intergenic regions for each individual cell.

#### Result processing at cell population

* Sum UMI counts 

All gene-UMI count matrix that belongs to the same experiment and species, as well as, cell-typeId, stageId, uberonId, strain and sex are retrieved and the UMI counts are summed across genomic features (gene level) and then the summed UMI count matrix is used to compute the CPM values of the cell population.


#### Post-processing: expression calls and rank computation

##### Expression calls

* Per individual cell
    
To define the call of expression of a gene in a cell as "present", we check whether its level of expression is over the background transcriptional noise in this correspondent cell. To estimate the background transcriptional noise in each individual cell, we use the level of expression of a set of intergenic regions defined in the RNA-Seq pipeline.

* Per cell population

Computation of the ratio calls for each gene that belongs to the same cell population (as well as, all assumptions described above).

To define the call of expression of a gene in a cell population as "present", we compute the density ratio of protein coding and intergenic region and then we check whether its level of expression is over the background transcriptional noise in this correspondent cell population type. 

The calls at the cell population level are based in confidence, this means present with low or high confidence.


##### Rank computation

To calculate the ranks of expression of genes per cell population we use the CPM values, this means, the output files from the Sum UMI counts (see below).