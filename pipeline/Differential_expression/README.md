# run differential expression analysis



* **Requirements**: having successfully run the step database creation AND insert species and taxonomy AND insert genes AND insert ontology AND insert stages AND insert RNA_Seq AND insert Affymetrix

* **Goal**: Create Differential expression files for RNA-Seq and Affymetrix and insert them in the Bgee database.

## Details

* This step Create Differential expression analysis files for RNA-Seq and Affymetrix and insert them in the Bgee database.
* It will backup previous release data
* It will run differential expression analysis for one stage and different anatomical entities.
* It will run differential expression analysis for one anatomical entity and different stages.
* It will insert above differential expression analysis into the bgee database by filling in deaAffymetrixProbesetSummary, deaRNASeqSummary, and differentialExpression tables.
* It will run sex differential expression analysis for RNA-Seq (same anatomical entity AND same stage).

## Data generation

* If it is the first time you execute this step in this pipeline run:
```
make clean
```

* Run Makefile:
```
make
```

## Details on stage/organ differential expression

**launch_diff_analysis_rna_seq.pl** and **launch_diff_analysis_affymetrix.pl**
* Retrieve information from the bgee database
* Remove "life stage" (UBERON:0000104) from analyses because it is the development stage root!
* Run run differential expression analysis for one stage and different anatomical entities AND for one anatomical entity and different stages
* Done only if at least 2 replicates for a given condition in one experiment
* Done only if at least 3 condition have replicates in one experiment. The more condition you have the more you will be able to detect tissue/stade differentially expressed genes 
(e.g a comparison of only liver and heart will not allow to determine that a gene is more expressed in liver than in other tisues but only that it is more expressed in liver than in heart. It is totally different if the comparison is between liver, heart, lung, kidney, ...)

**diff_analysis_rna_seq.R** and **diff_analysis_affymetrix.R**
* calcNormFactors function: 
	* normalize biological differences in RNA composition between samples (because the proportion of reads attributed to a given gene in a library depends on the expression properties of the whole sample rather than just the expression level of that gene)
	* more information in Robinson MD and Oshlack A. 2010. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol 11:R25.
* limma - voom : 
	* more information in Law C, Chen Y, and al. 2014. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15:R29.
	* apply normal-based microarray-like statistical methods to RNA-seq read counts
	* deal with the problem of large read counts having much larger standard deviations than small read counts
* voomWithQualityWeights function
 	* more information in Liu R, Holik AZ, and al. 2015. Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses. Nucleic Acids Res 43:e97
 	* down-weight the observations from more variable samples, thus retaining the maximum degrees of freedom whilst discounting noisy observations
* the counts matrix used as input by edgeR is a matrix of non logged TPMs. Then edgeR calculate the Fold Change (FC) of the differential expression. Finally, these FC values are log2 transformed by edgeR. The result of DE analysis are then expressed using logFC

## Details on sex differential expression

**launch_sex_diff_analysis_rna_seq.pl**
* Retrieve information from the bgee database
* Remove "life stage" (UBERON:0000104) from analyses because it is the development stage root!
* Done only if at least 2 replicates of both male and female for a given condition in one experiment
* Only biological replicates are used for the analyses. Technical replicates are removed using the annotation files (source_files/RNA-Seq/RNASeqLibrary.tsv AND source_files/RNA-Seq/RNASeqLibrary_worm.tsv)

**diff_sex_analysis_rna_seq.R**
* calcNormFactors function: 
	* normalize biological differences in RNA composition between samples (because the proportion of reads attributed to a given gene in a library depends on the expression properties of the whole sample rather than just the expression level of that gene)
	* more information in Robinson MD and Oshlack A. 2010. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol 11:R25.
* limma - voom : 
	* more information in Law C, Chen Y, and al. 2014. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15:R29.
	* apply normal-based microarray-like statistical methods to RNA-seq read counts
	* deal with the problem of large read counts having much larger standard deviations than small read counts
* voomWithQualityWeights function
 	* more information in Liu R, Holik AZ, and al. 2015. Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses. Nucleic Acids Res 43:e97
 	* down-weight the observations from more variable samples, thus retaining the maximum degrees of freedom whilst discounting noisy observations
* the counts matrix used as input by edgeR is a matrix of non logged TPMs. Then edgeR calculate the Fold Change (FC) of the differential expression. Finally, these FC values are log2 transformed by edgeR. The result of DE analysis are then expressed using logFC

## Error handling

* For RNA-Seq data, be sure that the value of both variables RNASEQALLRES and RNASEQABUNDANCEFILE (in the Makefile.common file) match the location of the input data
* Check the queries in the files:
    * launch_diff_analysis_affy.pl
    * launch_diff_analysis_rna_seq.pl
    * launch_sex_diff_analysis_rna_seq.pl

## Other notable Makefile targets

* generate species specific sex differential expression files for RNA-Seq
```
make generate_files_sex_diff_expression_rna_seq
```
