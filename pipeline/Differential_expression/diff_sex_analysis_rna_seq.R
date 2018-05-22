## Julien Wollbrett 01/2018, based on diff_analysis_rna_seq.R written by Marta Rosikiewicz and Julien Roux

## This script perform sex differential expression analysis on RNA-seq datasets using voom normalisation and Empirical Bayes Statistic (package limma)

##### arguments to provide
# target_file_path   - path to the file with target info (3 colums: 1 - organID, 2 - stageID, 3 - names of files with results for the samples)
# output_folder_path - path to the output folder
# input_folder_path  - path to the directory with all sample results for all libraries
# abundance_file_name- name of the file containing all abundance gene level
# speciesID          - id of the species
# R_log_file         - (optional) path to file to where R will write output

#### launching command:
# whole --args expression quoted with single quoting : ' '
# argument names are not quoted ex: target_file_path
# argument value (after =) quoted  with double quoting : " "

#### usage
# R CMD BATCH --no-save --no-restore '--args target_file_path="target_file_path" output_folder_path="output_folder_path" input_folder_path="input_folder_path" abundance_file_name="abundance_file_name" speciesID="speciesID" ' diff_analysis_RNA-seq.R R_log_file

#### example:
# R CMD BATCH --no-save --no-restore '--args target_file_path="../targets/sex/GSE30352___10090.target" output_folder_path="../processed_differential_bgee_v14/sex/GSE30352" input_folder_path="../all_results_bgee_v14" abundance_file_name="abundance_gene_level+new_tpm+new_fpkm+calls.tsv" speciesID="10090" ' diff_analysis_RNA-seq.R R_log_file

#### output files
# 4 columns with header, first field starting with '#' ("#gene_names","p_value","logFC","present_calls")
# output file name: speciesID, factor, nr of organs, nr of stages (separeted by "_")
# ".out' ending
# example: 10090_MA_0000072_MmusDO_0000040_6_1.out
#TODO change output file name: speciesID, experiment, organ, devStage (separeted by "_")

#### required packages:
library(edgeR)

# reading in arguments provided in command line
cmd_args = commandArgs(TRUE);

print (cmd_args)
if (length(cmd_args)==0){
  stop("no arguments provided\n")
} else {
  for(i in 1:length(cmd_args)){
    eval(parse(text=cmd_args[i]))
  }
}

# checking if all necessary arguments were passed in command line
command_arg <- c("target_file_path", "output_folder_path", "input_folder_path", "abundance_file_name", "speciesID")

for(c_arg in command_arg){
  if(!exists(c_arg)){
    stop(paste(c_arg, "  command line argument not provided\n"))
  } else {
    cat(paste(c_arg, ":\t", eval(c_arg), sep=""), "\n")
  }
}
# reading in file with target info ("target_file_path"), if file not exists script stops
if(file.exists(target_file_path)){
  target <- read.table(target_file_path, sep="\t", as.is=T)
} else {
  stop(paste("target file: ", target_file_path, " not exists\n"))
}

#################################################
## create an exprsMatrix from read count data ##
#################################################

sampleNames <- target[,4];

## reading in one file in order to get the number of genes
example_file <- read.table(file.path(input_folder_path, target[1,4], abundance_file_name), sep="\t", header=T, as.is=TRUE)
gene_nr <- nrow(example_file)
gene_names <- example_file[,1]
anatEntity_name <- target[1,1]
stage_name <- target[1,2]

## creating expression and call matrix
exprMatrix=matrix(ncol=length(sampleNames), nrow=gene_nr)
callMatrix=matrix(ncol=length(sampleNames), nrow=gene_nr)


## filling in expression and call matrix
for (i in 1:length(sampleNames)){
  ## first column used as row names
  expr <- read.table(file.path(input_folder_path, target[i,4], abundance_file_name), sep="\t", row.names=1, header=T, as.is=TRUE)
  ## matrix of TPMs
  exprMatrix[,i] = as.numeric(expr[,2])
  ## expression calls
  callMatrix[,i] = as.character(expr[,5])
}
rownames(exprMatrix) <- gene_names
colnames(exprMatrix) <- sampleNames
rownames(callMatrix) <- gene_names
colnames(callMatrix) <- sampleNames

## absent_calls - equal to 1 for genes that are called absent in all samples (RNA-seq libraries)
absent_calls <- apply(callMatrix, 1, function(x){sum(x == "absent") / length(x)})
## present_calls - true for genes with at least one present call
present_calls <- absent_calls != 1

## phenodata (phD) - "sex"
phD <- as.vector(target[,3])

#########################################################
cat ("factors for each sample:\n")
cat(paste(paste(target[,4], phD, sep="\t"), "\n", sep=""))
pheno_groups <- factor(phD)

## removing of genes that are absent in all samples
present_exprMatrix <- exprMatrix[present_calls, ]

## creating DGEList object from package edgeR
## Beware, the counts matrix is in fact a matrix of (non logged) TPMs
DGE_object <- DGEList(counts=present_exprMatrix, group=pheno_groups)
DGE_object <- calcNormFactors(DGE_object)

## voom normalisation - compute observational weights using TPMs values
## voomWithQualityWeights is used to down-weight low quality samples
design <- model.matrix(~ -1 + pheno_groups)
voom_object <- voomWithQualityWeights(DGE_object, design, plot=FALSE)

## Extract sample weights to export them later
samplew <- voom_object$sample.weights
names(samplew) <- colnames(voom_object$E)

####################
## limma analysis ##
####################

## creating contrast matrix
contrast.matrix <- makeContrasts(sexDiff = pheno_groupsmale - pheno_groupsfemale, levels=design)

## Empirical Bayes Analysis

fit <- lmFit(voom_object, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

## getting p-values
p_value_matrix <- fit2$p.value

## p-value adjusting with Benjamini-Hochberg procedure (FDR); adjustment for each column (sample) separately
adj_p_value_matrix <- apply(p_value_matrix, 2, p.adjust, method="BH")

#all_p_value_matrix - p-value for all genes in input files, p-value for genes absent in all conditions set to 1
all_p_value_matrix <- matrix(1, ncol=1, nrow=gene_nr)
all_p_value_matrix[present_calls, ] <- adj_p_value_matrix

#all_logFC_matrix - log2 fold change for all genes, set to 0 for genes absent in all conditions
all_logFC_matrix <- matrix(0, ncol=1, nrow=gene_nr)
all_logFC_matrix[present_calls, ] <- fit2$coefficients


############
## output ##
############

results <- cbind(gene_names, as.character(all_p_value_matrix[,1]), as.character(all_logFC_matrix[,1]), present_calls)
colnames(results) <- c('#gene_names', 'p_value', 'logFC', 'present_calls')

## write all results in folder specified in output_path
write.table(results, file=file.path(output_folder_path, paste(speciesID, "_sex_", gsub(":", "_",  anatEntity_name), "_", gsub(":", "_",  stage_name), ".out", sep="")), sep="\t", quote=F, col.names=TRUE, row.names=FALSE)


## save sample weigths: we want to be able to know if an sample had a lower quality and was down-weighted in teh analysis
write.table(samplew, file=file.path(output_folder_path, paste(speciesID, "_sex_", gsub(":", "_",  anatEntity_name), "_", gsub(":", "_",  stage_name), "_sampleWeights.tsv", sep="")), sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)

## Session information
print(sessionInfo())

