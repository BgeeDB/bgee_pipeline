## Julien Roux 05/05/09; modified Marta Rosikiewicz 07/14; modified Julien Roux 10/16

## This script perform differential expression analysis on array datasets using Empirical Bayes Statistic (package limma)

## Marta's modifications:
# - false discovery rate (fdr) correction separately for each comparison
# - fdr correction for all genes not only significantly changed
# - differential analysis only for genes that are called present at least in one sample
# - added third column in output file - "present call" which specify whether particular gene was called present in at least one sample from particular experiment
# - minor correction concerning signal values below 0
# - all paths to files and folders pass as command line arguments
# - number of organs and stages included in the name of output file

## Julien's modifications
# - use quality weights for samples in limma analysis. Array quality weights increase statistical power to detect true differential expression without increasing the false discovery rate, see Ritchie et al. 2006 BMC Bioinformatics
# - updated formating and readbility

#### arguments to provide
# target_file_path   - path to the file with target info (3 colums: 1 - organID, 2 - stageID, 3 - names of files with results for the samples ex.: GSM198412.CEL.out)
# output_folder_path - path to the output file
# input_folder_path  - path to the directory with all sample results for given experiment
# input_type         - one of the two terms: "MAS5", "Schuster"
# array_type         - id of the array platform
# R_log_file       - (optional) path to file to where R will write output

#### launching command:
# whole --args expression quoted with single quoting : ' '
# argument names are not quoted ex: target_file_path
# argument value (after =) quoted  with double quoting : " "

# R CMD BATCH --no-save --no-restore '--args target_file_path="target_file_path" output_folder_path="output_folder_path" input_folder_path="input_folder_path" input_type="input_type" array_type="array_type"' diff_analysis.R R_log_file

#### example:
# R CMD BATCH --no-save --no-restore '--args   target_file_path="../Affymetrix/bioconductor/targets/GSE1572___A-AFFY-1.target" output_folder_path="../Affymetrix/processed_differential/GSE1572" input_folder_path="../Affymetrix/processed_schuster/GSE1572" input_type="Schuster" array_type="A-AFFY-1"' diff_analysis.R GSE1572___A-AFFY-1.log

#### output file
# 4 columns with header, first field starting with '#' ("#probe_set_ID","p_value","logFC","present_calls")
# output file name: array_type, factor, nr of organs, nr of stages (separeted by "_")
# ".out' ending


#reading in arguments provided in command line
cmd_args <- commandArgs(TRUE);

print (cmd_args)
if (length(cmd_args)==0){stop("no arguments provided\n")
} else {
  for(i in 1:length(cmd_args)){
    eval(parse(text=cmd_args[i]))
  }
}

# checking if all necessairy arguments were passed in command line
command_arg <- c("target_file_path", "output_folder_path", "input_folder_path", "input_type", "array_type")

for(c_arg in command_arg){
  if (!exists(c_arg)){
    stop(paste(c_arg, "  command line argument not provided\n"))
  } else {
    cat(paste(c_arg, ":\t", eval(parse(text=c_arg)), sep=""), "\n")
  }
}
# reading in file with target info ("target_file_path"), if file not exists script stops
if (file.exists(target_file_path)){
  target <- read.table(target_file_path, sep="\t")
} else {
  stop(paste("target file: ", target_file_path, " not exists\n"))
}

## source("http://www.bioconductor.org/biocLite.R")
## biocLite("affy")
library(affy)

cat("Reading target file:\n", target_file_path, "\n")
target <- read.table(target_file_path, sep="\t")

#################################################
## create an exprs object from normalized data ##
#################################################

sampleNames <- target[,3];

## the format of the data file differs between mas5 and schuster
if (input_type == "MAS5") {
  pbsets <- read.table(file.path(input_folder_path, target[1,3]), sep="\t", header=TRUE, as.is=TRUE)[,1]
}
if (input_type == "Schuster") {
  pbsets <- read.table(file.path(input_folder_path, target[1,3]), sep="\t", header=FALSE, as.is=TRUE)[,1]
}

pbsets_nr <- length(pbsets)
exprMatrix <- matrix(ncol=length(sampleNames), nrow=length(pbsets))
callMatrix <- matrix(ncol=length(sampleNames), nrow=length(pbsets))

for (i in 1:length(sampleNames)){
  ## print(i)
  if (input_type == "MAS5") {
    expr_data <- read.table(file.path(input_folder_path, target[i,3]), sep="\t", row.names=1, header=TRUE, as.is=TRUE)
    expr <- expr_data[,2]
    expr[expr == "null"] <- 0
    expr <- as.numeric(expr)

    ## problem with  values equal or lower than 0, when using the log2 (Marta's correction:  values lower than 0 set to 0.1, calls for all probesets with values equal or lower than 0 are set to "absent")
    calls <- expr_data[,1] == "present"
    calls[expr <= 0] <- FALSE
    expr[expr <= 0] <- 0.1
    expr <- log2(expr);
  }
  if  (input_type == "Schuster") {
    expr_data <- read.table(file.path(input_folder_path, target[i,3]), sep="\t", row.names=1, header=FALSE, as.is=TRUE);
    expr <- as.numeric(expr_data[,1])
    calls <- expr_data[,2] == "P"
  }
  exprMatrix[,i] <- expr
  callMatrix[,i] <- calls
}
row.names(exprMatrix) <- pbsets
colnames(exprMatrix) <- sampleNames

row.names(callMatrix) <- pbsets
colnames(callMatrix) <- sampleNames

## absent_calls - equal to 1 for genes that are called absent in all samples (RNA-seq libraries)
absent_calls <- apply(callMatrix, 1, function(x){sum(x == FALSE) / length(x)})
## present_calls - true for genes with at least one present call
present_calls <- absent_calls != 1


phD <- matrix(nrow=length(sampleNames), ncol=1)
row.names(phD) <- sampleNames
colnames(phD) <- "organ_stage"
phD[,1] <- as.vector(paste(target[,1], target[,2], sep="_"))
## remove : in the name and organ IDs
phD[,1] <- gsub(":", "_", phD[,1])

metaData <- data.frame(labelDescription="organ_stage")
pD <- new("AnnotatedDataFrame",
          data=as.data.frame(phD),
          varMetadata=metaData
)

cat ("factors for each sample:\n")
cat(paste(paste(target[,3], phD, sep="\t"), "\n", sep=""))

## removing of genes that are absent in all samples
present_exprMatrix <- exprMatrix[present_calls, ]
data.norm <- new("ExpressionSet", phenoData=pD, exprs=present_exprMatrix)

####################
## limma analysis ##
####################
library(limma)
factors <- factor(pData(data.norm)$organ_stage)

## cell means model
design <- model.matrix(~ -1+factors)

## Calculate array weigths
arrayw <- arrayWeights(data.norm, design)

## Fit linear model
fit <- lmFit(object=data.norm, weights=arrayw, design=design)

## create a string for the mean in the comparison to the mean contrasts
global_mean <- "("
for (i in 1:length(unique(factors))){
  global_mean <- paste(global_mean, "factors", unique(factors)[i], "+", sep="")
}
global_mean <- substr(global_mean, 1, nchar(global_mean)-1)
global_mean <- paste(global_mean, ")/", length(unique(factors)), sep="")

## create the string of factors combinations for the contrasts
x <- vector()
for (i in 1:length(unique(factors))){
  x <- c(x, paste("factors", unique(factors)[i], "-", global_mean, sep=""))
}

contrast.matrix <- makeContrasts(contrasts=x, levels=design)
##replace the names of the columns because there is a problem when names are too long
colnames(contrast.matrix) <- paste(unique(factors), "-mean", sep="")

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

## correct for multiple testing on all concatenated p-values (all genes sig. for F-test foreach condition)
## getting p-values
p_value_matrix <- fit2$p.value

## p-value adjusting with Benjamini-Hochberg procedure (FDR); adjustment for each column (sample) separately
adj_p_value_matrix <- apply(p_value_matrix, 2, p.adjust, method="BH")

#all_p_value_matrix - p-value for all genes in input files, p-value for genes absent in all conditions set to 1
all_p_value_matrix <- matrix(1, ncol=length(unique(factors)), nrow=pbsets_nr)
all_p_value_matrix[present_calls, ] <- adj_p_value_matrix

#all_logFC_matrix - log2 fold change for all genes, set to 0 for genes absent in all conditions
all_logFC_matrix <- matrix(0, ncol=length(unique(factors)), nrow=pbsets_nr)
all_logFC_matrix[present_calls, ] <- fit2$coefficients

#defining the number of organs and stages
organ_nr <- length(unique(target[,1]))
stage_nr <- length(unique(target[,2]))

cat("\nnr of organs\t", organ_nr, "\n")
cat("\nnr of stages\t", stage_nr, "\n")

############
## output ##
############

for (i in 1:length(unique(factors))){

  #if logFC > 0: over-expressed; if coefficients < 0: under-expressed,
  results <- cbind(as.character(all_p_value_matrix[,i]), as.character(all_logFC_matrix[,i]), present_calls)
  results <- cbind(pbsets, results)
  colnames(results) <- c("#probe_set_ID", "p_value", "logFC", "present_calls")
  ## write all results to folder specified in output_folder_path
  write.table(results, file=file.path(output_folder_path, paste(array_type, "_", unique(factors)[i], "_", organ_nr, "_", stage_nr, ".out",  sep="")), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

## save array weigths: we want to be able to know if an array had a lower quality and was down-weighted in teh analysis
names(arrayw) <- colnames(data.norm)
write.table(arrayw, file=file.path(output_folder_path, "arrayWeights.tsv"), sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)

## Session information
print(sessionInfo())
