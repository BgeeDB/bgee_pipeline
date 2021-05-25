## Julien Roux 09/2015

# This script calculates edgeR TMM normalization factors for all libraries of an experiment, within one species, platform, libraryType and libraryOrientation
# It is inspired from script pipeline/Differential_expression/diff_analysis_RNAseq.R, which also uses TMM normalization as a first step
# It is launched by launch_calculate_TMM_factors.pl

##### arguments to provide #####
# target_file_path  - path to the file with info on libraries to use (3 colums: 1 - experimentId; 2 - library ID; 3 - path to file with read counts for library)
# output_file_path  - path to the output folder + name of the output file
# R_log_file        - (optional) path to file to where R will write output

## Usage:
## R CMD BATCH --no-save --no-restore '--args target_file_path="target_file_path" output_file_path="output_folder_path"' calculate_TMM_factors.R R_log_file

## Output: one TSV file for each target file

## Required packages
library(edgeR)

## Session info
print(sessionInfo())

## Reading in arguments provided in command line
cmd_args = commandArgs(TRUE)
print (cmd_args)
if (length(cmd_args)==0){
  stop("no arguments provided\n")
} else {
  for (i in 1:length(cmd_args)){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessairy arguments were passed in command line
command_arg <- c("target_file_path", "output_file_path")

for (c_arg in command_arg){
  if (!exists(c_arg)){
    stop(paste(c_arg, "  command line argument not provided\n"))
  } else {
    cat(paste(c_arg,":\t", eval(c_arg),sep=""),"\n")
  }
}
## reading in file with target info ("target_file_path"), if file not exists script stops
if (file.exists(target_file_path)){
	target <- read.table(target_file_path, sep="\t", as.is=T, h=T)
} else {
  stop(paste("target file: ", target_file_path," not exists\n"))
}

#################################################
## create an exprsMatrix from read count data ##
#################################################

## reading in one file in order to get the number of features (first line of target file)
example_file <- read.table(target$file[1], sep="\t", header=T, as.is=TRUE)
gene_nr <- nrow(example_file)
gene_ids <- example_file$gene_id

## creating expression
exprMatrix=matrix(ncol=length(target$rnaSeqLibraryId), nrow=gene_nr)

## filling in expression matrix
for (i in 1:length(target$rnaSeqLibraryId)){
  expr <- read.table(target$file[i], sep="\t", header=T, as.is=TRUE)
  exprMatrix[,i] = as.numeric(expr$counts)
}
rownames(exprMatrix) <- gene_ids
colnames(exprMatrix) <- target$rnaSeqLibraryId

##############################################
## create DGEList and calculate TMM factors ##
##############################################
DGE_object <- DGEList(counts=exprMatrix)
DGE_object <- calcNormFactors(DGE_object)

## Issue a warning if one TMM factor is > 100 (will not be able to insert it into database)
if ( any(DGE_object$samples$norm.factors >= 100) ) {
  warning( paste0("TMM factor for library ", DGE_object$samples[DGE_object$samples$norm.factors >= 100], " is bigger than 100. This will pose problem for insertion into database (need to update Decimal format).\n") )
}

############
## output ##
############
## Append the TMM factor to the target file
target <- cbind(target$rnaSeqExperimentId, target$rnaSeqLibraryId, DGE_object$samples$norm.factors)
colnames(target) <- c("rnaSeqExperimentId", "rnaSeqLibraryId", "tmmFactor")
## write results in folder specified in output_file_path
write.table(target, file=file.path(output_file_path), sep="\t", quote=F, col.names=TRUE, row.names=FALSE)
