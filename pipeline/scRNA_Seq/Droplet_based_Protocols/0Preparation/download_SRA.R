## SFC, 17 Sep 2019
## This script is used to download the raw data for target-based protocols (for the moment 10X)
## FOr each library we should have at least 3 files: R1, R2 and I (fastq.gz format)

## Usage:
## R CMD BATCH --no-save --no-restore '--args metadata_info="metadata_info_10X.txt" librariesDownloadedJura="librariesDownloadedJura.tsv" output="output_folder"' download_SRA.R download_SRA.Rout
## metadata_info --> File with all libraries that need to be download 
## librariesDownloadedJura --> File with information about the libraries that already exist in jura server and don't need to be downloaded
## output --> Path where should be saved the output results for each library

## libraries used
library(dplyr)
library(stringr)
library(data.table)
library(HelpersMG)
library(tools)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("metadata_info", "librariesDownloadedJura","output")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read EBI file. If file not exists, script stops
if( file.exists(metadata_info) ){
  readFile <- read.table(metadata_info, h=T, sep="\t")
} else {
  stop( paste("EBI file not found [", ENAfile, "]\n"))
}

## Read EBI file. If file not exists, script stops
if( file.exists(librariesDownloadedJura) ){
  libDownloadJura <- read.table(librariesDownloadedJura, h=T, sep="\t")
} else {
  stop( paste("Libraries from JURA file not found [", ENAfile, "]\n"))
}

###### Check libraries ######################################################################################################################
checkId <- data.frame(readFile$experiment_accession[!(readFile$experiment_accession %in% libDownloadJura$experiment_accession)])
colnames(checkId) <- "experiment_accession"
finalLibsToDown <- dplyr::filter(readFile, experiment_accession %in% checkId$experiment_accession)

##############################################################################################################################################
## remove NA and fastq.gz information that belongs to HCA (download with other script) or from EBI (download .bam files)
finalLibsToDown <- dplyr::filter(readFile, fastq_ftp != "NA" & fastq_ftp != "fastq.gz") 
for (library in  unique(finalLibsToDown$experiment_accession)) {
  
  ## create output for each library
  cat("treating the library: ", library, "\n")
  InfoFile <- paste0(output, "/", library)
  if (!dir.create(InfoFile)){
    dir.create(InfoFile)
  } else {
    print("File already exist.....")
  }
  
  ## select the SRR referent to the library to downlaod
  SRR_file <- finalLibsToDown$run_accession[finalLibsToDown$experiment_accession == library]
  setwd(InfoFile)
  
  system(sprintf('prefetch --type fastq %s', paste0(SRR_file)))
}
