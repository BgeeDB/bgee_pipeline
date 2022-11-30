## SFonsecaCosta, 17 Sep 2019

## This script is used to download the raw data for target-based protocols (for the moment 10X)
## For each library we should have at least 3 files: R1, R2 and I (fastq.gz format)

## Usage:
## R CMD BATCH --no-save --no-restore '--args metadata_info="metadata_info_10X.txt" librariesDownloadedJura="librariesDownloadedJura.tsv" output="output_folder"' download_SRA.R download_SRA.Rout
## metadata_info --> File with all libraries that need to be download
## librariesDownloadedJura --> File with information about the libraries that already exist in jura server and don't need to be downloaded
## output --> Path where should be saved the output results for each library

## libraries used
library(dplyr)
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
  stop( paste("metadata file not found [", metadata_info, "]\n"))
}

## Read EBI file. If file not exists, script stops
if( file.exists(librariesDownloadedJura) ){
  libDownloadJura <- read.table(librariesDownloadedJura, h=T, sep="\t")
} else {
  stop( paste("Libraries from JURA file not found [", librariesDownloadedJura, "]\n"))
}

###### Check libraries ######################################################################################################################
checkId <- data.frame(readFile$library_id[!(readFile$library_id %in% libDownloadJura$library_id)])
colnames(checkId) <- "library_id"
finalLibsToDown <- dplyr::filter(readFile, library_id %in% checkId$library_id)

##############################################################################################################################################
## remove NA and fastq.gz information that belongs to HCA (download with other script) or from EBI (download .bam files)
finalLibsToDown <- dplyr::filter(readFile, fastq_ftp != "NA" & fastq_ftp != "fastq.gz")

folder_experiments <- file.path(output, "EXPERIMENTS")
for (library in  unique(finalLibsToDown$library_id)) {

  ## create output for each library
  message("treating the library: ", library)
  species <- unique(finalLibsToDown$tax_id[finalLibsToDown$library_id == library])
  experiment <- unique(finalLibsToDown$experiment_id[finalLibsToDown$library_id == library])
  folder_species <- file.path(output, species)
  if (!file.exists(folder_species)) {
    dir.create(folder_species)
  }
  folder_library <- file.path(folder_species, library)
  if (!file.exists(folder_library)) {
    dir.create(folder_library)
  }
  ## select the SRR referent to the library to download
  SRR_file <- as.character(finalLibsToDown$run_accession[finalLibsToDown$library_id == library])
  setwd(folder_library)

  ## download run files using SRA Toolkit
  for(run in SRR_file) {
    system(sprintf('prefetch %s ', run))
    system(sprintf('fastq-dump --outdir %s --split-files %s', file.path(run, "FASTQ"),
      file.path(run, paste0(run, ".sra"))))
  }
    ## create a symlink following pattern */EXPERIMENTS/experiment_ID/library_id/
  file.symlink(from = folder_library, to = file.path(folder_experiments, experiment, library), overwrite = FALSE, recursive = TRUE)
}
