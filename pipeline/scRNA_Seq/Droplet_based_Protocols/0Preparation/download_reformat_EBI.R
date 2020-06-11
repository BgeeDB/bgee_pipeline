## SFC, 17 Sep 2019
## This script is used to download the raw data for target-based protocols (for the moment 10X)
## The files are bam file from EBI and converted to FASTQ.GZ in order to run all the analysis from the scratch

## Usage:
## R CMD BATCH --no-save --no-restore '--args metadata_info="metadata_info_10X.txt" librariesDownloadedJura="librariesDownloadedJura.tsv" output="output_folder" bamtofastq="bamtofastq"' download_reformat_EBI.R download_reformat_EBI.Rout
## metadata_info --> File with all libraries that need to be download 
## librariesDownloadedJura --> File with information about the libraries that already exist in jura server and don't need to be downloaded
## output --> Path where should be saved the output results for each library
## bamtofastq --> directory where is bamtofastq software from 10X

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
command_arg <- c("metadata_info", "librariesDownloadedJura","output", "bamtofastq")
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

libraryName <- select(finalLibsToDown, experiment_accession)
ftp <- as.character(finalLibsToDown$submitted_ftp)
ftp <- data.frame(unlist(strsplit(ftp, ";")))
colnames(ftp)<- "ftp_link"
ftp$ftp_link <- as.character(ftp$ftp_link)
ftp[is.na(ftp)] <- 0
ftp <- ftp %>% filter(!str_detect(ftp_link, ".bam.bai"))

## create final information of library and bam file path
generalInfo <- data.frame(libraryName, ftp)
## remove all libraries that don't have information to make the download from the bam files
generalInfo <- dplyr::filter(generalInfo, ftp_link != 0)
generalInfo <- dplyr::filter(generalInfo, ftp_link != "NA")
## remove libraries that have HCA in the information
generalInfo <- dplyr::filter(generalInfo, ftp_link != "HCA")


for (library in  unique(generalInfo$experiment_accession)) {
  
  ## create output for each library
  cat("treating the library: ", library, "\n")
  InfoFile <- paste0(output, "/", library)
  if (!dir.create(InfoFile)){
    dir.create(InfoFile)
  } else {
    print("File already exist.....")
  }
  
  ## select URL for the correspondent library
  ftpID <- generalInfo[generalInfo$experiment_accession %like% library,][,2]
  setwd(InfoFile)
  ## download file
  wget(c(paste0(ftpID)))
}

AllFilesDownload <- list.files(output,  pattern=".*bam$", full.names=T, recursive = TRUE)

## After download all the libraries convert the bam files to Fastq format
for (folder in AllFilesDownload) {
  #bamFile <- dir(paste0(output,library), pattern = "*.bam", full.names = TRUE, ignore.case = TRUE)
  output <- normalizePath(dirname(folder))
  output <- paste0(output,"/FASTQ")
  system(sprintf('%s %s %s', bamtofastq, paste0(folder), paste0(output)))
}

