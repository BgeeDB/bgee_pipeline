## SFonsecaCosta, September 2019

## This script is used to download the data that pass the minimum requirement of 100 cells and are present in the metadata_info (this means files that are in agreement with Bgee annotation and metadata from EBI)

## Usage:
## R CMD BATCH --no-save --no-restore '--args metadata_info="metadata_info.txt" librariesDownloadedJura="librariesDownloadedJura.tsv" output_folder="output_folder"' download_cleaning_data.R download_cleaning_data.Rout
## metadata_info -> File with all libraries that need to be download
## librariesDownloadedJura -> File with information about the libraries that already exist in jura server and don't need to be downloaded
## output_folder -> Path where should be saved the output results for each library

## libraries used
library(data.table)
library(HelpersMG)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("metadata_info", "librariesDownloadedJura", "output_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read metadata file from EBI. If file not exists, script stops
if( file.exists(metadata_info) ){
  readFile <- read.table(metadata_info, h=T, sep="\t")
} else {
  stop( paste("metadata file not found [", metadata_info, "]\n"))
}
## Read libraries_Downloaded_Jura_file. If file not exists, script stops
if( file.exists(librariesDownloadedJura) ){
  alreadyDownloadedLibraries <- read.table(librariesDownloadedJura, h=F, sep="\t")
  colnames(alreadyDownloadedLibraries) <- "experiment_accession"
} else {
  stop( paste("File with libraries already downloaded not found [", librariesDownloadedJura, "]\n"))
}

###### Check libraries #################################################################

checkId <- data.frame(readFile$experiment_accession[!(readFile$experiment_accession %in% alreadyDownloadedLibraries$experiment_accession)])
colnames(checkId) <- "experiment_accession"
finalLibsToDown <- dplyr::filter(readFile, experiment_accession %in% checkId$experiment_accession)

##########################################################################################
if (nrow(finalLibsToDown) == 0){

  cat("All the libraries present in metadata are already downloaded in the server")

} else {
  ftp <- as.character(finalLibsToDown$fastq_ftp)
  ftp <- data.frame(unlist(strsplit(ftp, ";")))
  colnames(ftp)<- "ftp_link"
  ftp$ftp_link <- as.character(ftp$ftp_link)
  ftp[is.na(ftp)] <- 0

  LibInfo <- c()
  for (i in 1:nrow(finalLibsToDown)) {

    rowinfo <- finalLibsToDown[i,]
    infoL <- rowinfo$experiment_accession

    if(rowinfo$library_layout == "PAIRED"){

      infoL <- as.character(rowinfo$experiment_accession)
      infoL1 <- as.character(rowinfo$experiment_accession)
      infoL2 <- rbind(infoL,infoL1)
      LibInfo <- rbind(LibInfo,infoL2)

    } else {
      libdup <- as.character(rowinfo$experiment_accession)
      LibInfo <- rbind(LibInfo, libdup)
    }
  }

  ## create final information of the library and FASZQ.gz file path
  generalInfo <- data.frame(LibInfo, ftp)

  for (library in  unique(generalInfo$LibInfo)) {

    ## create output for each library
    cat("treating the library: ", library, "\n")
    InfoFile <- file.path(output_folder, library)
    dir.create(InfoFile, showWarnings = FALSE)
    files<-list.files(path = InfoFile, pattern="*.fastq.gz$")

    ## control the downloading process if stops
    if (dir.exists(InfoFile) == TRUE && file_test("-f", paste0(InfoFile,"/",files)) == TRUE){

      cat("The folder for the library ", library, " was already created and the FASTQ.gz files downloaded.", "\n")

    } else {

      cat("The folder for the library ", library, " exist but .gz files are missing.", "\n")
      ## select URL for the correspondent library
      ftpID <- generalInfo[generalInfo$LibInfo %like% library,][,2]
      setwd(InfoFile)
      ## download data
      wget(c(paste0(ftpID)))

    }
  }
}



