## SFonsecaCosta 2019
## This script allow to verify the number of cells per cell-population regarding the: experimentID, species, cellTypeId, stageId, strain, uberonId and sex after the annotation process.
## Just cell-population that belongs to the same experimentID, species, cellTypeId, stageId, strain, uberonId and sex and have at least 50 cells will be keeped to continue in the pipeline.
## The output file generated (NEW_scRNASeqLibrary.tsv) will be used to download the data and to continue in the pipeline.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeqLibrary="scRNASeqLibrary.tsv" output_folder="output_folder"' pre_process_control_annotation.R pre_process_control_annotation.Rout
## scRNASeqLibrary --> File from manual annotation
## passScRNASeqLibrary --> Output file containing libraries that pass the threshold of minimum number of cells
## notPassScRNASeqLibrary --> Output file containing libraries that did not pass the threshold of minimum number of cells
## cellsThreshold --> minimum number of cells in a library

## Libraries
library(tools)
library(readr)
library(stringr)
library(dplyr)
library(forcats)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("scRNASeqLibrary", "passScRNASeqLibrary", "notPassScRNASeqLibrary", "cellsThreshold")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read scRNASeqLibrary file, if file not exists, script stops
if( file.exists(scRNASeqLibrary) ){
  annotation <- read.table(scRNASeqLibrary, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
  annotation <- dplyr::filter(annotation, annotation$protocolType == "Full-length")
  ## just replace NA in sex and strain per "missing information" to be a level in factor (levels())
  annotation$sex <- forcats::fct_explicit_na(annotation$sex)
  annotation$strain <- forcats::fct_explicit_na(annotation$strain)
} else {
  stop( paste("scRNASeqLibrary file not found [", scRNASeqLibrary, "]\n"))
}

#################################################################################
## create new output files or truncate existing ones
file.create(passScRNASeqLibrary)
file.create(notPassScRNASeqLibrary)

#TODO: maybe collect all pass/not_pass info in a DF and write the 2 files after (directly with the proper colnames)
for (species in unique(annotation$speciesId)) {
  message("Species:", species)
  for (experiment in unique(annotation$experimentId[annotation$speciesId == species])){
    message("Name of experiments that belongs to the species:", experiment)
    for (cellId in unique(annotation$cellTypeId[annotation$experimentId == experiment])){
      message("CellId info:", cellId)
      for (stageId in unique(annotation$stageId[annotation$cellTypeId == cellId])){
        message("StageId info:", stageId)
        for (strain in unique(annotation$strain[annotation$stageId == stageId])){
          message("Strain info:", strain)
          for (uberonId in unique(annotation$uberonId[annotation$strain == strain])){
            message("UberonId info:", uberonId)
            for (sex in unique(annotation$sex[annotation$uberonId == uberonId])){
              message("sex info:", sex)

              infoLib <- annotation$libraryId[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex]
              infoLib2 <- as.data.frame(infoLib); colnames(infoLib2) <- "libraryId"
              message("Libraries retrieved : ", length(infoLib))

              if(length(infoLib) >= cellsThreshold){
                extractInfo <- merge(infoLib2, annotation, by="libraryId")
                information_file <- data.frame(extractInfo)
                information_file <- information_file %>% filter(!str_detect(information_file$libraryId, "^#"))
                write.table(information_file, file = passScRNASeqLibrary, col.names = FALSE , row.names = FALSE ,append = TRUE, quote = FALSE, sep = "\t")
              } else {
                extractInfo <- merge(infoLib2, annotation, by="libraryId")
                information_file <- data.frame(extractInfo)
                information_file <- information_file %>% filter(!str_detect(information_file$libraryId, "^#"))
                write.table(information_file, file = notPassScRNASeqLibrary, col.names = FALSE, row.names = FALSE , append = TRUE, quote = FALSE, sep = "\t")
              }
            }
          }
        }
      }
    }
  }
}

## samples that passed
info1 <- file.info(passScRNASeqLibrary)
if (info1$size == "0"){
  message("File is empty")
} else {
  newFile <- read.table(passScRNASeqLibrary, h=F, sep="\t", quote="\\")
  colnames(newFile) <- colnames(annotation)
  write.table(newFile, file = passScRNASeqLibrary, col.names =TRUE , row.names = FALSE, quote = FALSE, sep = "\t")
}
## samples that should be excluded
info2 <- file.info(outfile_NOT_PASS)
if (info2$size == "0"){
  message("File is empty")
} else {
  excludeSamplesFile <- read.table(notPassScRNASeqLibrary, h=FALSE, sep="\t", quote="\\")
  colnames(excludeSamplesFile) <- colnames(annotation)
  write.table(excludeSamplesFile, file = notPassScRNASeqLibrary, col.names =TRUE , row.names = FALSE, quote = FALSE, sep = "\t")
}
