## SFonsecaCosta 2019
## This script allow to verify the number of cells per cell-population regarding the: experimentID, species, cellTypeId, stageId, strain, uberonId and sex after the annotation process.
## Just cell-population that belongs to the same experimentID, species, cellTypeId, stageId, strain, uberonId and sex and have at least 50 cells will be keeped to continue in the pipeline.
## The output file generated (passScRNASeqLibrary.tsv) will be used to download the data and to continue in the pipeline.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeqLibrary="scRNASeqLibrary.tsv" cellsThreshold="minimum number of cells" output_folder="output_folder"' pre_process_control_annotation.R pre_process_control_annotation.Rout
## scRNASeqLibrary --> File from manual annotation (scRNASeqFLLibrary.tsv)
## cellsThreshold --> minimum number of cells in a library
## output_folder --> folder where the generated file should be saved

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
command_arg <- c("scRNASeqLibrary", "cellsThreshold", "output_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read scRNASeqLibrary file, if file not exists, script stops
if( file.exists(scRNASeqLibrary) ){
  annotation <- read.table(scRNASeqLibrary, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
  ## just replace NA in sex and strain per "missing information" to be a level in factor (levels())
  annotation$sex <- forcats::fct_explicit_na(annotation$sex)
  annotation$strain <- forcats::fct_explicit_na(annotation$strain)
} else {
  stop( paste("scRNASeqLibrary file not found [", scRNASeqLibrary, "]\n"))
}

#################################################################################
## create new output files
pass_Libraries <- file.path(output_folder, "passScRNASeqLibrary.tsv")
if (file.exists(pass_Libraries)){
  message("File already exists and will be removed to create a new one to avoid overwritting!")
  file.remove(pass_Libraries)
  file.create(pass_Libraries)
} else {
  file.create(pass_Libraries)
}
notpass_Libraries <- file.path(output_folder,"notPassScRNASeqLibrary.tsv")
if (file.exists(notpass_Libraries)){
  message("File already exists and will be removed to create a new one to avoid overwritting!")
  file.remove(notpass_Libraries)
  file.create(notpass_Libraries)
} else {
  file.create(notpass_Libraries)
}

#TODO: maybe collect all pass/not_pass info in a DF and write the 2 files after (directly with the proper colnames)
for (species in unique(annotation$speciesId)) {
  message("Species:", species)
  for (experiment in unique(annotation$experimentId[annotation$speciesId == species])){
    message("Name of experiments that belongs to the species:", experiment)
    for (cellId in unique(annotation$cellTypeId_abInitio[annotation$experimentId == experiment])){
      message("CellId info:", cellId)
      for (stageId in unique(annotation$stageId[annotation$cellTypeId_abInitio == cellId])){
        message("StageId info:", stageId)
        for (strain in unique(annotation$strain[annotation$stageId == stageId])){
          message("Strain info:", strain)
          for (uberonId in unique(annotation$anatId[annotation$strain == strain])){
            message("UberonId info:", uberonId)
            for (sex in unique(annotation$sex[annotation$anatId == uberonId])){
              message("sex info:", sex)

              infoLib <- data.frame(annotation$libraryId[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId_abInitio == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$anatId == uberonId & annotation$sex == sex])
              colnames(infoLib) <- "libraryId"
              message("Libraries retrieved :", nrow(infoLib))
              
              if (nrow(infoLib) == 0){
                message("No Libraries retrieved!")
              } else if (nrow(infoLib) >= cellsThreshold){
                extractInfo <- data.frame(merge(infoLib, annotation, by="libraryId"))
                information_file <- extractInfo %>% filter(!str_detect(extractInfo$libraryId, "^#"))
                write.table(extractInfo, file = pass_Libraries, col.names = TRUE , row.names = FALSE ,append = TRUE, quote = FALSE, sep = "\t")
              } else {
                extractInfo <- data.frame(merge(infoLib, annotation, by="libraryId"))
                information_file <- extractInfo %>% filter(!str_detect(extractInfo$libraryId, "^#"))
                write.table(information_file, file = notpass_Libraries, col.names = TRUE, row.names = FALSE , append = TRUE, quote = FALSE, sep = "\t")
              }
            }
          }
        }
      }
    }
  }
}
