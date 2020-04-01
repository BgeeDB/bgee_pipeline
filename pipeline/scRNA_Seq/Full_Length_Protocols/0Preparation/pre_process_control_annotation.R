## SFonsecaCosta 2019
## This script allow to verify the number of cells per cell-population regarding the: experimentID, species, cellTypeId, stageId, strain, uberonId and sex after the annotation process.
## Just cell-population that belongs to the same experimentID, species, cellTypeId, stageId, strain, uberonId and sex and have at least 50 cells will be keeped to continue in the pipeline.
## The output file generated (NEW_scRNASeqLibrary.tsv) will be used to download the data and to continue in the pipeline.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeqLibrary="scRNASeqLibrary.tsv" output_folder="output_folder"' pre_process_control_annotation.R pre_process_control_annotation.Rout
## scRNASeqLibrary --> File from manual annotation
## output_folder --> Output folder where should be saved the new output file (after filtering cell types that not pass)

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
command_arg <- c("scRNASeqLibrary", "output_folder")
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
  pathAnnotation <- dirname(scRNASeqLibrary)
} else {
  stop( paste("scRNASeqLibrary file not found [", scRNASeqLibrary, "]\n"))
}

#################################################################################
## create output files
setwd(output_folder)
if (dir.exists(output_folder)){
  file.create("NEW_scRNASeqLibrary.tsv")
  file.create("NOT_PASS_scRNASeqLibrary.tsv")
  } else {
  print("Directoty not exists.....")
  }

infile <- file.path(pathAnnotation, "scRNASeqLibrary.tsv")
outfile_PASS <- file.path(output_folder, "NEW_scRNASeqLibrary.tsv")
outfile_NOT_PASS <- file.path(output_folder, "NOT_PASS_scRNASeqLibrary.tsv")

for (species in unique(annotation$speciesId)) {
  cat("Species:", species, "\n")
  for (experiment in unique(annotation$experimentId[annotation$speciesId == species])){
    cat("Name of experiments that belongs to the species:", experiment, "\n")
    for (cellId in unique(annotation$cellTypeId[annotation$experimentId == experiment])){
      cat("CellId info:", cellId, "\n")
      for (stageId in unique(annotation$stageId[annotation$cellTypeId == cellId])){
        cat("StageId info:", stageId, "\n")
        for (strain in unique(annotation$strain[annotation$stageId == stageId])){
          cat("Strain info:", strain, "\n")
          for (uberonId in unique(annotation$uberonId[annotation$strain == strain])){
            cat("UberonId info:", uberonId, "\n")
            for (sex in unique(annotation$sex[annotation$uberonId == uberonId])){
              cat("sex info:", sex, "\n")

        infoLib <- annotation$libraryId[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex]
        infoLib2 <- as.data.frame(infoLib); colnames(infoLib2) <- "libraryId"
        cat("Libraries retrieved : ", length(infoLib), "\n")

        if(length(infoLib) >= 50){
          extractInfo <- merge(infoLib2, annotation, by="libraryId")
          information_file <- data.frame(extractInfo)
          information_file <- information_file %>% filter(!str_detect(information_file$libraryId, "^#"))
          write.table(information_file, file = outfile_PASS, col.names = FALSE , row.names = FALSE ,append = TRUE, quote = FALSE, sep = "\t")
        } else {
          extractInfo <- merge(infoLib2, annotation, by="libraryId")
          information_file <- data.frame(extractInfo)
          information_file <- information_file %>% filter(!str_detect(information_file$libraryId, "^#"))
          write.table(information_file, file = outfile_NOT_PASS, col.names = FALSE, row.names = FALSE , append = TRUE, quote = FALSE, sep = "\t")
        }
            }
          }
        }
      }
    }
  }
}

## samples that passed
info1 <- file.info(outfile_PASS)
if (info1$size == "0"){
  cat("File is empty", "\n")
} else {
  newFile <- read.table(outfile_PASS, h=F, sep="\t", quote="\\")
  colnames(newFile) <- colnames(annotation)
  write.table(newFile, file = file.path(output_folder, "NEW_scRNASeqLibrary.tsv"), col.names =T , row.names = F, quote = FALSE, sep = "\t")
}
## samples that should be excluded
info2 <- file.info(outfile_NOT_PASS)
if (info2$size == "0"){
  cat("File is empty", "\n")
} else {
  excludeSamplesFile <- read.table(outfile_NOT_PASS, h=F, sep="\t", quote="\\")
  colnames(excludeSamplesFile) <- colnames(annotation)
  write.table(excludeSamplesFile, file = file.path(output_folder, "NOT_PASS_scRNASeqLibrary.tsv"), col.names =T , row.names = F, quote = FALSE, sep = "\t")
}
