## SFonsecaCosta, June 2019

## This script create the scrna_seq_sample_info file to run the scRNA-Seq pipeline.

## Usage:
## R CMD BATCH --no-save --no-restore '--args NEW_scRNASeqLibrary="NEW_scRNASeqLibrary.tsv" raw_cells_folder="raw_cells_folder" output_folder="output_folder"' prepare_scrna_seq_sample_info.R prepare_scrna_seq_sample_info.Rout
## NEW_scRNASeqLibrary --> Use the file NEW_scRNASeqLibrary (this means use just cell-types with more then 50 cells)
## raw_cells_folder --> Folder where is localized all raw libraries (this means all FASTQ.gz files per cell) of all experiments, species and cell-types
## output_folder --> Folder where should be saved the scrna_seq_sample_info.txt file

## Libraries
library(rjson)
library(plyr)
library(dplyr)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("NEW_scRNASeqLibrary", "raw_cells_folder", "output_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read annotation file. If file not exists, script stops
if( file.exists(NEW_scRNASeqLibrary) ){
  annotation <- read.table(NEW_scRNASeqLibrary, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
} else {
  stop( paste("NEW_scRNASeqLibrary file not found [", NEW_scRNASeqLibrary, "]\n"))
}

##################################################################### FUNCTION #############################################################################################
## run fastp and collect information
collectInformationFASTP <- function(raw_cells_folder, library){
  setwd(file.path(raw_cells_folder, library))
  rawFiles <- (list.files(path=file.path(raw_cells_folder, library), pattern = "*.gz"))
  nameRaw <- sub("\\..*", "", rawFiles)

  ## check if fastp has already run
  fastpJSON <- list.files(path=file.path(raw_cells_folder, library), pattern = "*.fastp.json.xz")

  if (isTRUE(file.exists(fastpJSON))){
    cat("The fastpJSON file already exist for this library ", library, "\n")
    libraryType <- ifelse(length(rawFiles) == 1, "SINGLE", "PAIRED")
    readJsonOutput <- fromJSON(file = (list.files(path=file.path(raw_cells_folder, library), pattern = "*.fastp.json.xz")))
    readLength <- readJsonOutput$summary$before_filtering$read1_mean_length
  } else {
    cat("Need to run fastp", "\n")
    if (length(rawFiles) == 1){
      cat("The library, ", library ," is SINGLE-end ", "\n")
      system(sprintf('%s -i %s -h %s -j %s', paste0("fastp"), file.path(raw_cells_folder, library, rawFiles), paste0(nameRaw, ".fastp.html"), paste0(nameRaw, ".fastp.json")))
      libraryType <- "SINGLE"
      ## collect readLength
      readJsonOutput <- fromJSON(file = (list.files(path=file.path(raw_cells_folder, library), pattern = "*.fastp.json$")))
      readLength <- readJsonOutput$summary$before_filtering$read1_mean_length
      system(sprintf('%s %s %s %s', paste0("xz"), paste0("-9"), paste0(nameRaw, ".fastp.html"), paste0(nameRaw, ".fastp.json")))
    } else {
      cat("The library, ", library ," is PAIRED-end ", "\n")
      read1 <- rawFiles[1]
      read2 <- rawFiles[2]
      nameFile <-  nameRaw <- sub("\\_.*", "", rawFiles)
      system(sprintf('%s -i %s -I %s -h %s -j %s', paste0("fastp"), file.path(raw_cells_folder, library, read1), file.path(raw_cells_folder, library, read2), paste0(unique(nameFile), ".fastp.html"), paste0(unique(nameFile), ".fastp.json")))
      libraryType <- "PAIRED"
      ## collect readLength
      readJsonOutput <- fromJSON(file = (list.files(path=file.path(raw_cells_folder, library), pattern = "*.fastp.json$")))
      readLength <- readJsonOutput$summary$before_filtering$read1_mean_length
      system(sprintf('%s %s %s %s', paste0("xz"), paste0("-9"), paste0(nameRaw, ".fastp.html"), paste0(nameRaw, ".fastp.json")))
    }
  }

  ## collect information per library
  information_file <- data.frame(library, libraryType, readLength)
  return(information_file)
  }

## Function to add species name
speciesName <- function(scrna_seq_sample_info){

  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 6239] <- "Caenorhabditis_elegans"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 7217] <- "Drosophila_ananassae"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 7227] <- "Drosophila_melanogaster"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 7230] <- "Drosophila_mojavensis"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 7237] <- "Drosophila_pseudoobscura"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 7240] <- "Drosophila_simulans"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 7244] <- "Drosophila_virilis"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 7245] <- "Drosophila_yakuba"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 7955] <- "Danio_rerio"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 8364] <- "Xenopus_tropicalis"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9031] <- "Gallus_gallus"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9258] <- "Ornithorhynchus_anatinus"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9365] <- "Erinaceus_europaeus"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9544] <- "Macaca_mulatta"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9593] <- "Gorilla_gorilla"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9597] <- "Pan_paniscus"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9598] <- "Pan_troglodytes"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9606] <- "Homo_sapiens"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9615] <- "Canis_lupus_familiaris"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9685] <- "Felis_catus"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9796] <- "Equus_caballus"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9823] <- "Sus_scrofa"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9913] <- "Bos_taurus"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 9986] <- "Oryctolagus_cuniculus"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 10090] <- "Mus_musculus"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 10116] <- "Rattus_norvegicus"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 10141] <- "Cavia_porcellus"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 13616] <- "Monodelphis_domestica"
  scrna_seq_sample_info$organism[scrna_seq_sample_info$speciesId == 28377] <- "Anolis_carolinensis"

  return(scrna_seq_sample_info)

}

##################################### OUTPUT #############################################################################################
## create an output file to collect information about each library....
InfoFile <- file.path(output_folder, "InfoFile.tsv")
if (!file.exists(InfoFile)){
  file.create(InfoFile)
  cat("libraryId\tlibraryType\treadLength\n",file = file.path(output_folder, "InfoFile.tsv"), sep = "\t")

} else {
  print("File already exist.....")
}

###################################### RUN PER lIBRARY ##################################################################################
for (libraryID in unique(annotation$libraryId)) {

  fastpInfo <- collectInformationFASTP(raw_cells_folder = raw_cells_folder, library = libraryID)

  write.table(fastpInfo, file = InfoFile, row.names = FALSE , col.names = FALSE , append = TRUE, quote = FALSE, sep = "\t")
}

### read final output file --> InfoFile
fileInfo <- file.path(output_folder,"InfoFile.tsv")
fileInfo <- read.table(fileInfo, header=TRUE, sep="\t")

## Create the scrna_seq_sample_info
scrna_seq_sample_info <- merge(annotation, fileInfo, by = "libraryId", incomparables = NaN)
scrna_seq_sample_info <- scrna_seq_sample_info %>% select("libraryId", "experimentId", "cellTypeName", "cellTypeId", "speciesId",
                                                          "platform", "protocol", "protocolType", "libraryType", "infoOrgan", "stageId",
                                                          "uberonId", "sex", "strain", "readLength")
scrna_seq_sample_info$organism <- "NaN"
finalTable <- speciesName(scrna_seq_sample_info = scrna_seq_sample_info)
write.table(finalTable,file = file.path(output_folder, "scrna_seq_sample_info.txt"),quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
## remove intermediary file
file.remove(file.path(output_folder, "InfoFile.tsv"))
