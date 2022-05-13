## SFonsecaCosta, June 2019

## This script create the scrna_seq_sample_info file to run the scRNA-Seq pipeline.

## Usage:
## R CMD BATCH --no-save --no-restore '--args NEW_scRNASeqLibrary="NEW_scRNASeqLibrary.tsv" raw_cells_folder="raw_cells_folder" output_sample_info_file="output_info_file"' prepare_scrna_seq_sample_info.R prepare_scrna_seq_sample_info.Rout
## NEW_scRNASeqLibrary --> Use the file NEW_scRNASeqLibrary (this means use just cell-types with more then 50 cells)
## raw_cells_folder --> Folder where is localized all raw libraries (this means all FASTQ.gz files per cell) of all experiments, species and cell-types
## output_sample_info_file --> path to the scrna_seq_sample_info.txt file containing information for each library

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
command_arg <- c("NEW_scRNASeqLibrary", "raw_cells_folder", "output_sample_info_file")
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
  fastp_dir <- file.path(raw_cells_folder, library)
  rawFiles <- (list.files(path = fastp_dir, pattern = "*.gz"))
  nameRaw <- sub("\\..*", "", rawFiles)

  ## check if fastp has already run
  fastpJSON <- list.files(path=fastp_dir, pattern = "*.fastp.json.xz")

  if (isTRUE(file.exists(file.path(fastp_dir, fastpJSON)))){
    cat("The fastpJSON file already exist for this library ", library, "\n")
    libraryType <- ifelse(length(rawFiles) == 1, "SINGLE", "PAIRED")
    readJsonOutput <- fromJSON(file = file.path(fastp_dir, list.files(path = fastp_dir, pattern = "*.fastp.json.xz")))
    readLength <- readJsonOutput$summary$before_filtering$read1_mean_length
  } else {
    cat("Need to run fastp", "\n")
    if (length(rawFiles) == 1){
      cat("The library, ", library ," is SINGLE-end ", "\n")
      system(sprintf('%s -i %s -h %s -j %s', paste0("fastp"), file.path(fastp_dir, rawFiles), paste0(nameRaw, ".fastp.html"), paste0(nameRaw, ".fastp.json")))
      libraryType <- "SINGLE"
      ## collect readLength
      readJsonOutput <- fromJSON(file = file.path(fastp_dir, list.files(path = fastp_dir, pattern = "*.fastp.json$")))
      readLength <- readJsonOutput$summary$before_filtering$read1_mean_length
      system(sprintf('%s %s %s %s', paste0("xz"), paste0("-9"), paste0(nameRaw, ".fastp.html"), paste0(nameRaw, ".fastp.json")))
    } else {
      cat("The library, ", library ," is PAIRED-end ", "\n")
      read1 <- rawFiles[1]
      read2 <- rawFiles[2]
      nameFile <-  nameRaw <- sub("\\_.*", "", rawFiles)
      system(sprintf('%s -i %s -I %s -h %s -j %s', paste0("fastp"), file.path(fastp_dir, read1), file.path(fastp_dir, read2), paste0(unique(nameFile), ".fastp.html"), paste0(unique(nameFile), ".fastp.json")))
      libraryType <- "PAIRED"
      ## collect readLength
      readJsonOutput <- fromJSON(file = file.path(fastp_dir, list.files(path = fatp_dir, pattern = "*.fastp.json$")))
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
## create a intermediary file to collect information about each library, or truncate already existing one....
#TODO: Should not create a tmp file but directly collect data in memory
tmpInfoFile <- paste0(output_sample_info_file, ".tmp")
file.create(tmpInfoFile)
cat("libraryId\tlibraryType\treadLength\n",file = tmpInfoFile, sep = "\t")

###################################### RUN PER lIBRARY ##################################################################################
for (libraryID in unique(annotation$libraryId)) {

  fastpInfo <- collectInformationFASTP(raw_cells_folder = raw_cells_folder, library = libraryID)

  write.table(fastpInfo, file = tmpInfoFile, row.names = FALSE , col.names = FALSE , append = TRUE, quote = FALSE, sep = "\t")
}

### read final output file --> InfoFile
fileInfo <- read.table(tmpInfoFile, header=TRUE, sep="\t")

## Create the scrna_seq_sample_info
scrna_seq_sample_info <- merge(annotation, fileInfo, by = "libraryId", incomparables = NaN)

scrna_seq_sample_info <- scrna_seq_sample_info %>% select("libraryId", "experimentId", "cellTypeName", "cellTypeId", "speciesId",
                                                          "platform", "protocol", "protocolType", "libraryType", "infoOrgan", "stageId",
                                                          "uberonId", "sex", "strain", "readLength")
scrna_seq_sample_info$organism <- "NaN"
#TODO: move "(missing)" values for columns sex and strains to "NA"
finalTable <- speciesName(scrna_seq_sample_info = scrna_seq_sample_info)
write.table(finalTable,file = output_sample_info_file, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
## remove intermediary file
file.remove(tmpInfoFile)
