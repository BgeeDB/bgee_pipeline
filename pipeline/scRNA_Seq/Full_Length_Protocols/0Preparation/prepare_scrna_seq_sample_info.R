## SFonsecaCosta, June 2019

## This script create the scrna_seq_sample_info file to run the scRNA-Seq pipeline.

## Usage:
## R CMD BATCH --no-save --no-restore '--args pass_annotationControl="passScRNASeqLibrary.tsv" raw_cells_folder="raw_cells_folder" output_folder="output_folder"' prepare_scrna_seq_sample_info.R prepare_scrna_seq_sample_info.Rout
## pass_annotationControl --> Use the file passScRNASeqLibrary.tsv (this means use just cell-types with more then 50 cells)
## raw_cells_folder --> Folder where is localized all raw libraries (this means all FASTQ.gz files per cell) of all experiments, species and cell-types
## output_folder --> Folder where the scrna_seq_sample_info.tsv should be saved

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
command_arg <- c("pass_annotationControl", "metadata_file", "raw_cells_folder", "output_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read annotation file. If file not exists, script stops
if( file.exists(pass_annotationControl) ){
  annotation <- read.table(pass_annotationControl, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
} else {
  stop( paste("pass_annotationControl file not found [", pass_annotationControl, "]\n"))
}

if( file.exists(metadata_file) ){
   metadata <- read.table(metadata_file, h=T, sep="\t", comment.char="")
} else {
  stop("metadata file not found [", metadata_file, "]")
}


##################################################################### FUNCTION #############################################################################################
## run fastp and collect information
collectInformationFASTP <- function(raw_cells_folder, annotation, library){
  species <- annotation$speciesId[annotation$libraryId == library]
  fastp_dir <- file.path(raw_cells_folder, species, library)
  rawFiles <- (list.files(path = fastp_dir, pattern = "*.gz"))
  nameRaw <- sub("\\..*", "", rawFiles)

  ## check if fastp has already run
  fastpJSON <- list.files(path=fastp_dir, pattern = "*.fastp.json.xz")

  if (isTRUE(file.exists(file.path(fastp_dir, fastpJSON)))){
    cat("The fastpJSON file already exist for this library ", library, "\n")
    libraryType <- ifelse(length(rawFiles) == 1, "SINGLE", "PAIRED")
    readJsonOutput <- fromJSON(file = file.path(fastp_dir, list.files(path = fastp_dir, pattern = "*.fastp.json.xz")))
    readLength <- readJsonOutput$read1_before_filtering$total_cycles
  } else {
    cat("Need to run fastp", "\n")
    if (length(rawFiles) == 1){
      cat("The library, ", library ," is SINGLE-end ", "\n")
      system(sprintf('%s -i %s -h %s -j %s', paste0("fastp"), file.path(fastp_dir, rawFiles), file.path(fastp_dir, paste0(nameRaw, ".fastp.html")), file.path(fastp_dir, paste0(nameRaw, ".fastp.json"))))
      libraryType <- "SINGLE"
      ## collect readLength
      readJsonOutput <- fromJSON(file = file.path(fastp_dir, list.files(path = fastp_dir, pattern = "*.fastp.json$")))
      readLength <- readJsonOutput$read1_before_filtering$total_cycles
      system(sprintf('%s %s %s %s', paste0("xz"), paste0("-9"), paste0(nameRaw, ".fastp.html"), paste0(nameRaw, ".fastp.json")))
    } else {
      cat("The library, ", library ," is PAIRED-end ", "\n")
      read1 <- rawFiles[1]
      read2 <- rawFiles[2]
      nameFile <-  nameRaw <- sub("\\_.*", "", rawFiles)
      system(sprintf('%s -i %s -I %s -h %s -j %s', paste0("fastp"), file.path(fastp_dir, read1), file.path(fastp_dir, read2), file.path(fastp_dir, paste0(unique(nameFile), ".fastp.html")),
          file.path(fastp_dir, paste0(unique(nameFile), ".fastp.json"))))
      libraryType <- "PAIRED"
      ## collect readLength
      readJsonOutput <- fromJSON(file = file.path(fastp_dir, list.files(path = fastp_dir, pattern = "*.fastp.json$")))
      readLength <- readJsonOutput$read1_before_filtering$total_cycles
      system(sprintf('%s %s %s %s', paste0("xz"), paste0("-9"), file.path(fastp_dir, paste0(nameRaw, ".fastp.html")), file.path(fastp_dir, paste0(nameRaw, ".fastp.json"))))
    }
  }

  ## collect information per library
  information_file <- data.frame(library, libraryType, readLength)
  return(information_file)
  }


## Function to add species name
addSpeciesName <- function(scrna_seq_sample_info, metadata){
  scrna_seq_sample_info$organism <- metadata$scientific_name[metadata$tax_id == scrna_seq_sample_info$speciesId]
  return(scrna_seq_sample_info)
}

##################################### OUTPUT #############################################################################################
## create a intermediary file to collect information about each library, or truncate already existing one....
#TODO: Should not create a tmp file but directly collect data in memory
tmpInfoFile <- paste0(output_folder, ".tmp")
file.create(tmpInfoFile)
cat("libraryId\tlibraryType\treadLength\n",file = tmpInfoFile, sep = "\t")

###################################### RUN PER lIBRARY ##################################################################################
for (libraryID in unique(annotation$libraryId)) {

  fastpInfo <- collectInformationFASTP(raw_cells_folder = raw_cells_folder, annotation = annotation, library = libraryID)

  write.table(fastpInfo, file = tmpInfoFile, row.names = FALSE , col.names = FALSE , append = TRUE, quote = FALSE, sep = "\t")
}

### read final output file --> InfoFile
fileInfo <- read.table(tmpInfoFile, header=TRUE, sep="\t")

## Create the scrna_seq_sample_info
scrna_seq_sample_info <- merge(annotation, fileInfo, by = "libraryId", incomparables = NaN)
scrna_seq_sample_info <- scrna_seq_sample_info %>% dplyr::select("libraryId", "experimentId", "cellTypeName_abInitio", "cellTypeId_abInitio", "speciesId",
                                                          "platform", "protocol", "protocolType", "libraryType", "infoOrgan", "stageId",
                                                          "anatId", "sex", "strain", "readLength", "speciesId", "genotype")
scrna_seq_sample_info$organism <- "NaN"
finalTable <- addSpeciesName(scrna_seq_sample_info = scrna_seq_sample_info, metadata = metadata)
message(finalTable)
write.table(finalTable,file = file.path( output_folder, "scrna_seq_sample_info.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
## remove intermediary file
file.remove(tmpInfoFile)
