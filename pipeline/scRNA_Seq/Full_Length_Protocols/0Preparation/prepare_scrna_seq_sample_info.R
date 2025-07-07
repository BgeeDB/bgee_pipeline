## SFonsecaCosta, June 2019

## This script create the scrna_seq_sample_info file to run the scRNA-Seq pipeline.

## Usage:
## R CMD BATCH --no-save --no-restore '--args filteredLibraryAnnotation="passScRNASeqLibrary.tsv" fastqDir="fastqDir" outputDir="outputDir"' prepare_scrna_seq_sample_info.R prepare_scrna_seq_sample_info.Rout
## filteredLibraryAnnotation --> list of raw libraries nnotation filtered to contain only full-length scRNA-Seq libraries
## fastqDir                  --> directory containing the fastq files
## metadataFile              --> file with the metadata information
## outputDir                 --> directory where the scrna_seq_sample_info.tsv should be saved

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
command_arg <- c("filteredLibraryAnnotation", "metadataFile", "fastqDir", "outputDir")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read annotation file. If file not exists, script stops
if( file.exists(filteredLibraryAnnotation) ){
  annotation <- read.table(filteredLibraryAnnotation, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
} else {
  stop( paste("filteredLibraryAnnotation file not found [", filteredLibraryAnnotation, "]\n"))
}

if( file.exists(metadataFile) ){
   metadata <- read.table(metadataFile, h=T, sep="\t", comment.char="")
} else {
  stop("metadata file not found [", metadataFile, "]")
}
#keep only library type and the library ID from the metadata file
metadata <- metadata %>% dplyr::select("library_id", "library_layout", "scientific_name")
colnames(metadata) <- c("libraryId", "libraryType", "scientific_name")

# then merge those info with the annotation file
scrna_seq_sample_info <- merge(annotation, metadata, by = "libraryId", incomparables = NaN)

##################################################################### FUNCTION #############################################################################################
## retrieve information from fastp
collectReadLengthFASTP <- function(fastq_dir, speciesId, library){
  fastq_library_dir <- file.path(fastq_dir, speciesId, library)
  if (!dir.exists(fastq_library_dir)) {
    cat(paste0("FASTQ files for library ", library, " are not available. The library will not be added to the sample info file.\n"))
    return(NA)
  }
  fastq_files <- (list.files(path = fastq_library_dir, pattern = "*.fastq.gz"))

  ## check if fastp has already run
  fastpJSON <- list.files(path=fastq_library_dir, pattern = "*.fastp.json.xz")
  if (length(fastpJSON) == 0) {
    stop("No fastp JSON for library ", library, ".")
    return(NA)
  }
  json_output <- fromJSON(file = file.path(fastq_library_dir, fastpJSON))
  read_length <- json_output$read1_before_filtering$total_cycles
  return(read_length)
  ## collect information per library
}

# Collect fastp info in memory
scrna_seq_sample_info$readLength <- NA
for (library_id in unique(annotation$libraryId)) {
  if (grepl("^#", library_id)) next  # Skip commented library_ids
  species_id <- annotation$speciesId[annotation$libraryId == library_id]
  scrna_seq_sample_info$readLength[scrna_seq_sample_info$libraryId == library_id] <- collectReadLengthFASTP(fastq_dir = fastqDir, speciesId = species_id, library = library_id)
}

## Create the scrna_seq_sample_info
scrna_seq_sample_info <- scrna_seq_sample_info %>% dplyr::select("libraryId", "experimentId", "cellTypeName", "cellTypeId", "speciesId","platform", "protocol", "protocolType", "libraryType", "infoOrgan", "stageId", "anatId", "sex", "strain", "readLength", "genotype", "scientific_name")

scrna_seq_sample_info$organism <- gsub(" ", "_", scrna_seq_sample_info$scientific_name)
scrna_seq_sample_info$scientific_name <- NULL
write.table(scrna_seq_sample_info, file = file.path( outputDir, "scrna_seq_sample_info.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

