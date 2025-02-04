## SFonsecaCosta, Sep 17 2019

## This script is used to retrieve the metadata for the target based
## protocols from SRA, to compare the annotation information (speciesId and protocol). Those metadata are mandatory
## to download fastq files.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeqExperiment="scRNASeqExperiment.tsv" scRNASeqTBLibrary="scRNASeqTBLibrary.tsv" output_folder="output_folder"' retrieve_metadata.R retrieve_metadata.Rout
## scRNASeqExperiment --> File with information about all experiments annotated
## scRNASeqTBLibrary  --> File with all libraries annotated by bgee
## metadata_file      --> Path to the location where the metadata file will be saved
## information_file   --> Path to the location where the information file will be generated

## libraries used
library(stringr)
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
command_arg <- c("scRNASeqExperiment","scRNASeqTBLibrary", "metadata_file", "information_file")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read experiments and library annotation files. Do not consider If file not exists, script stops
if( file.exists(scRNASeqExperiment) ){
  experiments <- read.table(scRNASeqExperiment, header=TRUE, sep="\t", quote = "\"")
  colnames(experiments)[1]<-"experimentId"
  experiments <- experiments %>% filter(!str_detect(experiments$experimentId, "^#"))
  
} else {
  stop( paste("The experiment file not found [", scRNASeqExperiment, "]\n"))
}

if( file.exists(scRNASeqTBLibrary) ){
  ##XXX could be necessary to add `quote = ""` depending on the final format of the file
  annotation <- read.table(scRNASeqTBLibrary, header=TRUE, sep="\t", comment.char = "", quote = "\"")
  colnames(annotation)[1]<-"libraryId"
  annotation <- annotation %>% filter(!str_detect(annotation$libraryId, "^#"))
} else {
  stop( paste("The library file not found [", scRNASeqTBLibrary, "]\n"))
}

###################################################################################################
## function to download metadata from SRA
SRA_metadata <- function(libraryID){
  
  PID <- as.character(libraryID)
  ena.url <- paste("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
                   PID,
                   "&result=read_run",
                   "&fields=sample_accession,experiment_accession,run_accession,",
                   "read_count,tax_id,scientific_name,",
                   "instrument_model,library_layout,fastq_ftp,submitted_ftp,",
                   "&download=TRUE",
                   sep="")
  readInfo <- read.table(url(ena.url), header=TRUE, sep="\t")
  return(readInfo)
}

# define header of SRA metadata files. Used to reorder columns before writing
metadata_file_header <- c("sample_accession", "experiment_id", "library_id", "run_accession",
                          "read_count", "tax_id", "scientific_name", "instrument_model",
                          "library_layout", "fastq_ftp", "submitted_ftp", "source")

## select 10X Genomics and 3' end target-based protocols to retrieve metadata
## We use a grep to detect 10X Genomics as the protocol name can contain the version (e.g 10X Genomics V3)
##TODO to remove once the filtering based on protocols is done separatly
selected_libraries <- as.data.frame(dplyr::filter(annotation, protocolType == "3'end" &
  grepl("10X Genomics", protocol, perl=T)))
selected_experiments <- as.data.frame(dplyr::filter(experiments, experimentId %in%
  unique(selected_libraries$experimentId)))

## extract metadata from SRA
metadata <- c()
metadata_with_mismatch <- c()

for (expId in unique(selected_experiments$experimentId)) {
  ## select source  of the experiment
  sourceId <- as.character(unique(selected_experiments$experimentSource[selected_experiments$experimentId == expId]))
  if(sourceId == "SRA" || sourceId == "EBI" || sourceId == "HCA"){
    for (libraryId in selected_libraries$libraryId[selected_libraries$experimentId == expId]) {
      library <- as.data.frame(selected_libraries[selected_libraries$libraryId == libraryId,])
      ## retrieve SRA metadata
      extractSRA <- SRA_metadata(libraryID = libraryId)
      extractSRA$source <- sourceId
      ## compare with Bgee annotation
      ## probably too stringent as it does not use a controlled vocabulary
      if (!identical(as.character(library$platform),as.character(unique(extractSRA$instrument_model)))) {
        warning("Mismatch platform for library ", library$libraryId, " expected", library$platform,
                " but was ", unique(extractSRA$instrument_model))
        metadata <- rbind(metadata, extractSRA)
      } else if (!identical(as.character(library$speciesId),as.character(unique(extractSRA$tax_id)))) {
        warning("Mismatch species for library ", library$libraryId, " expected", library$speciesId,
                " but was ", unique(extractSRA$tax_id))
        metadata_with_mismatch <- rbind(metadata_with_mismatch, extractSRA)
      } else {
        metadata <- rbind(metadata,extractSRA)
      }
    }
  }  else {
    stop("Unknown Source ", sourceId)
  }
}

## update name and order of columns
if (!is.null(metadata)) {
  names(metadata)[names(metadata) == 'experiment_accession'] <- 'library_id'
  metadata <- merge(metadata, selected_libraries[, c("libraryId","experimentId")],
    by.x="library_id", by.y="libraryId")
  names(metadata)[names(metadata) == 'experimentId'] <- 'experiment_id'
  metadata <- metadata[, metadata_file_header]
}

information <- merge(annotation, metadata[, c("library_id", "run_accession", "read_count", "tax_id",
  "scientific_name", "library_layout")], by.x="libraryId", by.y="library_id")

# write file with metadata from SRA
write.table(metadata, file = metadata_file, quote = FALSE, sep = "\t", col.names = TRUE,
  row.names = FALSE)
# write file merging annotation and some metadata from SRA
write.table(information, file = information_file, quote = TRUE, sep = "\t", col.names = TRUE,
  row.names = FALSE)

# update name and order columns of mismatch metadata in case there was some
if (!is.null(metadata_with_mismatch)) {
  names(metadata_with_mismatch)[names(metadata_with_mismatch) == 'experiment_accession'] <-
    'library_id'
  metadata_with_mismatch <- merge(metadata_with_mismatch, selected_libraries[,
    c("libraryId","experimentId")], by.x="library_id", by.y="libraryId")
  names(metadata_with_mismatch)[names(metadata_with_mismatch) == 'experimentId'] <- 'experiment_id'
  metadata_with_mismatch <- metadata_with_mismatch[, c(metadata_file_header)]
}
# the name of this file is not exported to Makefile.common as it just a control file
metadata_notmatch_file <- file.path(dirname(metadata_file),"metadata_notMatch_10X.tsv")
write.table(metadata_with_mismatch, file = metadata_notmatch_file, quote = FALSE, sep = "\t", 
  col.names = TRUE, row.names = FALSE)
