## SFonsecaCosta, Sep 17 2019

## This script is used to retrieve the metadata for the target based
## protocols from SRA source. And then to compare the annotation information (speciesId and protocol)
## for each library with metadata.

##XXX Is this script really useful??? Is it possible that libraries were not annotated to proper species? What does it mean? Bgee curation should be trusted

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeqExperiment="scRNASeqExperiment.tsv" scRNASeqTBLibrary="scRNASeqTBLibrary.tsv" output_folder="output_folder"' retrieve_metadata.R retrieve_metadata.Rout
## scRNASeqExperiment --> File with information about all experiments annotated
## scRNASeqTBLibrary --> File with all libraries annotated by bgee
## output_folder --> Folder where the output files should be saved

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
command_arg <- c("scRNASeqExperiment","scRNASeqTBLibrary", "output_file")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read experiments and library annotation files. Do not consider If file not exists, script stops
if( file.exists(scRNASeqExperiment) ){
  experiments <- read.table(scRNASeqExperiment, header=TRUE, sep="\t")
  colnames(experiments)[1]<-"experimentId"
  experiments <- experiments %>% filter(!str_detect(experiments$experimentId, "^#"))
  
} else {
  stop( paste("The experiment file not found [", scRNASeqExperiment, "]\n"))
}

if( file.exists(scRNASeqTBLibrary) ){
  ##XXX could be necessary to add `quote = ""` depending on the final format of the file
  annotation <- read.table(scRNASeqTBLibrary, header=TRUE, sep="\t", comment.char = "")
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
                   "&fields=experiment_accession,run_accession,",
                   "read_count,tax_id,scientific_name,",
                   "instrument_model,library_layout,fastq_ftp,submitted_ftp,",
                   "&download=TRUE",
                   sep="")
  readInfo <- read.table(url(ena.url), header=TRUE, sep="\t")
  return(readInfo)
}

# define header of SRA metadata files. Used to reorder columns before writing
metadata_file_header <- c("sample_accession","experiment_id", "library_id", "run_accession",
                          "read_count", "tax_id", "scientific_name", "instrument_model",
                          "library_layout", "fastq_ftp", "submitted_ftp", "source")

## select 10X Genomics and 3' end target-based protocols to retrieve metadata
## We use a grep to detect 10X Genomics as the protocol name can contain the version (e.g 10X Genomics V3)
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
  if(sourceId == "SRA"){
    for (libraryId in selected_libraries$libraryId[selected_libraries$experimentId == expId]) {
      library <- as.data.frame(selected_libraries[selected_libraries$libraryId == libraryId,])
      ## retrieve SRA metadata
      extractSRA <- SRA_metadata(libraryID = libraryId)
      extractSRA$source <- "SRA"
      ## compare with Bgee annotation
      ## probably too stringent as it does not use a controlled vocabulary
      if (!identical(as.character(library$platform),as.character(unique(extractSRA$instrument_model)))) {
        warning("Mismatch platform for library ", library$libraryId, " expected", library$platform,
                " but was ", unique(extractSRA$instrument_model))
        metadata_with_mismatch <- rbind(metadata_with_mismatch, extractSRA)
      } else if (!identical(as.character(library$speciesId),as.character(unique(extractSRA$tax_id)))) {
        warning("Mismatch species for library ", library$libraryId, " expected", library$speciesId,
                " but was ", unique(extractSRA$tax_id))
        metadata_with_mismatch <- rbind(metadata_with_mismatch, extractSRA)
      } else {
        metadata <- rbind(metadata,extractSRA)
      }
    }
  } else if (sourceId == "HCA") {
    ## add all HCA annotation because it is not possible to check the
    ## annotation since HCAExplore R package was deprecated
    message("HCA data! for experiment ", expId,
            ". Do not retrieve metadata!")
    hca_libraries <- annotation[annotation$experimentId == expId,]
    hca_metadata <- data.frame(NA, NA, hca_libraries$libraryId, NA, NA, hca_libraries$speciesId,
                               "Homo sapiens", hca_libraries$platform, NA, NA, NA, "HCA")
    # merge lines of the two df. Do not use rbind as column names are different 
    metadata <- as.data.frame(mapply(c, metadata, hca_metadata))
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
  metadata <- metadata[, c(metadata_file_header)]
}
# write file with metadata from SRA
write.table(metadata, file = output_file, quote = FALSE, sep = "\t", col.names = TRUE,
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
metadata_notmatch_file <- file.path(dirname(output_file),"metadata_notMatch_10X.tsv")
write.table(metadata_with_mismatch, file = metadata_notmatch_file, quote = FALSE, sep = "\t", 
  col.names = TRUE, row.names = FALSE)
