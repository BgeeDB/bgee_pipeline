## SFonsecaCosta, Sep 17 2019

## This script is used to retrieve the metadata for the target based protocols from SRA source.
## And then to compare the annotation information for each library with metadata.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeqExperiment="scRNASeqExperiment.tsv" scRNASeqTBLibrary="scRNASeqTBLibrary.tsv" output_folder="output_folder"' retrieve_metadata.R retrieve_metadata.Rout
## scRNASeqExperiment --> File with information about all experiments annotated
## scRNASeqTBLibrary --> File with all libraries annotated by bgee
## output_folder --> Folder where the output files should be saved

## libraries used
library(data.table)
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
command_arg <- c("scRNASeqExperiment","scRNASeqTBLibrary", "output_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read experiments and library annotation files. If file not exists, script stops
if( file.exists(scRNASeqExperiment) ){
  experiments <- fread(scRNASeqExperiment, h=T, sep="\t")
  colnames(experiments)[1]<-"experimentId"
} else {
  stop( paste("The experiment file not found [", scRNASeqExperiment, "]\n"))
}
if( file.exists(scRNASeqTBLibrary) ){
  annotation <- fread(scRNASeqTBLibrary, h=T, sep="\t")
  colnames(annotation)[1]<-"libraryId"
  annotation <- annotation %>% filter(!str_detect(annotation$libraryId, "^#"))
} else {
  stop( paste("The library file not found [", scRNASeqTBLibrary, "]\n"))
}
##########################################################################################################################
## function to download metadata from SRA
SRA_metadata <- function(libraryID){

  PID <- paste0(libraryID)
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

## generate output files 
metadata_file <- file.path(output_folder,"metadata_info_10X.tsv")
if (file.exists(metadata_file)){
  message("File already exists and will be removed to create a new one to avoid overwritting!")
  file.remove(metadata_file)
  file.create(metadata_file)
  cat("sample_accession\texperiment_accession\trun_accession\tread_count\ttax_id\tscientific_name\tinstrument_model\tlibrary_layout\tfastq_ftp\tsubmitted_ftp\n",file = metadata_file, sep = "\t")
} else {
  file.create(metadata_file)
  cat("sample_accession\texperiment_accession\trun_accession\tread_count\ttax_id\tscientific_name\tinstrument_model\tlibrary_layout\tfastq_ftp\tsubmitted_ftp\n",file = metadata_file, sep = "\t")
}
metadata_notmatch_file <- file.path(output_folder,"metadata_notMatch_10X.tsv")
if (file.exists(metadata_notmatch_file)){
  message("File already exists and will be removed to create a new one to avoid overwritting!")
  file.remove(metadata_notmatch_file)
  file.create(metadata_notmatch_file)
  cat("sample_accession\texperiment_accession\trun_accession\tread_count\ttax_id\tscientific_name\tinstrument_model\tlibrary_layout\tfastq_ftp\tsubmitted_ftp\n",file = metadata_notmatch_file, sep = "\t")
} else {
  file.create(metadata_notmatch_file)
  cat("sample_accession\texperiment_accession\trun_accession\tread_count\ttax_id\tscientific_name\tinstrument_model\tlibrary_layout\tfastq_ftp\tsubmitted_ftp\n",file = metadata_notmatch_file, sep = "\t")
}

## select just 10X target-based protocols to retrieve metadata
targetBased <- experiments[experiments$protocol %like% "10X Genomics", ]
targetBased_libraries <- dplyr::filter(annotation, protocolType == "3'end" & protocol == "10X Genomics")
targetBased_libraries <- targetBased_libraries[!grepl("#", targetBased_libraries$libraryId),]
experimentsIDs <- unique(targetBased_libraries$experimentId)
targetBased <- dplyr::filter(targetBased, experimentId %in% experimentsIDs)

## extract metadata from SRA
metadata <- c()
for (experiment in unique(targetBased$experimentId)) {
  ## select source  of the experiment
  sourceID <- as.character(unique(targetBased$experimentSource[targetBased$experimentId == experiment]))
  libraries_from_Exp <- annotation$libraryId[annotation$experimentId == experiment & annotation$protocol == "10X Genomics" & annotation$protocolType == "3'end"]
  
  if(sourceID == "SRA"){
    for (i in libraries_from_Exp) {
      extractSRA <- SRA_metadata(libraryID = i)
      metadata <- rbind(metadata,extractSRA)
    }
  } else if (sourceID == "HCA") {
    message("HCA data! Not retrieve metadata!")
  } else {
    message("Source ERROR!")
  }
}
## add HCA by default to metadata file directly from annotation, and not making comparison between the annotation and metadata. This is because HCAExplore R package was deprecated
hca_experiment <- experiments$experimentId[experiments$experimentSource == "HCA" ]
getLib <- annotation[annotation$experimentId == hca_experiment]
getLib <- data.frame(NA, getLib$libraryId, NA, NA, getLib$speciesId, "Homo sapiens", getLib$platform, NA, NA, NA)
colnames(getLib) <- colnames(metadata)
metadata <- rbind(metadata, getLib)
  
## compare information from annotation and metadata
## SRA: libraryID, plataform and speciesID
for(i in unique(targetBased_libraries$libraryId)) {

  annotationInfo <- targetBased_libraries[targetBased_libraries$libraryId == i,]
  metadataInfo <- metadata[metadata$experiment_accession == i,]

  compare_library <- identical(as.character(annotationInfo[['libraryId']]),as.character(metadataInfo[['experiment_accession']]))
  compare_machine <- identical(as.character(annotationInfo[['platform']]),as.character(metadataInfo[['instrument_model']]))
  compare_speciesID <- identical(as.character(annotationInfo[['speciesId']]),as.character(metadataInfo[['tax_id']]))
  
  if (compare_library == "TRUE" && compare_machine == "TRUE" && compare_speciesID == "TRUE"){
      message(as.character(annotationInfo$libraryId[1]), " complete match between annotation and metadata")

      ## export libraries that pass and will be downloaded
      write.table(metadataInfo, file = metadata_file, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)

    } else {
      message("For the library: ", as.character(annotationInfo$libraryId[1]), " the comparison library is: ", compare_library)
      message("For the library: ", as.character(annotationInfo$libraryId[1]), " the comparison platform is: ", compare_machine)
      message("For the library: ", as.character(annotationInfo$libraryId[1]), " the comparison species is: ", compare_speciesID)

      ## export libraries that will not be used to download
      write.table(metadataInfo, file = metadata_notmatch_file, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
    }
}
