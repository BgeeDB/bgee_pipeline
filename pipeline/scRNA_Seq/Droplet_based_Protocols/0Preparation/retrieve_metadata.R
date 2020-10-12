## SFonsecaCosta, Sep 17 2019

## This script is used to retrieve the metadata for the target based protocols from sources as SRA and HCA.
## And then compare the annotation information for each library with metadata specially with SRA source.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeqExperiment="scRNASeqExperiment.tsv" scRNASeqLibrary="scRNASeqLibrary.tsv" output="output_folder"' retrieve_metadata.R retrieve_metadata.Rout
## scRNASeqExperiment --> File with information about all experiments annotated
## scRNASeqLibrary --> File with all libraries annotated by bgee
## output --> Path where should be saved the output results after collect the metadata

## librarie used
library(HCAExplorer)
library(data.table)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("scRNASeqExperiment","scRNASeqLibrary", "output")
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
if( file.exists(scRNASeqLibrary) ){
  annotation <- fread(scRNASeqLibrary, h=T, sep="\t")
  colnames(annotation)[1]<-"libraryId"
} else {
  stop( paste("The library file not found [", scRNASeqLibrary, "]\n"))
}
##########################################################################################################################
## download metadata from SRA
SRA_metadata <- function(libraryID){

  PID <- paste0(libraryID)
  ena.url <- paste("http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=",
                   PID,
                   "&result=read_run",
                   "&fields=experiment_accession,run_accession,",
                   "read_count,tax_id,scientific_name,",
                   "instrument_model,library_layout,fastq_ftp,submitted_ftp,",
                   "&download=text",
                   sep="")
  readInfo <- read.table(url(ena.url), header=TRUE, sep="\t")
  return(readInfo)
}

## download metadata from HCA
HCA_metadata <- function(libraryID, experimentName){

  hca <- HCAExplorer(url = 'https://service.explore.data.humancellatlas.org')
  page1_HCA <- data.frame(hca@results)
  hca_2 <- nextResults(hca)
  page2_HCA <- data.frame(hca_2@results)

  ## detect experiment per experiment name
  metadataExperiment1 <- nrow(page1_HCA[page1_HCA$projects.projectTitle %like% paste0(experimentName),]) == 0
  metadataExperiment2 <- nrow(page2_HCA[page2_HCA$projects.projectTitle %like% paste0(experimentName),]) == 0

  if(metadataExperiment1 == "FALSE" & metadataExperiment2 == "TRUE"){
    metadataExperiment <- page1_HCA[page1_HCA$projects.projectTitle %like% paste0(experimentName),]

    experiment_accession <- libraryID
    run_accession <- NA
    read_count <- NA
    tax_id <- NA
    scientific_name <- metadataExperiment$donorOrganisms.genusSpecies
    instrument_model <- metadataExperiment$protocols.instrumentManufacturerModel
    instrument_model <- substr(instrument_model,1,regexpr(",",instrument_model)-1)
    library_layout <- metadataExperiment$protocols.pairedEnd
    fastq_ftp <- 	metadataExperiment$fileTypeSummaries.fileType
    submitted_ftp	<- "HCA"

    infoMeatadata <- c(experiment_accession, run_accession, read_count, tax_id,scientific_name, instrument_model, library_layout, fastq_ftp, submitted_ftp)

  } else if (metadataExperiment1 == "TRUE" & metadataExperiment2 == "FALSE") {
    metadataExperiment <- page2_HCA[page2_HCA$projects.projectTitle %like% paste0(experimentName),]

    experiment_accession <- libraryID
    run_accession <- NA
    read_count <- NA
    tax_id <- NA
    scientific_name <- metadataExperiment$donorOrganisms.genusSpecies
    instrument_model <- metadataExperiment$protocols.instrumentManufacturerModel
    instrument_model <- substr(instrument_model,1,regexpr(",",instrument_model)-1)
    library_layout <- metadataExperiment$protocols.pairedEnd
    fastq_ftp <- 	metadataExperiment$fileTypeSummaries.fileType
    submitted_ftp	<- "HCA"

    infoMeatadata <- c(experiment_accession, run_accession, read_count, tax_id,scientific_name, instrument_model, library_layout, fastq_ftp, submitted_ftp)

  } else {
    cat("The experiment title was not found in HCA metadata.", "\n", "Please check the annotation file.")
  }
  return(infoMeatadata)
}

## select just target-based protocols to retrieve metadata
targetBased <- experiments[experiments$protocol %like% "10X Genomics", ]
targetBased_libraries <- dplyr::filter(annotation, protocolType == "3'end" & protocol == "10X Genomics")
experimentsIDs <- unique(targetBased_libraries$experimentId)
targetBased <- dplyr::filter(targetBased, experimentId %in% experimentsIDs)

metadata <- c()
for (experiment in unique(targetBased$experimentId)) {
  ## select source  of the experiment
  sourceID <- as.character(unique(targetBased$experimentSource[targetBased$experimentId == experiment]))

  if(sourceID == "SRA"){
    libraries_from_Exp <- annotation$libraryId[annotation$experimentId == experiment & annotation$protocol == "10X Genomics" & annotation$protocolType == "3'end"]
    for (i in libraries_from_Exp) {
      extractSRA <- SRA_metadata(libraryID = i)
      metadata <- rbind(metadata,extractSRA)
    }
  } else if (sourceID == "HCA"){
    libraries_from_Exp <- annotation$libraryId[annotation$experimentId == experiment & annotation$protocol == "10X Genomics" & annotation$protocolType == "3'end"]
    experimentNAME <- as.character(unique(targetBased$experimentName[targetBased$experimentId == experiment]))
    for (i in libraries_from_Exp) {
      extractHCA <- HCA_metadata(libraryID = i, experimentName = experimentNAME)
      metadata <- rbind(metadata,extractHCA)
    }
  } else {
    cat("Source is not recognized!")
  }
}

## Create the output files to write the comparison between annotation and metada!
metadata_info <- paste0(output, "/","metadata_info_10X.txt")
if (!file.exists(metadata_info)){
  file.create(metadata_info)
  cat("experiment_accession\trun_accession\tread_count\ttax_id\tscientific_name\tinstrument_model\tlibrary_layout\tfastq_ftp\tsubmitted_ftp\n",file = file.path(output, "metadata_info_10X.txt"), sep = "\t")
} else {
  print("File already exist.....")
}
metadata_notMatch <- paste0(output, "/","metadata_notMatch_10X.txt")
if (!file.exists(metadata_notMatch)){
  file.create(metadata_notMatch)
  cat("experiment_accession\trun_accession\tread_count\ttax_id\tscientific_name\tinstrument_model\tlibrary_layout\tfastq_ftp\tsubmitted_ftp\n",file = file.path(output,"metadata_notMatch_10X.txt"), sep = "\t")
} else {
  print("File already exist.....")
}

## compare information from annotation and metadata
## SRA: libraryID, plataform and speciesID
## HCA: just plataform for the moment
for(i in 1:nrow(targetBased_libraries)) {

  annotationInfo <- targetBased_libraries[i,]
  metadataInfo <- metadata[i,]
  sourceInfo <- metadataInfo$submitted_ftp

  if (sourceInfo != "HCA" | is.na(sourceInfo) == TRUE){
    compare_library <- identical(as.character(annotationInfo[['libraryId']]),as.character(metadataInfo[['experiment_accession']]))
    compare_machine <- identical(as.character(annotationInfo[['platform']]),as.character(metadataInfo[['instrument_model']]))
    compare_speciesID <- identical(as.character(annotationInfo[['speciesId']]),as.character(metadataInfo[['tax_id']]))

    if (compare_library == "TRUE" && compare_machine == "TRUE" && compare_speciesID == "TRUE"){
      cat(as.character(annotationInfo$X.libraryId[1]), "complete match between annotation and metadata", "\n")

      ## export libraries that pass and will be downloaded
      write.table(metadataInfo, file = metadata_info, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)

    } else {
      cat("WARNING : ", "\n")
      cat("For the library: ", as.character(annotationInfo$libraryId[1]), "the comparison library is: ", compare_library, "\n")
      cat("For the library: ", as.character(annotationInfo$libraryId[1]), "the comparison platform is: ", compare_machine, "\n")
      cat("For the library: ", as.character(annotationInfo$libraryId[1]), "the comparison species is: ", compare_speciesID, "\n")

      ## export libraries that will not be used to download
      write.table(metadataInfo, file = metadata_notMatch, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
    }
  } else {
    ## The HCA metadata we just compare the instrument model (since we not retrieve experiment_accession and tax_id just species_name)
    compare_machine <- identical(as.character(annotationInfo[['platform']]),as.character(metadataInfo[['instrument_model']]))
    if (compare_machine == "TRUE"){
      cat(as.character(annotationInfo$X.libraryId[1]), "complete match between annotation and metadata", "\n")
      ## export libraries that pass and will be downloaded
      write.table(metadataInfo, file = metadata_info, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
    } else {
      cat("WARNING : ", "\n")
      cat("For the library: ", as.character(annotationInfo$libraryId[1]), "the comparison platform is: ", compare_machine, "\n")
      ## export libraries that will not be used to download
      write.table(metadataInfo, file = metadata_notMatch, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
    }
  }
}
