## SFonsecaCosta, September 2019
## This script is used to retrieve the metadata from EBI
## Then compare the annotation information from BGEE for each library with metadata from EBI.

## Usage:
## R CMD BATCH --no-save --no-restore '--args NEW_scRNASeqLibrary="NEW_scRNASeqLibrary.tsv" output_folder="output_folder"' retrieve_metadata.R retrieve_metadata.Rout
## NEW_scRNASeqLibrary -> File with all libraries annotated by bgee that pass the minimum requirement (>= 100 cells)
## output_folder -> Path where should be saved the output results after collect the metadata

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments exist....
command_arg <- c("NEW_scRNASeqLibrary", "output_folder")
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
  stop( paste("The annotation file not found [", NEW_scRNASeqLibrary, "]\n"))
}
##########################################################################################################################
## select just annotated protocols that are full length
fullLength <- data.frame(dplyr::filter(annotation, protocolType == "Full-length"))

# retrieve information from EBI metadata for each library annotated
libraryIDs <-fullLength[,1]
metadata <- c()
for (libraryID in unique(libraryIDs)) {

  PID <- paste0(libraryID)

  ena.url <- paste("https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=",
                   PID,
                   "&result=read_run",
                   "&fields=experiment_accession,run_accession,",
                   "read_count,tax_id,scientific_name,",
                   "instrument_model,library_layout,fastq_ftp,submitted_ftp,",
                   "&download=TRUE",
                   sep="")
  readInfo <- read.table(url(ena.url), header=TRUE, sep="\t")
  metadata <- rbind(metadata,readInfo)
}

## Create the output files to write the comparison between annotation and metadata
metadata_info <- file.path(output_folder, "metadata_info.txt")
if (!file.exists(metadata_info)){
  file.create(metadata_info)
  cat("sample_accession\texperiment_accession\trun_accession\tread_count\ttax_id\tscientific_name\tinstrument_model\tlibrary_layout\tfastq_ftp\tsubmitted_ftp\n",file = file.path(output_folder, "metadata_info.txt"), sep = "\t")
} else {
  print("File already exist.....")
}
metadata_notMatch <- file.path(output_folder, "metadata_notMatch.txt")
if (!file.exists(metadata_notMatch)){
  file.create(metadata_notMatch)
  cat("sample_accession\texperiment_accession\trun_accession\tread_count\ttax_id\tscientific_name\tinstrument_model\tlibrary_layout\tfastq_ftp\tsubmitted_ftp\n",file = file.path(output_folder, "metadata_notMatch.txt"), sep = "\t")
} else {
  print("File already exist.....")
}

## compare information between annotation and metadata (like: libraryID, plataform and speciesID)
for(i in 1:nrow(fullLength)) {

  annotationInfo <- fullLength[i,]
  metadataInfo <- metadata[i,]

  compare_library <- identical(as.character(annotationInfo[['libraryId']]),as.character(metadataInfo[['experiment_accession']]))
  compare_machine <- identical(as.character(annotationInfo[['platform']]),as.character(metadataInfo[['instrument_model']]))
  compare_speciesID <- identical(as.character(annotationInfo[['speciesId']]),as.character(metadataInfo[['tax_id']]))


  if (compare_library == "TRUE" && compare_machine == "TRUE" && compare_speciesID == "TRUE"){
    cat(as.character(annotationInfo$libraryId[1]), "complete match between annotation and metadata", "\n")

    ## export libraries that pass and will be downloded
    write.table(metadataInfo, file = metadata_info, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)

  } else {
    cat("WARNING : ", "\n")
    cat("For the library: ", as.character(annotationInfo$libraryId[1]), "the comparison library is: ", compare_library, "\n")
    cat("For the library: ", as.character(annotationInfo$libraryId[1]), "the comparison platform is: ", compare_machine, "\n")
    cat("For the library: ", as.character(annotationInfo$libraryId[1]), "the comparison species is: ", compare_speciesID, "\n")

    ## export libraries that will not be used to download
    write.table(metadataInfo, file = metadata_notMatch, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
  }
}
