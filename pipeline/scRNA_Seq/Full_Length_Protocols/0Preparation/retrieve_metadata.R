## SFonsecaCosta, September 2019
## This script is used to retrieve the metadata from EBI
## Then to compare the annotation information from BGEE for each library with metadata from EBI.

## Usage:
## R CMD BATCH --no-save --no-restore '--args pass_annotationControl="passScRNASeqLibrary.tsv" output_folder="output_folder"' retrieve_metadata.R retrieve_metadata.Rout
## pass_annotationControl -> File with all libraries annotated by bgee that pass the minimum requirement (>= 50 cells)
## output_folder --> folder where the output files should be saved

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments exist....
command_arg <- c("pass_annotationControl", "output_folder")
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
  stop( paste("The annotation file not found [", pass_annotationControl, "]\n"))
}
##########################################################################################################################
# retrieve information from EBI metadata for each library annotated
libraryIDs <-annotation[,1]
metadata <- c()
for (libraryID in unique(libraryIDs)) {

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
  metadata <- rbind(metadata,readInfo)
}
## What is called experiment_accession in ENA API is called library_id in our pipeline
names(metadata)[names(metadata) == 'experiment_accession'] <- 'library_id' 

## Create the output files to write the comparison between annotation and metadata
metadata_info <- file.path(output_folder,"metadata_info.tsv")
if (file.exists(metadata_info)){
  message("File already exists and will be removed to create a new one to avoid overwritting!")
  file.remove(metadata_info)
  file.create(metadata_info)
} else {
  file.create(metadata_info)
}
cat("sample_accession\texperiment_id\tlibrary_id\trun_accession\tread_count\ttax_id\tscientific_name\tinstrument_model\tlibrary_layout\tfastq_ftp\tsubmitted_ftp\n",file = metadata_info, sep = "\t")

metadata_info_not_match <- file.path(output_folder,"metadata_info_not_match.tsv")
if (file.exists(metadata_info_not_match)){
  message("File already exists and will be removed to create a new one to avoid overwritting!")
  file.remove(metadata_info_not_match)
  file.create(metadata_info_not_match)
} else {
  file.create(metadata_info_not_match)
}
cat("sample_accession\texperiment_id\tlibrary_id\trun_accession\tread_count\ttax_id\tscientific_name\tinstrument_model\tlibrary_layout\tfastq_ftp\tsubmitted_ftp\n",file = metadata_info_not_match, sep = "\t")

## compare information between annotation and metadata (like: libraryID, plataform and speciesID)
for(i in 1:nrow(annotation)) {

  annotationInfo <- annotation[i,]
  metadataInfo <- metadata[i,]

  # what we call library is called experiment accession at EBI
  compare_library <- identical(as.character(annotationInfo[['libraryId']]),as.character(metadataInfo[['library_id']]))
  compare_machine <- identical(as.character(annotationInfo[['platform']]),as.character(metadataInfo[['instrument_model']]))
  compare_speciesID <- identical(as.character(annotationInfo[['speciesId']]),as.character(metadataInfo[['tax_id']]))


  if (compare_library == "TRUE" && compare_machine == "TRUE" && compare_speciesID == "TRUE"){
    message(as.character(annotationInfo$libraryId[1]), "complete match between annotation and metadata")
    ## before writing the line to the file we first update the line to fit the experiment_id/library_id
    ## nomenclature used in the pipeline
    merge(metadataInfo, annotationInfo[, c("libraryId","experimentId")], by.x="library_id", by.y="libraryId")
    ## update column order
    metadataInfo <- metadataInfo[, c("sample_accession","experimentId","library_id","run_accession","read_count","tax_id","scientific_name","instrument_model","library_layout","fastq_ftp","submitted_ftp")]
    ## homogeneise use of snake case for column names
    names(metadata)[names(metadata) == 'experimentId'] <- 'experiment_id'
    ## export libraries that pass and will be downloaded
    write.table(metadataInfo, file = metadata_info, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)

  } else {
    warning("For the library: ", as.character(annotationInfo$libraryId[1]), "the comparison library is: ", compare_library)
    warning("For the library: ", as.character(annotationInfo$libraryId[1]), "the comparison platform is: ", compare_machine)
    warning("For the library: ", as.character(annotationInfo$libraryId[1]), "the comparison species is: ", compare_speciesID)

    ## export libraries that will not be used to download
    write.table(metadataInfo, file = metadata_info_not_match, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
  }
}

