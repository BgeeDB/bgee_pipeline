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

#create output dir if not already done
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

## Create the output files to write the comparison between annotation and metadata
metadata_info <- file.path(output_folder,"metadata_info.tsv")
file.create(metadata_info)
metadata_info_not_match <- file.path(output_folder,"metadata_info_not_match.tsv")
file.create(metadata_info_not_match)

#create two data.frame for libraries passing/not passing the verification
metadata_colnames <- c("sample_accession","experimentId","library_id","run_accession","read_count",
  "tax_id","scientific_name","instrument_model","library_layout","fastq_ftp","submitted_ftp");
passed <- data.frame(matrix(nrow = 0, ncol = length(metadata_colnames)))
colnames(passed) <- metadata_colnames
# add a column reason for exclusion to this file
not_passed <- data.frame(matrix(nrow = 0, ncol = length(metadata_colnames) + 1))
colnames(not_passed) <- c(metadata_colnames, "exclusion_reason")

for (row in seq(nrow(annotation))) {
  library <- annotation[row,]
  libraryID <- library$libraryId
  metadata_info <- tryCatch(
    {
      # retrieve information from EBI metadata for each library annotated
      ena.url <- paste("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
                   libraryID,
                   "&result=read_run",
                   "&fields=experiment_accession,run_accession,",
                   "read_count,tax_id,scientific_name,",
                   "instrument_model,library_layout,fastq_ftp,submitted_ftp,",
                   "&download=TRUE",
                   sep="")
      read.table(url(ena.url), header=TRUE, sep="\t")
    },
    error=function(cond) {
      error <- TRUE
      # Choose a return value in case of error
      return(NA)
    }
  )
  # check if there was an error during EBI metadata retrieval
  if (!is.data.frame(metadata_info) && is.na(metadata_info)) {
    metadata_info <- data.frame(matrix(nrow = 1, ncol = length(metadata_colnames) + 1))
    colnames(metadata_info) <- c(metadata_colnames, "exclusion_reason")
    metadata_info$library_id <- libraryID
    metadata_info$exclusion_reason <- "EBI URL error"
    not_passed <- rbind(not_passed, metadata_info)
  } else {
    ## What is called experiment_accession in ENA API is called library_id in our pipeline
    names(metadata_info)[names(metadata_info) == 'experiment_accession'] <- 'library_id' 
    compare_machine <- identical(as.character(library$platform),
      as.character(metadata_info$instrument_model))
    compare_speciesID <- identical(as.character(library[['speciesId']]),
      as.character(metadata_info[['tax_id']]))
    
    ## before writing the line to the file we first update the line to fit the
    ## experiment_id/library_id nomenclature used in the pipeline
    metadata_info <- merge(metadata_info, library[, c("libraryId","experimentId")],
      by.x="library_id", by.y="libraryId")
    ## update column order
    metadata_info <- metadata_info[, metadata_colnames]
    ## homogeneise use of snake case for column names
    names(metadata_info)[names(metadata_info) == 'experimentId'] <- 'experiment_id'      
    if (isTRUE(compare_machine) && isTRUE(compare_speciesID)) {
      ## export libraries that pass and will be downloaded
      passed <- rbind(passed, metadata_info)
    } else {
      if(isFALSE(compare_machine)) {
        metadata_info$exclusion_reason <- "protocol_mismatch"
      }
      if(isFALSE(compare_speciesID)) {
        metadata_info$exclusion_reason <- "species_mismatch"
      }
      not_passed <- rbind(not_passed, metadata_info)
    }
  }
}

# write metadata files
write.table(passed, metadata_info, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(not_passed, metadata_info_not_match, sep = "\t", col.names = TRUE, row.names = FALSE,
  quote = FALSE)

