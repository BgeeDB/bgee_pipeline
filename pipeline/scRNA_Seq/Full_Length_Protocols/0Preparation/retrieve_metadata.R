## SFonsecaCosta, September 2019
## This script is used to retrieve the metadata from EBI
## Then to compare the annotation information from BGEE for each library with metadata from EBI.

## parameters:
## filteredLibraryAnnotation  --> File with RNASeq FL libraries annotation
## metadataFile               --> file where metadata info are saved
## threads                    --> number of threads used to retrieve metadata

library(foreach)
library(doParallel)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments exist....
command_arg <- c("filteredLibraryAnnotation", "metadataFile", "threads")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read annotation file. If file not exists, script stops
if( file.exists(filteredLibraryAnnotation) ){
  annotation <- read.table(filteredLibraryAnnotation, h=TRUE, sep="\t", comment.char="")
} else {
  stop( paste("Library annotation file not found [", filteredLibraryAnnotation, "]\n"))
}

#create two data.frame for libraries passing/not passing the verification
metadata_colnames <- c("sample_accession","experimentId","library_id","run_accession","read_count",
  "tax_id","scientific_name","instrument_model","library_layout","fastq_ftp","submitted_ftp");
passed <- data.frame(matrix(nrow = 0, ncol = length(metadata_colnames)))
colnames(passed) <- metadata_colnames
# add a column reason for exclusion to this file
not_passed <- data.frame(matrix(nrow = 0, ncol = length(metadata_colnames) + 1))
colnames(not_passed) <- c(metadata_colnames, "exclusion_reason")

print("Starting parallelized loop to retrieve metadata from EBI")
registerDoParallel(cores = threads)

# Parallelized loop
results <- foreach(row = seq(nrow(annotation)), .combine = list, .multicombine = TRUE, .packages = c("utils")) %dopar% {
  library <- annotation[row,]
  libraryID <- library$libraryId
  metadata_info <- tryCatch(
    {
      # retrieve information from EBI metadata for each library annotated
      ena.url <- paste("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
                   libraryID,
                   "&result=read_run",
                   "&fields=sample_accession,experiment_accession,run_accession,",
                   "read_count,tax_id,scientific_name,",
                   "instrument_model,library_layout,fastq_ftp,submitted_ftp,",
                   "&download=TRUE",
                   sep="")
      read.table(url(ena.url), header=TRUE, sep="\t")
    },
    error=function(cond) {
      return(NA)
    }
  )
  # check if there was an error during EBI metadata retrieval
  if (!is.data.frame(metadata_info) && is.na(metadata_info)) {
    metadata_info <- data.frame(matrix(nrow = 1, ncol = length(metadata_colnames) + 1))
    colnames(metadata_info) <- c(metadata_colnames, "exclusion_reason")
    metadata_info$library_id <- libraryID
    metadata_info$exclusion_reason <- "EBI URL error"
    return(list(type = "not_passed", data = metadata_info))
  } else {
    ## What is called experiment_accession in ENA API is called library_id in our pipeline
    names(metadata_info)[names(metadata_info) == 'experiment_accession'] <- 'library_id' 
    compare_machine <- identical(as.character(library$platform),
      as.character(metadata_info$instrument_model))
    compare_speciesID <- identical(as.character(library[['speciesId']]),
      as.character(metadata_info[['tax_id']]))

    metadata_info <- merge(metadata_info, library[, c("libraryId","experimentId")],
      by.x="library_id", by.y="libraryId")
    metadata_info <- metadata_info[, metadata_colnames]
    names(metadata_info)[names(metadata_info) == 'experimentId'] <- 'experiment_id'      
    if (isTRUE(compare_machine) && isTRUE(compare_speciesID)) {
      return(list(type = "passed", data = metadata_info))
    } else {
      if(isFALSE(compare_machine)) {
        metadata_info$exclusion_reason <- "protocol_mismatch"
      }
      if(isFALSE(compare_speciesID)) {
        metadata_info$exclusion_reason <- "species_mismatch"
      }
      return(list(type = "not_passed", data = metadata_info))
    }
  }
}

passed <- do.call(rbind, lapply(results, function(x) {
  if (!is.null(x) && is.list(x) && !is.null(x$type) && x$type == "passed") {
    return(x$data)
  } else {
    return(NULL)
  }
}))

not_passed <- do.call(rbind, lapply(results, function(x) {
  if (!is.null(x) && is.list(x) && !is.null(x$type) && x$type == "not_passed") {
    return(x$data)
  } else {
    return(NULL)
  }
# write metadata files
write.table(passed, metadataFile, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
warning("Libraries for which metadata were not retrieved:\n", not_passed)

