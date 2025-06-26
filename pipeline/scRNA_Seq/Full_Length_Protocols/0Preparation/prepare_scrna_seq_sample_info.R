## SFonsecaCosta, June 2019

## This script create the scrna_seq_sample_info file to run the scRNA-Seq pipeline.

## Usage:
## R CMD BATCH --no-save --no-restore '--args filteredLibraryAnnotation="passScRNASeqLibrary.tsv" fastqDir="fastqDir" outputDir="outputDir"' prepare_scrna_seq_sample_info.R prepare_scrna_seq_sample_info.Rout
## filteredLibraryAnnotation --> list of raw libraries annotation filtered to contain only full-length scRNA-Seq libraries
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


##################################################################### FUNCTION #############################################################################################
## retrieve information from fastp
collectInformationFASTP <- function(fastq_dir, annotation, library){
  species <- annotation$speciesId[annotation$libraryId == library]
  fastq_library_dir <- file.path(fastq_dir, species, library)
  fastq_files <- (list.files(path = fastq_library_dir, pattern = "*.fastq.gz"))

  ## check if fastp has already run
  fastpJSON <- list.files(path=fastq_library_dir, pattern = "*.fastp.json.xz")

  if (isTRUE(file.exists(file.path(fastq_library_dir, fastpJSON)))){
    library_type <- ifelse(length(fastq_files) == 1, "SINGLE", "PAIRED")
    json_output <- fromJSON(file = file.path(fastq_library_dir, fastpJSON))
    read_length <- json_output$read1_before_filtering$total_cycles
  } else {
    stop("fastp has not been run for the library ", library)
  }

  ## collect information per library
  information_file <- data.frame(library, library_type, read_length)
  return(information_file)
  }


##################################### OUTPUT #############################################################################################
## create a intermediary file to collect information about each library, or truncate already existing one....
#TODO: Should not create a tmp file but directly collect data in memory
tmp_info_file <- paste0(outputDir, ".tmp")
file.create(tmp_info_file)
cat("libraryId\tlibraryType\treadLength\n",file = tmp_info_file, sep = "\t")

###################################### RUN PER lIBRARY ##################################################################################
for (library_id in unique(annotation$libraryId)) {

  fastp_info <- collectInformationFASTP(fastq_dir = fastqDir, annotation = annotation, library = library_id)

  write.table(fastp_info, file = tmp_info_file, row.names = FALSE , col.names = FALSE , append = TRUE, quote = FALSE, sep = "\t")
}

### read final output file --> InfoFile
file_info <- read.table(tmp_info_file, header=TRUE, sep="\t")

## Create the scrna_seq_sample_info
scrna_seq_sample_info <- merge(annotation, file_info, by = "libraryId", incomparables = NaN)
scrna_seq_sample_info <- scrna_seq_sample_info %>% dplyr::select("libraryId", "experimentId", "cellTypeName_abInitio", "cellTypeId_abInitio", "speciesId","platform", "protocol", "protocolType", "libraryType", "infoOrgan", "stageId", "anatId", "sex", "strain", "readLength", "speciesId", "genotype")
scrna_seq_sample_info$organism <- "NaN"
final_table <- metadata$scientific_name[metadata$tax_id == scrna_seq_sample_info$speciesId]
write.table(final_table,file = file.path( outputDir, "scrna_seq_sample_info.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
## remove intermediary file
file.remove(tmp_info_file)
