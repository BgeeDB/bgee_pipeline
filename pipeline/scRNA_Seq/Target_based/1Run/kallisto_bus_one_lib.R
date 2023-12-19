## SFonsecaCosta, 20 Sep 2019
## Modified May 2020

## This script is used to run Kallisto bus in one library

## Usage:
## R CMD BATCH --no-save --no-restore '--args metadata_file="metadata_info_10X.txt" libraryId="libraryId" speciesId="speciesId" fastqDir="fastqDir" gtfDir="gtfDir" scRNASeqInfoFile="scRNASeqInfoFile" kallisto_bus_results="folder_save_Results_from_kallisto_bus"' kallisto_bus_one_lib.R kallisto_bus_one_lib.Rout
## libraryId            --> ID of the library to process
## speciesId            --> ensembl ID of the species for which the library corresponds to
## fastqDir             --> Folder where all the libraries are located after download (the libaries are 
##                          organized by: species_ID/library_ID)
## gtfDir               --> Folder where is placed the informative files as transcriptomes index + gtf_all
## scRNASeqInfoFile     --> Path to the scRNA_Seq_info_TargetBased file
## kallisto_bus_results --> Path where should be saved the kallisto bus results

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

## checking if all necessary arguments were passed.
command_arg <- c("libraryId","speciesId", "fastqDir", "gtfDir", "scRNASeqInfoFile", "kallisto_bus_results")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

##########################################################################################################################################################
## first merge information from the metadata (like: reads information and SRR) and annotation file (from Bgee)
scRNASeqInfo <- read.table(scRNASeqInfoFile, sep = "\t", header = TRUE, quote = "\"")

species <- unique(as.character(scRNASeqInfo$scientific_name[scRNASeqInfo$libraryId == libraryId]))
message("Species:", species)
## collect species info
speciesName <- gsub(" ", "_", species)
transcriptomeIndexFiles <- list.files(gtfDir, pattern = paste0("^", speciesName, ".*transcriptome.idx$"))
print(transcriptomeIndexFiles)
#TODO: transcriptome index can potentially be compressed with xz.

if (! file.exists(file.path(kallisto_bus_results, libraryId))) {
  dir.create(file.path(kallisto_bus_results, libraryId))
}

## verify if library exist
if (file.exists(file.path(fastqDir, speciesId, libraryId))) {
  already_processed_file <- file.path(kallisto_bus_results, libraryId, "Done_kallisto_bus.txt");
  ## verify if kallisto_bus was already executed for the library
  if (file.exists(already_processed_file)){
    message("Kallisto was already executed for this library ",libraryId,"\n")
  } else {
    ## collect info about whitelist
    protocol <- unique(as.character(scRNASeqInfo$protocol[scRNASeqInfo$libraryId == libraryId]))
    whiteLInfo <- sub(pattern = "10X Genomics ", replacement = "", x = protocol)

    message("The whitelist used is:", whiteLInfo)

    ## verify if exist FASTQ (with or without recursive folder) or SRR folder with fastq.gz files for the library
    fastqLibDir <- file.path(fastqDir, speciesId, libraryId, "FASTQ")

    ## create directory for bus_output for each library
    busOutput <- file.path(kallisto_bus_results, libraryId)
    if (!dir.exists(busOutput)){
      dir.create(busOutput, recursive=TRUE)
    }
    detectFastqPath <- file.path(fastqDir, speciesId, libraryId)
    readBarcodes <- list.files(path=detectFastqPath, pattern = "*_R1.fastq.gz$", recursive = TRUE,
      full.names = TRUE, include.dirs = FALSE)
    readSeq <- list.files(path = detectFastqPath, pattern = "*_R2.fastq.gz$", recursive = TRUE,
      full.names = TRUE, include.dirs = FALSE)
    print(readBarcodes)
    print(readSeq)

    ## collect all fastq.gz files
    filesKallisto <- rbind(readBarcodes, readSeq)
    filesKallisto <- toString(filesKallisto)
    filesKallisto <- gsub(",", " " ,filesKallisto)
    message("Files to pass to Kallisto: ")
    message(filesKallisto)
    
    # select index depending on cell compartment used (cell or single nucléus)
    transcriptomeIndexFile <- transcriptomeIndexFiles[grep(pattern = "\\.transcriptome.idx", x = transcriptomeIndexFiles)]
    if (unique(as.character(scRNASeqInfo$RNAseqTags[scRNASeqInfo$libraryId == libraryId])) == "Sn-scRNA-seq") {
      transcriptomeIndexFile <- transcriptomeIndexFiles[grep(pattern = "single_nucleus_transcriptome.idx$",
        x = transcriptomeIndexFiles)]
    }

    #RUN Kallisto bus
    system(sprintf('kallisto bus -m -i %s -o %s -x %s -t 4 %s', file.path(gtfDir, transcriptomeIndexFile), paste0(busOutput), paste0("10x", whiteLInfo), paste0(filesKallisto)))
    ## control file
    file.create(already_processed_file)
  }
} else {
  message("Library not present in the folder ", file.path(fastqDir, libraryId), " to run Kallisto bus!")
}
