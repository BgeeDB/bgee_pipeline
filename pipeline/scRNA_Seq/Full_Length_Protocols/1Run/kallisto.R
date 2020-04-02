## SFonsecaCosta, Oct 2018
## This script is used to run kallisto for each library

## Usage:
## R CMD BATCH --no-save --no-restore '--args library_id="library_id" scrna_seq_sample_info="scrna_seq_sample_info.txt" raw_cells_folder="raw_cells_folder" infoFolder="infoFolder" output_folder="output_folder"' kallisto.R kallisto.Rout
## library_id --> library where the kallisto should be runned
## scrna_seq_sample_info --> annotation file that collect information about each cell/experiment
## raw_cells_folder --> folder where is located the raw data fastq.gz files per cell-type
## infoFolder --> folder where is localized the transcriptome index + gene2transcript + gene2biotype
## output_folder --> folder where should be saved the kallisto output

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}
## checking if all necessary arguments were passed....
command_arg <- c("library_id", "scrna_seq_sample_info", "raw_cells_folder", "infoFolder", "output_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}
## Read scrna_seq_sample_info file. If file not exists, script stops
if( file.exists(scrna_seq_sample_info) ){
  annotation <- read.table(scrna_seq_sample_info, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
} else {
  stop( paste("scrna_seq_sample_info file not found [", scrna_seq_sample_info, "]\n"))
}
############################################### FUNCTION ############################################################
## Function to run kallisto
kallisto <- function(library_id, raw_cells_folder, infoFolder, output_folder){

  ## collect info library
  libraryTypeInfo <- as.character(annotation$libraryType[annotation$libraryId == library_id])
  libraryReadLengthInfo <- annotation$readLength[annotation$libraryId == library_id]
  species <- as.character(annotation$organism[annotation$libraryId == library_id])

  ## collect path to transcriptome
  indexFile31 <- list.files(path=infoFolder, pattern = paste0("*.transcriptome.idx$"), full.names=T, recursive = TRUE)
  indexFile31 <- grep(species,indexFile31, value = TRUE)
  indexFile15 <- list.files(path=infoFolder, pattern = paste0("*.transcriptome_k15.idx$"), full.names=T, recursive = TRUE)
  indexFile15 <- grep(species,indexFile15, value = TRUE)

  ## path to fastq.gz files
  libTarget <- file.path(raw_cells_folder,library_id)
  fastqFile <- list.files(path=libTarget, pattern = ".*gz$", full.names=T, recursive = TRUE)

  ## create a folder for each library in the output directory and write pseudo-alignment
  ifelse(!dir.exists(file.path(output_folder, library_id)), dir.create(file.path(output_folder, library_id)), FALSE)
  outputLib <- file.path(output_folder, library_id)

  if (libraryTypeInfo == "SINGLE" && libraryReadLengthInfo >= 50){
    cat("The library ", library_id, " is SINGLE and read length > 50", "\n")
    system(sprintf('%s -i %s -o %s  --single -l 180 -s 20 %s',paste0("kallisto quant"), indexFile31, outputLib, fastqFile))
  } else if (libraryTypeInfo == "SINGLE" && libraryReadLengthInfo < 50){
    cat("The library ", library_id, " is SINGLE and read length < 50", "\n")
    system(sprintf('%s -i %s -o %s  --single -l 180 -s 20 %s',paste0("kallisto quant"), indexFile15, outputLib, fastqFile))
  } else if (libraryTypeInfo == "PAIRED" && libraryReadLengthInfo >= 50){
    cat("The library ", library_id, " is PAIRED and read length > 50", "\n")
    system(sprintf('%s -i %s -o %s %s %s',paste0("kallisto quant"), indexFile31, outputLib, fastqFile[1],fastqFile[2]))
  } else {
    cat("The library ", library_id, " is PAIRED and read length < 50", "\n")
    system(sprintf('%s -i %s -o %s %s %s',paste0("kallisto quant"), indexFile15, outputLib, fastqFile[1],fastqFile[2]))
  }

  ## control file in case server stops in the middle of the job.
  file.create(file.path(outputLib, "DONE.txt"))
}
############################################### APPLY PER LIBRARY ############################################################
## Apply for all libraries from different species
kallisto(library_id = library_id, raw_cells_folder = raw_cells_folder, infoFolder = infoFolder, output_folder = output_folder)

