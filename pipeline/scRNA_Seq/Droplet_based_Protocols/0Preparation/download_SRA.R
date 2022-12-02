## SFonsecaCosta, 17 Sep 2019

## This script is used to download the raw data for target-based protocols (for the moment 10X)
## For each library we should have at least 3 files: R1, R2 and I (fastq.gz format)

## /!\ Be aware of the hack used to rename fastq files at the end of this script. When running this
## script please first check if a better solution exists now.

## Usage:
## R CMD BATCH --no-save --no-restore '--args metadata_info="metadata_info_10X.txt" librariesDownloadedJura="librariesDownloadedJura.tsv" output="output_folder"' download_SRA.R download_SRA.Rout
## metadata_info --> File with all libraries that need to be download
## librariesDownloadedJura --> File with information about the libraries that already exist in jura server and don't need to be downloaded
## output --> Path where should be saved the output results for each library

## libraries used
library(dplyr)
library(HelpersMG)
library(tools)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("metadata_info", "librariesDownloadedJura","output")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read EBI file. If file not exists, script stops
if( file.exists(metadata_info) ){
  readFile <- read.table(metadata_info, h=T, sep="\t")
} else {
  stop( paste("metadata file not found [", metadata_info, "]\n"))
}

## Read EBI file. If file not exists, script stops
if( file.exists(librariesDownloadedJura) ){
  libDownloadJura <- read.table(librariesDownloadedJura, h=T, sep="\t")
} else {
  stop( paste("Libraries from JURA file not found [", librariesDownloadedJura, "]\n"))
}

## Check libraries
checkId <- data.frame(readFile$library_id[!(readFile$library_id %in% libDownloadJura$library_id)])
colnames(checkId) <- "library_id"
finalLibsToDown <- dplyr::filter(readFile, library_id %in% checkId$library_id)

## remove NA and fastq.gz information that belongs to HCA (download with other script) or from EBI (download .bam files)
finalLibsToDown <- dplyr::filter(readFile, fastq_ftp != "NA" & fastq_ftp != "fastq.gz")


####################### FUNCTIONS ###############################
check_seq_length <- function(fileName) {
  connection <- file(fileName,"r")
  first_line <- readLines(connection, n=1)
  close(connection)
  seq_length <- strsplit(x = regmatches(x = first_line, m = regexpr("length=[0-9]+$",first_line)),
    split = "=")[[1]][2]
  return(seq_length)
}
#################################################################


folder_experiments <- file.path(output, "EXPERIMENTS")
for (library in  unique(finalLibsToDown$library_id)) {

  ## create output for each library
  message("treating the library: ", library)
  species <- unique(finalLibsToDown$tax_id[finalLibsToDown$library_id == library])
  experiment <- unique(finalLibsToDown$experiment_id[finalLibsToDown$library_id == library])
  folder_species <- file.path(output, species)
  if (!file.exists(folder_species)) {
    dir.create(folder_species)
  }
  folder_library <- file.path(folder_species, library)
  if (!file.exists(folder_library)) {
    dir.create(folder_library)
  }
  ## select the SRR referent to the library to download
  run_accessions <- as.character(finalLibsToDown$run_accession[finalLibsToDown$library_id == library])
  setwd(folder_library)

  ## download run files using SRA Toolkit
  for(run in run_accessions) {
    system(sprintf('prefetch %s ', run))
    ## find the name as it can have .sra or .sralite as extension
    sra_file <- list.files(path = file.path(folder_library, run), pattern = run, full.names = TRUE)
    fastq_run_dir <- file.path(folder_library, run, "FASTQ")
    if (length(sra_file) == 1) {
      system(sprintf('fastq-dump --outdir %s --split-files %s', fastq_run_dir,
        sra_file))
    } else {
      stop("more than one sra file has been downloaded for the run ", run, ". Delete one of them",
        " to be able to generate the fastq files")
    }
    ## SRA does not provide precise names. There are no information about I1, R1 or R2 in the file names.
    ## It is not easy to find back which file contains which info.
    ## It looks like, for 10x v2 and 10x v3:
    ## 3 files ==> _1.fastq corresponds to I1. _2.fastq corresponds to R1. _3.fastq corresponds to R2.
    ## 2 files ==> _1.fastq corresponds to R1. _2.fastq corresponds to R2.
    ## This assumption is used in this script BUT in order to be more sure one check has been used based on
    ## expected sequence length: I1 < R1 < R2.
    ## As it is definitly not bullet proof (and because we did not find any other criteria), the script is
    ## stopped with an error message when this rule is not ok.
    ## once this rule has been validated, fastq file names is updated.
    fastq_files_names <- list.files(path = fastq_run_dir, pattern = ".fastq", full.names = TRUE)
    name_to_seq_length <- data.frame(matrix(NA, nrow = length(fastq_files_names), ncol = 2))
    ## we never expect to have 1 file or more than 3 files. Update implementation if it can happen
    if (length(fastq_files_names) == 1 || length(fastq_files_names) > 3) {
      stop("wrong number of files. Expected 2 or 3 but have ", length(fastq_files_names),
        " for run ", run)
    }
    for (file_name in fastq_files_names) {
      if (grepl(pattern = "_1\\.fastq$", file_name, perl = TRUE)) {
        name_to_seq_length[1,] <- c(file_name, check_seq_length(file_name))
      } else if (grepl(pattern = "_2\\.fastq$", file_name, perl = TRUE)) {
        name_to_seq_length[2,] <- c(file_name, check_seq_length(file_name))
      } else if (grepl(pattern = "_3\\.fastq$", file_name, perl = TRUE)) {
        name_to_seq_length[3,] <- c(file_name, check_seq_length(file_name))
      } else {
        stop("fastq file name ", file_name, " can not be transformed.")
      }
    }
    ## check that size sequences from I1 < R1 < R2
    if (is.unsorted(name_to_seq_length[,2]) ||
        name_to_seq_length[1,2] > name_to_seq_length[2,2]) {
      stop("the sequence size in SRA files does not match the assumption that",
        " file _1.fastq has sequence smaller than _2.fastq file which itself",
        " has smaller sequences than _3.fastq (if it exists). It is then not",
        " possible to rename SRA files for run ", run)
    }
    print(name_to_seq_length)
    ## if 3 files we consider that there are I1, R1 and R2
    ## if 2 files we consider that there are only R1 and R2 (no )
    with_I1_file = length(fastq_files_names) == 3
    print(with_I1_file)
    for (row_number in seq(nrow(name_to_seq_length))) {
      fastq_dir_path <- dirname(name_to_seq_length[row_number,1])
      new_file_name <- NULL
      old_file_name <- basename(name_to_seq_length[row_number,1])
      if (with_I1_file) {
        if (row_number == 1) {
	  new_file_name <- gsub("_1", "_I1", old_file_name)
        } else if (row_number == 2) {
	  new_file_name <- gsub("_2", "_R1", old_file_name);
        } else if (row_number == 3) {
          new_file_name <- gsub("_3", "_R2", old_file_name);
	}
      } else {
        if (row_number == 1) {
          new_file_name <- gsub("_1", "_R1", old_file_name);
          print(paste0("row 1 rename ", old_file_name, " to ", new_file_name))
        } else if (row_number == 2) {
	  new_file_name <- gsub("_2", "_R2", old_file_name);
          print(paste0("row 2 rename ", old_file_name, " to ", new_file_name))
        }
      }
      file.rename(from = file.path(fastq_dir_path, old_file_name), to = file.path(fastq_dir_path, new_file_name))
      system(sprintf("gzip %s", file.path(fastq_dir_path, new_file_name)))
    }
    # everything was done properly for this run. It is now possible to delete the .sra file
    file.remove(sra_file)
  }

  ## create a symlink following pattern */EXPERIMENTS/experiment_ID/library_id/
  file.symlink(from = folder_library, to = file.path(folder_experiments, experiment, library), overwrite = FALSE, recursive = TRUE)
}

