## Julien Wollbrett March 28, 2023
## This script creates a report containing number reads, number aligned, proportion aligned,
## aligned unique, and number target sequences
## It loops through all libraries

## Usage:
## R CMD BATCH --no-save --no-restore '--args metadata_file="/path/to/bgeecall_input.tsv" kallisto_dir=$(RNASEQ_CLUSTER_BGEECALL_OUTPUT) fastq_dir=$(FASTQ_DIR) output_file=$(TARGET_BASED_STATS)' rna_seq_calls_plot.R rna_seq_calls_plot.Rout
## metadata_file              - path to file containing Bgee annotation for all libraries to process
## kallisto_dir               - path to folder where BgeeCall wrote the calls.
## fastq_dir                  - path to the parent dir containing all fastq
## output_file                - path to output file containing report of kallisto and read size for all samples

## load required packages
library("rjson")

## Session info
print(sessionInfo())

## reading in arguments provided in command line
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed in command line
command_arg <- c("metadata_file", "kallisto_dir", "fastq_dir", "output_file")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

##load data
libraries <- read.table(metadata_file, sep = "\t", header = TRUE, comment.char = "")
## in this file each row corrspond to one run. As we want info at library level, we have to remove
## duplicated library IDs
libraries <- libraries[!duplicated(libraries$library_id),]

## init variables for calls report
kallisto_info_file   <- "run_info.json"
reads_R_stat_suffix  <- ".R.stat"
kallisto_info_all_samples <- data.frame(matrix(nrow = 0, ncol = 11))
kallisto_info_report_columns <- c("libraryId", "reads", "min_read_size", "max_read_size", "number_aligned", "number_unique_aligned", 
  "prop_aligned", "prop_unique_aligned", "number_targets", "start_time", "kallisto_version")
names(kallisto_info_all_samples) <- kallisto_info_report_columns

## check that calls have been generated for all libraries
libraries_wo_calls <- 0
lib_dirs <- list.dirs(path = kallisto_dir, full.names = FALSE, recursive = FALSE)
message(nrow(libraries), " libraries in the bgeecall info file")
for(line in seq(nrow(libraries))) {
  library_id <- as.character(libraries$library_id[line])
  species_id <- as.character(libraries$tax_id[line])
  if(dir.exists(file.path(kallisto_dir, library_id))) {
  ##if(library_id %in% lib_dirs) {

    # retrieve info for kallisto report
    kallisto_info_file_path <- file.path(kallisto_dir, library_id, kallisto_info_file)
    library_r_stat_files <- list.files(path = file.path(fastq_dir, species_id, library_id), 
      pattern = reads_R_stat_suffix, full.names = TRUE, recursive = TRUE)
    if(file.exists(kallisto_info_file_path)) {
      if ((length(library_r_stat_files) > 0) ) {
        json_info <- fromJSON(file = kallisto_info_file_path)
        min_reads <- Inf
        max_reads <- 0
        for (r_stat_file in library_r_stat_files) {
          r_stats <- read.table(file = r_stat_file, header = TRUE, sep = "\t", comment.char = "")
          if (r_stats[1,1] < min_reads) {
            min_reads <- r_stats[1,1]
          }
          if (r_stats[1,2] > max_reads) {
            max_reads <- r_stats[1,2]
          }
        }
        kallisto_info <- data.frame(library_id, json_info$n_processed, r_stats[1,1], r_stats[1,2], 
          json_info$n_pseudoaligned, json_info$n_unique, json_info$p_pseudoaligned, json_info$p_unique, 
          json_info$n_targets, json_info$start_time, json_info$kallisto_version)
        names(kallisto_info) <- kallisto_info_report_columns
        kallisto_info_all_samples <- rbind(kallisto_info_all_samples, kallisto_info)
        } else {
          warning(library_id, " : R.stat file not generated")
          libraries_wo_calls <- libraries_wo_calls + 1
        }
    } else {
      warning(library_id, " : kallisto info file was not generated")
      libraries_wo_calls <- libraries_wo_calls + 1
    }
  } else {
    warning(library_id, " : library directory not created ", file.path(kallisto_dir, library_id))
    libraries_wo_calls <- libraries_wo_calls + 1
  }
}

if (libraries_wo_calls > 0 ) {
  warning("Calls or R.stat file were not generated for ", libraries_wo_calls, " libraries.")
}

# save kallisto report
write.table(kallisto_info_all_samples, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE)

