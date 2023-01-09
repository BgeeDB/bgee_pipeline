## SFonsecaCosta, May 2020

## This script is used to create a transcript_to_gene file from the gtf_all file that will be used later in bustools.
## Also a file gene_to_biotype is created for each species.
## NOTE: if the index contain the intergenic regions the file transcript_to_gene should also contain all information: genic + intergenic

## Usage:
## R CMD BATCH --no-save --no-restore '--args folder_gtf="folder_gtf"' generateInfo.R generateInfo.Rout
## folder_gtf --> Folder where are all gtf files (gtf_all with intergenic) are for the different species

## libraries used
library(stringr)
library(dplyr)
# used to read gtf file and extract transcript/gene infos
library(rtracklayer)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("folder_gtf", "metadata_info_file")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

#read metadata file
metadata_info <- read.table(metadata_info_file, sep = "\t", header = TRUE)
# retrieve scientific name of species for which single cell libraries exist
species_to_use <- str_replace_all(unique(metadata_info$scientific_name), " ", "_")
all_files <- list.files(folder_gtf, pattern="*gtf_all$", full.names=TRUE, recursive = TRUE)
for (i in all_files) {
  gtf_file_name <- basename(i)
  # check if current gtf_all file corresponds to a species for which single 
  # cell libraries exist. If not move to next gtf_all file
  if (is.na(match(1,pmatch(species_to_use, gtf_file_name)))) {
    next
  }
  message("generate info using ", i)
  species_id <- data.frame(str_split(gtf_file_name, "\\."))[1,1]

  gtf_info <- import(i, format = "gtf")

  # extract all info from gtf. gene_name and gene_biotype can have NA values
  extracted_gtf_info <- as.data.frame(mcols(gtf_info)[,c("gene_id","gene_name","transcript_id", "gene_biotype")])
  
  # generate output files  
  write.table(extracted_gtf_info[,c("transcript_id", "gene_id")], file = file.path(folder_gtf, paste0(species_id, "_transcript_to_gene_with_intergenic.tsv")), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(extracted_gtf_info[,c("gene_id", "gene_biotype")], file = file.path(folder_gtf, paste0(species_id, "_gene_to_biotype_with_intergenic.tsv")), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(extracted_gtf_info[,c("gene_id", "gene_name")], file = file.path(folder_gtf, paste0(species_id, "_gene_to_geneName_with_intergenic.tsv")), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}
