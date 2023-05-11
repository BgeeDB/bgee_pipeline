# Julien Wollbrett
# created: 5/05/2023

## From GTF file downloaded from Ensembl, this script created a new GTF file containing nascent (not matured) transcripts: 
## - retrieve the start of the 1st exon and the stop of the last exon of the transcript
## - create a new transcript id by adding the prefix "_nascent" to the initial transcript ID e.g FBtr0306592_nascent
##   for the transcript FBtr0306592
## - create a gtf file with each line corresponding to a nascent transcript
## - also create a tsv file mapping the new transcript ID to the proper gene ID
##
## This gtf file is used only for single nucleus rnaseq libraries. It allows to map reads to not spliced intron or exon/intron regions

## Invoking:
# R CMD BATCH --no-save --no-restore --slave '--args gtf_dir = "gtf_dir" metadata_file = "metadata_file"' generate_GTF_nascent.R Rout_path

## Arguments to provide:
# "metadata_file"   - full path to the metadata file created by the pipeline
# "gtf_dir"         - full path to the folder containing GTF files

library(R.utils)

cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("gtf_dir", "metadata_file")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

######################## Functions #############################

## Function for obtaining the part of the annotation field from gtf file
## Input: whole field, already splitted, for example: gene_id "FBgn0264003"; gene_name "mir-5613"; gene_source "FlyBase"; gene_biotype "pre_miRNA"
## Output: could be for example gene_id: FBgn0264003
get_annot_value <- function(split_annotation, field_name){
  ## find the right field
  field_all <- split_annotation[grep(field_name, split_annotation, fixed = T)];

  ## split the field
  field_value <- strsplit(field_all, ' ', fixed = T)[[1]][2];

  ## remove the last ';' if necessairy
  field_value <- sub(';', '', field_value,fixed = T)

  return(field_value)
}

remove_annot_value <- function(split_annotation, field_name){
  ## find the right field
  return(split_annotation[grep(field_name, split_annotation, fixed = T, invert = TRUE)])
}

add_annot_value <- function(split_annotation, field_to_add, add_after_field){
  return(append(split_annotation, field_to_add, after = grep(add_after_field, split_annotation)))
}

update_transcript_id <- function(split_annotation){
  transcript_id_position <- grep("transcript_id", split_annotation)
  split_annotation[transcript_id_position] <- paste0(split_annotation[transcript_id_position], "_nascent")
  return(split_annotation)
}

###################################################################

# read metadata file to retrieve all species for which single cell target base libraries are annotated
metadata <- read.table(metadata_file, header = TRUE, sep = "\t", comment.char = "")
speciesNames <- gsub(" ", "_", unique(metadata$scientific_name))

for (speciesName in speciesNames) {

  gtf_file <- list.files(path = gtf_dir, pattern = paste0(speciesName, ".*gtf$"), full.names = TRUE)
  ## there is potentially several gtf files already created.
  ## All gtf files created from the original Ensembl one have a prefix to specify what they correspond to
  ## e.g "wo_intergenic" or "nascent". In order to select the proper one we filter on the already known file
  ## names. It is definitely not bullet proof but the simplest solution to implement
  gtf_gz_file <- list.files(path = gtf_dir, pattern = paste0(speciesName, ".*gtf.gz$"), full.names = TRUE)
  gtf_gz_file <- gtf_gz_file[!grepl(pattern = "wo_intergenic.gtf.gz$|nascent.gtf.gz$", x = gtf_gz_file)]
  ## in case the gtf file is gzipped
  if (isTRUE(gtf_gz_file)) {
    gunzip(gtf_gz_file)
    gtf_file <- list.files(path = gtf_dir, pattern = paste0(speciesName, ".*gtf$"), full.names = TRUE)
  }
  gtf_file <- gtf_file[!grepl(pattern = "wo_intergenic.gtf$|nascent.gtf$", x = gtf_file)]
  
  ## reading in gtf file (gzipped, no need to uncompress)
  cat("Reading GTF file...\n")
  gene_gtf <- as.matrix(read.table(file = gtf_file, sep = "\t", strip.white = TRUE, as.is = TRUE,
    colClasses = "character", comment.char = '#'))

  ## selecting transcript lines
  gene_gtf_transcript <- gene_gtf[gene_gtf[,3] == "transcript",]

  #transform "transcript" to "exon"
  exon_gtf_nascent <- gene_gtf_transcript
  exon_gtf_nascent[,3] <- "exon"

  ## splitting the annotation field using single space "; " as a pattern
  split_annotation_list <- strsplit(exon_gtf_nascent[,9], "; ",  fixed  = T)
  #remove exon_id
  field_all <- sapply(split_annotation_list, function(x){ remove_annot_value(x, 'exon_id') })
  field_all <- sapply(field_all, function(x){ add_annot_value(x, 'exon_number 1', 'transcript_id') })
  field_all <- sapply(field_all, function(x){ update_transcript_id(x) })

  exon_gtf_nascent[,9] <- sapply(field_all, function(x){ paste(x,  collapse = "; ") })


  ## getting the vector of the gene IDs (1 for every exon)
  gene_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_id') })
  ## getting the vector of the transcript IDs (1 for every exon)
  transcript_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'transcript_id') })

  ## getting the table with mappings between gene IDs and transcript (for tximport)
  tx2gene_nascent_ids <- as.data.frame(unique(cbind(transcript_ids, gene_ids), MARGIN = 1))
  colnames(tx2gene_nascent_ids) <- c("TXNAME", "GENEID")
  tx2gene_nascent_ids$TXNAME <- paste0(tx2gene_nascent_ids$TXNAME, "_nascent")

  ## Output:
  ## GTF file with both genic exons and intergenic regions
  nascent_gtf_file <- file.path(gtf_dir, gsub(pattern = ".gtf", replacement = ".nascent.gtf",
    x = basename(gtf_file)))

  message("write file : ", nascent_gtf_file)
  write.table(x = exon_gtf_nascent,
              file = nascent_gtf_file,
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)

  ## Table of mapping between transcript_id and gene_id
  tx2gene_file <- file.path(gtf_dir, gsub(pattern = ".gtf", replacement = ".tx2gene",
    x = basename(gtf_file)))
  if (!file.exists(tx2gene_file)) {
    if(file.exists(paste0(tx2gene_file, ".xz"))) {
      system(sprintf("unxz %s", paste0(tx2gene_file, ".xz")))
    } else {
      stop("transcript to gene file [", tx2gene_file, "] not found.")
    }
  }
  tx2gene_ids <- read.table(file = tx2gene_file, sep = "\t", header = TRUE)
  tx2gene_single_nucleus_ids <- rbind(tx2gene_ids, tx2gene_nascent_ids)

  single_nucleus_tx2gene_file <- file.path(gtf_dir, gsub(pattern = ".gtf",
    replacement = ".single_nucleus.tx2gene", x = basename(gtf_file)))

  message("write file : ", single_nucleus_tx2gene_file)
  write.table(x = tx2gene_single_nucleus_ids,
              file = single_nucleus_tx2gene_file,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)

  # now that all files are created for the species we can gzip the original gtf file
  gzip(gtf_file)
}

