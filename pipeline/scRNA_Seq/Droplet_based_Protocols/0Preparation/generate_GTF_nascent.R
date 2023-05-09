# Julien Wollbrett
# created: 5/05/2023

## From GTF file downloaded from Ensembl, this script created a new GTF file containing nascent (not matured) transcripts: 
## - retrieve the start of the 1st exon and the stop of the last exon of the transcript
## - create a new transcript id by adding the prefix "_nascent" to the initial transcript ID e.g FBtr0306592_nascent
##   for the transcript FBtr0306592
## - create a gtf file with each line corresponding to a nascent transcript
## - also create a tsv file mapping the new transcript ID to the proper gene ID
##
## This gtf file is used only for single nucléus rnaseq libraries. It allows to map reads to not spliced intron or exon/intron regions

## Invoking:
# R CMD BATCH --no-save --no-restore --slave '--args gene_gtf_path = "gene_gtf_path" output_gtf_path = "output_gtf_path"' generate_GTF_nascent.R Rout_path

## Arguments to provide:
# "gene_gtf_path"   - full path to gtf file containing gene, transcript, and exon tags
# "output_gtf_dir" - full path to output folder
# Rout_path - path to  output .Rout file (optional)

cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("gene_gtf_path","output_dir")
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

## reading in gtf file (gzipped, no need to uncompress)
cat("Reading GTF file...\n")
gene_gtf <- as.matrix(read.table(file = gzfile(gene_gtf_path, "r"), sep = "\t", strip.white = TRUE, as.is = TRUE, colClasses = "character", comment.char = '#'))
## Header lines starting with # are not read

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
tx2gene_ids <- unique(cbind(transcript_ids, gene_ids), MARGIN = 1)

## Output:
## GTF file with both genic exons and intergenic regions
nascent_gtf_file <- file.path(output_dir, gsub(pattern = ".gtf", replacement = "_nascent.gtf", x = basename(gene_gtf_path)))

message("write file : ", nascent_gtf_file)
write.table(x = exon_gtf_nascent,
            file = nascent_gtf_file,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

## Table of mapping between transcript_id and gene_id
nascent_tx2gene_file <- file.path(output_dir, gsub(pattern = ".gtf", replacement = "_nascent_tx2gene.tsv", x = basename(gene_gtf_path)))

message("write file : ", nascent_tx2gene_file)
write.table(x = tx2gene_ids,
            file = nascent_tx2gene_file,
            sep = "\t",
            row.names = FALSE,
            col.names = c("TXNAME", "GENEID"),
            quote = FALSE)
