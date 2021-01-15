## SFonsecaCosta, May 2020

## This script is used to create a transcript_to_gene file from the gtf_all file that will be used later in bustools.
## Also a file gene_to_biotype is created for each species.
## NOTE: if the index contain the intergenic regions the file transcript_to_gene should also contain all information: genic + intergenic

## Usage:
## R CMD BATCH --no-save --no-restore '--args folder_gtf="folder_gtf"' generateInfo.R generateInfo.Rout
## folder_gtf --> Folder where are all gtf files (gtf_all with intergenic) are for the different species

## libraries used
library(data.table)
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

## checking if all necessary arguments were passed....
command_arg <- c("folder_gtf")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

AllFiles <- list.files(folder_gtf, pattern="*gtf_all", full.names=T, recursive = TRUE)

for (i in AllFiles) {

  speciesID <- basename(i)
  speciesID <- data.frame(str_split(speciesID, "\\."))[1,1]

  readFile <- fread(i)

  ## extract information: gene_id, transcript_id, biotype_id
  system(sprintf('grep %s %s | sed -e "s/^.*gene_id //; s/;.*$//" > %s',paste0('gene_id '), paste0(i), paste0(folder_gtf,"/gene_id_", speciesID, ".tsv")))
  system(sprintf('grep %s %s | sed -e "s/^.*transcript_id //; s/;.*$//" > %s',paste0('transcript_id '), paste0(i), paste0(folder_gtf,"/transcript_id_", speciesID, ".tsv")))
  system(sprintf('grep %s %s | sed -e "s/^.*gene_biotype //; s/;.*$//" > %s',paste0('gene_biotype '), paste0(i), paste0(folder_gtf,"/biotype_id_", speciesID, ".tsv")))

  gene_id <- fread(paste0(folder_gtf,"/gene_id_", speciesID, ".tsv"), header = FALSE)
  transcript_id <- fread(paste0(folder_gtf,"/transcript_id_", speciesID, ".tsv"), header = FALSE)
  biotype <- fread(paste0(folder_gtf,"/biotype_id_", speciesID, ".tsv"), header = FALSE)

  ## collect transcript_to_gene
  transcript_to_gene <- cbind(transcript_id, gene_id)
  colnames(transcript_to_gene) <- c("transcript", "gene")
  transcript_to_gene <- transcript_to_gene[!duplicated(transcript_to_gene[,c('transcript')]),]
  ## collect gene_to_biotype
  gene_to_biotype <- seq(max(nrow(gene_id), nrow(biotype)))
  gene_to_biotype <- data.frame(gene_id[gene_to_biotype], biotype[gene_to_biotype])
  colnames(gene_to_biotype) <- c("gene", "biotype")
  gene_to_biotype <- gene_to_biotype[!duplicated(gene_to_biotype[,c('gene')]),]

  write.table(transcript_to_gene, file = paste0(folder_gtf, "/transcript_to_gene_with_intergenic_", speciesID,".tsv"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(gene_to_biotype, file = paste0(folder_gtf, "/gene_to_biotype_with_intergenic_", speciesID,".tsv"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

  ## remove intermediary files
  file.remove(paste0(folder_gtf, "/gene_id_", speciesID, ".tsv"))
  file.remove(paste0(folder_gtf, "/transcript_id_", speciesID, ".tsv"))
  file.remove(paste0(folder_gtf, "/biotype_id_", speciesID, ".tsv"))

}
