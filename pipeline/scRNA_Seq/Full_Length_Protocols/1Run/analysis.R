## SFonsecaCosta, Oct 2018
## This script is adapted from Julien Roux (Bgee RNA-Seq pipeline mainly from rna_seq_analysis.R)

## Usage:
## R CMD BATCH --no-save --no-restore '--args scrna_seq_sample_info="scrna_seq_sample_info.txt" cells_folder="cells_folder" infoFolder="infoFolder"' analysis.R analysis.Rout
## scrna_seq_sample_info --> annotation file that collect information about each cell/experiment
## cells_folder --> where we is located all the libraries/cells after Kallisto (treated data)
## infoFolder --> folder where is localized the transcriptome index + gene2transcript + gene2biotype

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}
## checking if all necessary arguments were passed....
command_arg <- c("scrna_seq_sample_info", "cells_folder", "infoFolder")
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
## Function to retrieve informative files
infoFiles <- function(cells_folder, libraryId, species){

  gene2transcript <- list.files(path=infoFolder, pattern = paste0("*_without_intergenic_gene2transcript.tsv$"), full.names=T, recursive = TRUE)
  gene2transcript_file <- grep(species,gene2transcript, value = TRUE)
  gene2transcript <- read.table(gene2transcript_file, h=F, sep="\t")
  names(gene2transcript) <- c("gene_id", "transcript_id")

  gene2biotype <- list.files(path=infoFolder, pattern = paste0("*_without_intergenic_gene2biotype.tsv$"), full.names=T, recursive = TRUE)
  gene2biotype_file <- grep(species,gene2biotype, value = TRUE)
  gene2biotype <- read.table(gene2biotype_file, h=F, sep="\t")
  names(gene2biotype) <- c("gene_id", "biotype")
  gene2biotype <- gene2biotype[order(gene2biotype$gene_id), ] ## order by gene Id

  ## reading abundance file from kallisto
  abundanceKallisto <- read.table(file.path(cells_folder, libraryId, "abundance.tsv"), h=T, sep="\t")
  return(list(gene2transcript,gene2biotype,abundanceKallisto))

}

## Function to normalize data (TPM)
countToTpm <- function(counts, effLen){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

## Function to recalculate TPM and Sum TranscriptID to Gene level
analysisData <- function(abundance, gene2transcript, gene2biotype){
  ## Add gene ids to the kallisto output. Effectively, this removes all intergenic regions,
  genic_count <- merge(abundance, gene2transcript, by.x=1, by.y=2)[, c(1,6,2,3,4,5)]
  names(genic_count)[2] <- "gene_id"
  genic_count <- merge(genic_count, gene2biotype, by.x=2, by.y=1)[, c(2,1,7,3,4,5,6)]
  genic_count <- genic_count[order(genic_count$target_id), ]
  ## Add gene id column to kallisto's output (genic regions only)
  abundance$gene_id <- rep(NA, times=length(abundance$target_id))
  abundance <- abundance[order(abundance$target_id),] ## sort by transcript_id
  ## Check transcripts are matching
  abundance$gene_id[abundance$target_id %in% genic_count$target_id] <- as.character(genic_count$gene_id)
  ## Add region type
  abundance$type <- rep("intergenic", times=length(abundance$target_id))
  abundance$type[abundance$target_id %in% genic_count$target_id] <- "genic"
  ## Add biotype for genic regions
  abundance$biotype <- rep(NA, times=length(abundance$target_id))
  abundance$biotype[abundance$target_id %in% genic_count$target_id] <- as.character(genic_count$biotype)
  ## reorder columns and rows before exporting
  abundance <- abundance[order(abundance$gene_id), c(1, 6, 2:5,7, 8)]

  ## Recalculate  using only genic regions!
  genic_count$tpm <- countToTpm(genic_count$est_counts, genic_count$eff_length)
  ## reorder columns and rows before exporting
  genic_count <- genic_count[order(genic_count$gene_id), c(1:2, 4:7, 3)]

  ## Gene-level expression
  ## Sum TPMs and counts of all transcripts of each gene
  gene_count <- aggregate(genic_count[,5:6], list(genic_count$gene_id), sum)
  names(gene_count)[1] <- "gene_id"
  gene_count$biotype <- gene2biotype$biotype
  ## Gene-level expression for kallisto's output:
  ## Nothing is summed for intergenic regions
  kallisto_gene_count <- abundance[abundance$type == "intergenic", c(1, 5, 6, 7, 8)]
  names(kallisto_gene_count)[1] <- "gene_id"
  ## For genic regions sum read counts, TPMs
  temp <- abundance[abundance$type == "genic", c(2, 5, 6)]
  temp <- aggregate(temp[,2:3], list(temp$gene_id), sum)
  names(temp)[1] <- "gene_id"
  temp$type <-  rep("genic", times=length(temp$gene_id))
  temp$biotype <- gene2biotype$biotype
  ## Make final table with both genic and intergenic regions
  kallisto_gene_count <- rbind(temp, kallisto_gene_count)
  ## return genic and intergenic regions at transcript level (abundance), gene level with intergenic (kallisto_gene_count) and just genic level (gene_count)
  return(list(abundance,kallisto_gene_count,gene_count))
}

############################################### APPLY PER LIBRARY ############################################################
## Apply for all libraries from different species
for(species in unique(annotation$organism)){
  for(libraryId in annotation$libraryId[annotation$organism == species]){

    cat("Species ", species ," and library ", libraryId, "\n")

    ## collect information
    cat("Collect information .... ", "\n")
    speciesID <- gsub(" ", "_", species)
    collectInfo <- infoFiles(cells_folder = cells_folder, libraryId = libraryId, species = speciesID)
    ## run rna_seq_analysis
    cat("Start analysis .... ", "\n")
    resultAnalysis <- analysisData(abundance = as.data.frame(collectInfo[3]), gene2transcript = as.data.frame(collectInfo[1]), gene2biotype = as.data.frame(collectInfo[2]))

    cat("Writing results ... ", "\n")
    ## result transcript level with intergenic information
    write.table(as.data.frame(resultAnalysis[1]), file = file.path(cells_folder, libraryId, "abundance+transcriptLevel+intergenic.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    ## result gene level with intergenic information
    write.table(as.data.frame(resultAnalysis[2]), file = file.path(cells_folder, libraryId, "abundance+geneLevel+intergenic.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    ## result gene level just genic regions re-calculate values TPM
    write.table(as.data.frame(resultAnalysis[3]), file = file.path(cells_folder, libraryId, "abundance+geneLevel+newTPM.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  }
}
