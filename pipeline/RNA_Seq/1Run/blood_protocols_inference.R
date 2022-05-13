## SFonsecaCosta Mars 2020

## This script is used to re-annotate blood libraries from RNA-Seq data 
## Mainly to separate libraries from two distinct protocols: with globin depletion or without globin depletion

## Usage:
## R CMD BATCH --no-save --no-restore '--args RNASeqLibrary="RNASeqLibrary.tsv" globin_file="globin_file.tsv" kallisto_count_folder="kallisto_count_folder" output="output"' blood_protocols_inference.R blood_protocols_inference.Rout
## RNASeqLibrary -> annotation file from Bgee
## globin_file -> file with all sub-units of hemoglobins per specie
## kallisto_count_folder -> path to kallisto result file
## output --> folder where the plots should be saved to verify a posteriori if needed

## Libraries used
library(data.table)
library(dplyr)
library(tools)
library(reshape2)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}
## checking if all necessary arguments were passed....
command_arg <- c("RNASeqLibrary", "globin_file", "kallisto_count_folder", "output")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}
#####################################################################################################################
## Reading input files: RNASeqLibrary and globin_file
if( file.exists(RNASeqLibrary) ){
  annotation <- read.table(RNASeqLibrary, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
} else {
  stop( paste("RNASeqLibrary.tsv file not found [", RNASeqLibrary, "]\n"))
}
if( file.exists(globin_file) ){
  globinIDs <- read.table(globin_file, h=T, sep="\t", comment.char="")
} else {
  stop( paste("globin_file.tsv file not found [", globin_file, "]\n"))
}
#####################################################################################################################
## function to collect information per library
infoPerLibraries <- function(libraries, globinUnits){
  
  ## retrieve information from abundance file at gene level for each library
  filePathLib <- file.path(kallisto_count_folder, libraries)
  AllFiles <- list.files(filePathLib,  pattern="abundance_gene_level\\+fpkm\\+intergenic.tsv$", full.names=T, recursive = TRUE)
  
  ## collect info to each globin subunit regarding the referent species
  collectInfo <- c()
  for (lib in AllFiles) {
    libraryId <- gsub(".*//\\s*|/abundance_gene_level\\+fpkm\\+intergenic.tsv.*", "", lib)
    librariesRead <- fread(lib)
    selectGlobins <- dplyr::filter(librariesRead, gene_id %in% globinUnits)
    selectGlobins <- selectGlobins[c(1,3)]
    ## collect libraryId + gene_id + TPM
    allInfo <- data.frame(libraryId, selectGlobins)
    collectInfo <- rbind(collectInfo,allInfo)
  }
  ## reshape the data
  reshapeTable <- reshape(collectInfo, direction='wide',idvar='libraryId', timevar='gene_id')
  reshapeTable$sum <- rowSums(reshapeTable[2:ncol(reshapeTable)])
  
  return(reshapeTable)
}

## apply cut-off and plot the data per species
cutoff_plot <- function(species, table_per_species, output){
  ## define threshold cut-off
  wgd <- 5e+05
  gd <- 2.5e+05
  table_per_species$globin_reduction_inference <- ifelse(table_per_species$sum < gd, "globin depletion", ifelse(table_per_species$sum >= wgd, "without globin depletion", "unsure"))
  
  ## plot per species
  pdf(file =  paste0(output, "/Sum_globin_subunits_per_library_",species,".pdf"), width = 8, height = 8)
  par(mgp=c(2,1,0),mar=c(3,4,3,2))
  table_per_species <- table_per_species[order(table_per_species$sum),] 
  plot(table_per_species$sum, pch=20, main=paste0("Globin subunits of the specie ",species),col=ifelse(table_per_species$sum < gd, "purple", ifelse(table_per_species$sum >= wgd, "darkblue", "gray")), 
       xlab="Number of libraries", ylab=expression(paste(sum(TPM), " globin subunits")), ylim=c(0,max(table_per_species$sum)), type = "o")
  abline(h=wgd, col="black", lty=2)
  abline(h=gd, col="black", lty=2)
  dev.off()
  
  ##export info per library
  allInfoCombined <- table_per_species[c(1,length(table_per_species))]
  return(allInfoCombined)
}

allInfoPerSpecies <- c()
for (species in unique(annotation$speciesId)) {
  ## specify just blood samples
  bloodAnnotation <- dplyr::filter(annotation, speciesId == species & uberonId == "UBERON:0000178")
  globinUnits <- globinIDs$ensembl_ID[globinIDs$speciesId == species]
  
  if((length(bloodAnnotation) != "0") & (length(globinUnits) != "0") ){
    libraries <- subset(annotation, speciesId %in% species)[,1]
    ## collect info
    subunitsInfo <- infoPerLibraries(libraries = libraries, globinUnits = globinUnits)
    ## apply cut-off and plot data
    cutoff <- cutoff_plot(species = species, table_per_species = subunitsInfo, output = output)
    allInfoPerSpecies <- rbind(allInfoPerSpecies, cutoff)
  } else {
    cat("Species without blood libraries or information about hemoglobin subunits.")
  }
}

## re-write annotation file with information about globin protocol
annotation <- merge(annotation, allInfoPerSpecies, by="libraryId", incomparables=NaN, all = TRUE)
write.table(annotation, file = file.path(output,"RNASeqLibrary_PostProcessing.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
