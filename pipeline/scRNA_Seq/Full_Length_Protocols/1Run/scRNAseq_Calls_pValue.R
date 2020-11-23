## SFonsecaCosta, Oct 2020
## This script is used to calculate pValue threoretical for each gene ID using the reference intergenic regions 
## to perform the calls a pValue_cutoff should be specified
## This is applied per libraryId (means per individual cell)

## Usage:
## R CMD BATCH --no-save --no-restore '--args NEW_scRNASeq_sample_info="NEW_scRNASeq_sample_info" cells_folder="cells_folder" refIntergenic_folder="refIntergenic_folder" desired_pValue_cutoff="desired_pValue_cutoff" output_folder="output_folder"' scRNAseq_Calls_pValue.R scRNAseq_Calls_pValue.Rout
## NEW_scRNASeq_sample_info --> file with info on mapped libraries
## cells_folder --> where we is located all the libraries/cells after Kallisto (treated data)
## refIntergenic_folder --> folder where is located the reference intergenic regions of all species (.fa files)
## desired_pValue_cutoff --> desired p-value cutoff to call present genes
## output_folder --> folder where to export the plots and information about the libraries

## libraries used
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Biostrings)

## reading in arguments provided in command line
cmd_args = commandArgs(TRUE);

print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}
## checking if all necessary arguments were passed in command line
command_arg <- c("NEW_scRNASeq_sample_info", "cells_folder", "refIntergenic_folder", "desired_pValue_cutoff" ,"output_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}
## Read NEW_scRNASeq_sample_info file. If file not exists, script stops
if( file.exists(NEW_scRNASeq_sample_info) ){
  annotation <- read.table(NEW_scRNASeq_sample_info, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
} else {
  stop( paste("NEW_scRNASeq_sample_info file not found [", NEW_scRNASeq_sample_info, "]\n"))
}

###################################################################################################################
## Provide the reference intergenic regions = TRUE
refIntergenic <- function(counts, folder_refIntergenic, speciesID){
  
  referenceIntergenic <- paste0(folder_refIntergenic, "/", speciesID, "_intergenic.fa")
  referenceIntergenic <- readDNAStringSet(referenceIntergenic)
  seq_name <- names(referenceIntergenic)
  seq_name <- gsub( " .*$", "", seq_name )
  ##seq_name <- gsub( "_", "-", seq_name)
  seqNamesFinal <- as.data.frame(seq_name)
  colnames(seqNamesFinal) <- "gene_id"
  seqNamesFinal$refIntergenic <- "TRUE"
  ## seqNamesFinal with same size that counts Table (here everything different refIntergenic is FALSE)
  collectIDs <- as.data.frame(counts$gene_id)
  colnames(collectIDs) <- "gene_id"
  seqNamesFinal <- merge(collectIDs, seqNamesFinal, by = "gene_id", all.x = TRUE)
  seqNamesFinal$refIntergenic <- ifelse(is.na(seqNamesFinal$refIntergenic) == TRUE, "FALSE", seqNamesFinal$refIntergenic)
  
  return(seqNamesFinal)
}

## function to calculate pValue from the theoretical data
theoretical_pValue <- function(counts, referenceIntergenic){
  
  selected_Ref_Intergenic <- dplyr::filter(referenceIntergenic, refIntergenic == "TRUE")
  selected_Ref_Intergenic <- merge(counts, selected_Ref_Intergenic, by = "gene_id")
  selected_Ref_Intergenic <- dplyr::filter(selected_Ref_Intergenic, type == "intergenic" & tpm > 0)
  
  ## select genic and intergenic region from the library
  regions <- dplyr::filter(counts, tpm > 0)
  ## calculate z-score for each ID using the reference intergenic
  regions$zScore <- (log2(regions$tpm) - mean(log2(selected_Ref_Intergenic$tpm))) / sd(log2(selected_Ref_Intergenic$tpm))
  ## calculate p-values for each ID
  regions$pValue <- pnorm(regions$zScore, lower.tail = FALSE)
  
  return(list(regions, 2^(mean(log2(selected_Ref_Intergenic$tpm))), 2^(sd(log2(selected_Ref_Intergenic$tpm)))))
}

### Collect information after cut-off
cutoff_info <- function(library_id, counts, desired_pValue_cutoff, meanRefIntergenic, sdRefIntergenic){
  
  ## collect stats using pValue_cutoff
  genic_region_present <- nrow(dplyr::filter(counts, type == "genic" & calls_pValue == "present"))
  proportion_genic_present <- (nrow(dplyr::filter(counts, type == "genic" & calls_pValue == "present"))/nrow(dplyr::filter(counts, type == "genic")))*100
  coding_region_present  <- nrow(dplyr::filter(counts, biotype == "protein_coding" & calls_pValue == "present"))
  proportion_coding_present <- (nrow(dplyr::filter(counts, biotype == "protein_coding" & calls_pValue == "present"))/nrow(dplyr::filter(counts, biotype == "protein_coding")))*100
  intergenic_region_present <- nrow(dplyr::filter(counts, type == "intergenic" & calls_pValue == "present"))
  proportion_intergenic_present <- (nrow(dplyr::filter(counts, type == "intergenic" & calls_pValue == "present"))/nrow(dplyr::filter(counts, type == "intergenic")))*100
  
  TPM_Threshold <- min(counts$tpm[counts$type == "genic" & counts$calls_pValue == "present"])
  ## Export cutoff_info_file
  collectInfo <- c(library_id,TPM_Threshold, sum(counts$type == "genic"), genic_region_present,proportion_genic_present, 
                   sum(counts$biotype %in% "protein_coding"),coding_region_present,proportion_coding_present, 
                   sum(counts$type %in% "intergenic"),intergenic_region_present,proportion_intergenic_present, 
                   desired_pValue_cutoff, meanRefIntergenic, sdRefIntergenic)
  names(collectInfo) <- c("libraryId","TPM_Threshold","Genic","Genic_region_present","Proportion_genic_present", 
                          "Protein_coding","Coding_region_present","Proportion_coding_present", 
                          "Intergenic","Intergenic_region_present","Proportion_intergenic_present", 
                          "pValue_cutoff", "meanRefIntergenic", "sdRefIntergenic")
  return(collectInfo)
}

## export the plot
plotData <- function(libraryId, kallisto_gene_counts, refIntergenic, TMP_threshold){
  ## export distribution
  dens <- density(log2(kallisto_gene_counts$tpm+1e-6), na.rm=T)
  dens_genic <- density(log2(kallisto_gene_counts$tpm[kallisto_gene_counts$type == "genic"]+1e-6), na.rm=T)
  dens_genic$y <- dens_genic$y * nrow(dplyr::filter(kallisto_gene_counts, type == "genic")) / length(kallisto_gene_counts$tpm)
  dens_coding <- density(log2(kallisto_gene_counts$tpm[kallisto_gene_counts$biotype == "protein_coding"]+1e-6), na.rm=T)
  dens_coding$y <- dens_coding$y * nrow(dplyr::filter(kallisto_gene_counts, biotype == "protein_coding")) / length(kallisto_gene_counts$tpm)
  dens_intergenic <- density(log2(kallisto_gene_counts$tpm[kallisto_gene_counts$type == "intergenic"]+1e-6), na.rm=T)
  dens_intergenic$y <- dens_intergenic$y * nrow(dplyr::filter(kallisto_gene_counts, type == "intergenic")) / length(kallisto_gene_counts$tpm)
  refIntergenic <- merge(dplyr::filter(kallisto_gene_counts, type == "intergenic"), referenceIntergenic, by = "gene_id")
  refIntergenic <- dplyr::filter(refIntergenic, refIntergenic == "TRUE")
  dens_Refintergenic <- density(log2(refIntergenic$tpm+1e-6), na.rm=T)
  dens_Refintergenic$y <- dens_Refintergenic$y * nrow(refIntergenic) / length(kallisto_gene_counts$tpm)
  
  pdf(file.path(output_folder, libraryId, paste0("Distribution_", libraryId, ".pdf")), width=10, height=6)
  par(mfrow=c(1,2))
  plot(dens, lwd=2, main=paste0(libraryId), xlab="Log2(TPM)")
  lines(dens_genic,col="red", lwd=2)
  lines(dens_coding,col="indianred", lwd=2)
  lines(dens_intergenic,col="darkblue", lwd=2)
  lines(dens_Refintergenic,col="cyan", lwd=2)
  abline(v=TMP_threshold, col="gray", lty=2, lwd=2)
  legend("topright", legend = c("All", "Genic", "PC", "Int", "RefInt", "cutoff"), col=c("black", "red" ,"indianred", "darkblue", "cyan", "gray"),lty=c(1,1,1,1,1,2), bty = "n")
  
  ## export frequency of pValue
  typeGenic <- as.numeric(kallisto_gene_counts$pValue[kallisto_gene_counts$type == "genic"])
  hist(na.omit(typeGenic), main=paste0(libraryId), xlab="pValue", xlim=c(0,1))
  mtext(paste0("Genic region"))
  abline(v=desired_pValue_cutoff, col="red", lwd=2)
  dev.off()
}

## export info of all libraries
All_samples <- paste0(output_folder, "/All_samples_stats_FL.tsv")
if (!file.exists(All_samples)){
  file.create(All_samples)
  cat("libraryId\tTPM_Threshold\tGenic\tGenic_region_present\tProportion_genic_present\tProtein_coding\tCoding_region_present\tProportion_coding_present\tIntergenic\tIntergenic_region_present\tProportion_intergenic_present\tpValue_cutoff\tmeanRefIntergenic\tsdRefIntergenic\tspecies\torganism\n",file = file.path(output_folder,"All_samples_stats_FL.tsv"), sep = "\t")
} else {
  print("File already exist.....")
}

collectSamples <- data.frame()
for(species in unique(annotation$speciesId)){
  print(species)
  for(libraryId in annotation$libraryId[annotation$speciesId == species]){
    
    cat("Doing the library: ",libraryId, "\n")
    ## genic + intergenic
    kallisto_gene_counts <- read.table(file.path(cells_folder, libraryId, "abundance+geneLevel+intergenic.tsv"), h=T, sep="\t")
    ## retrieve reference intergenic
    referenceIntergenic <- refIntergenic(counts = kallisto_gene_counts, folder_refIntergenic=refIntergenic_folder, speciesID=species)
    
    ## pValue calculation and call to kallisto_gene_counts
    calculatePvalue <- theoretical_pValue(counts = kallisto_gene_counts, referenceIntergenic = referenceIntergenic)
    
    pValueCalculation <- calculatePvalue[[1]]
    pValueCalculation$calls_pValue <- ifelse(pValueCalculation$pValue <= as.numeric(desired_pValue_cutoff), "present", "absent" )
    ## add info also about genesID where TPM = 0 and were not used for the pValue calculation (but are important for the final stats)
    regionZero <- kallisto_gene_counts[ !kallisto_gene_counts$gene_id %in% pValueCalculation$gene_id, ]
    regionZero$zScore <- "NA"; regionZero$pValue <- "NA"; regionZero$calls_pValue <-"absent"
    
    allData <- rbind(pValueCalculation, regionZero)
    genicRegion <- dplyr::filter(allData, type=="genic")
    genicRegion <- genicRegion[order(genicRegion$gene_id),]
    kallisto_gene_counts <- rbind(genicRegion, dplyr::filter(allData, type == "intergenic"))
    
    ## Export cutoff information file + new files with calls
    cutoff_info_file <- cutoff_info(libraryId, kallisto_gene_counts, desired_pValue_cutoff = as.numeric(desired_pValue_cutoff), meanRefIntergenic = calculatePvalue[[2]], sdRefIntergenic = calculatePvalue[[3]])
    TMP_threshold <- log2(as.numeric(cutoff_info_file[2]))
    
    ## export data
    plotData(libraryId=libraryId, kallisto_gene_counts=kallisto_gene_counts, refIntergenic=referenceIntergenic, TMP_threshold=TMP_threshold)
    
    pathExport <- file.path(cells_folder, libraryId)
    write.table(kallisto_gene_counts,file = file.path(pathExport, "abundance+geneLevel+intergenic+calls.tsv"),quote=F, sep = "\t", col.names=T, row.names=F)
    write.table(t(t(cutoff_info_file)),file = file.path(pathExport, "cutoff_info_file.tsv"),quote=F, sep = "\t", col.names=F, row.names=T)
    
    ## add this to big data frame with all samples
    this_sample <- as.data.frame(t(cutoff_info_file), stringsAsFactors=F)
    this_sample[, "species"]  <- species
    this_sample[, "organism"] <- as.character(unique(annotation$organism[annotation$speciesId == species]))
    collectSamples <- rbind(collectSamples, this_sample)
  }
  }

write.table(collectSamples, file = file.path(output_folder, "All_samples_stats_FL.tsv"),col.names =F , row.names = F, append = T,quote = FALSE, sep = "\t")
## final plot per species
pdf(file.path(output_folder, paste0("All_libraries_information.pdf")), width=16, height=6) 
g1 <- ggplot(collectSamples, aes(x=organism, y=as.numeric(Proportion_genic_present))) + 
  geom_boxplot(notch=TRUE)+ylim(0,100)+xlab(" ")+ylab("% Genic Present") 
g2 <- ggplot(collectSamples, aes(x=organism, y=as.numeric(Proportion_coding_present))) + 
  geom_boxplot(notch=TRUE)+ylim(0,100)+xlab(" ")+ylab("% Protein Coding Present")
g3 <- ggplot(collectSamples, aes(x=organism, y=as.numeric(Proportion_intergenic_present))) + 
  geom_boxplot(notch=TRUE)+ylim(0,100)+xlab(" ")+ylab("% Intergenic Present")
g4 <- ggplot(collectSamples, aes(x=organism, y=as.numeric(meanRefIntergenic))) + 
  geom_boxplot(notch=TRUE)+ylim(0,1)+xlab(" ")+ylab("Mean Reference Intergenic")
g5 <- ggplot(collectSamples, aes(x=organism, y=as.numeric(sdRefIntergenic))) + 
  geom_boxplot(notch=TRUE)+ylim(0,10)+xlab(" ")+ylab("SD Reference Intergenic")
g6 <- ggplot(collectSamples, aes(x=organism, y=as.numeric(pValue_cutoff))) + 
  geom_boxplot(notch=TRUE)+ylim(0,1)+xlab(" ")+ylab("pValue_cutoff")
grid.arrange(g1,g2,g3,g4,g5,g6, nrow = 2, ncol = 3)
dev.off()
