## SFonsecaCosta, October 2020

## This script is used make the calls of present and absent genes per library/cell-type population using pValue theoretical.
## This means: sum the UMI that belongs to the same cell-type population, normalize CPM and then call present and absent genes based on the pValue_theoretical cut-off.
## NOTE: Here we just use library/cell population that have more than 50 cells per library/cell-type pop and pass the bimodality test.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNASeq_Info.txt" bimodalityFile="bimodality_targetBased.txt" folder_data="folder_data" folder_refIntergenic="folder_refIntergenic" desired_pValue_cutoff="desired_pValue_cutoff" output_folder="output_folder"' Calls_cell_pop_per_library.R Calls_cell_pop_per_library.Rout
## scRNASeq_Info --> File that results from annotation and metadata (libraries downloaded and with extra information as readlength or SRR) 
## bimodalityFile --> File with information about cell populations that not pass the bimodality test
## folder_data --> Folder where are all the libraries after cell identification
## folder_refIntergenic --> Folder where is located the reference intergenic files for each species
## desired_pValue_cutoff --> desired pValue cutoff to call present genes
## output_folder --> Folder where we should save the results 

## libraries used
library(Biostrings)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggExtra)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("scRNASeq_Info", "bimodalityFile", "folder_data", "folder_refIntergenic" ,"desired_pValue_cutoff", "output_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read scRNASeq_Info file. If file not exists, script stops
if( file.exists(scRNASeq_Info) ){
  scRNASeqAnnotation <- read.table(scRNASeq_Info, h=T, sep="\t", comment.char="", quote = "")
} else {
  stop( paste("The scRNASeq information file was not found [", scRNASeq_Info, "]\n"))
}
## Read bimodality file. If file not exists, script stops
if( file.exists(bimodalityFile) ){
  bimodalityInfo <- read.table(bimodalityFile, h=T, sep="\t", comment.char="", quote = "")
} else {
  stop( paste("The cell information file was not found [", bimodalityFile, "]\n"))
}
##########################################################################################################################################################
## Sum the UMI of all barcodes and then compute the CPM normalization
sumUMICellPop <- function(folder_data, library, cellPop){
  
  cellPop <- fread(file.path(folder_data,library, paste0("Raw_Counts_",cellPop,".tsv")))
  cellPop$sumUMI <- rowSums(cellPop[ ,2:(length(cellPop)-4)])
  cellPop$CPM <- cellPop$sumUMI / sum(cellPop$sumUMI) * 1e6
  
  ## export cell pop info table
  cellPop <- data.frame(cellPop$gene_id, cellPop$sumUMI, cellPop$CPM, cellPop$type, cellPop$biotype)
  colnames(cellPop) <- c("gene_id", "sumUMI", "CPM", "type", "biotype")
  ## just re-order
  cellPopGenic <- data.frame(dplyr::filter(cellPop, type == "genic"))
  cellPopGenic <- cellPopGenic[order(cellPopGenic$gene_id),]
  cellPop <- rbind(cellPopGenic, dplyr::filter(cellPop, type == "intergenic"))
  return(cellPop)
}

## Provide the reference intergenic regions = TRUE
refIntergenic <- function(counts, folder_refIntergenic, speciesID){
  
  referenceIntergenic <- file.path(folder_refIntergenic, paste0(speciesID, "_intergenic.fa.gz"))
  if(!file.exists(referenceIntergenic)) {
    referenceIntergenic <- file.path(folder_refIntergenic, paste0(speciesID, "intergenic.fa"))
    if(!file.exists(referenceIntergenic)) {
      stop("reference intergenic file ", referenceIntergenic, "not available (compressed or not compressed")
    }
  }
  referenceIntergenic <- readDNAStringSet(referenceIntergenic)
  seq_name <- names(referenceIntergenic)
  seq_name <- gsub( " .*$", "", seq_name )
  seq_name <- gsub( "_", "-", seq_name)
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
theoretical_pValue <- function(counts, refrenceIntergenic){
  
  ## select all the intergenic region from the library
  intergenicRegionsLibrary <- dplyr::filter(counts, type == "intergenic")
  ## select just the true intergenic
  intergenicRegions <- dplyr::filter(refrenceIntergenic, refIntergenic == "TRUE")
  
  ## keep just information about reference intergenic region detected in the counts file to the calculation
  selected_Ref_Intergenic <- merge(intergenicRegionsLibrary, intergenicRegions, by="gene_id")
  
  ## select values with CPM > 0 (because we will use log2 scale)
  selected_Ref_Intergenic <- dplyr::filter(selected_Ref_Intergenic, CPM > 0 & type == "intergenic")
  
  ## select genic and intergenic region from the library with CPM > 0
  regions <- dplyr::filter(counts, CPM > 0)
  ## calculate z-score for each gene_id using the reference intergenic
  regions$zScore <- (log2(regions$CPM) - mean(log2(selected_Ref_Intergenic$CPM))) / sd(log2(selected_Ref_Intergenic$CPM))
  ## calculate p-values for each gene_id
  regions$pValue <- pnorm(regions$zScore, lower.tail = FALSE)
  
  return(list(regions, 2^(mean(log2(selected_Ref_Intergenic$CPM))), 2^(sd(log2(selected_Ref_Intergenic$CPM)))))
}

cutoff_info <- function(library_id, cellTypeName, cellTypeId, counts, desired_pValue_cutoff, meanRefIntergenic, sdRefIntergenic){
  
  ## collect stats using pValue_cutoff
  genic_region_present <- nrow(dplyr::filter(counts, type == "genic" & calls_pValue == "present"))
  proportion_genic_present <- (nrow(dplyr::filter(counts, type == "genic" & calls_pValue == "present"))/nrow(dplyr::filter(counts, type == "genic")))*100
  coding_region_present  <- nrow(dplyr::filter(counts, biotype == "protein_coding" & calls_pValue == "present"))
  proportion_coding_present <- (nrow(dplyr::filter(counts, biotype == "protein_coding" & calls_pValue == "present"))/nrow(dplyr::filter(counts, biotype == "protein_coding")))*100
  intergenic_region_present <- nrow(dplyr::filter(counts, type == "intergenic" & calls_pValue == "present"))
  proportion_intergenic_present <- (nrow(dplyr::filter(counts, type == "intergenic" & calls_pValue == "present"))/nrow(dplyr::filter(counts, type == "intergenic")))*100
  
  CPM_Threshold <- min(counts$CPM[counts$type == "genic" & counts$calls_pValue == "present"])
  ## Export cutoff_info_file
  collectInfo <- c(library_id, cellTypeName, cellTypeId, CPM_Threshold, sum(counts$type == "genic"), genic_region_present,proportion_genic_present, 
                   sum(counts$biotype %in% "protein_coding"),coding_region_present,proportion_coding_present, 
                   sum(counts$type %in% "intergenic"),intergenic_region_present,proportion_intergenic_present, 
                   desired_pValue_cutoff, meanRefIntergenic, sdRefIntergenic)
  names(collectInfo) <- c("libraryId", "cellTypeName", "cellTypeId", "CPM_Threshold", "Genic","Genic_region_present","Proportion_genic_present", 
                          "Protein_coding","Coding_region_present","Proportion_coding_present", 
                          "Intergenic","Intergenic_region_present","Proportion_intergenic_present", 
                          "pValue_cutoff", "meanRefIntergenic", "sdRefIntergenic")
  return(collectInfo)
}


## export the plot
plotData <- function(libraryId, cellPopName, counts, refIntergenic, CPM_threshold){
  ## export distribution
  dens <- density(log2(counts$CPM+1e-6), na.rm=T)
  dens_genic <- density(log2(counts$CPM[counts$type == "genic"]+1e-6), na.rm=T)
  dens_genic$y <- dens_genic$y * nrow(dplyr::filter(counts, type == "genic")) / length(counts$CPM)
  dens_coding <- density(log2(counts$CPM[counts$biotype == "protein_coding"]+1e-6), na.rm=T)
  dens_coding$y <- dens_coding$y * nrow(dplyr::filter(counts, biotype == "protein_coding")) / length(counts$CPM)
  dens_intergenic <- density(log2(counts$CPM[counts$type == "intergenic"]+1e-6), na.rm=T)
  dens_intergenic$y <- dens_intergenic$y * nrow(dplyr::filter(counts, type == "intergenic")) / length(counts$CPM)
  refIntergenic <- merge(dplyr::filter(counts, type == "intergenic"), referenceIntergenic, by = "gene_id")
  refIntergenic <- as.data.table(dplyr::filter(refIntergenic, refIntergenic == "TRUE"))
  dens_Refintergenic <- density(log2(refIntergenic$CPM+1e-6), na.rm=T)
  dens_Refintergenic$y <- dens_Refintergenic$y * nrow(refIntergenic) / length(counts$CPM)
  if(!dir.exists(file.path(output_folder, libraryId))) {
    dir.create(file.path(output_folder,libraryId))
  }
  pdf(file.path(output_folder, libraryId, paste0("Distribution_", libraryId, "_",cellPopName, ".pdf")), width=10, height=6) 
  par(mfrow=c(1,2))
  plot(dens, lwd=2, main=paste0(libraryId), xlab="Log2(CPM)")
  mtext(paste0(cellPopName))
  lines(dens_genic,col="red", lwd=2)
  lines(dens_coding,col="indianred", lwd=2)
  lines(dens_intergenic,col="darkblue", lwd=2)
  lines(dens_Refintergenic,col="cyan", lwd=2)
  abline(v=CPM_threshold, col="gray", lty=2, lwd=2)
  legend("topright", legend = c("All", "Genic" ,"PC", "Int", "RefInt", "cutoff"), col=c("black", "red", "indianred", "darkblue", "cyan", "gray"),lty=c(1,1,1,1,1,2), lwd=2, bty = "n")
  
  ## export frequency of pValue for all genic region
  genicRegion <- as.numeric(counts$pValue[counts$type == "genic"])
  hist(na.omit(genicRegion), main=paste0(libraryId), xlab="pValue", xlim=c(0,1))
  mtext(paste0("Genic region_",cellPopName))
  abline(v=desired_pValue_cutoff, col="red", lwd=2)
  dev.off()
}

## export info stats of all libraries/cell-population
if(!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}
All_samples <- file.path(output_folder, "All_cellPopulation_stats_10X.tsv")
if (!file.exists(All_samples)){
  file.create(All_samples)
  cat("libraryId\tcellTypeName\tCPM_Threshold\tGenic\tGenic_region_present\tProportion_genic_present\tProtein_coding","Coding_region_present\tProportion_coding_present\tIntergenic\tIntergenic_region_present\tProportion_intergenic_present\tpValue_cutoff\tmeanRefIntergenic\tsdRefIntergenic\tspecies\torganism\n",file = file.path(output_folder,"All_cellPopulation_stats_10X.tsv"), sep = "\t")
} else {
  print("File already exist.....")
}

collectSamplesStats <- data.frame()
## loop thought all libraries
for (libraryid in unique(scRNASeqAnnotation$libraryId)) {
  
    message("Library: ", libraryid)
  
    ## collect all cell populations that belongs to the library
    AllCellPop <- list.files(path = file.path(folder_data, libraryid), pattern = "^Raw_Counts_")
    collectCellNames <- str_remove(AllCellPop, "Raw_Counts_")
    collectCellNames <- str_remove(collectCellNames, ".tsv")
    ## remove from the analysis cell-type populations that not pass the bimodality quality control (means < 50 cells or are not bimodal distribution)
    cellsInfo <- bimodalityInfo$Cell_Name[bimodalityInfo$library == libraryid]
    
    ## select just cell pop that pass the quality control 
    AllCellPop <- setdiff(collectCellNames,cellsInfo)
    speciesID <- scRNASeqAnnotation$speciesId[scRNASeqAnnotation$libraryId == libraryid]
   
    message("AllCellPop: ", AllCellPop)
    message("speciesID: ", speciesID)
 
    for (cellPop in AllCellPop) {
      
      message("Doing cell population: ", cellPop)
      
        ## collect cellName and cellId
        cellName <- sub("\\_.*", "", cellPop)
        cellID <- sub(".*_", "", cellPop)

        ## collect the sumUMI + normalization for the target cellPop
        cellPop_normalized <- sumUMICellPop(folder_data = folder_data, library = libraryid, cellPop = cellPop)
        
        ## Information about reference intergenic
        referenceIntergenic <- refIntergenic(counts = cellPop_normalized, folder_refIntergenic = folder_refIntergenic, speciesID = speciesID)
        
        ## calls with pValue theoretical
        calculatePvalues <- theoretical_pValue(counts = cellPop_normalized, refrenceIntergenic = referenceIntergenic)
        calculationInfo <- calculatePvalues[[1]]
        calculationInfo$calls_pValue <- ifelse(calculationInfo$pValue <= as.numeric(desired_pValue_cutoff), "present", "absent" )
        
        ## add info also about genesID where CPM = 0 and were not used for the pValue calculation (but are important for the final stats)
        regionZero <- cellPop_normalized[ !cellPop_normalized$gene_id %in% calculationInfo$gene_id, ]
        regionZero$zScore <- "NA"; regionZero$pValue <- "NA"; regionZero$calls_pValue <-"absent"
        
        allData <- rbind(calculationInfo, regionZero)
        genicRegion <- dplyr::filter(allData, type=="genic")
        genicRegion <- genicRegion[order(genicRegion$gene_id),]
        finalData <- rbind(genicRegion, dplyr::filter(allData, type == "intergenic"))
        finalData$cellTypeName <- cellName
        finalData$cellTypeId <- gsub("-",":",cellID)
        
        ## collect just genic region and re-calculate CPM
        finalData_genic <- dplyr::filter(finalData, type == "genic")
        finalData_genic$CPM <- finalData_genic$sumUMI / sum(finalData_genic$sumUMI) * 1e6
        
        ## Export cutoff information file + new files with calls
        cutoff_info_file <- cutoff_info(libraryid, cellTypeName = cellName, cellTypeId = gsub("-",":",cellID), counts = finalData, desired_pValue_cutoff = as.numeric(desired_pValue_cutoff), meanRefIntergenic = calculatePvalues[[2]], sdRefIntergenic = calculatePvalues[[3]])
        CPM_threshold <- log2(as.numeric(cutoff_info_file[3]))
        
        ## export data
        plotData(libraryId=libraryid, cellPopName =cellPop, counts=finalData, refIntergenic=referenceIntergenic, CPM_threshold=CPM_threshold)
        
        pathExport <- file.path(output_folder, libraryid)
        write.table(finalData,file = file.path(pathExport, paste0("Calls_cellPop_",libraryid, "_",cellPop,"_genic+intergenic.tsv")),quote=F, sep = "\t", col.names=T, row.names=F)
        write.table(finalData_genic,file = file.path(pathExport, paste0("Calls_cellPop_",libraryid, "_",cellPop,"_genic.tsv")),quote=F, sep = "\t", col.names=T, row.names=F)
        write.table(t(t(cutoff_info_file)),file = file.path(pathExport, paste0("cutoff_info_file_",libraryid, "_",cellPop,".tsv")),quote=F, sep = "\t", col.names=F, row.names=T)
        
        ## add this to big data frame with all samples information
        this_sample <- as.data.frame(t(cutoff_info_file), stringsAsFactors=F)
        this_sample[, "species"]  <- speciesID
        this_sample[, "organism"] <- as.character(unique(scRNASeqAnnotation$scientific_name[scRNASeqAnnotation$speciesId == speciesID]))
        collectSamplesStats <- rbind(collectSamplesStats, this_sample)
    }
}
write.table(collectSamplesStats, file = file.path(output_folder, "All_cellPopulation_stats_10X.tsv"),col.names =F , row.names = F, append = T,quote = FALSE, sep = "\t")

## final plot per species
pdf(file.path(output_folder, paste0("All_libraries_stats_information_10X.pdf")), width=16, height=6) 
g1 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(Proportion_genic_present))) + 
  geom_boxplot()+ylim(0,100)+xlab(" ")+ylab("% Genic Present")
g2 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(Proportion_coding_present))) + 
  geom_boxplot(notch=TRUE)+ylim(0,100)+xlab(" ")+ylab("% Protein Coding Present")
g3 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(Proportion_intergenic_present))) + 
  geom_boxplot(notch=TRUE)+ylim(0,100)+xlab(" ")+ylab("% Intergenic Present")
g4 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(meanRefIntergenic))) + 
  geom_boxplot(notch=TRUE)+ylim(0,1)+xlab(" ")+ylab("Mean Reference Intergenic")
g5 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(sdRefIntergenic))) + 
  geom_boxplot(notch=TRUE)+ylim(0,10)+xlab(" ")+ylab("SD Reference Intergenic")
g6 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(pValue_cutoff))) + 
  geom_boxplot(notch=TRUE)+ylim(0,1)+xlab(" ")+ylab("pValue_cutoff")
grid.arrange(g1,g2,g3,g4,g5,g6, nrow = 2, ncol = 3)
dev.off()
