## SFonsecaCosta, October 2020

## This script is used make the calls of present and absent genes per library/cell-type
## population using pValue theoretical.
## This means: sum the UMI that belongs to the same cell-type population, normalize CPM and then
## call present and absent genes based on the pValue_theoretical cut-off.

## NOTE: Since Bgee 15.1 do not filter anymore library/cell population that have more than 50 cells
## per library/cell-type pop and pass the bimodality test.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNASeq_Info.txt" folder_data="folder_data" folder_gtf="folder_gtf" desired_pValue_cutoff="desired_pValue_cutoff" output_folder="output_folder"' Calls_cell_pop_per_library.R Calls_cell_pop_per_library.Rout
## scRNASeq_Info --> File that results from annotation and metadata (libraries downloaded and with extra information as readlength or SRR) 
## folder_data --> Folder where are all the libraries after cell identification
## folder_gtf --> Folder where is located the reference intergenic files for each species
## desired_pValue_cutoff --> desired pValue cutoff to call present genes
## output_folder --> Folder where we should save the results

## libraries used
library(Biostrings)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(ggExtra)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){
  stop("no arguments provided\n")
} else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("scRNASeq_Info", "folder_data", "folder_gtf" ,"desired_pValue_cutoff",
  "output_folder")
for (c_arg in command_arg) {
  if (!exists(c_arg)) {
    stop(paste(c_arg,"command line argument not provided\n"))
  }
}

## Read scRNASeq_Info file. If file not exists, script stops
if (file.exists(scRNASeq_Info)) {
  scRNASeqAnnotation <- read.table(scRNASeq_Info, h=T, sep="\t", comment.char="", quote = "\"")
} else {
  stop(paste("The scRNASeq information file was not found [", scRNASeq_Info, "]\n"))
}
##########################################################################################################################################################
## Provide the reference intergenic regions = TRUE
##TODO Why bother retrieving ref. intergenic??? the ampping/quantification should already be donne
## on the reference intergenic.
##TODO: check that only ref intergenic were used for each step of the pipeline( e.g gene annotation, kallisto index, ....)
refIntergenic <- function(counts, folderGtf, speciesName){
  gene2biotype_file <- list.files(path = folderGtf, pattern = paste0(speciesName, ".*gene2biotype"), full.names = TRUE)
  gene2biotype <- read.table(gene2biotype_file, sep = "\t", header = TRUE)
  gene2biotype$refIntergenic <- ifelse(is.na(gene2biotype$biotype), FALSE, TRUE)
  gene2biotype$gene_id <- ifelse(is.na(gene2biotype$biotype), gsub( "_", "-", gene2biotype$id),
    gene2biotype$id)
  ref_intergenic <- gene2biotype[,c("id","refIntergenic")]
  return(ref_intergenic)
}

## Sum the UMI of all barcodes and then compute the CPM normalization
sumUMICellPop <- function(rawCountFile) {
  cellPop <- read.table(rawCountFile, header = TRUE, sep = "\t")
  # sum counts for all barcodes.
  ## it is possible to have only one barcode per celltype. It is then not possible to
  ## use the function rowSums.
  ## In that case the sumUMI should be the value of the number of UMI for the only
  ## available column.
  ##TODO find a most elegant way to detect number of barcode per celltype
  ##     maybe check the file barcode to UMI and count from that file.
  if (ncol(cellPop) == 6) {
    cellPop$sumUMI <- cellPop[ ,2]
  } else {
    cellPop$sumUMI <- rowSums(cellPop[ ,2:(length(cellPop)-4)])
  }
  cellPop$CPM <- cellPop$sumUMI / sum(cellPop$sumUMI) * 1e6
  
  ## export cell pop info table
  cellPop <- data.frame(cellPop$gene_id, cellPop$sumUMI, cellPop$CPM, cellPop$type,
                        cellPop$biotype)
  colnames(cellPop) <- c("gene_id", "sumUMI", "CPM", "type", "biotype")
  ## just re-order
  cellPopGenic <- data.frame(dplyr::filter(cellPop, type == "genic"))
  cellPopGenic <- cellPopGenic[order(cellPopGenic$gene_id),]
  cellPop <- rbind(cellPopGenic, dplyr::filter(cellPop, type == "intergenic"))
  return(cellPop)
}

## function to calculate pValue from the theoretical data
theoretical_pValue <- function(counts, refrenceIntergenic){
  ## select all the intergenic region from the library
  ##TODO: once again here it should only contain reference intergenic. To check...
  intergenicRegionsLibrary <- dplyr::filter(counts, type == "intergenic")
  intergenicRegions <- dplyr::filter(refrenceIntergenic, refIntergenic == "TRUE")
  ## keep just information about reference intergenic region detected in the counts file to the
  ## calculation
  selectedRefIntergenic <- merge(intergenicRegionsLibrary, intergenicRegions, by.x="gene_id", by.y="id")
  ## select values with CPM > 0 (because we will use log2 scale)
  selectedRefIntergenic <- dplyr::filter(selectedRefIntergenic, CPM > 0 & type == "intergenic")
  ## select genic and intergenic region from the library with CPM > 0
  regions <- dplyr::filter(counts, CPM > 0)
  ## calculate z-score for each gene_id using the reference intergenic
  regions$zScore <- (log2(regions$CPM) - mean(log2(selectedRefIntergenic$CPM))) / 
    sd(log2(selectedRefIntergenic$CPM))
  ## calculate p-values for each gene_id
  regions$pValue <- pnorm(regions$zScore, lower.tail = FALSE)
  return(list(regions, 2^(mean(log2(selectedRefIntergenic$CPM))), 
    2^(sd(log2(selectedRefIntergenic$CPM)))))
}

cutoff_info <- function(library_id, cellTypeId, counts, desired_pValue_cutoff,
  meanRefIntergenic, sdRefIntergenic){
  ## collect stats using pValue_cutoff
  genic_region_present <- nrow(dplyr::filter(counts, type == "genic" & calls_pValue == "present"))
  proportion_genic_present <- (nrow(dplyr::filter(counts, type == "genic" &
    calls_pValue == "present"))/nrow(dplyr::filter(counts, type == "genic")))*100
  coding_region_present  <- nrow(dplyr::filter(counts, biotype == "protein_coding" &
    calls_pValue == "present"))
  proportion_coding_present <- (nrow(dplyr::filter(counts, biotype == "protein_coding" &
    calls_pValue == "present"))/nrow(dplyr::filter(counts, biotype == "protein_coding")))*100
  intergenic_region_present <- nrow(dplyr::filter(counts, type == "intergenic" &
    calls_pValue == "present"))
  proportion_intergenic_present <- (nrow(dplyr::filter(counts, type == "intergenic" &
    calls_pValue == "present"))/nrow(dplyr::filter(counts, type == "intergenic")))*100
  
  CPM_Threshold <- min(counts$CPM[counts$type == "genic" & counts$calls_pValue == "present"])
  ## Export cutoff_info_file
  collectInfo <- c(library_id, cellTypeId, CPM_Threshold, sum(counts$type == "genic"),
    genic_region_present,proportion_genic_present, sum(counts$biotype %in% "protein_coding"),
    coding_region_present,proportion_coding_present, sum(counts$type %in% "intergenic"),
    intergenic_region_present,proportion_intergenic_present, desired_pValue_cutoff,
    meanRefIntergenic, sdRefIntergenic)
  names(collectInfo) <- c("libraryId", "cellTypeId", "CPM_Threshold", "Genic",
                          "Genic_region_present","Proportion_genic_present",
                          "Protein_coding","Coding_region_present","Proportion_coding_present",
                          "Intergenic","Intergenic_region_present","Proportion_intergenic_present",
                          "pValue_cutoff", "meanRefIntergenic", "sdRefIntergenic")
  return(collectInfo)
}


## export the plot
plotData <- function(libraryId, cellTypeId, counts, refIntergenic, CPM_threshold){
  ## export distribution
  dens <- density(log2(counts$CPM+1e-6), na.rm=T)
  dens_genic <- density(log2(counts$CPM[counts$type == "genic"]+1e-6), na.rm=T)
  dens_genic$y <- dens_genic$y * nrow(dplyr::filter(counts, type == "genic")) / length(counts$CPM)
  dens_coding <- density(log2(counts$CPM[counts$biotype == "protein_coding"]+1e-6), na.rm=T)
  dens_coding$y <- dens_coding$y * nrow(dplyr::filter(counts, biotype == "protein_coding")) /
    length(counts$CPM)
  dens_intergenic <- density(log2(counts$CPM[counts$type == "intergenic"]+1e-6), na.rm=T)
  dens_intergenic$y <- dens_intergenic$y * nrow(dplyr::filter(counts, type == "intergenic")) /
    length(counts$CPM)
  refIntergenic <- merge(dplyr::filter(counts, type == "intergenic"), referenceIntergenic,
    by = "gene_id")
  refIntergenic <- as.data.frame(dplyr::filter(refIntergenic, refIntergenic == "TRUE"))
  dens_Refintergenic <- density(log2(refIntergenic$CPM+1e-6), na.rm=T)
  dens_Refintergenic$y <- dens_Refintergenic$y * nrow(refIntergenic) / length(counts$CPM)
  if(!dir.exists(file.path(output_folder, libraryId))) {
    dir.create(file.path(output_folder,libraryId))
  }
  pdf(file.path(output_folder, libraryId, paste0("Distribution_", libraryId, "_",cellTypeId,
    ".pdf")), width=10, height=6) 
  par(mfrow=c(1,2))
  plot(dens, lwd=2, main=paste0(libraryId), xlab="Log2(CPM)")
  mtext(paste0(cellTypeId))
  lines(dens_genic,col="red", lwd=2)
  lines(dens_coding,col="indianred", lwd=2)
  lines(dens_intergenic,col="darkblue", lwd=2)
  lines(dens_Refintergenic,col="cyan", lwd=2)
  abline(v=CPM_threshold, col="gray", lty=2, lwd=2)
  legend("topright", legend = c("All", "Genic" ,"PC", "Int", "RefInt", "cutoff"),
    col=c("black", "red", "indianred", "darkblue", "cyan", "gray"),lty=c(1,1,1,1,1,2),
    lwd=2, bty = "n")
  
  ## export frequency of pValue for all genic region
  genicRegion <- as.numeric(counts$pValue[counts$type == "genic"])
  hist(na.omit(genicRegion), main=paste0(libraryId), xlab="pValue", xlim=c(0,1))
  mtext(paste0("Genic region_",cellTypeId))
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
  cat("libraryId\tcellTypeId\tCPM_Threshold\tGenic\tGenic_region_present\tProportion_genic_present\tProtein_coding","Coding_region_present\tProportion_coding_present\tIntergenic\tIntergenic_region_present\tProportion_intergenic_present\tpValue_cutoff\tmeanRefIntergenic\tsdRefIntergenic\tspecies\torganism\n",file = file.path(output_folder,"All_cellPopulation_stats_10X.tsv"), sep = "\t")
} else {
  print("File already exist.....")
}

collectSamplesStats <- data.frame()
## loop thought all libraries
for (libraryId in unique(scRNASeqAnnotation$libraryId)) {
  
  if (file.exists(file.path(output_folder, libraryId))) {
    message("calls directory for library ", libraryId, " already processed.")
  } else {
    message("Library: ", libraryId)
    ## collect all cell populations that belongs to the library
    rawCountFiles <- list.files(path = file.path(folder_data, libraryId), pattern = "^Raw_Counts_", full.names = TRUE)
    if (length(rawCountFiles) > 0) {
      dir.create(file.path(output_folder, libraryId))
    }
    species_name <- gsub(" ", "_", unique(scRNASeqAnnotation$scientific_name[scRNASeqAnnotation$libraryId == libraryId]))
    speciesId <- unique(scRNASeqAnnotation$speciesId[scRNASeqAnnotation$libraryId == libraryId])
    for (rawCountFile in rawCountFiles) {
      message(rawCountFile)
      cellPop <- str_remove(basename(rawCountFile), "Raw_Counts_")
      cellPop <- str_remove(cellPop, ".tsv")
      cellTypeId <- gsub("_", ":", cellPop)
      
      message("Process celltype: ", cellTypeId)

      ## collect the sumUMI + normalization for the target cellPop
      cellPop_normalized <- sumUMICellPop(rawCountFile = rawCountFile)
      ## Information about reference intergenic
      referenceIntergenic <- refIntergenic(counts = cellPop_normalized, folderGtf = folder_gtf,
        speciesName = species_name)
      ## calls with pValue theoretical
      calculatePvalues <- theoretical_pValue(counts = cellPop_normalized,
        refrenceIntergenic = referenceIntergenic)
      calculationInfo <- calculatePvalues[[1]]
      calculationInfo$calls_pValue <- ifelse(calculationInfo$pValue <=
        as.numeric(desired_pValue_cutoff), "present", "absent" )
      
      ## add info also about genesID where CPM = 0 and were not used for the pValue calculation
      ## (but are important for the final stats)
      ##TODO: this part will have to be removed if we consider genes with CPM = 0 for pValue calculation.
      ##      Otherwise, remove this todo
      regionZero <- cellPop_normalized[!cellPop_normalized$gene_id %in% calculationInfo$gene_id,]
      regionZero$zScore <- "NA"; regionZero$pValue <- "NA"; regionZero$calls_pValue <-"absent"
      
      allData <- rbind(calculationInfo, regionZero)
      #TODO Actually, why bother reordering??
      #What a creepy reordering.
      #genicRegion <- dplyr::filter(allData, type=="genic")
      #genicRegion <- genicRegion[order(genicRegion$gene_id),]
      #finalData <- rbind(genicRegion, dplyr::filter(allData, type == "intergenic"))
      #finalData$cellTypeName <- cellName
      finalData <- allData[order(allData$type, allData$gene_id), ]
      finalData$cellTypeId <- cellTypeId
      
      ## collect just genic region and re-calculate CPM
      finalData_genic <- dplyr::filter(finalData, type == "genic")
      finalData_genic$CPM <- finalData_genic$sumUMI / sum(finalData_genic$sumUMI) * 1e6
      
      ## Export cutoff information file + new files with calls
      cutoff_info_file <- cutoff_info(libraryId, cellTypeId = cellTypeId, counts = finalData,
        desired_pValue_cutoff = as.numeric(desired_pValue_cutoff),
        meanRefIntergenic = calculatePvalues[[2]], sdRefIntergenic = calculatePvalues[[3]])
      CPM_threshold <- log2(as.numeric(cutoff_info_file[3]))
      
      ## export data
      ## for some libraries it is not possible to create the plot. Use a try/catch to be able
      ## to generate the plot for all possible libraries
      ##TODO: check why it is not possible to generate the plot
      tryCatch(
        expr = {
          plotData(libraryId = libraryId, cellTypeId = cellPop, counts = finalData,
            refIntergenic = referenceIntergenic, CPM_threshold = CPM_threshold)
          },
          error = function(e) {
            warning("did not manage to create plot for library ", libraryId)
          }
      )
      
      pathExport <- file.path(output_folder, libraryId)
      write.table(finalData,file = file.path(pathExport, paste0("Calls_cellPop_",libraryId,
        "_",cellPop,"_genic+intergenic.tsv")),quote=FALSE, sep = "\t", col.names=TRUE,
        row.names=FALSE)
      write.table(finalData_genic,file = file.path(pathExport, paste0("Calls_cellPop_",
        libraryId, "_",cellPop,"_genic.tsv")),quote=FALSE, sep = "\t", col.names=TRUE,
        row.names=FALSE)
      write.table(t(t(cutoff_info_file)),file = file.path(pathExport,
        paste0("cutoff_info_file_",libraryId, "_",cellPop,".tsv")),quote=FALSE, sep = "\t",
      col.names=FALSE, row.names=TRUE)
      
      ## add this to big data frame with all samples information
      this_sample <- as.data.frame(t(cutoff_info_file), stringsAsFactors=F)
      this_sample[, "species"]  <- speciesId
      this_sample[, "organism"] <- as.character(unique(
        scRNASeqAnnotation$scientific_name[scRNASeqAnnotation$speciesId == speciesId]))
      collectSamplesStats <- rbind(collectSamplesStats, this_sample)
    }
  }
}
write.table(collectSamplesStats, file = file.path(output_folder, "All_cellPopulation_stats_10X.tsv"),
  col.names = FALSE, row.names = FALSE, append = TRUE, quote = FALSE, sep = "\t")

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
