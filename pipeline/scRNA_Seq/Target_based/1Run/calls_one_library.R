## SFonsecaCosta, October 2020

## This script is used make the calls of present and absent genes per library/cell-type
## population using pValue theoretical.
## This means: sum the UMI that belongs to the same cell-type population, normalize CPM and then
## call present and absent genes based on the pValue_theoretical cut-off.

## NOTE: Since Bgee 15.1 do not filter anymore library/cell population that have more than 50 cells
## per library/cell-type pop and pass the bimodality test.

## Starting from Bgee 15.2 this script process one library. It is an other script that manages all libraries
## by creating jobs on a slurm cluster. Each job correspond to one library

## Usage:
## R CMD BATCH --no-save --no-restore '--args libraryId="libraryId" speciesId="speciesId" speciesName="speciesName" celltypeFolder="celltypeFolder" refIntergenicFolder="refIntergenicFolder" pValueCutoff="pValueCutoff" callsOutputFolder="callsOutputFolder"' calls_one_library.R calls_one_library.Rout
## libraryId                --> ID of the library
## speciesId                --> NCBI taxon ID of the species (e.g 9606)
## speciesName              --> species name (e.g Homo sapiens)
## celltypeFolder           --> Folder containing celltype quantification for that library
## refIntergenicFolder      --> Folder containing all ref intergenic
## pValueCutoff             --> desired pValue cutoff to call present genes
## callsOutputFolder        --> Folder where we should save the results

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
command_arg <- c("libraryId", "speciesId", "speciesName" ,"celltypeFolder",
  "refIntergenicFolder", "pValueCutoff", "callsOutputFolder")
for (c_arg in command_arg) {
  if (!exists(c_arg)) {
    stop(paste(c_arg,"command line argument not provided\n"))
  }
}

##########################################################################################################################################################
## Provide the reference intergenic regions = TRUE
##TODO Why bother retrieving ref. intergenic??? the mapping/quantification should already be done
## on the reference intergenic.
##TODO: check that only ref intergenic were used for each step of the pipeline( e.g gene annotation, kallisto index, ....)
refIntergenicIds <- function(refIntergenicFolder, speciesId){
  intergenicCoordinates <- read.table(file = file.path(refIntergenicFolder,
    paste0(speciesId,"_coordinates.tsv"), sep = "\t", header = T)
  intergenicIds <- paste0(intergenic$chr, "-", intergenic$start,"-",intergenic$end)
  return(intergenicIds)
}

## Sum the UMI of all barcodes and then compute the CPM normalization
sumUMICellPop <- function(rawCountFile, refIntergenicIds) {
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
  #TODO: filtering to remove for Bgee 16.0 as only ref. intergenic should be used
  # filter to keep ref. intergenic only
  cellPop <- dplyr::filter(cellPop, (type == "genic") ||
    (type == "intergenic" & ! gene_id %in% refIntergenicIds))$gene_id
  ## just re-order, why bother reordering that way? couldn't we simply order on gene_id???
  cellPopGenic <- data.frame(dplyr::filter(cellPop, type == "genic"))
  cellPopGenic <- cellPopGenic[order(cellPopGenic$gene_id),]
  cellPop <- rbind(cellPopGenic, dplyr::filter(cellPop, type == "intergenic"))
  return(cellPop)
}

## function to calculate pValue from the theoretical data
theoretical_pValue <- function(counts){
  ##TODO: once again here it should only contain reference intergenic. To check...
  ## select values with CPM > 0 (because we will use log2 scale)
  selectedRefIntergenic <- dplyr::filter(counts, CPM > 0 & type == "intergenic")
  ## select genic and ref. intergenic region from the library with CPM > 0
  regions <- dplyr::filter(counts, CPM > 0)
  ## calculate z-score for each gene_id using the reference intergenic
  regions$zScore <- (log2(regions$CPM) - mean(log2(selectedRefIntergenic$CPM))) / 
    sd(log2(selectedRefIntergenic$CPM))
  ## calculate p-values for each gene_id
  regions$pValue <- pnorm(regions$zScore, lower.tail = FALSE)
  return(list(regions, 2^(mean(log2(selectedRefIntergenic$CPM))), 
    2^(sd(log2(selectedRefIntergenic$CPM)))))
}

cutoff_info <- function(library_id, cellTypeId, counts, pValueCutoff,
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
    intergenic_region_present,proportion_intergenic_present, pValueCutoff,
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
  pdf(file.path(libraryOutputFolder, paste0("Distribution_", libraryId, "_",cellTypeId,
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
  abline(v=pValueCutoff, col="red", lwd=2)
  dev.off()
}

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

##################### Main part of the script ##########################

libraryOutputFolder <- file.path(callsOutputFolder, libraryId)
libStatsFile <- file.path(libraryOutputFolder, paste0(libraryId,"_stats.tsv"))

if (file.exists(file.path(libraryOutputFolder, libStatsFile))) {
  message("calls directory for library ", libraryId, " already processed.")
  stop_quietly()
}

## collect all cell population files that belongs to the library
rawCountFiles <- list.files(path = file.path(celltypeFolder, libraryId), pattern = "^Raw_Counts_", full.names = TRUE)
if (length(rawCountFiles) > 0) {
  dir.create(libraryOutputFolder, recursive = TRUE)
} else {
  message("no celltype raw count file for library ", libraryId, ". No calls will be generated for that library")
  stop_quietly()
}

## Information about reference intergenic
referenceIntergenicIds <- refIntergenicIds(refIntergenicFolder, speciesId)
allCelltypeInfo <- file.path(libraryOutputFolder, "All_cellPopulation_stats_10X.tsv")

collectSamplesStats <- ()

for (rawCountFile in rawCountFiles) {
  message(rawCountFile)
  cellPop <- str_remove(basename(rawCountFile), "Raw_Counts_")
  cellPop <- str_remove(cellPop, ".tsv")
  cellTypeId <- gsub("_", ":", cellPop)
  
  message("Process celltype: ", cellTypeId)
  ## collect the sumUMI + normalization for the target cellPop
  cellPop_normalized <- sumUMICellPop(rawCountFile = rawCountFile)

  ## calculate pValue of presence/absence of expression
  calculatePvalues <- theoretical_pValue(counts = cellPop_normalized)
  calculationInfo <- calculatePvalues[[1]]
  calculationInfo$calls_pValue <- ifelse(calculationInfo$pValue <=
    as.numeric(pValueCutoff), "present", "absent" )
  
  ## add info also about genesID where CPM = 0 and were not used for the pValue calculation
  ## (but are important for the final stats)
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
    pValueCutoff = as.numeric(pValueCutoff),
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
  
  write.table(finalData,file = file.path(libraryOutputFolder, paste0("Calls_cellPop_",libraryId,
    "_",cellPop,"_genic+intergenic.tsv")),quote=FALSE, sep = "\t", col.names=TRUE,
    row.names=FALSE)
  write.table(finalData_genic,file = file.path(libraryOutputFolder, paste0("Calls_cellPop_",
    libraryId, "_",cellPop,"_genic.tsv")),quote=FALSE, sep = "\t", col.names=TRUE,
    row.names=FALSE)
  write.table(t(t(cutoff_info_file)),file = file.path(libraryOutputFolder,
    paste0("cutoff_info_file_",libraryId, "_",cellPop,".tsv")),quote=FALSE, sep = "\t",
  col.names=FALSE, row.names=TRUE)
  
  ## add this to big data frame with all samples information
  this_sample <- as.data.frame(t(cutoff_info_file), stringsAsFactors=F)
  this_sample[, "species"]  <- speciesId
  this_sample[, "organism"] <- as.character(unique(
    scRNASeqAnnotation$scientific_name[scRNASeqAnnotation$speciesId == speciesId]))
  collectSamplesStats <- rbind(collectSamplesStats, this_sample)
}

## export info stats of all libraries/cell-population
file.create(libStatsFile)
cat("libraryId\tcellTypeId\tCPM_Threshold\tGenic\tGenic_region_present\tProportion_genic_present\tProtein_coding\t",
  "Coding_region_present\tProportion_coding_present\tIntergenic\tIntergenic_region_present\tProportion_intergenic_present\t",
  "pValue_cutoff\tmeanRefIntergenic\tsdRefIntergenic\tspecies\torganism\n",
  file = libStatsFile, sep = "\t")
} else {
  print("File already exist.....")
}
write.table(collectSamplesStats, file = libStatsFile, col.names = FALSE, row.names = FALSE,
  append = TRUE, quote = FALSE, sep = "\t")

## not anymore generated in this script. Need to create an other script to generate those plots
## final plot per species
#pdf(file.path(output_folder, paste0("All_libraries_stats_information_10X.pdf")), width=16, height=6) 
#g1 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(Proportion_genic_present))) + 
#  geom_boxplot()+ylim(0,100)+xlab(" ")+ylab("% Genic Present")
#g2 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(Proportion_coding_present))) + 
#  geom_boxplot(notch=TRUE)+ylim(0,100)+xlab(" ")+ylab("% Protein Coding Present")
#g3 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(Proportion_intergenic_present))) + 
#  geom_boxplot(notch=TRUE)+ylim(0,100)+xlab(" ")+ylab("% Intergenic Present")
#g4 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(meanRefIntergenic))) + 
#  geom_boxplot(notch=TRUE)+ylim(0,1)+xlab(" ")+ylab("Mean Reference Intergenic")
#g5 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(sdRefIntergenic))) + 
#  geom_boxplot(notch=TRUE)+ylim(0,10)+xlab(" ")+ylab("SD Reference Intergenic")
#g6 <- ggplot(collectSamplesStats, aes(x=organism, y=as.numeric(pValue_cutoff))) + 
#  geom_boxplot(notch=TRUE)+ylim(0,1)+xlab(" ")+ylab("pValue_cutoff")
#grid.arrange(g1,g2,g3,g4,g5,g6, nrow = 2, ncol = 3)
#dev.off()
