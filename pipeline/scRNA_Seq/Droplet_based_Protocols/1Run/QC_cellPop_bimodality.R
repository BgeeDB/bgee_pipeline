## SFonsecaCosta, 2020

## This script is used to check the bimodality of the cell population per library.
## After cell-type identification we test if the cell-population per library follow a bimodal distribution,
## by excluding some noisy genes from the population

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNA_Seq_info_TargetBased.txt" folder_data="folder_data" output_folder="output_folder"' QC_cellPop_bimodality.R QC_cellPop_bimodality.Rout
## scRNASeq_Info --> File that results from annotation and metadata (libraries downloaded and with extra information as SRR)
## folder_data --> Folder where are all the libraries with correspondent cell-type population identification (Raw and normalized files)
## output_folder --> Folder where the results should be saved (normally same that folder_data)

## libraries used
library(data.table)
library(stringr)
library(dplyr)
library(mclust)
library(LaplacesDemon)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed.
command_arg <- c("scRNASeq_Info", "folder_data", "output_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read scRNA-Seq info file. If file not exists, script stops
if( file.exists(scRNASeq_Info) ){
  scRNA_annotation_file <- fread(scRNASeq_Info)
} else {
  stop( paste("The annotation file not found [", scRNASeq_Info, "]\n"))
}
###############################################################################################
## function to make the deconvolution of the protein coding density
deconvolution <- function(UMIgeneID){
  ## select just protein_coding genes to the Mclust
  proteinCoding <- UMIgeneID[UMIgeneID$biotype == "protein_coding", ]
  decov = densityMclust(proteinCoding$ratio)
  proteinCoding$classification <- decov$classification
  return(proteinCoding)
}

## plot the data for the cell population if pass the QC
plotData <- function(libraryID, allInformation, deconvolutionInfo, classification, cutoff, cellPopName, modesInfo){
  
  pdf(file = paste0(output_folder, "/", libraryID, "/QC_bimodality_",cellPopName, ".pdf"))
  layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE))
  ## plot info about UMI counts and genes detected per individual cell.
  plot(allInformation$genes, allInformation$UMIcounts, xlab="Number of genes", ylab="UMI counts", main=paste0("cell pop ",cellPopName), pch=20)
  mtext(paste0(sizeData, " Cells"))
  hist(allInformation$genes, xlab="genes", main=paste0(cellPopName))
  hist(allInformation$UMIcounts, xlab="UMI counts", main=paste0(cellPopName))
  
  ## plot deconvolution curves
  plot(density((deconvolutionInfo$ratio)), lwd=3, main="Density protein coding", xlab="ratio")
  for (i in classification) {
    densityPlot <- density((deconvolutionInfo$ratio[deconvolutionInfo$classification == i]))
    densityPlot$y <- densityPlot$y * length(deconvolutionInfo$ratio[deconvolutionInfo$classification == i]) / length(deconvolutionInfo$ratio)
    lines(densityPlot, col="indianred", lwd=2, lty=2)
    ## Print gaussian number on plot
    text(densityPlot$x[densityPlot$y == max(densityPlot$y)], 0.005, labels = i, col="indianred", font=3)
  }
  abline(v=cutoff, lwd=1, lty=2, col="gray")
  legend("topright", legend = c("Protein_coding", "Deconvolution","cutoff genes"), col=c("black", "indianred", "gray"), pch=20, bty = "n")
  
  
  ## plot density after filtering low detected genes across the cell population
  pc_afterFiltering <- deconvolutionInfo[deconvolutionInfo$ratio >= cutoff, ]
  densityPC <- density(pc_afterFiltering$ratio)
  plot(densityPC,col="darkblue",lwd=2, main="Density of protein coding", xlab="ratio", xlim=c(min(densityPC$x), max(densityPC$x)))
  mtext(paste0("Noisy genes < ", round(cutoff, digits = 3)))
  legend("topright", legend = c(paste0("Bimodality = ", is.bimodal(pc_afterFiltering$ratio, min.size=0.1)), paste0("Mode_1 = ", round(modesInfo[1], digits = 2)), paste0("Mode_2 = ", round(modesInfo[2], digits = 2))), bty = "n")
  dev.off()
}

bimodality_targetBased <- file.path(output_folder, "bimodality_targetBased.txt")
if (!file.exists(bimodality_targetBased)){
  file.create(bimodality_targetBased)
  cat("library\texperimentID\tCell_Name_ID\tcomment\n",file = bimodality_targetBased, sep = "\t")
} else {
  print("File already exist.....")
}

## apply for each library/cellpop
for (libraryID in unique(scRNA_annotation_file$libraryId)) {
  
  ## verify if the code already run for this library (check if the bimodality_DONE file already exist)
  bimodalityDone <- file.exists(file.path(folder_data, libraryID, "bimodality_DONE.txt"))
  
  if (bimodalityDone == TRUE){
    message("Bimodality done for this library: ", libraryID)
  } else {
    
    message("Treating library: ", libraryID)
    
    ## select all cell population for the library
    AllCellPop <- list.files(path = file.path(folder_data, libraryID), pattern = "^Raw_Counts_")
    experimentID <- scRNA_annotation_file$experimentId[scRNA_annotation_file$libraryId == libraryID]
    
    for (cellPop in AllCellPop) {
      
      cellPopName <- str_remove(cellPop, "Raw_Counts_")
      cellPopName <- str_remove(cellPopName, ".tsv")
      
      message("Doing: ", cellPopName)
      
      cellpop <- fread(file.path(folder_data,libraryID,cellPop))
      ## remove info about gene_name, biotype, type
      sizeData <- length(cellpop)-5
      colectInfo <- cellpop %>% dplyr::select("gene_id", "biotype", "type", "cellTypeName", "cellTypeId")
      
      ## just make QC bimodality if cell population have at least 50 cells
      if (sizeData >= 50){
        
        genicRegion <- dplyr::filter(cellpop, type == "genic")
        genicRegion <- genicRegion[, -c("biotype", "type", "cellTypeName", "cellTypeId")]
        
        ## UMI counts and genes detected per cell using all genic region
        UMIcounts <- as.data.frame(colSums((genicRegion[,2:ncol(genicRegion)])))
        colnames(UMIcounts) <- "UMIcounts"
        UMIgenes <- as.data.frame(apply(genicRegion[,2:length(genicRegion)],2,function(x)sum(x != 0)))
        colnames(UMIgenes) <- "genes"
        allInformation <- cbind(UMIcounts, UMIgenes)
        allInformation$cells <- rownames(allInformation)
        
        ## verify in how many genes we have UMI higher 0 across the cell population
        UMIgeneID <- as.data.frame(apply(genicRegion[,2:length(genicRegion)],1,function(x)sum(x != 0)))
        colnames(UMIgeneID) <- "genes"
        UMIgeneID$ratio <- UMIgeneID$genes/sizeData
        UMIgeneID$gene_id <- genicRegion$gene_id
        UMIgeneID <- merge(UMIgeneID, colectInfo, by = "gene_id")
        ## remove genes never detected in the cell population = 0
        UMIgeneID <- UMIgeneID[UMIgeneID$ratio > 0, ]
        
        ## deconvolution of protein coding density
        deconvolutionInfo <- deconvolution(UMIgeneID = UMIgeneID)
        classification <- sort(unique(deconvolutionInfo$classification))
        ## verify the minimum amount of genes that can be classified as noisy genes (this means detected in few cells)
        for (i in classification) {
          
          cutoff <- min((deconvolutionInfo$ratio[deconvolutionInfo$classification == i]))
          proteinCoding <- deconvolutionInfo[deconvolutionInfo$ratio >= cutoff, ]
          bimodalityCalculation <- is.bimodal(proteinCoding$ratio, min.size=0.1)
          modesInfo <- Modes(proteinCoding$ratio, min.size=0.1)
          
          if (bimodalityCalculation == TRUE){
            break
          } 
        }
        
        ## verify if 1 mode > 0.5 and 2 mode < 0.8 in the ratio (in order to use the bimodality distribution to validate genes that are always detected in the cell-population)
        mode_1 <- modesInfo$modes[1] > 0.5
        mode_2 <- modesInfo$modes[2] < 0.8
        
        if (bimodalityCalculation == TRUE & mode_1 == FALSE & mode_2 == FALSE){
          message("The cell population ", cellPopName, " from the library ", libraryID ," follow a bimodal distribution.")
          ## retrieve modes info
          modesInfo <- modesInfo$modes
          ## plot data 
          plotData(libraryID = libraryID, allInformation = allInformation, deconvolutionInfo = deconvolutionInfo, classification = classification, cutoff = cutoff, cellPopName = cellPopName, modesInfo = modesInfo)
        } else if (bimodalityCalculation == TRUE & mode_1 == TRUE | bimodalityCalculation == TRUE & mode_2 == TRUE){
          message("The cell population ", cellPopName, " from the library ", libraryID ," not present confidence enough to detect a set of genes that should be always detected for the cell-type.")
          ## add to the excluded file libraries/cell pop
          infoCollected <- data.frame(libraryID, experimentID, cellPopName, "mode_1 or mode_2 not fit the minimum requirement")
          write.table(infoCollected, file = bimodality_targetBased, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
        } else {
          message("The cell population ", cellPopName, " from the library ", libraryID ," is not bimodal.")
          ## add to the excluded file libraries/cell pop
          infoCollected <- data.frame(libraryID, experimentID, cellPopName, "not bimodal")
          write.table(infoCollected, file = bimodality_targetBased, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
        }
      ## collect information for libraries/cell pop not pass the minimum requirement 50 cells or are not bimodal population
      } else {
        message("The cell population ", cellPopName, " from the library ", libraryID ," have < 50 cells.")
        ## add to the excluded file libraries/cell pop
        infoCollected <- data.frame(libraryID, experimentID, cellPopName, "< 50 cells")
        write.table(infoCollected, file = bimodality_targetBased, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
      }
    }
    ## control file in case densityMclust stops in some library because of (.Machine$double.xmax)
    file.create(file.path(folder_data, libraryID, "bimodality_DONE.txt"))
  }
}
