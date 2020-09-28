## SFonsecaCosta, Jully 2020

## This script is used separately from Call_PresentGenes_indivCell_cellPop.R to make the code more clean and to run independently. 
## Here for each condiction we compute the sum of the raw UMI counts across cells at gene level and then we compute the correspondent CPM for each gene
## The output files exported are: raw+normalization just for genic regions (used later in the rank computation) and raw+normalization for genic + intergenic region.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNASeq_Info.txt" InformationAllLibraries="InformationAllLibraries.txt" folder_data="folder_data" output="output"' sum_raw_UMI.R sum_raw_UMI.Rout
## scRNASeq_Info --> File that results from annotation and metadata (libraries downloaded and with extra information as readlength or SRR) 
## InformationAllLibraries --> File with information about each cell after barcode and gene markers annotation (per library contain total number of cells and cell Name)
## folder_data --> Folder where are all the libraries after cell identification
## output --> Folder where we should save the results 

## libraries used
library(dplyr)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed.
command_arg <- c("scRNASeq_Info", "InformationAllLibraries", "folder_data", "output")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read scRNASeq_Info file. If file not exists, script stops.
if( file.exists(scRNASeq_Info) ){
  scRNASeqAnnotation <- read.table(scRNASeq_Info, h=T, sep="\t", comment.char="", quote = "")
} else {
  stop( paste("The scRNASeq information file was not found [", scRNASeq_Info, "]\n"))
}
## Read InformationAllLibraries file. If file not exists, script stops.
if( file.exists(InformationAllLibraries) ){
  cellInfo <- read.table(InformationAllLibraries, h=T, sep="\t", comment.char="", quote = "")
} else {
  stop( paste("The cell information file was not found [", InformationAllLibraries, "]\n"))
}
##########################################################################################################################################################
## write info file about cell population
cellPopInfo <- paste0(output, "/sum_scRNASeq_sample_info.txt")
if (!file.exists(cellPopInfo)){
  file.create(cellPopInfo)
  cat("libraryId\texperimentId\tuberonId\tcellTypeId\tstageId\tsex\tstrain\tspeciesId\n",file = paste0(output, "/sum_scRNASeq_sample_info.txt"), sep = "\t")
} else {
  print("File already exist.....")
}

sumUMI <- function(finalTable, output){
  
  finalTable[is.na(finalTable)] <- 0
  ## remove rows with zeros in all cells
  finalTable <- finalTable[which(rowSums(finalTable[,4:ncol(finalTable)]) > 0), ]
  ## add UMI sum column
  finalTable$sum_raw_UMI <- rowSums(finalTable[,4:ncol(finalTable)])
  
  ## Raw UMI counts + normalized counts using just the genic region
  justGenic <- dplyr::filter(finalTable, type == "genic")
  justGenic$CPM <- justGenic$sum_raw_UMI / sum(justGenic$sum_raw_UMI) * 1e6
  justGenic <- justGenic[c("gene_id", "biotype", "type", "sum_raw_UMI", "CPM")]
   
  ## Raw UMI counts + normalized counts with genic and intergenic regions
  allRegions <- finalTable
  allRegions$CPM <- allRegions$sum_raw_UMI / sum(allRegions$sum_raw_UMI) * 1e6
  allRegions <- allRegions[c("gene_id", "biotype", "type", "sum_raw_UMI", "CPM")]
  
   return(list(justGenic, allRegions))
}

## loop through all data
for (species in unique(scRNASeqAnnotation$speciesId)) {
  cat("Species:", species, "\n")
  for (experiment in unique(scRNASeqAnnotation$experimentId[scRNASeqAnnotation$speciesId == species])){
    cat("Name of experiments that belongs to the species:", experiment, "\n")
    for (uberonId in unique(scRNASeqAnnotation$uberonId[ scRNASeqAnnotation$speciesId == species & scRNASeqAnnotation$experimentId == experiment])) {
      cat("Uberon info:", uberonId, "\n")
      for (cellId in unique(scRNASeqAnnotation$cellTypeId[scRNASeqAnnotation$speciesId == species & scRNASeqAnnotation$experimentId == experiment & scRNASeqAnnotation$uberonId == uberonId])){
        cat("CellId info:", cellId, "\n")
        for (stageId in unique(scRNASeqAnnotation$stageId[scRNASeqAnnotation$speciesId == species & scRNASeqAnnotation$uberonId == uberonId & scRNASeqAnnotation$cellTypeId == cellId & scRNASeqAnnotation$experimentId == experiment])){
          cat("StageId info:", stageId, "\n")
          for (sex in unique(scRNASeqAnnotation$sex[scRNASeqAnnotation$speciesId == species & scRNASeqAnnotation$experimentId == experiment & scRNASeqAnnotation$uberonId == uberonId & scRNASeqAnnotation$cellTypeId == cellId & scRNASeqAnnotation$stageId == stageId])){
            cat("Sex info:", sex, "\n")
            sex <- paste0(sex)
            for (strain in unique(scRNASeqAnnotation$strain[scRNASeqAnnotation$speciesId == species & scRNASeqAnnotation$experimentId == experiment & scRNASeqAnnotation$uberonId == uberonId & scRNASeqAnnotation$cellTypeId == cellId & scRNASeqAnnotation$stageId == stageId & scRNASeqAnnotation$sex == sex])){
              cat("Strain info:", strain, "\n")
              strain <- paste0(strain)
              
              
              if (sex == "NA" & strain == "NA"){
                librariesInfo <- scRNASeqAnnotation$libraryId[scRNASeqAnnotation$speciesId == species & scRNASeqAnnotation$experimentId == experiment & scRNASeqAnnotation$uberonId == uberonId & scRNASeqAnnotation$cellTypeId == cellId & scRNASeqAnnotation$stageId == stageId]
              } else if (sex != "NA" & strain == "NA") {
                librariesInfo <- scRNASeqAnnotation$libraryId[scRNASeqAnnotation$speciesId == species & scRNASeqAnnotation$experimentId == experiment & scRNASeqAnnotation$uberonId == uberonId & scRNASeqAnnotation$cellTypeId == cellId & scRNASeqAnnotation$stageId == stageId & scRNASeqAnnotation$sex == sex]
              } else {
                librariesInfo <- scRNASeqAnnotation$libraryId[scRNASeqAnnotation$speciesId == species & scRNASeqAnnotation$experimentId == experiment & scRNASeqAnnotation$uberonId == uberonId & scRNASeqAnnotation$cellTypeId == cellId & scRNASeqAnnotation$stageId == stageId & scRNASeqAnnotation$sex == sex & scRNASeqAnnotation$strain == strain]
              }
              
              ## collect cell-type information from info file
              cell_Name <- unique(cellInfo$Cell_Name[cellInfo$experimentID == experiment])
              
              for (cell in cell_Name) {
                file <- as.data.frame(paste0(folder_data,  librariesInfo, "/Raw_Counts_", cell, ".tsv"))
                colnames(file) <- "path"
                
                fileTRUE <- c()
                ## select just files that exist in the correspondent library
                for (i in file$path) {
                  verify <- file.exists(i)
                  if (verify == TRUE){
                    fileTRUE <- rbind(fileTRUE, i)
                  } else {
                    cat("The file ", i , "not exist in the correspondent library","\n",
                        "The pricipal cause is mainly the barcode or gene markers not identify any cell in this library for this cell type.", "\n")
                  }
                }
                
                cat("Number of libraries where the cell ", cell, " is detected :", length(fileTRUE), "\n")
                ## collect info to export as file name
                fileName <- paste0(cell , "_", experiment,"_", uberonId,"_", cellId,"_", stageId,"_", sex,"_", strain,"_", species)
                fileName <- gsub(":","-",fileName)
                fileName <- gsub("/","-",fileName)
                
                ## create a folder with name of the cells summed
                outFolder <- paste0(output, "/", fileName)
                if (!dir.exists(outFolder)){
                  dir.create(outFolder)
                } else {
                  print("Folder already exist.....")
                }
                
                ## read all files for the correspondent cell-type
                if (length(fileTRUE) == 0){
                  
                  cat("Libraries not contain this cell type with this correspondent condictions!", "\n")
                  
                } else if (length(fileTRUE) == 1){
                  finalTable <- read.table(fileTRUE, header = TRUE, sep = "\t")
                  finalTableIdentifiers <- finalTable[ ,c("gene_id", "biotype", "type")]
                  finalTable$gene_id <- finalTable$biotype <- finalTable$type <- NULL
                  finalTable <- data.frame(finalTableIdentifiers, finalTable)
                  
                  ## collect raw info + normalization
                  sumUMI_information <- sumUMI(finalTable = finalTable, output = outFolder)
                  ## write tables
                  write.table(sumUMI_information[[1]], file = paste0(outFolder, "/SumRaw+Normalization_Genic_" , fileName, ".tsv"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
                  write.table(sumUMI_information[[2]], file = paste0(outFolder, "/SumRaw+Normalization_Genic+Intergenic_" , fileName, ".tsv"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
                  ## write information to the info file about cell population
                  infoCellPop <- c(fileName, cell, experiment, uberonId, cellId, stageId, sex, strain, species)
                  write.table(t(infoCellPop), file = cellPopInfo, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
                } else {
                  All_libs <- lapply(fileTRUE, read.delim, stringsAsFactors=FALSE)
                  ## merge by geneID from all samples that have the same cell-type
                  finalTable <- Reduce(function(...) merge(..., by = c("gene_id", "biotype","type"), all=TRUE), All_libs)
                  
                  ## collect raw info + normalization
                  sumUMI_information <- sumUMI(finalTable = finalTable, output = outFolder)
                  ## write tables
                  write.table(sumUMI_information[[1]], file = paste0(outFolder, "/SumRaw+Normalization_Genic_" , fileName, ".tsv"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
                  write.table(sumUMI_information[[2]], file = paste0(outFolder, "/SumRaw+Normalization_Genic+Intergenic_" , fileName, ".tsv"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
                  ## write information to the info file about cell population
                  infoCellPop <- c(fileName, cell, experiment, uberonId, cellId, stageId, sex, strain, species)
                  write.table(t(infoCellPop), file = cellPopInfo, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
                }
                
              }
            }
          }
        }
      }
    }
  }
}
