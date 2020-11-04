## SFonsecaCosta, October 2020

## This script is used to generate the file library/cell-type population pValue theoretical.
## This means: sum the UMI that belongs to the same cell-type population, normalize CPM and then call present genes based on the pValue_theoretical cut-off.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNASeq_Info.txt" InformationAllLibraries="InformationAllLibraries.txt" folder_data="folder_data" folder_refIntergenic="folder_refIntergenic" desired_r_cutoff="desired_r_cutoff" desired_pValue_cutoff="desired_pValue_cutoff" output="output"' Call_PresentGenes_indivCell_cellPop.R Call_PresentGenes_indivCell_cellPop.Rout
## scRNASeq_Info --> File that results from annotation and metadata (libraries downloaded and with extra information as readlength or SRR) 
## InformationAllLibraries --> File with information about each cell after barcode and gene markers annotation (per library contain total number of cells and cell Name)
## folder_data --> Folder where are all the libraries after cell identification
## folder_refIntergenic --> Folder where is located the reference intergenic files for each species
## desired_r_cutoff --> proportion of intergenic allowed (both individual cell and cell population to define the ratio cutoff)
## desired_pValue_cutoff --> desired pValue cutoff to call present genes
## output --> Folder where we should save the results 

## libraries used
library(Biostrings)
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
command_arg <- c("scRNASeq_Info", "InformationAllLibraries", "folder_data", "folder_refIntergenic" ,"desired_r_cutoff", "desired_pValue_cutoff", "output")
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
## Read InformationAllLibraries file. If file not exists, script stops
if( file.exists(InformationAllLibraries) ){
  cellInfo <- read.table(InformationAllLibraries, h=T, sep="\t", comment.char="", quote = "")
} else {
  stop( paste("The cell information file was not found [", InformationAllLibraries, "]\n"))
}
##########################################################################################################################################################
## Sum the UMI of all barcodes and then compute the CPM normalization
sumUMICellPop <- function(folder_data, library, cellPop){
  
  cellPop <- fread(file.path(folder_data,library, cellPop))
  cellPop$sumUMI <- rowSums(cellPop[ ,2:(length(cellPop)-2)])
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
  
  referenceIntergenic <- paste0(folder_refIntergenic, "/", speciesID, "_intergenic.fa")
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

## Perform calls using the ratio (old approach)
calculate_r <- function(counts, seqNamesFinal, ratioValue){

 selected_coding <- counts$biotype %in% "protein_coding" 
 selected_intergenic <- (counts$type %in% "intergenic" & seqNamesFinal$refIntergenic == "TRUE") 
    
 summed_intergenic <- sapply(unique(sort(counts$CPM[selected_coding])), function(x){
    return( sum(counts$CPM[selected_intergenic] >= x) )
  })
  summed_coding <- c(0, cumsum(rle(sort(counts$CPM[selected_coding]))$lengths))
  summed_coding <- summed_coding[-(length(summed_coding))]
  summed_coding <- sum(selected_coding) - summed_coding
  
  r <- ( summed_intergenic / sum(selected_intergenic) ) / ( summed_coding / sum(selected_coding) )
  percent <- (1-ratioValue)*100
  
  if (sum(r < ratioValue) == 0){
    CPM_cutoff <- sort(unique(counts$CPM[selected_coding]))[which(r == min(r))[1]]
    r_cutoff <- min(r)
    cat(paste0("There is no CPM cutoff for which " , percent,"%", " of the expressed genes would be coding.", "\n",
               "CPM cutoff is fixed at the first value with maximum coding/intergenic ratio.", "\n",
               "r=", r_cutoff, "at CPM=", CPM_cutoff,"\n"))
  } else {
    CPM_cutoff <- sort(unique(counts$CPM[selected_coding]))[which(r < ratioValue)[1]]
    r_cutoff <- ratioValue
    cat(paste0("CPM cutoff for which " , percent,"%", " of the expressed genes are be coding found at CPM=", CPM_cutoff,"\n"))
  }
  return(list(CPM_cutoff, r_cutoff))
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
  
  ## select genic region from the library
  genicRegions <- dplyr::filter(counts, CPM > 0 & type == "genic")
  ## calculate z-score for each gene_id using the reference intergenic
  genicRegions$zScore <- (log2(genicRegions$CPM) - mean(log2(selected_Ref_Intergenic$CPM))) / sd(log2(selected_Ref_Intergenic$CPM))
  ## calculate p-values for each gene_id
  genicRegions$pValue <- pnorm(genicRegions$zScore, lower.tail = FALSE)
  return(genicRegions)
}


## loop thought all libraries
for (library in unique(scRNASeqAnnotation$libraryId)) {
  
  ## collect all cell populations that belongs to the library
  AllCellPop <- list.files(path = file.path(folder_data,library), pattern = "^Raw_Counts_")
  speciesID <- scRNASeqAnnotation$speciesId[scRNASeqAnnotation$libraryId == library]
  
  for (cellPop in AllCellPop) {
    
    cellPopName <- str_remove(cellPop, "Raw_Counts_")
    cellPopName <- str_remove(cellPopName, ".tsv")
    ## verify if cell pop information
    numberCells <- cellInfo$cells_afterFiltering[cellInfo$library == library & cellInfo$Cell_Name == cellPopName]
    
    ## just use  ell-populations with at least 10 cells in the library
    if(numberCells >= 10){
      
    ## collect the sumUMI + normalization for the target cellPop
    cellPop_normalized <- sumUMICellPop(folder_data = folder_data, library = library, cellPop = cellPop)
    
    ## Information about reference intergenic
    referenceIntergenic <- refIntergenic(counts = cellPop_normalized, folder_refIntergenic = folder_refIntergenic, speciesID = speciesID)
      
    ## calls with Bgee ratio
    callRatio <- calculate_r(counts = cellPop_normalized , seqNamesFinal = referenceIntergenic, ratioValue = desired_r_cutoff)
    cellPop_normalized$call <- ifelse(cellPop_normalized$CPM >= callRatio[[1]], "present", "-")
    
    ## calls with pValue theoretical
    calculatePvalues <- theoretical_pValue(counts = cellPop_normalized, refrenceIntergenic = referenceIntergenic)
    calculatePvalues <- calculatePvalues %>% dplyr::select(gene_id, zScore, pValue)
    ## merge all information for genic region 
    callsCellPop_Final <- merge(dplyr::filter(cellPop_normalized, type == "genic"), calculatePvalues, by = "gene_id", all = TRUE)
    callsCellPop_Final$calls_pValue <- ifelse(callsCellPop_Final$pValue <= as.numeric(desired_pValue_cutoff), "present", "-")
    ## add intergenic region to the final matrix
    justIntergenic <- dplyr::filter(cellPop_normalized, type == "intergenic")
    justIntergenic$zScore <- "NA"; justIntergenic$pValue <- "NA"; justIntergenic$calls_pValue <-"NA" 
    
    callsCellPop_Final <- rbind(callsCellPop_Final, justIntergenic) 
    
    ## write the final data_frame per cell population
    fileName <- paste0("Calls_cell_pop_", cellPopName, "_",library,".tsv")
    ifelse(!dir.exists(file.path(output, library)), dir.create(file.path(output, library)), FALSE)
    write.table(callsCellPop_Final, file = file.path(output, library, fileName), quote=F, sep = "\t", col.names=T, row.names=F)
    
    }else {
      cat("The cell population will not be used.")
    }
 }
}
