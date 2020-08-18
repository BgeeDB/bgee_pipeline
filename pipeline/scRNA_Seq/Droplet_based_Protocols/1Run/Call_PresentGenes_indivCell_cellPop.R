## SFC, 21 Oct 2019
## Updated June 2020

## This script is used to do the call of present genes at individual cell and to do the calls at the cell population level.
## the aggregation of cells is done by cells that belongs to the same experiment, uberonId, cellTypeId, stageId, strain and sex and then the calls is done
## by a fraction of cells that have the gene classified as present.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNASeq_Info.txt" InformationAllLibraries="InformationAllLibraries.txt" folder_data="folder_data" desired_r_cutoff="desired_r_cutoff" output="output"' Call_PresentGenes_indivCell_cellPop.R Call_PresentGenes_indivCell_cellPop.Rout
## scRNASeq_Info --> File that results from annotation and metadata (libraries downloaded and with extra information as readlength or SRR) 
## InformationAllLibraries --> File with information about each cell after barcode and gene markers annotation (per library contain total number of cells and cell Name)
## folder_data --> Folder where are all the libraries after cell identification
## desired_r_cutoff --> proportion of intergenic allowed (both individual cell and cell population to define the ratio cutoff)
## output --> Folder where we should save the results 

## libraries used
library(ggplot2)
library(gridExtra)
library(gridExtra)
library(data.table)
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
command_arg <- c("scRNASeq_Info", "InformationAllLibraries", "folder_data", "desired_r_cutoff", "output")
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
### Creat file to save information per individual cell
presentCell <- paste0(output, "/Present_Info_individual_Cell.tsv")
if (!file.exists(presentCell)){
  file.create(presentCell)
  cat("barcode\texperiment\tcellName\tCPM_cutoff\tNumber_genic\tproportionGenic\tNumber_intergenic\tproportionIntergenic\tNumber_proteinCoding\tproportionproteinCoding\tratio\tspecies\n",file = paste0(output, "/Present_Info_individual_Cell.tsv"), sep = "\t")
} else {
  print("File already exist.....")
}

### Creat file to save information per cell population
presentCellPopulation <- paste0(output, "/Present_Info_CellPopulation.tsv")
if (!file.exists(presentCellPopulation)){
  file.create(presentCellPopulation)
  cat("experiment\tcellName\tgenic_all\tproportion_genic_all\tintergenic_all\tproportion_intergenic_all\tproteinCoding_all\tproportion_PC_all\tgenic_density\tproportion_genic_density\tintergenic_density\tproportion_intergenic_density\tproteinCoding_density\tproportion_PC_density\tgenic_ratio\tproportion_genic_ratio\tintergenic_ratio\tproportion_intergenic_ratio\tproteinCoding_ratio\tproportion_PC_ratio\tcutoff_density\tcutoff_ratio\tspeciesId\n", file = paste0(output, "/Present_Info_CellPopulation.tsv"), sep = "\t")
} else {
  print("File already exist.....")
}

## export information per individual cell
exportInfo_Cell <- function(cellCalls){
  
  genic <- nrow(dplyr::filter(cellCalls, type == "genic" & cellCalls[,ncol(cellCalls)] == "present"))
  proportionGenic <- genic / nrow(dplyr::filter(cellCalls, type == "genic"))*100
  
  intergenic <- nrow(dplyr::filter(cellCalls, type == "intergenic" & cellCalls[,ncol(cellCalls)] == "present"))
  proportionIntergenic <- intergenic / nrow(dplyr::filter(cellCalls, type == "intergenic"))*100
  
  proteinCoding <- nrow(dplyr::filter(cellCalls, biotype == "protein_coding" & cellCalls[,ncol(cellCalls)] == "present"))
  proportionproteinCoding <- proteinCoding / nrow(dplyr::filter(cellCalls, biotype == "protein_coding"))*100
  
  collectInfoData <- data.frame(genic, proportionGenic, intergenic, proportionIntergenic, proteinCoding, proportionproteinCoding)
  return(unlist(collectInfoData, use.names=FALSE))
  
}

## export information per cell population
exportInfo_CellPop <- function(collectInfoAllCells){
  
  ## info genes present after density cutoff (include density cutoff + ratio cutoff)
  genic_all <- nrow(dplyr::filter(collectInfoAllCells, type == "genic" & call_cutoff_density == "present")) 
  proportion_genic_all <- genic_all / nrow(dplyr::filter(collectInfoAllCells, type == "genic"))*100
  intergenic_all <- nrow(dplyr::filter(collectInfoAllCells, type == "intergenic" & call_cutoff_density == "present")) 
  proportion_intergenic_all <- intergenic_all / nrow(dplyr::filter(collectInfoAllCells, type == "intergenic"))*100
  proteinCoding_all <- nrow(dplyr::filter(collectInfoAllCells, biotype == "protein_coding" & call_cutoff_density == "present")) 
  proportion_PC_all <- proteinCoding_all / nrow(dplyr::filter(collectInfoAllCells, biotype == "protein_coding"))*100
  
  ## info genes present between density cutoff and ratio cutoff
  genic_density <- nrow(dplyr::filter(collectInfoAllCells, type == "genic" & call_cutoff_density == "present" & call_cutoff_ratio == "-")) 
  proportion_genic_density <- genic_density / nrow(dplyr::filter(collectInfoAllCells, type == "genic"))*100
  intergenic_density <- nrow(dplyr::filter(collectInfoAllCells, type == "intergenic" & call_cutoff_density == "present" & call_cutoff_ratio == "-")) 
  proportion_intergenic_density <- intergenic_density / nrow(dplyr::filter(collectInfoAllCells, type == "intergenic"))*100
  proteinCoding_density <- nrow(dplyr::filter(collectInfoAllCells, biotype == "protein_coding" & call_cutoff_density == "present" & call_cutoff_ratio == "-")) 
  proportion_PC_density <- proteinCoding_density/ nrow(dplyr::filter(collectInfoAllCells, biotype == "protein_coding"))*100
  
  ## info genes present just in ratio cutoff
  genic_ratio <- nrow(dplyr::filter(collectInfoAllCells, type == "genic" & call_cutoff_ratio == "present")) 
  proportion_genic_ratio <- genic_ratio / nrow(dplyr::filter(collectInfoAllCells, type == "genic"))*100
  intergenic_ratio <- nrow(dplyr::filter(collectInfoAllCells, type == "intergenic" & call_cutoff_ratio == "present")) 
  proportion_intergenic_ratio <- intergenic_ratio / nrow(dplyr::filter(collectInfoAllCells, type == "intergenic"))*100
  proteinCoding_ratio <- nrow(dplyr::filter(collectInfoAllCells, biotype == "protein_coding" & call_cutoff_ratio == "present")) 
  proportion_PC_ratio <- proteinCoding_ratio / nrow(dplyr::filter(collectInfoAllCells, biotype == "protein_coding"))*100

  collectInfoData <- data.frame(genic_all, proportion_genic_all, intergenic_all, proportion_intergenic_all, proteinCoding_all, proportion_PC_all, 
                                genic_density, proportion_genic_density, intergenic_density, proportion_intergenic_density, proteinCoding_density, proportion_PC_density, 
                                genic_ratio, proportion_genic_ratio, intergenic_ratio, proportion_intergenic_ratio, proteinCoding_ratio, proportion_PC_ratio)
  return(unlist(collectInfoData, use.names=FALSE))
}

## function to call present genes at cell population level depending on the cuoff applied
plotCellPop <- function(collectInfoAllCells, desired_r_cutoff, cell_Name, sizeData){
  
  selected_coding <- collectInfoAllCells$biotype %in% "protein_coding"
  selected_intergenic <- collectInfoAllCells$type %in% "intergenic"
  
  ## intergenic 
  summed_intergenic <- sapply(unique(sort(collectInfoAllCells$ratio[selected_coding])), function(x){
    return( sum(collectInfoAllCells$ratio[selected_intergenic] >= x) )
  })
  
  ## protein coding
  summed_coding <- c(0, cumsum(rle(sort(collectInfoAllCells$ratio[selected_coding]))$lengths))
  summed_coding <- summed_coding[-(length(summed_coding))]
  summed_coding <- sum(selected_coding) - summed_coding
  
  ## Now we can calculate r
  r <- ( summed_intergenic / sum(selected_intergenic) ) /
    ( summed_coding / sum(selected_coding) )
  
  percent <- (1-desired_r_cutoff)*100
  
  ## Select the minimal value of TPM for which the ratio of genes and intergenic regions is equal to 0.05 or lower (first test if at least 1 TPM value has this property):
  if (sum(r < desired_r_cutoff) == 0){
    ratio_cutoff <- sort(unique(collectInfoAllCells$ratio[selected_coding]))[which(r == min(r))[1]]
    r_cutoff <- min(r)
    cat(paste0("    There ratio cutoff for which " , percent,"%", " of the expressed genes would be coding. Ratio cutoff is fixed at the first value with maximum coding/intergenic ratio. r=", r_cutoff, " at ratio=", ratio_cutoff,"\n"))
  } else {
    ratio_cutoff <- sort(unique(collectInfoAllCells$ratio[selected_coding]))[which(r < desired_r_cutoff)[1]]
    r_cutoff <- desired_r_cutoff
    cat(paste0("    The ratio cutoff for which " , percent,"%", " of the expressed genes are be coding found at ratio=", ratio_cutoff,"\n"))
  }

  ProteinCoding_density <- as.data.frame(collectInfoAllCells$ratio[collectInfoAllCells$biotype == "protein_coding"])
  colnames(ProteinCoding_density) <- 'ratio'
  ProteinCoding_density <- data.frame(region=rep(c('protein_coding'),each=length(ProteinCoding_density)), ProteinCoding_density)
  
  intergenic_density <- as.data.frame(collectInfoAllCells$ratio[collectInfoAllCells$type == "intergenic"])
  colnames(intergenic_density) <- 'ratio'
  intergenic_density <- data.frame(region=rep(c('intergenic'),each=length(intergenic_density)), intergenic_density)
  
  finalTableInfo <- rbind(ProteinCoding_density, intergenic_density)
  data_ggplot <- reshape2::melt(finalTableInfo)

  p <- ggplot(data_ggplot) + geom_density(aes(x = log(value+0.000001), colour = region)) 
  p <- ggplot_build(p) 
  
  extractValues <- p$data[[1]]
  intergenicGroup <- dplyr::filter(extractValues, group == "1")
  codingGroup <- dplyr::filter(extractValues, group == "2")
  
  generalTab <- data.frame(intergenicGroup$x,intergenicGroup$density, codingGroup$density)
  generalTab$Int_less_Pc <- generalTab$intergenicGroup.density < generalTab$codingGroup.density
  
  first_cutoff <- generalTab$intergenicGroup.x[generalTab$Int_less_Pc == "TRUE" & generalTab$intergenicGroup.x >= log(1/sizeData)][1]
  
  return(list(p, exp(first_cutoff), ratio_cutoff))
}


### function to call present genes at individual cell and to export final file at cell population level
callPresentGenes <- function(finalTable, desired_r_cutoff, experiment, cell_Name, species){
  
  collectInfoAllCells <- data.frame()
  for (i in 4:ncol(finalTable)) {
    
    barcode <- names(finalTable)[i]
    infoType <- finalTable[,1:3]
    
    selected_coding <- finalTable$biotype %in% "protein_coding"
    selected_intergenic <- finalTable$type %in% "intergenic"
    
    summed_intergenic <- sapply(unique(sort(finalTable[,i][selected_coding])), function(x){
      return(sum(finalTable[,i][selected_intergenic] >= x) )})
    
    summed_coding <- c(0, cumsum(rle(sort(finalTable[,i][selected_coding]))$lengths))
    summed_coding <- summed_coding[-(length(summed_coding))]
    summed_coding <- sum(selected_coding) - summed_coding
    
    
    r <- ( summed_intergenic / sum(selected_intergenic) ) /  ( summed_coding / sum(selected_coding) )
    percent <- (1-desired_r_cutoff)*100
    
    if (sum(r < desired_r_cutoff) == 0){
      CPM_cutoff <- sort(unique(finalTable[,i][selected_coding]))[which(r == min(r))[1]]
      r_cutoff <- min(r)
      cat(paste0(" There is no CPM cutoff for which " , percent,"%", " of the expressed genes would be coding. CPM cutoff is fixed at the first value with maximum coding/intergenic ratio. r=", r_cutoff, " at CPM=", CPM_cutoff,"\n"))
    } else {
      ## else the sum of (r < desired_r_cutoff) is ≠ 0 the TPM cutoff will be the first row where the value of r is < then the desired_r_cutoff
      CPM_cutoff <- sort(unique(finalTable[,i][selected_coding]))[which(r < desired_r_cutoff)[1]]
      r_cutoff <- desired_r_cutoff
      cat(paste0(" CPM cutoff for which " , percent,"%", " of the expressed genes are be coding found at CPM=", CPM_cutoff,"\n"))
    }
    
    ## collect calls per cell
    cellCalls <- as.data.frame(ifelse(finalTable[,i]>=CPM_cutoff, "present", "-"))
    cellCalls <- cbind(infoType, cellCalls)
    colnames(cellCalls)[4] <- barcode
   
    ### collect information about each individual cell (barcode)
    infoProportion <- exportInfo_Cell(cellCalls = cellCalls)
    finalInfoProportion <- as.data.frame(t(c(barcode, experiment, cell_Name, CPM_cutoff, infoProportion, desired_r_cutoff, species)))
    write.table(finalInfoProportion, file = file.path(presentCell), append = TRUE, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    if(ncol(collectInfoAllCells) == 0){
      collectInfoAllCells = cellCalls
    } else{
      Id <- barcode
      cellCalls <- data.frame(cellCalls[,4])
      colnames(cellCalls) <- Id
      collectInfoAllCells <- cbind(collectInfoAllCells, cellCalls)
    }
  }
  
  ## collect information per cell population
  collectInfoAllCells$cellPresent <- rowSums(collectInfoAllCells == "present")
  collectInfoAllCells$ratio <- collectInfoAllCells$cellPresent / sizeData
  
  ### export plot per cell population with cut-offs
  cellPop <- plotCellPop(collectInfoAllCells = collectInfoAllCells, desired_r_cutoff = desired_r_cutoff, cell_Name = cell_Name, sizeData = sizeData)

  ## define cutoffs
  collectInfoAllCells$call_cutoff_density <- ifelse(collectInfoAllCells$ratio >= cellPop[[2]], "present", "-")
  collectInfoAllCells$call_cutoff_ratio <- ifelse(collectInfoAllCells$ratio >= cellPop[[3]], "present", "-")
  collectInfoAllCells$confidence <- ifelse(collectInfoAllCells$ratio < cellPop[[2]], "-", ifelse(collectInfoAllCells$ratio >= cellPop[[2]] & collectInfoAllCells$ratio < cellPop[[3]], "low", "high"))
  
  ## write information per cell population
  cellPopulationInfo <-  exportInfo_CellPop(collectInfoAllCells = collectInfoAllCells)
    
  ## plot distribution with cutoffs
  p1 <- cellPop[[1]]$plot +
    scale_color_manual(values=c("#00BFC4", "#F8766D")) +
    labs(title = "Density distribution of protein coding vs intergenic present", x="Log(ratio)",y="Density") +
    geom_vline(aes(xintercept=log(cellPop[[2]]), linetype="dotted"), color="black", size=0.5) + 
    geom_vline(aes(xintercept=log(cellPop[[3]]), linetype="solid"), color = "black", size=0.5) +
    scale_linetype_manual(name = "Cut-Off", values = c("dotted", "solid"), labels = c(paste0("Cut-off density (", round(cellPop[[2]]*sizeData, digits = 0), " cell/s )"), paste0("Cut-off ratio (", round(cellPop[[3]]*sizeData, digits = 0) , " cells)")))
  ## plot percentage of regions for this celltype
  infocellpop <- as.data.frame(cellPopulationInfo)
  rownames(infocellpop) <- c("genic_all", "proportion_genic_all", "intergenic_all", "proportion_intergenic_all", "proteinCoding_all", "proportion_PC_all", 
                                    "genic_density", "proportion_genic_density", "intergenic_density", "proportion_intergenic_density", "proteinCoding_density", "proportion_PC_density", 
                                    "genic_ratio", "proportion_genic_ratio", "intergenic_ratio", "proportion_intergenic_ratio", "proteinCoding_ratio", "proportion_PC_ratio")
 infocellpop$region <- rownames(infocellpop)
 subsetInfo <- subset(infocellpop,region %in% c("proportion_PC_all" , "proportion_PC_density", "proportion_PC_ratio"))
  
 p2 <- ggplot(data = subsetInfo, aes(x = region, y = cellPopulationInfo)) +
   geom_point() +
   geom_hline(yintercept=c(70,80), linetype="dashed", color = "red") +
   labs(title = "Proportion of coding genes present", x=" ",y="%") +
   coord_cartesian(ylim = c(0, 100))
 
  pdf(file.path(output, paste0("Cell_Population_",cell_Name,".pdf")), width = 14, height = 7)
  grid.arrange(p1,p2,ncol = 2, nrow = 1)
  dev.off()
  
  ## export info about the cell at population level
  finalInfoCellPop <- as.data.frame(t(c(experiment, cell_Name, cellPopulationInfo ,cellPop[[2]], cellPop[[3]], species)))
  write.table(finalInfoCellPop, file = file.path(presentCellPopulation), append = TRUE, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  return(collectInfoAllCells)
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
          
           ## collect cell-type information from info file (and ignore cell-types not identified)
          cell_Name <- unique(cellInfo$Cell_Name[cellInfo$experimentID == experiment & cellInfo$Cell_Name != "No match between annotation and analysis"])

             for (cell in cell_Name) {
               file <- as.data.frame(paste0(folder_data,  librariesInfo, "/Normalized_Counts_", cell, ".tsv"))
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
               ## read all files for the correspondent cell type
               if (length(fileTRUE) == 0){
                 
                 cat("Libraries not contain this cell type with this correspondent condictions!", "\n")
                 
               } else if (length(fileTRUE) == 1){
                 All_libs <- read.table(fileTRUE, header = TRUE, sep = "\t")
                 ## finalTableis just in this case one library
                 finalTable <- All_libs
                 ## re-order table in this case
                 finalTableIdentifiers <- finalTable[ ,c("gene_id", "biotype", "type")]
                 finalTable$gene_id <- finalTable$biotype <- finalTable$type <- NULL
                 finalTable <- data.frame(finalTableIdentifiers, finalTable)
                 finalTable[is.na(finalTable)] <- 0
                 ## remove rows with zeros in all cells
                 finalTable <- finalTable[which(rowSums(finalTable[,4:ncol(finalTable)]) > 0), ]
                 sizeData <- ncol(finalTable)-3
                 
                 ## collect information per cell/ cell population and plots
                 calls_Inv_cellPop <- callPresentGenes(finalTable = finalTable, desired_r_cutoff = as.numeric(desired_r_cutoff), experiment = experiment, cell_Name = cell, species = species)
                 ## adj output name
                 fileName <- paste0("calls_InvidCell_cellPopulation_", cell , "_", experiment,"_", uberonId,"_", cellId,"_", stageId,"_", sex,"_", strain,"_", species,".tsv")
                 fileName <- gsub(":","-",fileName)
                 fileName <- gsub("/","-",fileName)
                 ## write information of the calls
                 write.table(calls_Inv_cellPop, file = paste0(output, "/" ,fileName), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
                 
               } else {
                 All_libs <- lapply(fileTRUE, read.delim, stringsAsFactors=FALSE)
                 ## merge by geneID from all samples that have the same cell-type
                 finalTable <- Reduce(function(...) merge(..., by = c("gene_id", "biotype","type"), all=TRUE), All_libs)
                 finalTable[is.na(finalTable)] <- 0
                 ## remove rows with zeros in all cells
                 finalTable <- finalTable[which(rowSums(finalTable[,4:ncol(finalTable)]) > 0), ]
                 sizeData <- ncol(finalTable)-3
                 
                 ## collect information per cell/ cell population and plots
                 calls_Inv_cellPop <- callPresentGenes(finalTable = finalTable, desired_r_cutoff = as.numeric(desired_r_cutoff), experiment = experiment, cell_Name = cell, species = species)
                 ## adj output name
                 fileName <- paste0("calls_InvidCell_cellPopulation_", cell , "_", experiment,"_", uberonId,"_", cellId,"_", stageId,"_", sex,"_", strain,"_", species,".tsv")
                 fileName <- gsub(":","-",fileName)
                 fileName <- gsub("/","-",fileName)
                 ## write information of the calls
                 write.table(calls_Inv_cellPop, file = paste0(output, "/" ,fileName), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
                 
               }
             }
          
         }
        }
       }
      }
    }
  }
}

## final plots: per individual cell per species + plot per cell population per species
## individual cells
pdf(file.path(output, paste0("Plot_distribution_individual_cells.pdf")), width = 14, height = 7)
indCell <- fread(file.path(output, "Present_info_individual_Cell.tsv"))
proportion_all_samples <- indCell[,c(6,10,8,12)]
proportion_all_samples$species <- as.factor(proportion_all_samples$species)
proportion_all_samples <- reshape2::melt(proportion_all_samples, id=c("species"))

g1 <- ggplot(proportion_all_samples, aes(species,value, fill = species)) +
  geom_boxplot() + stat_smooth() +
  coord_cartesian(ylim = c(0, 100)) + facet_wrap(~variable) +
  geom_hline(yintercept=c(70,80), linetype="dashed", color = "black") +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab("%")
g1
dev.off()

## cell population
pdf(file.path(output, paste0("Plot_distribution_cell_population.pdf")), width = 16, height = 14)
cellPopulation <- fread(file.path(output, "Present_Info_CellPopulation.tsv"))
proportion_all_samples <- cellPopulation %>% select(experiment, cellName, proportion_PC_all,
                                                    proportion_PC_density, proportion_PC_ratio, speciesId)
proportion_all_samples$speciesId <- as.factor(proportion_all_samples$speciesId)
proportion_all_samples <- reshape2::melt(proportion_all_samples, id=c("speciesId","experiment", "cellName"))

g1 <- ggplot(proportion_all_samples, aes(x = speciesId, y = value, fill = variable, color = variable)) +
  geom_boxplot() + coord_cartesian(ylim = c(0, 100)) +
  geom_hline(yintercept=c(70,80), linetype="dashed", color = "black") +
  labs(title = paste0("Cell population - Fraction cells (sum calls)"), x="Organism",y="% protein-coding genes present")

g2 <- ggplot(proportion_all_samples, aes(x = experiment, y = value, fill = cellName, color = cellName)) +
  geom_point() + stat_smooth() +
  coord_cartesian(ylim = c(0, 100)) + facet_wrap(~variable) +
  geom_hline(yintercept=c(70,80), linetype="dashed", color = "black") +
  theme(legend.position="bottom") +
  labs(title = paste0("Cell population - Fraction cells (sum calls)"), x="Organism",y="% protein-coding genes present")

grid.arrange(g1,g2,ncol = 1, nrow = 2)
dev.off()

