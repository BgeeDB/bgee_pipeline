## SFC, 10 Oct 2019

## Updated June 2022
## manage barcode per experimentID and remove everything linked with marker genes

## This script is used to do:
## the initial quality control (filtering) based on the Knee plot,
## to map cell types using the information from barcode annotation file,
## to provide a variability info file (this means info about the proportion of cell-types per cluster)
## to provide a list of gene markers per cell-type after data analysis (validated with gene markers from the annotation file)

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNA_Seq_info_target.txt" kallisto_bus_results="kallisto_bus_results" folderSupport="folder_with_info_files" infoFolder="folder_with_barcodes_files" output="output"' QC_CelltypeIdentification.R QC_CelltypeIdentification.Rout
## scRNASeq_Info --> File that results from annotation and metadata (libraries downloaded and with extra information as readlength or SRR)
## kallisto_bus_results --> Folder where are all the libraries (after process busfile)
## folderSupport --> Folder where is saved the files: gene_to_biotype_with_intergenic_ per species ID
## infoFolder --> Folder where we have the files corresponding to barcodes per experiment ID
## output --> Folder where we should save the results

## libraries used

##TODO: remove these packages as they are not used
#library(magrittr)
#library(lattice)
#library(tibble)
#library(devtools)

library(Matrix)
library(BUSpaRse)
library(DropletUtils)
library(ggplot2)
library(gridExtra)
library(Seurat)
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
command_arg <- c("scRNASeq_Info", "kallisto_bus_results", "folderSupport", "infoFolder", "output")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read scRNASeq_Info file. If file not exists, script stops
if( file.exists(scRNASeq_Info) ){
  scRNASeqAnnotation <- read.table(scRNASeq_Info, header = TRUE, sep = "\t", quote = "\"")
} else {
  stop( paste("The scRNASeq_Info file not found [", scRNASeq_Info, "]\n"))
}

##########################################################################################################################################################
## function to filter the barcodes based on the inflection point
singlecellKnee <- function(sparseMatrix, libraryID){
  ## rank barcodes to do the knee
  bc_rank <- barcodeRanks(sparseMatrix)
  kneePlot <- qplot(bc_rank$total, bc_rank$rank, geom = "line") +
    geom_vline(xintercept = metadata(bc_rank)$knee, color = "blue", linetype = 2) +
    geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype = 2) +
    annotate("text", y = 10000, x = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
             label = c("knee", "inflection"), color = c("blue", "green")) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = "Barcode rank", x = "Total UMI count")
 
  ## Filtering barcodes based on inflection point.
  m_filtered <- sparseMatrix[, tot_counts > metadata(bc_rank)$inflection]
  ## Filtering all genes that are zero in all barcodes!
  ##m_filtered <- m_filtered[tot_genes > 0, ]
  
  pdf(file = paste0(output, libraryID, "/kneePlot.pdf"), width = 16, height = 10)
  grid.arrange(kneePlot,ncol = 1, nrow = 1)
  dev.off()
  
  return(m_filtered)
}

## function to create the seurat object and normalize the data (each cell CPM)
seuratObject <- function(m_filtered){
  object <- CreateSeuratObject(counts = m_filtered, project = "scRNA", min.cells = 0, min.features = 0)
  objectNormalized <- NormalizeData(object, normalization.method = "RC", scale.factor = 1e6)
  return(objectNormalized)
}

## function to target cells after knee filtering based on the barcode information (provide clustering information)
targetCells <- function(objectNormalized, barcodeIDs, biotypeInfo, libraryID){
  
  ## select cells based on the barcode ID's
  if (nrow(barcodeIDs) != 0){
    
    ##create data.frame mapping celltypeId to corresponding files
    ## this file will have 3 columns
    ## celltypeId       countType(raw or normalized)     filePath
    mappingCellTypeToFiles <- data.frame()

    onlyBarcode <- barcodeIDs[,1]
    onlyBarcode <- unlist(onlyBarcode, use.names=FALSE)
    ## subset barcodes
    myData <- subset(objectNormalized,  cells = onlyBarcode)
    
    ## verify if all are present (this is dependent of the Knee plot - quality control)
    presentSubset <- as.data.frame(colnames(myData))
    colnames(presentSubset) <- "barcode"
    #dim(presentSubset)
    
    ## change metadata: add info about the cell type for each barcode
    ## select from barcode file: barcode ID (column 1) and cellTypeName (column 9)
    barcode_cell <- barcodeIDs %>% dplyr::select(barcode, cellTypeName, cellTypeId)
    barcode_cell <- barcode_cell[ barcode_cell$barcode %in% presentSubset$barcode, ]
    
    barcode_cellName <- barcode_cell$cellTypeName
    barcode_cellName <- unlist(barcode_cellName, use.names=FALSE)
    barcode_cellId <- barcode_cell$cellTypeId
    barcode_cellId <- unlist(barcode_cellId, use.names=FALSE)
    myData$cell_type <- barcode_cellName
    myData$cell_id <- barcode_cellId
    #head(myData[[]])
    ## remove unsigned cells (this avoid the exportation of the file of unassigned cells)
    myData <- subset(myData,  cell_type != "Unassigned")
    
    set.seed(42)
    myData <- FindVariableFeatures(myData, selection.method = "vst", nfeatures = 2000)
    ## Identify the 20 most highly variable genes
    top20 <- head(VariableFeatures(myData), 20)
    ## scale the data
    all.genes <- rownames(myData)
    myData <- ScaleData(myData, features = all.genes)
    ## Run PCA
    myData <- RunPCA(myData, features = VariableFeatures(object = myData))
    myData <- FindNeighbors(myData, dims = 1:20)
    myData <- FindClusters(myData, resolution = 1)
    ## Run UMAP
    myData <- RunUMAP(myData, dims = 1:20)
    umapPlot <- DimPlot(myData, reduction = "umap", group.by='seurat_clusters')
    pdf(file.path(output, libraryID, "UMAP_seurat_cluster.pdf"))
    print(umapPlot)
    dev.off()
    umapPlot <- DimPlot(myData, reduction = "umap", group.by='cell_type')
    ## save UMAP information with cell_type
    pdf(file.path(output, libraryID, "UMAP_cellType_cluster.pdf"))
    print(umapPlot)
    dev.off()
    umapPlot <- DimPlot(myData, reduction = "umap", group.by='cell_id')
    ## save UMAP information with cell_ontology
    pdf(file.path(output, libraryID, "UMAP_cellID_cluster.pdf"))
    print(umapPlot)
    dev.off()
    
    ## collect raw UMI counts and normalized data counts for each cell
    finalRaw <- data.frame(myData@assays$RNA@counts)
    finalCPM <- data.frame(myData@assays$RNA@data)
    
    ## write in the output the info per cell type that is present in the library
    infoCollected <- c()
    for (cell in unique(myData$cell_type)) {
      for (cellid in unique(myData$cell_id[myData$cell_type == cell])) {
        message("cell : ", cellid)
        ## split information
        barcodesID <- colnames(myData)[myData$cell_type == cell & myData$cell_id == cellId] 
        ## export raw counts to each cell type
        rawCountsCell <- finalRaw[(names(finalRaw) %in% barcodesID)]
        rawCountsCell <- cbind(names = rownames(rawCountsCell), rawCountsCell)
        colnames(rawCountsCell)[1] <- "gene_id"
        ## add biotype info to raw counts
        collectBiotypeRaw <- merge(rawCountsCell, biotypeInfo, by = "gene_id", all.x = TRUE)
        ## add type info
        collectBiotypeRaw$type <- ifelse(is.na(collectBiotypeRaw$biotype), "intergenic", "genic")
        collectBiotypeRaw$cellTypeName <- cell
        collectBiotypeRaw$cellTypeId <- cellId
        
        ## export normalized counts to each cell type
        normalizedCountsCell <- finalCPM[(names(finalCPM) %in% barcodesID)]
        normalizedCountsCell <- cbind(names = rownames(normalizedCountsCell), normalizedCountsCell)
        colnames(normalizedCountsCell)[1] <- "gene_id"
        ## add biotype info to normalized counts
        collectBiotypeNorm <- merge(normalizedCountsCell, biotypeInfo, by = "gene_id", all.x = TRUE)
        ## add type info
        collectBiotypeNorm$type <- ifelse(is.na(collectBiotypeNorm$biotype), "intergenic", "genic")
        collectBiotypeNorm$cellTypeName <- cell
        ## write output information to integrate in Bgee
        rawCountFilePath <- file.path(output, libraryID, paste0("Raw_Counts_", cellIdModified,
          ".tsv"))
        normalizedCountFilePath <- file.path(output, libraryID, paste0("Normalized_Counts_",
          cellIdModified, ".tsv"))
        write.table(collectBiotypeRaw, file = rawCountFilePath, sep="\t", row.names = FALSE, 
          quote = FALSE)
        write.table(collectBiotypeNorm, file = normalizedCountFilePath, sep="\t", row.names = FALSE,
          quote = FALSE)

        ## add file info in mapping data.frame
        mappingCellTypeToFiles <- rbind(mappingCellTypeToFiles, as.data.frame(c>(cellId, "raw",
          rawCountFilePath)))
        mappingCellTypeToFiles <- rbind(mappingCellTypeToFiles, as.data.frame(c>(cellId,
          "normalized", normalizedCountFilePath)))
        
        ## info per cell (-5 because of geneId, biotype, type columns and cellName and cellID)
        cellInfo <- c(nrow(collectBiotypeNorm), ncol(collectBiotypeNorm)-5, cell, cellId)
        infoCollected <- rbind(infoCollected, cellInfo)
      }
    }
    infoCollected <- as.data.frame(infoCollected)
    colnames(infoCollected) <- c("Number_genes", "Number_cells", "Cell_Name", "Cell_Ontology")

    colnames(mappingCellTypeToFiles) <- c("cellTypeId", "countType", "filePath")
    write.table(mappingCellTypeToFiles, file.path(output, libraryId, "mappingCellTypeToFiles.tsv"),
      sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

    return(list(myData[[]],infoCollected))
  } else {
    message("Library", libraryID, "not take in consideration for posterior analysis.")
    message("Barcode not provided!")
  }
}

#############################################################################################################################
## collect information for all libraries
globalInfoLibraries <- file.path(output, "InformationAllLibraries.txt")
file.create(globalInfoLibraries)

fileHeader <- c("library", "experimentID", "Initial_UMI_barcode",
  "Initial_tot_genes", "genes_afterFiltering", "cells_afterFiltering",
  "cellTypeName", "cellTypeId")

allInfo <- data.frame()
## For each library that belongs to each experiment do:
for (libraryID in unique(scRNASeqAnnotation$libraryId)) {
  
  ## verify if library exist
  if (file.exists(file.path(kallisto_bus_results, libraryID, "gene_counts"))){

    message("Treating library ", libraryID)
    
    ## info about species and experiment
    speciesName <- unique(scRNASeqAnnotation$scientific_name[scRNASeqAnnotation$libraryId == libraryID])
    speciesName <- gsub(" ", "_", speciesName)
    experimentID <- unique(scRNASeqAnnotation$experimentId[scRNASeqAnnotation$libraryId == libraryID])

    path2Files <- file.path(kallisto_bus_results, libraryID, "gene_counts/")
    
    ## create folder per library
    dir.create(file.path(output, libraryID))
    ## Create a sparseMatrix
    sparseMatrix <- read_count_output(path2Files, "gene", tcc = FALSE)
    message("created sparse matrix")
    ## export global information
    ## barcodes detected (cell per column) and genes detected (rows)
    ## How many UMIs per barcode (cell)
    tot_counts <- Matrix::colSums(sparseMatrix)
    message("created total counts variable")
    ## How many genes detected
    tot_genes <- rowSums(sparseMatrix)
    message("created total genes variable")
    ## filtering cells based on the knee plot
    knee <- singlecellKnee(sparseMatrix = sparseMatrix, libraryID = libraryID)
    message("created kneeplot")
    ## perform the seurat object and get raw and normalized counts
    object <- seuratObject(m_filtered = knee)
  
    ## read biotype information
    biotypeInfoFile <- file.path(folderSupport, paste0(speciesName, "_gene_to_biotype_with_intergenic.tsv"))
    # the file is potentially compressed. Have to uncompress it first
    if (!file.exists(biotypeInfoFile)) {
      if (file.exists(paste0(biotypeInfoFile, ".xz"))) {
        system(sprintf("unxz %s", paste0(biotypeInfoFile,".xz")))
      } else {
        stop("gene to biotype file [", biotypeInfoFile, "] does not exist.")
      }
    }
    biotypeInfo <- read.table(biotypeInfoFile)
    ## hack needed for Bgee 15.1 as the gene to biotype file contain redundant genes
    ##TODO: update script generating gene to biotype
    biotypeInfo <- unique(biotypeInfo)
    colnames(biotypeInfo) <- c("gene_id", "biotype")
    ## get barcode information per experiment
    barcodeExperiment <- read.table(file.path(infoFolder, paste0("scRNASeq_barcode_", experimentID,".tsv")), header = TRUE,
      sep = "\t")
    barcodeLibrary <- dplyr::filter(barcodeExperiment, library == libraryID)
    
    ## link barcode to cell ID and export information
    finalInformationCells <- targetCells(objectNormalized = object, barcodeIDs = barcodeLibrary, biotypeInfo = biotypeInfo, libraryID = libraryID)
    infoCells <- finalInformationCells[[2]]
    
    ## verify if function targetCells is NULL
    if (is.null(finalInformationCells)){
      message("The library", libraryID, " not contain barcodes.")
      message("Any raw or normalized count is exported.")
      message("The library will not be present in the informative files: InformationAllLibraries.txt, variabilityClusters.txt and markerGenes_Validated.txt file")
    } else {
      ## write all information for the library
      info <- data.frame(libraryID, experimentID, length(tot_counts) , length(tot_genes), infoCells$Number_genes, infoCells$Number_cells, infoCells$Cell_Name, infoCells$Cell_Ontology)
      allInfo <- rbind(allInfo, info)
    }
  } else {
      warning("The library ", libraryID, " has not been porcessed by bustools")
  }
}
colnames(allInfo) <- fileHeader
write.table(allinfo, file = globalInfoLibraries, quote = FALSE, sep = "\t", col.names = TRUE,
  row.names = FALSE)
