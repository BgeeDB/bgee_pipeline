## SFC, 10 Oct 2019
## Updated June 2020

## This script is used to do the initial quality control (filtering) based on the Knee plot
## also to map cell types using the information from barcode annotation OR gene markers annotation.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNA_Seq_info_target.txt" folder_data="folder_data" infoFolder="infoFolder" output="output"' QC_CellTypeIdentification.R QC_CellTypeIdentification.Rout
## scRNASeq_Info --> File that results from annotation and metadata (libraries downloaded and with extra information as readlength or SRR) 
## folder_data --> Folder where are all the libraries (after process busfile)
## infoFolder --> Folder where we have the files corresponding to barcodes and gene markers annotation
## output --> Folder where we should save the results 

## libraries used
library(devtools)
library(BUSpaRse)
library(DropletUtils)
library(ggplot2)
library(magrittr)
library(data.table)
library(Seurat)
library(Matrix)
library(lattice)
library(UniprotR)
library(biomaRt)
library(gridExtra)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("scRNASeq_Info","folder_data", "infoFolder", "output")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read scRNASeq_Info file. If file not exists, script stops
if( file.exists(scRNASeq_Info) ){
  scRNASeqAnnotation <- fread(scRNASeq_Info, h=T, sep="\t")
} else {
  stop( paste("The scRNASeq_Info file not found [", scRNASeq_Info, "]\n"))
}
## read the barcode file, stop script if the files not exist!
if( file.exists(infoFolder)){
  barcodesFile <- paste0(infoFolder, "/scRNASeq_barcode.tsv")
  barcodes <- fread(barcodesFile, header = TRUE, sep="\t")
} else {
  stop( paste("The info folder was not found in [", infoFolder, "]\n"))
}
## read the markers file, stop script if the files not exist!
if( file.exists(infoFolder)){
  markersFile <- paste0(infoFolder, "/scRNASeq_markers_10X.tsv")
  markers <- fread(markersFile, header = TRUE, sep="\t")
} else {
  stop( paste("The markers file not found in [", infoFolder, "]\n"))
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
  
  ## Filtering barcodes based on inflection point!
  m_filtered <- sparseMatrix[, tot_counts > metadata(bc_rank)$inflection]
  ## Filtering all genes that are zero in all barcodes!
  m_filtered <- m_filtered[tot_genes > 0, ]
  
  pdf(file = paste0(output, libraryID, "/kneePlot.pdf"), width = 16, height = 10)   
  grid.arrange(kneePlot,ncol = 1, nrow = 1)
  dev.off()
  
  return(m_filtered)
}

## function to create the seurat object and normalize the data (each cell CPM)
seurtObject <- function(m_filtered){
  object <- CreateSeuratObject(counts = m_filtered, project = "scRNA", min.cells = 0, min.features = 0)
  objectNormalized <- NormalizeData(object, normalization.method = "RC", scale.factor = 1e6)
  return(objectNormalized)
}

## function to target cells after knee filtering based on the barcode or gene markers information.
targetCells <- function(objectNormalized, barcodeIDs, geneMarkers, biotypeInfo , geneNameFile, libraryID){
  
  ## verify if barcode annotation exist
  lenghtBarcodeFilter <- nrow(barcodeIDs)
  
  ## select cells based on the barcode ID's
  if (lenghtBarcodeFilter != 0){
    
    onlyBarcode <- barcodeIDs[,1]
    onlyBarcode <- unlist(onlyBarcode, use.names=FALSE)
    ## subset barcodes
    myData <- subset(objectNormalized,  cells = onlyBarcode)
    
    ## verify if all are present (this is dependent of the Knee plot - quality control)
    presentSubset <- as.data.frame(colnames(myData))
    colnames(presentSubset) <- "barcode"
    #dim(presentSubset)
    
    ## change metadata: add info about the cell type for each barcode
    ## select from barcode file: barcode ID (column 1) and cell_type_harmonization (column 10)
    barcode_cellName <- barcodeIDs[,c(1,10)] 
    barcode_cellName <- merge(presentSubset, barcode_cellName, by ="barcode")
    barcode_cellName <- barcode_cellName$cell_type_harmonization
    barcode_cellName <- unlist(barcode_cellName, use.names=FALSE)
    myData$cell_type <- barcode_cellName
    #head(myData[[]])
    ## remove unsigned cells (this avoid the exportation of the file of unassigned cells)
    myData <- subset(myData,  cell_type != "Unassigned")
    
    set.seed(42)
    myData <- FindVariableFeatures(myData, selection.method = "vst", nfeatures = 2000)
    ## Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(myData), 10)
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
    pdf(paste0(output, "/", libraryID, "/UMAP_seurat_cluster.pdf"))
    print(umapPlot)
    dev.off()
    umapPlot <- DimPlot(myData, reduction = "umap", group.by='cell_type')
    ## save UMAP information with cell_type
    pdf(paste0(output, "/", libraryID, "/UMAP_cellType_cluster.pdf"))
    print(umapPlot)
    dev.off()
    
    ## collect raw UMI counts and normalized data counts for each cell
    finalRaw <- data.frame(myData@assays$RNA@counts)
    finalCPM <- data.frame(myData@assays$RNA@data)
    
    ## write in the output the info per cell type that is present in the library
    infoCollected <- c()
    for (cell in unique(myData$cell_type)) {
      ## split information
      barcodesID <- colnames(myData)[myData$cell_type == cell]
      ## export raw counts to each cell type
      rawCountsCell <- finalRaw[(names(finalRaw) %in% barcodesID)]
      rawCountsCell <- setDT(rawCountsCell, keep.rownames = TRUE)[]
      colnames(rawCountsCell)[1] <- "gene_id"
      ## add biotype info to raw counts
      collectBiotypeRaw <- merge(rawCountsCell, biotypeInfo, by = "gene_id", all.x = TRUE)
      ## add type info
      collectBiotypeRaw$type <- ifelse(is.na(collectBiotypeRaw$biotype), "intergenic", "genic")
      
      ## export normalized counts to each cell type
      normalizedCountsCell <- finalCPM[(names(finalCPM) %in% barcodesID)]
      normalizedCountsCell <- setDT(normalizedCountsCell, keep.rownames = TRUE)[]
      colnames(normalizedCountsCell)[1] <- "gene_id"
      ## add biotype info to normalized counts
      collectBiotypeNorm <- merge(normalizedCountsCell, biotypeInfo, by = "gene_id", all.x = TRUE)
      ## add type info
      collectBiotypeNorm$type <- ifelse(is.na(collectBiotypeNorm$biotype), "intergenic", "genic")
      
      ## write output information to integrate in Bgee
      write.table(collectBiotypeRaw, file = paste0(output, "/", libraryID, "/Raw_Counts_", cell, ".tsv"), sep="\t", row.names = FALSE, quote = FALSE)
      write.table(collectBiotypeNorm, file = paste0(output, "/", libraryID, "/Normalized_Counts_", cell, ".tsv"), sep="\t", row.names = FALSE, quote = FALSE)
      ## info per cell (-3 because of geneId, biotype and type columns)
      cellInfo <- c(nrow(collectBiotypeNorm), ncol(collectBiotypeNorm)-3, cell)
      infoCollected <- rbind(infoCollected, cellInfo)
    }
    infoCollected <- as.data.frame(infoCollected)
    colnames(infoCollected) <- c("Number_genes", "Number_cells", "Cell_Name")
    return(infoCollected)
  } else {
    
   cat("In this library ", libraryID, " the cell identification will be done based on gene markers", "\n")
    
    myData <- FindVariableFeatures(objectNormalized, selection.method = "vst", nfeatures = 2000)
    ## Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(myData), 10)
    ## scale the data
    all.genes <- rownames(myData)
    myData <- ScaleData(myData, features = all.genes)
    ## Run PCA 
    myData <- RunPCA(myData, features = VariableFeatures(object = myData))
    myData <- FindNeighbors(myData, dims = 1:20)
    myData <- FindClusters(myData, resolution = 1)
    ## Run UMAP
    myData <- RunUMAP(myData, dims = 1:20)
    umapPlot <- DimPlot(myData, reduction = "umap", label = TRUE, label.size = 4)
    ## save UMAP information
    pdf(paste0(output, "/", libraryID, "/UMAP_seurat_cluster.pdf"))
    print(umapPlot)
    dev.off()
    
    ## find the marker genes between clusters
    markersLibrary <- FindAllMarkers(myData, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
    
    getInfo <- merge(markersLibrary, geneNameFile, by = "gene")
    names(getInfo)[length(names(getInfo))]<-"marker_gene_name_for_cell_type" 
    getInfo <- getInfo[,c(1,8)]
    
    ## verify if this markers genes from the analysis exist in annotation file
    finalMarker <- merge(geneMarkerLibrary, getInfo, by = "marker_gene_name_for_cell_type")
    
    if (nrow(finalMarker) != 0){
      
      genesID <- finalMarker$gene
      markersPlot <- VlnPlot(myData, features = genesID)
      pdf(paste0(output, "/", libraryID, "/markers_genes_validated.pdf"), height = 15, width = 20)
      print(markersPlot)
      dev.off()
      
      ## collect information about the cluster
      collectClusterInfo <- markersLibrary[markersLibrary$gene %in% genesID, ]
      finalMarker <- data.frame(finalMarker$gene, finalMarker$marker_gene_name_for_cell_type, finalMarker$cell_type)
      colnames(finalMarker) <- c("gene", "gene_name", "cell_type")
      collectClusterInfo <- merge(collectClusterInfo, finalMarker, by="gene")
      
      ## add extra info in the metadata
      myData$cellType <- "Unassigned"
      
      ## attribute cell type to metadata
      for (i in 1:nrow(collectClusterInfo)) {
        ## collect cluster and cell type
        cluster <- collectClusterInfo[i, "cluster"]
        cellType <- collectClusterInfo[i, "cell_type"]
        ## add to cellType column
        myData$cellType[myData$seurat_clusters == cluster] <- cellType
      }
      
      ## plot UMAP after markers identify the cell-type
      umapPlot <- DimPlot(myData, reduction = "umap", group.by='cellType')
      pdf(paste0(output, "/", libraryID, "/UMAP_cellType_cluster.pdf"))
      print(umapPlot)
      dev.off()
      
      ## collect raw UMI counts and normalized data counts
      finalRaw <- data.frame(myData@assays$RNA@counts)
      finalCPM <- data.frame(myData@assays$RNA@data)
      
      ## write output per cell type that is present in the library
      infoCollected <- c()
      for (cell in unique(myData$cellType)) {
        
        if (cell == "Unassigned"){
          cat("Barcodes don't have cellType information", "\n")
        } else {
          
          ## split information
          barcodesID <- colnames(myData)[myData$cellType == cell]
          ## export raw counts to each cell type
          rawCountsCell <- finalRaw[(names(finalRaw) %in% barcodesID)]
          rawCountsCell <- setDT(rawCountsCell, keep.rownames = TRUE)[]
          colnames(rawCountsCell)[1] <- "gene_id"
          ## add biotype info to raw counts
          collectBiotypeRaw <- merge(rawCountsCell, biotypeInfo, by = "gene_id", all.x = TRUE)
          ## add type info
          collectBiotypeRaw$type <- ifelse(is.na(collectBiotypeRaw$biotype), "intergenic", "genic")
          
          ## export normalized counts to each cell type
          normalizedCountsCell <- finalCPM[(names(finalCPM) %in% barcodesID)]
          normalizedCountsCell <- setDT(normalizedCountsCell, keep.rownames = TRUE)[]
          colnames(normalizedCountsCell)[1] <- "gene_id"
          ## add biotype info to normalized counts
          collectBiotypeNorm <- merge(normalizedCountsCell, biotypeInfo, by = "gene_id", all.x = TRUE)
          ## add type info
          collectBiotypeNorm$type <- ifelse(is.na(collectBiotypeNorm$biotype), "intergenic", "genic")
          
          ## write output information to integrate in Bgee
          write.table(collectBiotypeRaw, file = paste0(output, "/", libraryID, "/Raw_Counts_", cell, ".tsv"), sep="\t", row.names = FALSE, quote = FALSE)
          write.table(collectBiotypeNorm, file = paste0(output, "/", libraryID, "/Normalized_Counts_", cell, ".tsv"), sep="\t", row.names = FALSE, quote = FALSE)
          ## info per cell (-3 because of geneId, biotype and type columns)
          cellInfo <- c(nrow(collectBiotypeNorm), ncol(collectBiotypeNorm)-3, cell)
          infoCollected <- rbind(infoCollected, cellInfo)
        }
      }
      infoCollected <- as.data.frame(infoCollected)
      colnames(infoCollected) <- c("Number_genes", "Number_cells", "Cell_Name")
      return(infoCollected)
    } else {
      cat("For this library ", libraryID, " any gene marker from the analysis match the gene markers from annotation", "\n",
          "Library not take in consideration for posterior analysis", "\n")
      ## count unique cells annotated -1 (barcodes without cellType annotation)
      infoCollected <- c(nrow(myData), ncol(myData), paste0("No match between annotation and analysis"))
      infoCollected <- data.frame(matrix(infoCollected,1,3))
      colnames(infoCollected) <- c("Number_genes", "Number_cells", "Cell_Name")
      return(infoCollected)
    }
  }
}

#############################################################################################################################
## collect information for all libraries
globalInfoLibraries <- paste0(output, "/InformationAllLibraries.txt")
if (!file.exists(globalInfoLibraries)){
  file.create(globalInfoLibraries)
  cat("library\texperimentID\tInitial_UMI_barcode\tInitial_tot_genes\tgenes_afterFiltering\tcells_afterFiltering\tCell_Name\n",file = paste0(output, "/InformationAllLibraries.txt"), sep = "\t")
} else {
  print("File already exist.....")
}

## For each library that belongs to each experiment do:
for (libraryID in scRNASeqAnnotation$libraryId) {
  
  ## verify if library exist
  pathLib <- file.exists(file.path(folder_data, libraryID))
  
  if (pathLib == TRUE){
    
    cat("Treating library ", libraryID, "\n")
    
    ## info about species and experiment
    speciesID <- scRNASeqAnnotation$scientific_name[scRNASeqAnnotation$libraryId == libraryID]
    speciesID <- gsub(" ", "_", speciesID)
    experimentID <- scRNASeqAnnotation$experimentId[scRNASeqAnnotation$libraryId == libraryID]
    
    path2Files <- paste0(folder_data, libraryID, "/busOutput/gene_counts")
    path2ListFiles <- list.files(path2Files)
    
    ## create folder per library
    dir.create(file.path(output, libraryID))
    
    ## Create a sparseMatrix 
    sparseMatrix <- read_count_output(path2Files, "gene", tcc = FALSE)
    
    ## export global information
    ## barcodes detected (cell per column) and genes detected (rows)
    ## How many UMIs per barcode (cell)
    tot_counts <- Matrix::colSums(sparseMatrix)
    ## How many genes detected
    tot_genes <- rowSums(sparseMatrix)
    
    ## filtering cells based on the knee plot
    knee <- singlecellKnee(sparseMatrix = sparseMatrix, libraryID = libraryID)
    
    ## perform the seurat object and get raw and normalized counts
    object <- seurtObject(m_filtered = knee)
    
    ## read biotype information
    biotypeInfo <- fread(paste0(infoFolder, "/gene_to_biotype_with_intergenic_", speciesID,".tsv"))
    colnames(biotypeInfo) <- c("gene_id", "biotype")
    ## get barcode information for the library
    barcodeLibrary <- dplyr::filter(barcodes, library == libraryID)
    ## get gene markers information
    geneMarkerLibrary <- dplyr::filter(markers, experimentId == experimentID)
    ## read info that link gene_id with gene_name
    geneNameFile <- fread(paste0(infoFolder, "/gene_to_geneName_with_intergenic_", speciesID,".tsv"))
    
    ## link barcode to cell ID or link cell ID based on gene marker and export information
    finalInformationCells <- targetCells(objectNormalized = object, barcodeIDs = barcodeLibrary, geneMarkers = geneMarkerLibrary, biotypeInfo = biotypeInfo, geneNameFile = geneNameFile, libraryID = libraryID)
    
    ## write all information for the library
    allinfo <- data.frame(libraryID, experimentID, length(tot_counts) , length(tot_genes), finalInformationCells$Number_genes, finalInformationCells$Number_cells, finalInformationCells$Cell_Name)
    write.table(allinfo, file = globalInfoLibraries, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
  
  } else {
      cat("This library ", libraryID, " not exist in the directory.", "\n")
    } 
  }
