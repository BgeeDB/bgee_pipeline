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
library(devtools)
library(BUSpaRse)
library(DropletUtils)
library(ggplot2)
library(magrittr)
library(data.table)
library(Matrix)
library(lattice)
library(gridExtra)
library(dplyr)
library(Seurat)
library(tibble)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("scRNASeq_Info","kallisto_bus_results", "folderSupport","infoFolder", "output")
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
seurtObject <- function(m_filtered){
  object <- CreateSeuratObject(counts = m_filtered, project = "scRNA", min.cells = 0, min.features = 0)
  objectNormalized <- NormalizeData(object, normalization.method = "RC", scale.factor = 1e6)
  return(objectNormalized)
}

## function to target cells after knee filtering based on the barcode information (provide clustering information)
targetCells <- function(objectNormalized, barcodeIDs, biotypeInfo, libraryID){
  
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
    barcode_cell <- barcodeIDs %>% dplyr::select(barcode, cell_type_harmonization, cellTypeId)
    barcode_cell <- barcode_cell[ barcode_cell$barcode %in% presentSubset$barcode, ]
    
    barcode_cellName <- barcode_cell$cell_type_harmonization
    barcode_cellName <- unlist(barcode_cellName, use.names=FALSE)
    barcode_cellid <- barcode_cell$cellTypeId
    barcode_cellid <- unlist(barcode_cellid, use.names=FALSE)
    myData$cell_type <- barcode_cellName
    myData$cell_id <- barcode_cellid
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
    pdf(paste0(output, "/", libraryID, "/UMAP_seurat_cluster.pdf"))
    print(umapPlot)
    dev.off()
    umapPlot <- DimPlot(myData, reduction = "umap", group.by='cell_type')
    ## save UMAP information with cell_type
    pdf(paste0(output, "/", libraryID, "/UMAP_cellType_cluster.pdf"))
    print(umapPlot)
    dev.off()
    umapPlot <- DimPlot(myData, reduction = "umap", group.by='cell_id')
    ## save UMAP information with cell_ontology
    pdf(paste0(output, "/", libraryID, "/UMAP_cellID_cluster.pdf"))
    print(umapPlot)
    dev.off()
    
    ## collect raw UMI counts and normalized data counts for each cell
    finalRaw <- data.frame(myData@assays$RNA@counts)
    finalCPM <- data.frame(myData@assays$RNA@data)
    
    ## write in the output the info per cell type that is present in the library
    infoCollected <- c()
    for (cell in unique(myData$cell_type)) {
      for (cellid in unique(myData$cell_id[myData$cell_type == cell])) {
        ## split information
        barcodesID <- colnames(myData)[myData$cell_type == cell & myData$cell_id == cellid] 
        ## export raw counts to each cell type
        rawCountsCell <- finalRaw[(names(finalRaw) %in% barcodesID)]
        rawCountsCell <- setDT(rawCountsCell, keep.rownames = TRUE)[]
        colnames(rawCountsCell)[1] <- "gene_id"
        ## add biotype info to raw counts
        collectBiotypeRaw <- merge(rawCountsCell, biotypeInfo, by = "gene_id", all.x = TRUE)
        ## add type info
        collectBiotypeRaw$type <- ifelse(is.na(collectBiotypeRaw$biotype), "intergenic", "genic")
        collectBiotypeRaw$cellTypeName <- cell
        collectBiotypeRaw$cellTypeId <- cellid
        
        ## export normalized counts to each cell type
        normalizedCountsCell <- finalCPM[(names(finalCPM) %in% barcodesID)]
        normalizedCountsCell <- setDT(normalizedCountsCell, keep.rownames = TRUE)[]
        colnames(normalizedCountsCell)[1] <- "gene_id"
        ## add biotype info to normalized counts
        collectBiotypeNorm <- merge(normalizedCountsCell, biotypeInfo, by = "gene_id", all.x = TRUE)
        ## add type info
        collectBiotypeNorm$type <- ifelse(is.na(collectBiotypeNorm$biotype), "intergenic", "genic")
        collectBiotypeNorm$cellTypeName <- cell
        collectBiotypeNorm$cellTypeId <- cellid
        cellid <- gsub(":","-",cellid)
        ## write output information to integrate in Bgee
        write.table(collectBiotypeRaw, file = paste0(output, "/", libraryID, "/Raw_Counts_", cell, "_", cellid, ".tsv"), sep="\t", row.names = FALSE, quote = FALSE)
        write.table(collectBiotypeNorm, file = paste0(output, "/", libraryID, "/Normalized_Counts_", cell, "_", cellid,  ".tsv"), sep="\t", row.names = FALSE, quote = FALSE)
        ## info per cell (-5 because of geneId, biotype, type columns and cellName and cellID)
        cellInfo <- c(nrow(collectBiotypeNorm), ncol(collectBiotypeNorm)-5, cell, cellid)
        infoCollected <- rbind(infoCollected, cellInfo)
      }
    }
    infoCollected <- as.data.frame(infoCollected)
    colnames(infoCollected) <- c("Number_genes", "Number_cells", "Cell_Name", "Cell_Ontology")
    return(list(myData[[]],infoCollected))
  } else {
    message("Library", libraryID, "not take in consideration for posterior analysis.")
    message("Barcode not provided!")
  }
}

#############################################################################################################################
## collect information for all libraries
globalInfoLibraries <- paste0(output, "/InformationAllLibraries.txt")
if (file.exists(globalInfoLibraries)){
  message("File already exists and will be removed to create a new one to avoid overwritting!")
  file.remove(globalInfoLibraries)
  file.create(globalInfoLibraries)
  cat("library\texperimentID\tInitial_UMI_barcode\tInitial_tot_genes\tgenes_afterFiltering\tcells_afterFiltering\tcellTypeName\tcellTypeId\n",file = globalInfoLibraries, sep = "\t")
} else {
  file.create(globalInfoLibraries)
  cat("library\texperimentID\tInitial_UMI_barcode\tInitial_tot_genes\tgenes_afterFiltering\tcells_afterFiltering\tcellTypeName\tcellTypeId\n",file = globalInfoLibraries, sep = "\t")
}

## For each library that belongs to each experiment do:
for (libraryID in scRNASeqAnnotation$libraryId) {
  
  ## verify if library exist
  if (file.exists(file.path(kallisto_bus_results, libraryID))){

    message("Treating library ", libraryID)
    
    ## info about species and experiment
    speciesID <- unique(scRNASeqAnnotation$scientific_name[scRNASeqAnnotation$libraryId == libraryID])
    speciesID <- gsub(" ", "_", speciesID)
    experimentID <- unique(scRNASeqAnnotation$experimentId[scRNASeqAnnotation$libraryId == libraryID])

    path2Files <- file.path(kallisto_bus_results, libraryID, "gene_counts/")
    
    ## create folder per library
    if(!dir.exists(file.path(output,libraryID))) {
      dir.create(file.path(output, libraryID))
    }
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
    biotypeInfo <- fread(file.path(folderSupport, paste0(speciesID, "_gene_to_biotype_with_intergenic.tsv")))
    colnames(biotypeInfo) <- c("gene_id", "biotype")
    ## get barcode information per experiment
    barcodeExperiment <- fread(file.path(infoFolder, paste0("scRNASeq_barcode_", experimentID,".tsv")))
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
      allinfo <- data.frame(libraryID, experimentID, length(tot_counts) , length(tot_genes), infoCells$Number_genes, infoCells$Number_cells, infoCells$Cell_Name, infoCells$Cell_Ontology)
      write.table(allinfo, file = globalInfoLibraries, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
    }
  } else {
      warning("This library ", libraryID, " not exist in the directory.")
  }
}

