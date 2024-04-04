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
command_arg <- c("libraryId", "scRNASeq_Info", "kallisto_bus_results", "folderSupport", "infoFolder", "output")
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
## We decided not to filter data in our pipeline but rather to trust author
## cluster/celltype analysis. Then, we do not anymore filter cells based on the kneeplot.
## We still generate the kneeplot and provide number of cells that would have been
## filtered if we were using this QC.
singlecellKnee <- function(sparseMatrix, libraryID, tot_counts){
  ## rank barcodes to do the knee
  bc_rank <- NULL
  tryCatch(
    {
      bc_rank <- barcodeRanks(sparseMatrix)
    },
    error=function(e){
      if (grepl(pattern = "insufficient unique points for computing knee/inflection points", x = e)) {
        warning("Not enough reads mapped to the transcriptome. This low number of reads did not allow to to compute the knee/inflection points.",
		" This library will not be processed and will not be inserted in the database.")
      }
      stop("Not able to rank barcodes. The error message was : ", e)
    }
  )
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
  
  pdf(file = file.path(output, libraryID, "kneePlot.pdf"), width = 16, height = 10)
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
targetCells <- function(objectNormalized, objectNormalized_filtered, barcodeInfo, biotypeInfo, libraryID){
  
  ## select cells based on the barcode ID's
  if (nrow(barcodeInfo) != 0){
    message("library ID ", unique(barcodeInfo$library), " has ", nrow(barcodeInfo), " barcodes with annotation")
    barcode_ids_having_celltype_annot <- barcodeInfo$barcode
    ## subset barcodes for which celltype are provided
    seurat_object_with_celltype <- subset(objectNormalized,  cells = barcode_ids_having_celltype_annot)
    seurat_object_with_celltype_filtered <- subset(objectNormalized_filtered,  cells = barcode_ids_having_celltype_annot)
    message(ncol(seurat_object_with_celltype) - ncol(seurat_object_with_celltype_filtered), " cells out of ",
        ncol(seurat_object_with_celltype), " with celltype annotation would have been removed by filtering ",
        "on the inflection point of the kneeplot.")
    seurat_object_with_celltype$internal_cluster_id <- barcodeInfo$internal_cluster_id
    seurat_object_with_celltype$celltype_id <- barcodeInfo$cellTypeId
    seurat_object_with_celltype$celltype_name <- barcodeInfo$cellTypeName
    seurat_object_with_celltype$celltype_author_annotation <- barcodeInfo$cell_type
    #head(myData[[]])
    ## remove unsigned cells (this avoid the exportation of the file of unassigned cells)
    seurat_object_with_celltype <- subset(seurat_object_with_celltype,  celltype_name != "Unassigned")
    
    set.seed(42)

    seurat_object_with_celltype <- FindVariableFeatures(seurat_object_with_celltype, selection.method = "vst", nfeatures = 2000)
    ## Identify the 20 most highly variable genes
    message("20 most highly variable genes : ")
    head(VariableFeatures(seurat_object_with_celltype), 20)
    ## scale the data
    seurat_object_with_celltype <- ScaleData(seurat_object_with_celltype, features = rownames(seurat_object_with_celltype))
    ## Run PCA
    ## catch error specific to the version of Seurat we use when less than 50 cells per library. Updating Seurat should remove this error.
    ## For now we continue processing libraries even if it contains less than 50 cells but throw a warning and do not generate the PCA and UMAP
    tryCatch(
      {
        seurat_object_with_celltype <- RunPCA(seurat_object_with_celltype, features = VariableFeatures(object = seurat_object_with_celltype))
        seurat_object_with_celltype <- FindNeighbors(seurat_object_with_celltype, dims = 1:20)
        seurat_object_with_celltype <- FindClusters(seurat_object_with_celltype, resolution = 1)
        ## Run UMAP
        seurat_object_with_celltype <- RunUMAP(seurat_object_with_celltype, dims = 1:20)
        umapPlot <- DimPlot(seurat_object_with_celltype, reduction = "umap", group.by='seurat_clusters')
        pdf(file.path(output, libraryID, "UMAP_seurat_cluster.pdf"))
        print(umapPlot)
        dev.off()
        # use the internal cluster ID to cluster as it is not possible to use the author free text celltype
        # information for some dataset (e.g FCA)
        umapPlot <- DimPlot(seurat_object_with_celltype, reduction = "umap", group.by='internal_cluster_id')
        ## save UMAP information with cell_type
        pdf(file.path(output, libraryID, "UMAP_cellType_author_annotation_cluster.pdf"))
        print(umapPlot)
        dev.off()
        umapPlot <- DimPlot(seurat_object_with_celltype, reduction = "umap", group.by='celltype_id')
        ## save UMAP information with cell_ontology
        pdf(file.path(output, libraryID, "UMAP_cellID_cluster.pdf"))
        print(umapPlot)
        dev.off()
      },
      error=function(e){
        warning("an error occured while generating the PCA or the UMAP. No UMAP file will be created")
        if (grepl(pattern = "max\\(nu, nv\\) must be strictly less than min\\(nrow\\(A\\), ncol\\(A\\)\\)", x = e)) {
          message("Seurat can not generate a PCA for less than 50 cells")
	}
	message("the error message is : ", e)
      }
    )
    
    ## collect raw UMI counts and normalized data counts for each cell
    finalRaw <- data.frame(seurat_object_with_celltype@assays$RNA@counts)
    finalCPM <- data.frame(seurat_object_with_celltype@assays$RNA@data)
    
    ## write in the output the info per internal cluster ID
    ##
    infoCollected <- c()
    for (internalClusterId in unique(seurat_object_with_celltype$internal_cluster_id)) {
      cellId = unique(as.character(seurat_object_with_celltype$celltype_id[seurat_object_with_celltype$internal_cluster_id == internalClusterId]))
      cellName = unique(as.character(seurat_object_with_celltype$celltype_name[seurat_object_with_celltype$internal_cluster_id == internalClusterId]))
      celltypeAuthorAnnotation = unique(as.character(seurat_object_with_celltype$celltype_author_annotation[seurat_object_with_celltype$internal_cluster_id == internalClusterId]))
      if (is.na(cellId)) {
        stop("No cell-type ID for library ", libraryID, " and internal cluster ID ", internalClusterId)
      }
      if(length(cellId) > 1){
        stop("More than one cell-type ID for library ", libraryID, ", internal cluster ID ", internalClusterId, " celltypeIds : ",cellId)
      }
      message("cluster info[internalClusterId: ", internalClusterId, ", cellId: ", cellId, ", celltypeAuthorAnnotation= ", celltypeAuthorAnnotation, "]")
      ## split information
      barcodesID <- colnames(seurat_object_with_celltype)[seurat_object_with_celltype$internal_cluster_id == internalClusterId] 
      ## export raw counts to each cell type
      rawCountsCell <- finalRaw[(names(finalRaw) %in% barcodesID)]
      rawCountsCell <- cbind(names = rownames(rawCountsCell), rawCountsCell)
      colnames(rawCountsCell)[1] <- "gene_id"
      ## add biotype info to raw counts
      collectBiotypeRaw <- merge(rawCountsCell, biotypeInfo, by = "gene_id", all.x = TRUE)
      collectBiotypeRaw$internalClusterId <- internalClusterId
      collectBiotypeRaw$cellTypeId <- cellId
        
      ## export normalized counts to each internal cluster ID
      normalizedCountsCell <- finalCPM[(names(finalCPM) %in% barcodesID)]
      normalizedCountsCell <- cbind(names = rownames(normalizedCountsCell), normalizedCountsCell)
      colnames(normalizedCountsCell)[1] <- "gene_id"
      ## add biotype info to normalized counts
      collectBiotypeNorm <- merge(normalizedCountsCell, biotypeInfo, by = "gene_id", all.x = TRUE)
      ## add type info
      collectBiotypeNorm$type <- ifelse(is.na(collectBiotypeNorm$biotype), "intergenic", "genic")
      collectBiotypeNorm$cellTypeName <- cellName
      ## write output information to integrate in Bgee
      rawCountFilePath <- file.path(output, libraryID, paste0("Raw_Counts_", internalClusterId,
          ".tsv"))
      normalizedCountFilePath <- file.path(output, libraryID, paste0("Normalized_Counts_",
        internalClusterId, ".tsv"))
      write.table(collectBiotypeRaw, file = rawCountFilePath, sep="\t", row.names = FALSE, 
        quote = FALSE)
      write.table(collectBiotypeNorm, file = normalizedCountFilePath, sep="\t", row.names = FALSE,
        quote = FALSE)
 
      ## info per cell (-5 because of geneId, biotype, type columns and cellName and cellID)
      cellInfo <- c(nrow(collectBiotypeNorm), ncol(collectBiotypeNorm)-5, celltypeAuthorAnnotation, cellId)
      infoCollected <- rbind(infoCollected, cellInfo)
    }
    infoCollected <- as.data.frame(infoCollected)
    colnames(infoCollected) <- c("Number_genes", "Number_cells", "Celltype_author_annotation", "Celltype_id")

    return(list(seurat_object_with_celltype[[]],infoCollected))
  } else {
    message("Library", libraryID, "not take in consideration for posterior analysis.")
    message("Barcode not provided!")
  }
}

#############################################################################################################################
## collect information for all libraries
if (!dir.exists(output)) {
  dir.create(output)
}

## verify if library exist
if (file.exists(file.path(kallisto_bus_results, libraryId, "gene_counts"))){

  message("Treating library ", libraryId)
  
  ## info about species and experiment
  speciesName <- unique(scRNASeqAnnotation$scientific_name[scRNASeqAnnotation$library_id == libraryId])
  speciesName <- gsub(" ", "_", speciesName)
  experimentID <- unique(scRNASeqAnnotation$experiment_id[scRNASeqAnnotation$library_id == libraryId])

  path2Files <- file.path(kallisto_bus_results, libraryId, "gene_counts/")

  if(file.exists(file.path(output, libraryId, "done"))) {
    message("celltype identification already done for library : ", libraryId)
  } else {
  
    ## create folder per library
    if(!dir.exists(file.path(output, libraryId))) {
      dir.create(file.path(output, libraryId))
    }

    ## Create a sparseMatrix
    sparseMatrix <- read_count_output(path2Files, "gene", tcc = FALSE)
    message("created sparse matrix")
    ## export global information
    ## barcodes detected (cell per column) and genes detected (rows)
    ## How many UMIs per barcode (cell)
    tot_counts <- Matrix::colSums(sparseMatrix)
    ## How many genes detected
    tot_genes <- rowSums(sparseMatrix)
    ## filter cells based on inflection point of the kneeplot.
    ## this is not anymore used to filter the cells but just to generate a kneeplot
    ## and provide filtering info in case we want to check proportion of cells with cell-type annotation
    ## that would have been removed
    sparseMatrix_filtered <- singlecellKnee(sparseMatrix = sparseMatrix, libraryID = libraryId, tot_counts = tot_counts)
    ## perform the seurat object and get raw and normalized counts
    object_seurat <- seuratObject(m_filtered = sparseMatrix)
    object_seurat_filtered <- seuratObject(m_filtered = sparseMatrix_filtered)
    ## read biotype information
    biotypeInfoFile <- list.files(path = folderSupport, pattern = paste0("^", speciesName, ".*.gene2biotype$"),
      full.names = TRUE)
    message("biotype info file : ", biotypeInfoFile)
    if (length(biotypeInfoFile) != 1) {
      stop("No gene2biotype file found for the species ", speciesName)
    }
    biotypeInfo <- read.table(biotypeInfoFile)
    colnames(biotypeInfo) <- c("gene_id", "biotype", "type")
    ## get barcode information per experiment
    barcodeExperiment <- read.table(file.path(infoFolder,
      paste0("scRNASeq_barcode_", experimentID,".tsv")), header = TRUE,
      sep = "\t", quote = "\"")
    message("barcode experiment file has ", nrow(barcodeExperiment), " rows")
    barcodeLibrary <- dplyr::filter(barcodeExperiment, library == libraryId)
    ## link barcode to cell ID and export information
    finalInformationCells <- targetCells(objectNormalized = object_seurat, objectNormalized_filtered = object_seurat_filtered,
      barcodeInfo = barcodeLibrary, biotypeInfo = biotypeInfo, libraryID = libraryId)
    infoCells <- finalInformationCells[[2]]
    ## verify if function targetCells is NULL
    if (is.null(finalInformationCells)){
      message("The library", libraryId, " does not contain any barcode.")
      message("Any raw or normalized count is exported.")
      message("The library will not be present in the informative files: ",
        "InformationAllLibraries.txt, variabilityClusters.txt and markerGenes_Validated.txt file")
    } else {
      ## write all information for the library
      info <- data.frame(libraryId, experimentID, length(tot_counts) , length(tot_genes),
        infoCells$Number_genes, infoCells$Number_cells, as.character(infoCells$Celltype_author_annotation), as.character(infoCells$Celltype_id))
      colnames(info) <- c("library", "experimentID", "Initial_UMI_barcode",
        "Initial_tot_genes", "genes_afterFiltering", "cells_afterFiltering",
        "cellTypeAuthorAnnotation", "cellTypeId")

      # write file containing knee plot filtering stats
      write.table(info, file = file.path(output, libraryId, "kneePlotFilteringInfo.txt"),
        quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
    }
    # write an empty file to detect already processed libraries
    file.create(file.path(output, libraryId, "done"))
  }
} else {
    stop("The library ", libraryId, " has not been porcessed by bustools")
}

