## SFC, 10 Oct 2019
## Updated June 2020

## This script is used to do:
## the initial quality control (filtering) based on the Knee plot,
## to map cell types using the information from barcode annotation file,
## to provide a variability info file (this means info about the proportion of cell-types per cluster)
## to provide a list of gene markers per cell-type after data analysis (validated with gene markers from the annotation file)

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNA_Seq_info_target.txt" folder_data="folder_data" kallisto_bus_results="kallisto_bus_results" folderSupport="folderSupport" infoFolder="infoFolder" output="output"' QC_CelltypeIdentification.R QC_CelltypeIdentification.Rout
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
command_arg <- c("scRNASeq_Info","folder_data", "kallisto_bus_results", "folderSupport", "infoFolder", "output")
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
    
    ## find the marker genes between clusters
    markersLibrary <- FindAllMarkers(myData, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
    
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
    return(list(myData[[]],markersLibrary,infoCollected))
  } else {
    message("Library", libraryID, "not take in consideration for posterior analysis.")
    message("Barcode not provided!")
  }
}

## Function to export information of % cell-types per cluster and to provide 
## gene markers per cell-type (using markers found by the analysis and validated using the annotation file)
qualityControl <- function(finalInformationCells, geneMarkerLibrary, geneNameFile, speciesID, ensembl_Uniprot_file){
  
  #### 1) Variability in clusters
  
  myData <- finalInformationCells[[1]]
  ## collect variability of different cells per cluster, as well as, the correspondent %proportion of each cell type in the library
  variability_cluster <- c()
  global_proportion <- c()
  
  for (cluster in unique(myData$seurat_clusters)) {
    
    barcodesCluster <- myData[which(myData$seurat_clusters == cluster),]
    infoCells <- as.data.frame(table(barcodesCluster$cell_type, barcodesCluster$cell_id))
    infoCells <- dplyr::filter(infoCells, Freq != 0)
    
    for (i in seq(nrow(infoCells))) {
      
      cellName <- infoCells$Var1[i]
      cellId <- infoCells$Var2[i]
      numberCells <- infoCells$Freq[i]
      
      calculateVariabilityCluster <- numberCells/sum(infoCells$Freq)
      infoVariabilityCluster <- (t(c(cluster, calculateVariabilityCluster, as.character(cellName), as.character(cellId))))
      variability_cluster <- rbind(variability_cluster, infoVariabilityCluster)
      
      calculateGlobalProportion <- (numberCells/nrow(myData))*100
      infoGlobalVariability <- (t(c(cluster, calculateGlobalProportion, as.character(cellName), as.character(cellId))))
      global_proportion <- rbind(global_proportion, infoGlobalVariability)
    }
  }
  ## general Info about variability in the cluster (NOTE: this is also influenced by the parameters used in the targetCells function)
  ## so this works as informative file in case we need to go back to the annotation part
  variabilityFinalInfo <- data.frame(variability_cluster, global_proportion[,2])
  colnames(variabilityFinalInfo) <- c("cluster", "Variability_cluster", "Cell_Name", "cellTypeId" , "Global_Proportion_library")
  variabilityFinalInfo$Classification_Variability_cluster <- " "
  
  ## add a warning in the file in case some cluster have one cell-type with less 80% of occupancy in the cluster (this means high variability of cells in the same cluster)
  for (i in unique(variabilityFinalInfo$cluster)) {
    ## select rows in data.frame referent to the cluster
    rowID <- which(variabilityFinalInfo$cluster == i)
    checkCluster <- as.numeric(as.character(variabilityFinalInfo$Variability_cluster[variabilityFinalInfo$cluster == i]))
    check <- any(checkCluster >= 0.80)
    if (check == "TRUE"){
      variabilityFinalInfo$Classification_Variability_cluster[rowID]<- " "
    } else {
      variabilityFinalInfo$Classification_Variability_cluster[rowID]<- "high"
    }
  }
  
  #### 2) Gene markers validation
  
  ## gene markers found by the Bgee analysis
  markersLibrary <- finalInformationCells[[2]]
  ## add cell-type name detected in the cluster
  markersLibrary$cell_names_cluster <- " "
  for (i in unique(markersLibrary$cluster)) {
    ## select rows in data.frame referent to the cluster
    rowID <- which(markersLibrary$cluster == i)
    cellType <- unique(myData$cell_type[myData$seurat_clusters == i])
    cellType <- paste(cellType,collapse=",")
    markersLibrary$cell_names_cluster[rowID]<- cellType
  }
  ## provide gene name to gene IDs
  markersLibrary <- merge(markersLibrary, geneNameFile, by = "gene")
  colnames(markersLibrary)[1] <- "markerGene_ID_UniProt_Ensembl"
  
  ## compare and validate using the gene markers from the annotation
  if (nrow(geneMarkerLibrary) == 0){
    message("The experiment not provide info about gene markers in the annotation file.")
  } else {
    ## select from annotation
    subset_annotation_markers <- as.data.table(geneMarkerLibrary %>% dplyr::select(experimentId, uberonId, uberonName, cell_type, cell_type_harmonization, cellTypeId, cellTypeName, markerGene_ID_UniProt_Ensembl))
    ## because we subset, some rows can be duplicated (Explanation: what makes the row unique is the source column in the annotation, where we report from where the info comes, example cluster figure).
    subset_annotation_markers <- subset_annotation_markers[!duplicated(subset_annotation_markers), ]
    
    ## see if we have info from UNIPROT or ensembl. 
    isUniprot <- subset_annotation_markers[grepl('ENS', subset_annotation_markers$markerGene_ID_UniProt_Ensembl)]
    
    if (nrow(isUniprot) == 0){
      ## If info cames from UNIPROT use informative file for the referent species to match geneID_TO_UniprotID
      annotLookup <- ensembl_Uniprot_file
      
      subset_annotation_markers <- merge(subset_annotation_markers, annotLookup, by ="markerGene_ID_UniProt_Ensembl")
      colnames(subset_annotation_markers)[1] <- "UNIPROT_ID"
      colnames(subset_annotation_markers)[9] <- "markerGene_ID_UniProt_Ensembl"
      
      ## verify if this markers genes from the analysis exist in annotation file
      finalMarker <- merge(subset_annotation_markers, markersLibrary, by = "markerGene_ID_UniProt_Ensembl")
      ## just keep marker genes if p_val_adj is statistically significant
      geneMarkersValidatedFromBgee <- finalMarker[finalMarker$p_val_adj <= 0.05, ]
      
      ## verify if marker match between the cell-type of the annotation and one of the cell types identified in the analysis-cluster
      ## compare cell_type_hamonization column (annotation) and cell_names_cluster (analysis)
      geneMarkersValidatedFromBgee$cell_names_cluster <- strsplit(geneMarkersValidatedFromBgee$cell_names_cluster, ",")
      for (i in seq(length(geneMarkersValidatedFromBgee$cell_names_cluster))) {
        harmonization <- geneMarkersValidatedFromBgee$cell_type_harmonization[i]
        cluster <- geneMarkersValidatedFromBgee$cell_names_cluster[[i]]
        geneMarkersValidatedFromBgee$match[i] <- (harmonization %in% cluster)
      }
      ## keep everything that is in agreement between annotation and analysis
      geneMarkersValidatedFromBgee <- geneMarkersValidatedFromBgee[geneMarkersValidatedFromBgee$match == TRUE, ]
      if(nrow(geneMarkersValidatedFromBgee) == 0){
        namesCol <- names(geneMarkersValidatedFromBgee)
        ## if none marker is found in annotation this is "-" row
        geneMarkersValidatedFromBgee <- as.data.frame(t(rep(c("-"),times=ncol(geneMarkersValidatedFromBgee))))
        colnames(geneMarkersValidatedFromBgee) <- namesCol
      } else {
        geneMarkersValidatedFromBgee <- geneMarkersValidatedFromBgee
      }
      
    } else {
      ## if annotation provide ensembl id we don't need to retrieve from UNIPROT (do directly)
      finalMarker <- merge(subset_annotation_markers, markersLibrary, by = "markerGene_ID_UniProt_Ensembl")
      geneMarkersValidatedFromBgee <- finalMarker[finalMarker$p_val_adj <= 0.05, ]
      ## compare cell_type_hamonization column (annotation) and cell_names_cluster (analysis)
      geneMarkersValidatedFromBgee$cell_names_cluster <- strsplit(geneMarkersValidatedFromBgee$cell_names_cluster, ",")
      for (i in seq(length(geneMarkersValidatedFromBgee$cell_names_cluster))) {
        harmonization <- geneMarkersValidatedFromBgee$cell_type_harmonization[i]
        cluster <- geneMarkersValidatedFromBgee$cell_names_cluster[[i]]
        geneMarkersValidatedFromBgee$match[i] <- (harmonization %in% cluster)
      }
      ## keep everything that is agreement between annotation and analysis
      geneMarkersValidatedFromBgee <- geneMarkersValidatedFromBgee[geneMarkersValidatedFromBgee$match == TRUE, ]
      ## add empty column just to match the output file (since we have ensembl ID we don't have info for UNIPROT)
      geneMarkersValidatedFromBgee <- add_column(geneMarkersValidatedFromBgee, UNIPROT_ID = "-", .after = "markerGene_ID_UniProt_Ensembl")
    }
  }
  return(list(variabilityFinalInfo, geneMarkersValidatedFromBgee))
}

#############################################################################################################################
## collect information for all libraries
##XXX what happens if file already exist? add new lines to already existing one? Should maybe empty it and create a new one or at least stop the script
globalInfoLibraries <- paste0(output, "/InformationAllLibraries.txt")
if (!file.exists(globalInfoLibraries)){
  file.create(globalInfoLibraries)
  cat("library\texperimentID\tInitial_UMI_barcode\tInitial_tot_genes\tgenes_afterFiltering\tcells_afterFiltering\tcellTypeName\tcellTypeId\n",file = paste0(output, "/InformationAllLibraries.txt"), sep = "\t")
} else {
  print("File already exist.....")
}
## variability per cluster
variabilityCluster <- paste0(output, "/variabilityClusters.txt")
if (!file.exists(variabilityCluster)){
  file.create(variabilityCluster)
  cat("libraryID\texperimentID\tcluster\tVariability_cluster\tcellTypeName\tcellTypeId\tGlobal_Proportion_library\tClassification_Variability_cluster\n",file = paste0(output, "/variabilityClusters.txt"), sep = "\t")
} else {
  print("File already exist.....")
}
## marker genes validated
markerGenes <- paste0(output, "/markerGenes_Validated.txt")
if (!file.exists(markerGenes)){
  file.create(markerGenes)
  cat("libraryID\tmarkerGene_ID_UniProt_Ensembl\tUNIPROT_ID\texperimentId\tuberonId\tuberonName\tcell_type\tcell_type_harmonization\tcellTypeId\tcellTypeName\tp_val\tavg_logFC\tpct.1\tpct.2\tp_val_adj\tcluster\tcell_names_cluster\tgene_name\tmatch\n",file = paste0(output, "/markerGenes_Validated.txt"), sep = "\t")
} else {
  print("File already exist.....")
}

## For each library that belongs to each experiment do:
for (libraryID in scRNASeqAnnotation$libraryId) {
  
  ## verify if library exist
  if (file.exists(file.path(folder_data, libraryID))){

    message("Treating library ", libraryID)
    
    ## info about species and experiment
    speciesID <- scRNASeqAnnotation$scientific_name[scRNASeqAnnotation$libraryId == libraryID]
    speciesID <- gsub(" ", "_", speciesID)
    experimentID <- scRNASeqAnnotation$experimentId[scRNASeqAnnotation$libraryId == libraryID]

    path2Files <- file.path(kallisto_bus_results, libraryID, "gene_counts")

    path2ListFiles <- list.files(path2Files)
    
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
    biotypeInfo <- fread(paste0(folderSupport, "/gene_to_biotype_with_intergenic_", speciesID,".tsv"))
    colnames(biotypeInfo) <- c("gene_id", "biotype")
    ## get barcode information for the library
    barcodeLibrary <- dplyr::filter(barcodes, library == libraryID)
    ## get gene markers information
    geneMarkerLibrary <- dplyr::filter(markers, experimentId == experimentID)
    ## read info that link gene_id with gene_name
    geneNameFile <- read.table(paste0(folderSupport, "/gene_to_geneName_with_intergenic_", speciesID,".tsv"))
    colnames(geneNameFile) <- c("gene", "gene_name")
    
    ## link barcode to cell ID and export information
    finalInformationCells <- targetCells(objectNormalized = object, barcodeIDs = barcodeLibrary, biotypeInfo = biotypeInfo, libraryID = libraryID)
    infoCells <- finalInformationCells[[3]]
    
    ## verify if function targetCells is NULL
    if (is.null(finalInformationCells)){
      message("The library", libraryID, " not contain barcodes.")
      message("Any raw or normalized count is exported.")
      message("The library will not be present in the informative files: InformationAllLibraries.txt, variabilityClusters.txt and markerGenes_Validated.txt file")
    } else {
      ## write all information for the library
      allinfo <- data.frame(libraryID, experimentID, length(tot_counts) , length(tot_genes), infoCells$Number_genes, infoCells$Number_cells, infoCells$Cell_Name, infoCells$Cell_Ontology)
      write.table(allinfo, file = globalInfoLibraries, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
      
      ## collect file that link geneID_TO_Uniprot information for the species
      ensembl_Uniprot <- read.table(paste0(infoFolder, "/ensembl_Uniprot_", speciesID,".tsv"), header = TRUE, sep="\t")
      ## export information about variability of cells per cluster and export information about validated gene markers after data analysis
      infoQC <- qualityControl(finalInformationCells = finalInformationCells, geneMarkerLibrary = geneMarkerLibrary, geneNameFile = geneNameFile, speciesID = speciesID, ensembl_Uniprot_file = ensembl_Uniprot)
      variabilityClusterInfo <- data.frame(libraryID, experimentID, infoQC[[1]])
      markersInfoValidated <- data.frame(libraryID, infoQC[[2]])
      
      write.table(variabilityClusterInfo, file = variabilityCluster, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
      fwrite(markersInfoValidated, file = markerGenes, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
    }
  } else {
      warning("This library ", libraryID, " not exist in the directory.")
  }
}

