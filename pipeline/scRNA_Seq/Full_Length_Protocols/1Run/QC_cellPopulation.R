## SFonsecaCosta, Sept 2018

### This script is used to do the quality control of the data.
### In this script we work with abundance_gene_level+fpkm+intergenic.tsv output file from the previous step (kallisto_analysis.R)

## 1 Step --> Create a matrix with all cells that belong to the same cell-type, same experiment and same species.
## 2 Step --> Verify if each cell-type per experiment and species follow a bimodal distribution
## 3 Step --> generate output file with experiments that pass the requirement

## Usage:
## R CMD BATCH --no-save --no-restore '--args scrna_seq_sample_info="scrna_seq_sample_info.tsv" cells_folder="cells_folder" output_folder="output_folder" plot="no"' QC_cellPopulation.R QC_cellPopulation.Rout
## scrna_seq_sample_info --> collect information about each cell/experiment/species
## cells_folder --> where we is located the data (raw and treated data)
## output_folder --> where we save the output's (files and plots)
## plot --> if should plot the graphic results of each experiment (by default is no)

## Libraries Used
library(ggplot2)
library(gridExtra)
library(gghighlight)
library(LaplacesDemon)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}
## checking if all necessary arguments were passed....
command_arg <- c("scrna_seq_sample_info", "cells_folder", "output_folder", "plot")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}
## Read scrna_seq_sample_info file. If file not exists, script stops
if( file.exists(scrna_seq_sample_info) ){
  annotation <- read.table(scrna_seq_sample_info, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
} else {
  stop( paste("scrna_seq_sample_info file not found [", scrna_seq_sample_info, "]\n"))
}

## Create output files
if (dir.exists(output_folder)){
 
   file.create("Modality_Cell_type_per_experiment.tsv")
  cat("experiment\tcellId\tspecies\tmodeDistribution\n",file = file.path(output_folder, "Modality_Cell_type_per_experiment.tsv"), sep = "\t")
  file.create("NEW_scRNASeq_sample_info.tsv")
  cat("libraryId\texperimentId\tcellTypeName\tcellTypeId\tspeciesId\tplatform\tlibraryType\tprotocol\tprotocolType\treadLength\tinfoOrgan\torganism\n",file = file.path(output_folder,"NEW_scRNASeq_sample_info.tsv"), sep = "\t")	
  file.create("Discard_scRNASeq_sample_info.tsv")
  cat("libraryId\texperimentId\tcellTypeName\tcellTypeId\tspeciesId\tplatform\tlibraryType\tprotocol\tprotocolType\treadLength\tinfoOrgan\torganism\n",file = file.path(output_folder, "Discard_scRNASeq_sample_info.tsv"), sep = "\t")	
  
} else {
  print("Directoty not exists.....")
}
ModalFile <- file.path(output_folder, "Modality_Cell_type_per_experiment.tsv")


## Check bimodality + plot distribution of the protein coding for each cell-type
checkDataPlot <- function(matrixData, output, plot){
  
  ### Just protein coding genes
  proteinCoding <- dplyr::filter(matrixData, biotype == "protein_coding")
  ## Count TPM different zero for each cell
  cellsTPM <- as.data.frame(colSums(proteinCoding[,4:length(proteinCoding)] != 0))
  colnames(cellsTPM) <- "Sum_cells"
  cellsTPM$seq <- seq(1, nrow(cellsTPM), by=1)
  
  ### Verify in how many cells the gene have a TPM higher then zero
  codingGenes <- as.data.frame(apply(proteinCoding[,4:length(proteinCoding)],1,function(x)sum(x != 0)))
  colnames(codingGenes) <- "genes"
  codingGenes$seq <- seq(1, nrow(codingGenes), by=1)
  
  ## Create final table
  finalTable <- data.frame(proteinCoding[,1], codingGenes)
  colnames(finalTable) <- c("gene", "number_of_cells_geneDetected", "seq")
  finalTable <- as.data.frame(finalTable[order(finalTable$number_of_cells_geneDetected),])
  finalTable$seq <- seq(1, nrow(finalTable), by=1)
  
  number_cells <- length(matrixData[,4:length(matrixData)])
  cell_25 <- number_cells*0.25
  cell_100 <- number_cells*1
  all_cells <- sum(finalTable$number_of_cells_geneDetected >= cell_100)
  
  ##### PLOTING #####
  if (plot == "yes"){
    p1 <- ggplot(cellsTPM, aes(x=seq,y=Sum_cells)) +
      geom_point() + 
      labs(title = "Protein Coding (Normalized values - TPM)", x="Number of cells",y="Number of genes")
    p2 <- ggplot(codingGenes, aes(x = genes)) + 
      geom_density(alpha=1, color="#E69F00", fill="#E69F00") + 
      labs(title = "Distribution protein coding", x="Number of cells",y="Density")
    p3 <- ggplot(finalTable) +
      geom_point(aes(x=seq,y=number_of_cells_geneDetected, color=ifelse(finalTable$number_of_cells_geneDetected < cell_25, "black", ifelse(finalTable$number_of_cells_geneDetected == cell_100,"purple", "gray")))) + 
      labs(title = "Proportion of genes expressed over cells", x="Number of genes",y="Number of cells") +
      scale_color_manual(name="Proportion of cells",values=c("black", "gray","purple"),labels=c("=< 25%","25% < gene < 100%",paste0("all cells (", all_cells,")"))) +
      theme(legend.position=c(0.85, 0.1))
    
    pdf(file = file.path(output_folder, paste0("Plot","_", cellId, "_",experiment ,".pdf")), width = 12, height = 12)   
    grid.arrange(p1,p2,p3,ncol = 1, nrow = 3)
    dev.off()
   } else {
    cat("Analysis done without reporting graphics.", "\n")
  }
  
  ## test modality of the data
  if ( is.unimodal(codingGenes$genes) == TRUE){
    cat("The cell-type ",cellId," from the experiment:", experiment,"is unimodal for the protein coding genes!", "\n", "Warning: this experiment can be removed...", "\n")
    modeIs <- paste0("unimodal")
  } else if (is.bimodal(codingGenes$genes) == FALSE){
    multimodal <- is.multimodal(codingGenes$genes) 
    trimodal <- is.trimodal(codingGenes$genes)
    cat("The cell-type ",cellId," from the experiment:", experiment,"is trimodal", trimodal, "or multimodal:", multimodal, "\n")
    modeIs <- paste0("multimodal")
  } else {
    cat("The cell-type ",cellId," from the experiment:", experiment, "is bimodal!", "\n")
    modeIs <- paste0("bimodal")
  }
  collectInfoFile <- c(experiment, cellId, species, modeIs)
  ## Export all information
  write.table(t(collectInfoFile), file = ModalFile, col.names = FALSE , row.names = FALSE ,append = TRUE, quote = FALSE, sep = "\t")
}

## Function to split scrna_seq_sample_info into 2 files: NEW_scRNASeq_info_file_final.tsv and Discard_scRNASeq_sample_info.tsv
splitInfoFiles <- function(ModalityFIle){
  
  for (species in unique(ModalityFIle$species)) {
    cat("Species:", species, "\n")
    for (experiment in unique(ModalityFIle$experiment[ModalityFIle$species == species])){
      cat("Name of Experiments:", experiment, "\n")
      for (cellId in unique(ModalityFIle$cellId[ModalityFIle$experiment == experiment])){
        cat("Name of cells:", cellId, "\n")
        modality <- as.character(ModalityFIle$modeDistribution[ModalityFIle$cellId == cellId]) 
        cat("Modality:", modality, "\n")
        
          if (modality == "bimodal"){
          cat("Cell-type is bimodal!", "\n")
          infoLib <- data.frame(annotation$libraryId[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeName == cellId])
          colnames(infoLib) <- "libraries"
          passQC <- annotation[annotation$libraryId %in% infoLib$libraries,]
          write.table(passQC, file = file.path(output_folder, "NEW_scRNASeq_sample_info.tsv"), col.names = FALSE, row.names = FALSE , append = TRUE, quote = FALSE, sep = "\t")
        } else {
          cat("Discard cell-type/experiment from the scRNASeq_info_file", "\n")
          infoLib <- data.frame(annotation$libraryId[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeName == cellId])
          colnames(infoLib) <- "libraries"
          addDiscardToDiscardFile <- annotation[annotation$libraryId %in% infoLib$libraries,]
          write.table(addDiscardToDiscardFile, file = file.path(output_folder, "Discard_scRNASeq_sample_info.tsv"), col.names = FALSE, row.names = FALSE , append = TRUE, quote = FALSE, sep = "\t")
        } 
      }
    }
  }
}

## Run the QC per cell-Type/experiment/species
for (species in unique(annotation$speciesId)) {
  cat("Species:", species, "\n")
  for (experiment in unique(annotation$experimentId[annotation$speciesId == species])){
    cat("Name of Experiments:", experiment, "\n")
    for (cellId in unique(annotation$cellTypeName[annotation$experimentId == experiment])){
      cat("Name of cells:", cellId, "\n")
      
      ### 1 STEP --> create matrix
      infoLib <- annotation$libraryId[annotation$experimentId == experiment & annotation$cellTypeName == cellId]
      file <- file.path(cells_folder, infoLib, "abundance_gene_level+fpkm+intergenic.tsv")
      cat("Number of cells:", length(file), "\n")
      All_libs <- lapply(file, read.delim)
      DATA <- do.call("cbind", All_libs)
      ## select colummns with gene_id, type and biotype
      Info_table <- DATA[,c(1,5,6)]
      colnames(Info_table) <- c("gene_id","type", "biotype")
      # select all columns with tpm info referent to each cell
      tpm <- DATA[, grepl("tpm", names( DATA))]
      mergeCells <- data.frame(Info_table, tpm)
      
      ### 2 STEP --> Check bimodality and plot cell distribution per experiment/species
      bimodality <- checkDataPlot(matrixData = mergeCells, output = output_folder, plot = plot)
      
      ### 3 Step --> generate output with experiments that pass or not the QC
      readModalityFile <- read.table(file.path(output_folder, "Modality_Cell_type_per_experiment.tsv"), header = TRUE, sep="\t")
      QC <- splitInfoFiles(ModalityFIle = readModalityFile)
    }
  }
}
