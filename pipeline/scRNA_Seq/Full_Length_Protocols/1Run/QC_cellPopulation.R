## SFonsecaCosta, Sept 2018

### This script is used to do the quality control of the data.
### In this script we work with abundance+geneLevel+intergenic.tsv output file from the previous step (kallisto followed by the analysis)

## 1 Step --> Create a matrix with all cells that belong to the same: cellTypeId, stageId, strain, uberonId, sex that belongs to same experiment and same species.
## 2 Step --> Verify if each of the conditions selected in the 1 Step follow a bimodal distribution
## 3 Step --> generate output file with experiments that pass the requirement

## Usage:
## R CMD BATCH --no-save --no-restore '--args scrna_seq_sample_info="scrna_seq_sample_info.tsv" cells_folder="cells_folder" sample_info_pass="sample_info_pass" sample_info_discarded="sample_info_discarded" modality_info="modality_info" calls_file_name="calls_file_name" plot="no"' QC_cellPopulation.R QC_cellPopulation.Rout
## scrna_seq_sample_info    --> collect information about each cell/experiment/species
## cells_folder             --> where we is located all the libraries/cells after Kallisto (treated data)
## sample_info_pass         --> path to file containing sample info that passed the quality control
## sample_info_discarded    --> path to file containing sample info that did not pass the quality control
## modality_info            --> path to file describing modality distribution per library
## calls_file_name          --> name of files containing tpm abundance at gene level
## plot                     --> if should plot the graphic results of each experiment (by default is no)

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
command_arg <- c("scrna_seq_sample_info", "cells_folder", "sample_info_pass", "sample_info_discarded", "modality_info", "calls_file_name", "plot")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

############################################### FUNCTIONS ############################################################
## Check bimodality + plot distribution of the protein coding for each cell-type
checkDataPlot <- function(matrixData, modalityInfoFile, plot){
  
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
    # create plots in same folder than modalityInfoFile
    output_folder <- dirname(modalityInfoFile)
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
    
    pdf(file = file.path(output_folder, paste0("Plot","_", cellId, "_", stageId, "_", uberonId, "_", strain, "_", sex,"_", experiment , "_", species,".pdf")), width = 12, height = 12)
    grid.arrange(p1,p2,p3,ncol = 1, nrow = 3)
    dev.off()
  } else {
    message("Analysis done without reporting graphics.")
  }
  
  ## test modality of the data
  if ( is.unimodal(codingGenes$genes) == TRUE){
    message("The cell-type ",cellId," from the experiment:", experiment,"is unimodal for the protein coding genes!", "\n", "Warning: this experiment can be removed...")
    modeIs <- paste0("unimodal")
  } else if (is.bimodal(codingGenes$genes) == FALSE){
    multimodal <- is.multimodal(codingGenes$genes)
    trimodal <- is.trimodal(codingGenes$genes)
    message("The cell-type ",cellId," from the experiment:", experiment,"is trimodal", trimodal, "or multimodal:", multimodal)
    modeIs <- paste0("multimodal")
  } else {
    message("The cell-type ",cellId," from the experiment:", experiment, "is bimodal!")
    modeIs <- paste0("bimodal")
  }
  collectInfoFile <- c(experiment, cellId, species, stageId, uberonId, sex, strain, modeIs)
  ## Export all information
  write.table(t(collectInfoFile), file = modalityInfoFile, col.names = FALSE , row.names = FALSE ,append = TRUE, quote = FALSE, sep = "\t")
}

## Function to split scrna_seq_sample_info into 2 files: NEW_scRNASeq_info_file_final.tsv and Discard_scRNASeq_sample_info.tsv
splitInfoFiles <- function(modality_info, sample_info, sample_info_pass, sample_info_discarded){
  uniq_cond <- unique(modality_info[,c("speciesId", "experimentId", "cellTypeId", "stageId", "strain", "uberonId", "sex")])

  for (i in seq_len(nrow(all_uniq_cond))){

    message("Species:", uniq_cond$speciesId[i], " - experiment: ", uniq_cond$experimentId[i], " - cellType: ", uniq_cond$cellTypeId[i], 
      " - stage: ", uniq_cond$stageId[i], " - strain:", uniq_cond$strain[i], " - uberonId: ", uniq_cond$uberonId[i], 
      " - sex: ", uniq_cond$sex[i])
                 
    modality <- as.character(modality_info$modeDistribution[modality_info$speciesId == uniq_cond$speciesId[i] & 
      modality_info$experiment == uniq_cond$experimentId[i] & modality_info$cellTypeId == uniq_cond$cellTypeId[i] & 
      modality_info$stageId == uniq_cond$stageId[i] & modality_info$strain == uniq_cond$strain[i] & 
      modality_info$uberonId == uniq_cond$uberonId[i] & modality_info$sex == uniq_cond$sex[i]])
  
    message("Modality:", modality)
                
    if (length(modality != 0) && modality == "bimodal"){
      message("Cell-type is bimodal!")
      infoLib <- data.frame(annotation$libraryId[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex])
      colnames(infoLib) <- "libraries"
      passQC <- annotation[annotation$libraryId %in% infoLib$libraries,]
      write.table(passQC, file = sample_info_pass, col.names = FALSE, row.names = FALSE , append = TRUE, quote = FALSE, sep = "\t")
    } else if (length(modality != 0) && modality != "bimodal") {
      message("Discard cell-type/experiment from the scRNASeq_info_file")
      infoLib <- data.frame(annotation$libraryId[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex])
      colnames(infoLib) <- "libraries"
      addDiscardToDiscardFile <- annotation[annotation$libraryId %in% infoLib$libraries,]
      write.table(addDiscardToDiscardFile, file = sample_info_discarded, col.names = FALSE, row.names = FALSE , append = TRUE, quote = FALSE, sep = "\t")
    } else {
      message("Not take in consideration this combination to split scRNA-Seq info file.")
    }

  }
}



## Read scrna_seq_sample_info file. If file not exists, script stops
if( file.exists(scrna_seq_sample_info) ){
  annotation <- read.table(scrna_seq_sample_info, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
  ## remove special characters in strain
  annotation$strain <- gsub('\\/', '_', annotation$strain)
  annotation$strain <- gsub('\\s+', '_', annotation$strain)
} else {
  stop( paste("scrna_seq_sample_info file not found [", scrna_seq_sample_info, "]\n"))
}

## Create output files and check directories already exist
if (!dir.exists(dirname(modality_info))){
  dir.create(dirname(modality_info), recursive=TRUE)
}
file.create(modality_info)
cat("experiment\tcellTypeId\tspeciesId\tstageId\tuberonId\tsex\tstrain\tmodeDistribution\n",file = modality_info, sep = "\t")
 
if (!dir.exists(dirname(sample_info_pass))){
  dir.create(dirname(sample_info_pass))
}  
file.create(sample_info_pass)
cat("libraryId\texperimentId\tcellTypeName\tcellTypeId\tspeciesId\tplatform\tprotocol\tprotocolType\tlibraryType\tinfoOrgan\tstageId\tuberonId\tsex\tstrain\treadLength\torganism\n",file = sample_info_pass, sep = "\t")

if (!dir.exists(dirname(sample_info_discarded))){
  dir.create(dirname(sample_info_discarded))
}
file.create(sample_info_discarded)
cat("libraryId\texperimentId\tcellTypeName\tcellTypeId\tspeciesId\tplatform\tprotocol\tprotocolType\tlibraryType\tinfoOrgan\tstageId\tuberonId\tsex\tstrain\treadLength\torganism\n",file = sample_info_discarded, sep = "\t")

## Run the QC per condition
uniq_cond <- unique(annotation[,c("speciesId", "experimentId", "cellTypeId", "stageId", "strain", "uberonId", "sex")])

for (i in seq_len(nrow(all_uniq_cond))){

  message("Species:", uniq_cond$speciesId[i], " - experiment: ", uniq_cond$experimentId[i], " - cellType: ", uniq_cond$cellTypeId[i], 
    " - stage: ", uniq_cond$stageId[i], " - strain:", uniq_cond$strain[i], " - uberonId: ", uniq_cond$uberonId[i], " - sex: ", uniq_cond$sex[i])

  infoLib <- annotation$libraryId[annotation$speciesId == uniq_cond$speciesId[i] & annotation$experimentId == uniq_cond$experimentId[i] & 
    annotation$cellTypeId == uniq_cond$cellTypeId[i] & annotation$stageId == uniq_cond$stageId[i] & annotation$strain == uniq_cond$strain[i] & 
    annotation$uberonId == uniq_cond$uberonId[i] & annotation$sex == uniq_cond$sex[i]]

  file <- file.path(cells_folder, infoLib, calls_file_name)
  message("Number of cells:", length(infoLib))

  if (length(file) != 0){
    All_libs <- lapply(file, read.delim)
    DATA <- do.call("cbind", All_libs)
    ## select colummns with gene_id, type and biotype
    Info_table <- DATA[,c("id","type","biotype")]
    colnames(Info_table) <- c("gene_id","type", "biotype")
    # select all columns with tpm info referent to each cell
    tpm <- DATA[, grepl("abundance", names( DATA))]
    mergeCells <- data.frame(Info_table, tpm)
    ### 2 STEP --> Check bimodality and plot cell distribution per experiment/species
    ## save plot in same directory than file describing modality 
    bimodality <- checkDataPlot(matrixData = mergeCells, modalityInfoFile = modality_info, plot = plot)
  } else {
    message("Do not take in consideration the combination [", as.character(uniq_cond[i]), "]to analyse the QC (0 cells).")
  }
}

### 3 Step --> generate output with experiments that pass or not the QC
readModalityFile <- read.table(modality_info, header = TRUE, sep="\t")
QC <- splitInfoFiles(modality_info = readModalityFile, sample_info = annotation, sample_info_pass = sample_info_pass, sample_info_discarded = sample_info_discarded)
