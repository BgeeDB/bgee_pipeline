## SFonsecaCosta, Oct 2018

## This script is used after call present genes per individual library/cell
## We sum the calls across the cell population for each gene and compute the ratio for genic and intergenic regions + calculate the cut-offs based on density and proportion of intergenic for the cell population

## Usage:
## R CMD BATCH --no-save --no-restore '--args NEW_scRNASeq_sample_info="NEW_scRNASeq_sample_info.tsv" cells_folder="cells_folder" refIntergenicFolder="refIntergenic_folder" output_folder="output_folder" ratioValue="ratioValue"' Sum_Calls_cellPopulation.R Sum_Calls_cellPopulation.Rout
## NEW_scRNASeq_sample_info --> info file with all libraries/cell
## cells_folder --> where we is located all the libraries/cells after Kallisto (treated data)
## refIntergenicFolder --> Folder where is the reference Intergenic regions per species
## output_folder --> output where should be written the Sum_Calls_cellPopulation
## ratioValue --> proportion of intergenic allowed for the ratio cut-off

## Libraries used
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(Biostrings)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed in command line
command_arg <- c("NEW_scRNASeq_sample_info", "cells_folder", "refIntergenicFolder", "output_folder", "ratioValue")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read NEW_scRNASeq_sample_info file. If file not exists, script stops
if( file.exists(NEW_scRNASeq_sample_info) ){
  annotation <- read.table(NEW_scRNASeq_sample_info, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
  ## remove special characters in strain
  annotation$strain <- gsub('\\/', '_', annotation$strain)
  annotation$strain <- gsub('\\s+', '_', annotation$strain)
} else {
  stop( paste("NEW_scRNASeq_sample_info file not found [", NEW_scRNASeq_sample_info, "]\n"))
}

##########################################################################################################################################
## Collect information after the Sum_Calls per cell population
stats_data <- function(finalMatrix){
  
  ## genic present
  genic <- nrow(dplyr::filter(finalMatrix, type == "genic"))
  genic_present <- nrow(dplyr::filter(finalMatrix, type == "genic" & cellsPresent != 0))
  ## protein-coding present
  protein_coding <- nrow(dplyr::filter(finalMatrix, biotype == "protein_coding"))
  protein_coding_present <- nrow(dplyr::filter(finalMatrix, biotype == "protein_coding" & cellsPresent != 0))
  proportion_codingPresent <- (protein_coding_present/protein_coding)*100
  ## intergenic present
  intergenic <- nrow(dplyr::filter(finalMatrix, type == "intergenic"))
  intergenic_present <- nrow(dplyr::filter(finalMatrix, type == "intergenic" & cellsPresent != 0))
  proportion_intergenicPresent <- (intergenic_present/intergenic)*100
  
  ## Collect info for protein coding
  # Number genes (protein coding) that are always absent in all cells
  total_pc_zeros <- nrow(dplyr::filter(finalMatrix, stats == "0" & biotype == "protein_coding"))
  # Number genes (protein coding) that are always present in all cells
  total_pc_100 <- nrow(dplyr::filter(finalMatrix, stats == "1" & biotype == "protein_coding"))
  # Number genes (protein coding) that are in at least 50% of the cells
  total_pc_50 <- nrow(dplyr::filter(finalMatrix, stats >= "0.5" & biotype == "protein_coding"))
  # Number genes (protein coding) that are less than 50% of the cells
  total_pc_less_50 <- nrow(dplyr::filter(finalMatrix, stats < "0.5" & biotype == "protein_coding"))
  
  ## Collect info for Intergenic region
  # Number intergenic regions that are always absent in all cells
  total_intergenic_zeros <- nrow(dplyr::filter(finalMatrix, stats == "0" & type == "intergenic"))
  # Number intergenic regions that are always present in all cells
  total_intergenic_all <- nrow(dplyr::filter(finalMatrix, stats == "1" & type == "intergenic"))
  # Number intergenic regions that are in at least 50% of the cells
  total_intergenic_50 <- nrow(dplyr::filter(finalMatrix, stats >= "0.5" & type == "intergenic"))
  # Number intergenic regions that are less than 50% of the cells
  total_intergenic_less_50 <- nrow(dplyr::filter(finalMatrix, stats < "0.5" & type == "intergenic"))
  
  ## Collect info after cut-off density
  ## genic present after cutoff
  genic_present_cutoff <- nrow(dplyr::filter(finalMatrix, type == "genic" & call_CutOff_density == "present"))
  ## protein-coding present after cutoff
  protein_coding_present_cutoff <- nrow(dplyr::filter(finalMatrix, biotype == "protein_coding" & call_CutOff_density == "present"))
  proportion_codingPresent_cutoff <- (protein_coding_present_cutoff/protein_coding)*100
  ## intergenic present after cutoff
  intergenic_present_cutoff <- nrow(dplyr::filter(finalMatrix, type == "intergenic" & call_CutOff_density == "present"))
  proportion_intergenicPresent_cutoff <- (intergenic_present_cutoff/intergenic)*100
  
  ## Collect info after cut-off ratio
  ## genic present after cutoff
  genic_present_cutoff_ratio <- nrow(dplyr::filter(finalMatrix, type == "genic" & call_CutOff_ratio == "present"))
  ## protein-coding present after cutoff
  protein_coding_present_cutoff_ratio <- nrow(dplyr::filter(finalMatrix, biotype == "protein_coding" & call_CutOff_ratio == "present"))
  proportion_codingPresent_cutoff_ratio <- (protein_coding_present_cutoff_ratio/protein_coding)*100
  ## intergenic present after cutoff
  intergenic_present_cutoff_ratio <- nrow(dplyr::filter(finalMatrix, type == "intergenic" & call_CutOff_ratio == "present"))
  proportion_intergenicPresent_cutoff_ratio <- (intergenic_present_cutoff_ratio/intergenic)*100
  
  ## Export informative file
  information_file <- c(experiment, cellId, genic, genic_present, protein_coding, protein_coding_present, proportion_codingPresent, intergenic, intergenic_present,proportion_intergenicPresent,
                        total_pc_zeros, total_pc_100, total_pc_50, total_pc_less_50,
                        total_intergenic_zeros, total_intergenic_all, total_intergenic_50, total_intergenic_less_50,
                        genic_present_cutoff, protein_coding_present_cutoff, proportion_codingPresent_cutoff, intergenic_present_cutoff, proportion_intergenicPresent_cutoff,
                        genic_present_cutoff_ratio, protein_coding_present_cutoff_ratio, proportion_codingPresent_cutoff_ratio, intergenic_present_cutoff_ratio, proportion_intergenicPresent_cutoff_ratio,
                        species, organism)
  
  return(information_file)
}

## Provide the reference intergenic regions = TRUE
refIntergenic <- function(counts, folder_refIntergenic, speciesID){
  
  referenceIntergenic <- paste0(folder_refIntergenic, "/", speciesID, "_intergenic.fa")
  referenceIntergenic <- readDNAStringSet(referenceIntergenic)
  seq_name <- names(referenceIntergenic)
  seq_name <- gsub( " .*$", "", seq_name )
  ##seq_name <- gsub( "_", "-", seq_name)
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

######################################################################################################################
## Cut-off function (just select protein-coding and reference intergenic region) + Plot and export info per cell-type
cutoff_used <- function(finalMatrix, output_folder, ratioValue, folder_refIntergenic){
  
  selected_coding <- finalMatrix$biotype %in% "protein_coding"
  selected_intergenic <- (finalMatrix$type %in% "intergenic" & referenceIntergenic$refIntergenic == "TRUE")
  
  ## intergenic
  summed_intergenic <- sapply(unique(sort(finalMatrix$stats[selected_coding])), function(x){
    return( sum(finalMatrix$stats[selected_intergenic] >= x) )
  })
  
  ## protein coding
  summed_coding <- c(0, cumsum(rle(sort(finalMatrix$stats[selected_coding]))$lengths))
  summed_coding <- summed_coding[-(length(summed_coding))]
  summed_coding <- sum(selected_coding) - summed_coding
  
  ## Now we can calculate r
  r <- ( summed_intergenic / sum(selected_intergenic) ) /
    ( summed_coding / sum(selected_coding) )
  
  percent <- (1-ratioValue)*100
  
  ## Select the minimal value of ratio_cutoff for which the ratio of genes and intergenic regions is equal to 0.05 or lower
  if (sum(r < ratioValue) == 0){
    ratio_cutoff <- sort(unique(finalMatrix$stats[selected_coding]))[which(r == min(r))[1]]
    r_cutoff <- min(r)
    cat(paste0("    There ratio cutoff for which " , percent,"%", " of the expressed genes would be coding. Ratio cutoff is fixed at the first value with maximum coding/intergenic ratio. r=", r_cutoff, " at ratio=", ratio_cutoff,"\n"))
  } else {
    ratio_cutoff <- sort(unique(finalMatrix$stats[selected_coding]))[which(r < ratioValue)[1]]
    r_cutoff <- ratioValue
    cat(paste0("    The ratio cutoff for which " , percent,"%", " of the expressed genes are be coding found at ratio=", ratio_cutoff,"\n"))
  }
  
  ##  Cut-off using densities
  total_protein_coding <- dplyr::filter(finalMatrix, biotype == "protein_coding")
  total_protein_coding <- total_protein_coding %>% select(gene_id, type, biotype, stats)
  
  total_intergenic <- dplyr::filter(finalMatrix, type == "intergenic" & referenceIntergenic$refIntergenic == "TRUE")
  total_intergenic <- total_intergenic %>% select(gene_id, type, biotype, stats)
  
  FinalTable <- rbind(total_protein_coding, total_intergenic)
  reOrderTable <- melt(FinalTable)
  
  p <- ggplot(reOrderTable) + geom_density(aes(x = log(value+0.000001), colour = type))
  p <- ggplot_build(p)
  
  extractValues <- p$data[[1]]
  intergenicGroup <- dplyr::filter(extractValues, group == "1")
  codingGroup <- dplyr::filter(extractValues, group == "2")
  
  generalTab <- data.frame(intergenicGroup$x,intergenicGroup$density, codingGroup$density)
  generalTab$Int_less_Pc <- generalTab$intergenicGroup.density < generalTab$codingGroup.density
  
  first_cutoff <- generalTab$intergenicGroup.x[generalTab$Int_less_Pc != "TRUE" & generalTab$intergenicGroup.x >= log(1/sizeMatrix)][1]
  
  ## collect results --> without cutoff
  finalMatrix$calls_without_CutOff <- finalMatrix$stats > 0
  finalMatrix$calls_without_CutOff <- gsub('TRUE', 'present', finalMatrix$calls_without_CutOff)
  finalMatrix$calls_without_CutOff <- gsub('FALSE', '-', finalMatrix$calls_without_CutOff)
  
  ## collect results --> cutoff density
  finalMatrix$call_CutOff_density <- finalMatrix$stats >= exp(first_cutoff)
  finalMatrix$call_CutOff_density <- gsub('TRUE', 'present', finalMatrix$call_CutOff_density)
  finalMatrix$call_CutOff_density <- gsub('FALSE', '-', finalMatrix$call_CutOff_density)
  
  ## collect results --> cutoff ratio
  finalMatrix$call_CutOff_ratio <- finalMatrix$stats >= ratio_cutoff
  finalMatrix$call_CutOff_ratio <- gsub('TRUE', 'present', finalMatrix$call_CutOff_ratio)
  finalMatrix$call_CutOff_ratio <- gsub('FALSE', '-', finalMatrix$call_CutOff_ratio)
  
  ## add confidence column
  finalMatrix$confidence <- ifelse(finalMatrix$stats < exp(first_cutoff), "-", ifelse(finalMatrix$stats >= exp(first_cutoff) & finalMatrix$stats < ratio_cutoff, "Low", "High"))
  
  ## collect information about protein coding after cutoff
  finalMatrixInfo <- finalMatrix %>% select(gene_id, type, biotype, calls_without_CutOff, call_CutOff_density, call_CutOff_ratio)
  propotioncalls_without_CutOff <- (nrow(dplyr::filter(finalMatrixInfo, biotype == "protein_coding" & calls_without_CutOff == "present")) / nrow(total_protein_coding))*100
  propotionCalls_cutoffDensity <- (nrow(dplyr::filter(finalMatrixInfo, biotype == "protein_coding" & call_CutOff_density == "present")) / nrow(total_protein_coding))*100
  propotionCalls_cutoffRatio <- (nrow(dplyr::filter(finalMatrixInfo, biotype == "protein_coding" & call_CutOff_ratio == "present")) / nrow(total_protein_coding))*100
  
  InfoProportionCoding <- t(data.frame(propotioncalls_without_CutOff,propotionCalls_cutoffDensity,propotionCalls_cutoffRatio))
  colnames(InfoProportionCoding) <- "proportionCoding";
  InfoProportionCoding <- as.data.frame(InfoProportionCoding)
  InfoProportionCoding$ID <- c("without_cutoff", "cutoff_density","cutoff_ratio")
  
  ## plot the results
  pdf(file = file.path(output_folder, paste0("Genes_present_per_cell_ratio", "_", experiment, "_", cellId ,"_", stageId, "_", strain, "_", uberonId, "_", sex, "_",species ,".pdf")), width = 14, height = 7)
  p1 <- p$plot +
    labs(title = "Density distribution of protein coding vs reference intergenic", x="Log(ratio)",y="Density") +
    scale_color_discrete(name = "Region", labels = c("Protein coding", "Intergenic")) +
    geom_vline(aes(xintercept=first_cutoff, linetype="dotted"), color="black", size=0.5) +
    geom_vline(aes(xintercept=log(ratio_cutoff), linetype="solid"), color = "black", size=0.5) +
    scale_linetype_manual(name = "Cut-Off", values = c("dotted", "solid"), labels = c(paste0("Cut-off density (", round(exp(first_cutoff)*sizeMatrix, digits = 0), " cell/s)"), paste0("Cut-off ratio (", round(ratio_cutoff*sizeMatrix, digits = 0) , " cell/s)")))
  
  p2 <- ggplot(data = InfoProportionCoding, aes(x = ID, y = proportionCoding)) +
    geom_boxplot() +
    geom_hline(yintercept=c(70,80), linetype="dashed", color = "red") +
    labs(title = "Proportion of coding genes present", x=" ",y="%") +
    coord_cartesian(ylim = c(0, 100))
  grid.arrange(p1,p2,ncol = 2, nrow = 1)
  dev.off()
  
  ## write informative matrix for each independent experiment
  write.table(finalMatrix, file = file.path(output_folder, paste0("Sum_Call_", experiment, "_", cellId ,"_", stageId, "_", strain, "_", uberonId, "_", sex, "_", species ,".tsv")), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  return(finalMatrix)
}

################################################################################
statsFile <- file.path(output_folder, "Stats_SumFile.tsv")
if (!file.exists(statsFile)){
  file.create(statsFile)
  cat("experimentId\tcellID\tgenic\tgenic_present\tprotein_coding\tprotein_coding_present\tproportion_codingPresent\tintergenic\tintergenic_present\tproportion_intergenicPresent\tProteinCoding_genes_always_absent\tProteinCoding_genes_always_present\tProteinCoding_genes_present_at_least_50%_cells\tProteinCoding_genes_present_less_50%_cells\tIntergenic_region_always_absent\tIntergenic_region_always_present\tIntergenic_region_present_at_least_50%_cells\tIntergenic_region_present_less_50%_cells\tgenic_present_cutoff\tprotein_coding_present_cutoff\tproportion_codingPresent_cutoff\tintergenic_present_cutoff\tproportion_intergenicPresent_cutoff\tgenic_present_cutoff_ratio\tprotein_coding_present_cutoff_ratio\tproportion_codingPresent_cutoff_ratio\tintergenic_present_cutoff_ratio\tproportion_intergenicPresent_cutoff_ratio\tspeciesId\torganism\n", file = file.path(output_folder, "Stats_SumFile.tsv"), sep = "\t")
} else {
  print("File already exist.....")
}

for (species in unique(annotation$speciesId)) {
  cat("Species:", species, "\n")
  for (experiment in unique(annotation$experimentId[annotation$speciesId == species])){
    cat("Name of experiments that belongs to the species:", experiment, "\n")
    for (cellId in unique(annotation$cellTypeId[annotation$experimentId == experiment])){
      cat("CellId info:", cellId, "\n")
      for (stageId in unique(annotation$stageId[annotation$cellTypeId == cellId])){
        cat("StageId info:", stageId, "\n")
        for (strain in unique(annotation$strain[annotation$stageId == stageId])){
          cat("Strain info:", strain, "\n")
          for (uberonId in unique(annotation$uberonId[annotation$strain == strain])){
            cat("UberonId info:", uberonId, "\n")
            for (sex in unique(annotation$sex[annotation$uberonId == uberonId])){
              cat("sex info:", sex, "\n")
              
              organism <- as.character(unique(annotation$organism[annotation$speciesId == species]))
              infoLib <- annotation$libraryId[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex]
              file <- file.path(cells_folder, infoLib)
              AllFiles <- list.files(file,  pattern="abundance_gene_level\\+fpkm\\+intergenic\\+calls.tsv$", full.names=T, recursive = TRUE)
              print(AllFiles)
              AllFiles <- lapply(AllFiles, read.delim)
              
              DATA <- do.call("cbind", AllFiles)
              ## select colummns with gene_id, type and biotype
              info_table <- DATA[,c(1,5,6)]
              colnames(info_table) <- c("gene_id","type", "biotype")
              call <- DATA[, grepl("call", names( DATA))]
              ## Create a full matrix with all information
              merge_all_info <- data.frame(info_table, call)
              
              ## replace string present/- (note we use "-" because we not call absent genes in individual libraries) by number 1/0
              finalMatrix <- as.matrix(merge_all_info)
              finalMatrix <- gsub("present", "1", finalMatrix)
              finalMatrix <- gsub("-", "0", finalMatrix)
              finalMatrix <- as.data.frame(finalMatrix)
              ## extract size of the experiment
              sizeMatrix <- length(4:ncol(finalMatrix))
              ## How many times each gene is called present per cell
              finalMatrix$cellsPresent <- rowSums(finalMatrix[ ,4:ncol(finalMatrix)] != 0)
              ## Proportion of present regarding the number of cells (ratio)
              finalMatrix$stats <- (finalMatrix$cellsPresent)/sizeMatrix
              
              ## Information about reference intergenic
              referenceIntergenic <- refIntergenic(counts = finalMatrix, folder_refIntergenic = refIntergenicFolder, speciesID = species)
              
              ## calls 
              calls <- cutoff_used(finalMatrix = finalMatrix, output_folder = output_folder, ratioValue=ratioValue, folder_refIntergenic=referenceIntergenic)
              
              ## provide stats
              statistics_experiment <- stats_data(finalMatrix=calls)
              this_sample <- as.data.frame(t(statistics_experiment), stringsAsFactors=F)
              ## write informative matrix for each independent experiment
              write.table(this_sample, file = statsFile,col.names =F , row.names = F,append = T,quote = FALSE, sep = "\t")
            }
          }
        }
      }
    }
  }
}

### Final plot --> per cellId and organism
plotData <- function(dataFile){
  dataFile <- read.table(dataFile, header=TRUE, sep="\t")
  proportion_all_samples <- dataFile %>% select(experimentId, cellID, proportion_codingPresent,
                                                proportion_codingPresent_cutoff, proportion_codingPresent_cutoff_ratio, organism)
  names(proportion_all_samples)[names(proportion_all_samples) == 'proportion_codingPresent_cutoff'] <- 'proportion_codingPresent_cutoff_density'
  proportion_all_samples <- melt(proportion_all_samples)
  g1 <- ggplot(proportion_all_samples, aes(x = organism, y = value, fill = variable, color = variable)) +
    geom_point() + coord_cartesian(ylim = c(0, 100)) +
    geom_hline(yintercept=c(70,80), linetype="dashed", color = "black") +
    labs(title = paste0("Sum Calls - Intergenic"), x="Organism",y="% protein-coding genes present")
  
  g2 <- ggplot(proportion_all_samples, aes(x = experimentId, y = value, fill = cellID, color = cellID)) +
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 100)) + facet_wrap(~variable) +
    geom_hline(yintercept=c(70,80), linetype="dashed", color = "black") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    labs(title = paste0("Sum Calls - Intergenic"), x="Organism",y="% protein-coding genes present")
  
  pdf(file = file.path(output_folder, "Summary_all_experiments_without_and_with_cutoff.pdf"), width = 16, height = 10)
  grid.arrange(g1,g2,ncol = 1, nrow = 2)
  dev.off()
}
plotData(statsFile)
