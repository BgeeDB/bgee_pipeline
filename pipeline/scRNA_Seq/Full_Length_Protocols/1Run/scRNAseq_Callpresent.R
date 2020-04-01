## SFonsecaCosta, Nov 2018
## Adaptation from Julien Roux script (Bgee RNA-Seq pipeline - rna_seq_presence_absence.r)
## Apply the cutoff to call present genes for each library that correspond to each individual cell

## Usage:
## R CMD BATCH --no-save --no-restore '--args NEW_scRNASeq_sample_info="NEW_scRNASeq_sample_info" cells_folder="cells_folder" sum_species="sum_species" gaussian_choice="gaussian_choice" ratioValue="ratioValue" output_folder="output_folder"' scRNAseq_Callpresent.R scRNAseq_Callpresent.Rout
## NEW_scRNASeq_sample_info --> file with info on mapped libraries
## cells_folder --> where we is located all the libraries/cells after Kallisto (treated data)
## sum_species --> where is localized sum per species
## gaussian_choice --> gaussian choice
## ratioValue --> ratio used to perform the cutoff (should be between 0 and 1, where 1 menas 100%)
## output_folder --> folder where to export the plots, summed data, and classification of coding/intergenic regions

## libraries used
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)

## reading in arguments provided in command line
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}
## checking if all necessary arguments were passed in command line
command_arg <- c("NEW_scRNASeq_sample_info", "cells_folder", "sum_species", "gaussian_choice", "ratioValue", "output_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}
## Read NEW_scRNASeq_sample_info file. If file not exists, script stops
if( file.exists(NEW_scRNASeq_sample_info) ){
  annotation <- read.table(NEW_scRNASeq_sample_info, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
} else {
  stop( paste("NEW_scRNASeq_sample_info file not found [", NEW_scRNASeq_sample_info, "]\n"))
}
## File indicating gaussian_choice
if( file.exists(gaussian_choice) ){
  gaussian <- read.table(gaussian_choice, h=T, sep="\t")
} else {
  warning( paste("File with gaussian choices not found [", gaussian_choice, "]\n"))
}

###################################################################################################################
### Calculate r for the cut-off
calculate_r <- function(counts, selected_coding, selected_intergenic, ratioValue){
  summed_intergenic <- sapply(unique(sort(counts$tpm[selected_coding])), function(x){
    return( sum(counts$tpm[selected_intergenic] >= x) )
  })
  summed_coding <- c(0, cumsum(rle(sort(counts$tpm[selected_coding]))$lengths))
  summed_coding <- summed_coding[-(length(summed_coding))]
  summed_coding <- sum(selected_coding) - summed_coding

  r <- ( summed_intergenic / sum(selected_intergenic) ) / ( summed_coding / sum(selected_coding) )
  percent <- (1-ratioValue)*100

  if (sum(r < ratioValue) == 0){
    TPM_cutoff <- sort(unique(counts$tpm[selected_coding]))[which(r == min(r))[1]]
    r_cutoff <- min(r)
    cat(paste0("There is no TPM cutoff for which " , percent,"%", " of the expressed genes would be coding.", "\n",
    "TPM cutoff is fixed at the first value with maximum coding/intergenic ratio.", "\n",
    "r=", r_cutoff, "at TPM=", TPM_cutoff,"\n"))
  } else {
    TPM_cutoff <- sort(unique(counts$tpm[selected_coding]))[which(r < ratioValue)[1]]
    r_cutoff <- ratioValue
    cat(paste0("TPM cutoff for which " , percent,"%", " of the expressed genes are be coding found at TPM=", TPM_cutoff,"\n"))
  }
  return(list(TPM_cutoff, r_cutoff))
}

### Collect MaxIntergenic region
collectmaxIntergenic <- function(gaussian, sum_by_species){
  collect_intergenic <- paste0("intergenic_", gaussian$selectedGaussianIntergenic[gaussian$speciesId == species])

  if (gaussian$selectionSideIntergenic[gaussian$speciesId == species] == "Left"){
    max_intergenic <- max(sum_by_species$tpm[sum_by_species$classification %in% collect_intergenic])
  } else if (gaussian$selectionSideIntergenic[gaussian$speciesId == species] == "Right") {
    max_intergenic <- min(sum_by_species$tpm[sum_by_species$classification %in% collect_intergenic])
  }
  return(max_intergenic)
}

### Calculate FPKM and TPM cutOff (FPKM_cutoff = use genic+intergenic; FPKM_final_cutoff = just genic; TPM_final_cutoff = just genic)
calculateFinalTPM <- function(counts, gene_counts, TPM_cutoff){
  FPKM_cutoff <-  TPM_cutoff * na.omit(counts$fpkm / counts$tpm)[1]
  FPKM_final_cutoff <- FPKM_cutoff * sum(counts$est_counts)/sum(gene_counts$est_counts)
  TPM_final_cutoff <- FPKM_final_cutoff * na.omit(gene_counts$tpm / gene_counts$fpkm)[1]
  return(list(FPKM_cutoff,FPKM_final_cutoff, TPM_final_cutoff))
}

### Collect information after cut-off
cutoff_info <- function(library_id, counts, max_intergenic, TPM_cutoff, TPM_final_cutoff, r_cutoff){
  proportion_genic_present <- (nrow(dplyr::filter(counts, type == "genic" & call == "present"))/nrow(dplyr::filter(counts, type == "genic")))*100
  genic_region_present <- nrow(dplyr::filter(counts, type == "genic" & call == "present"))
  proportion_coding_present <- (nrow(dplyr::filter(counts, biotype == "protein_coding" & call == "present"))/nrow(dplyr::filter(counts, biotype == "protein_coding")))*100
  coding_region_present  <- nrow(dplyr::filter(counts, biotype == "protein_coding" & call == "present"))
  proportion_intergenic_present <- (nrow(dplyr::filter(counts, type == "intergenic" & call == "present"))/nrow(dplyr::filter(counts, type == "intergenic")))*100
  intergenic_region_present <- nrow(dplyr::filter(counts, type == "intergenic" & call == "present"))
  ## Export cutoff_info_file
  collectInfo <- c(library_id,TPM_cutoff,proportion_genic_present, genic_region_present,sum(counts$type == "genic"),
                   proportion_coding_present, coding_region_present, sum(counts$biotype %in% "protein_coding"),
                   proportion_intergenic_present, intergenic_region_present, sum(counts$type == "intergenic"),
                   r_cutoff)
  names(collectInfo) <- c("libraryId","cutoffTPM", "proportionGenicPresent", "numberGenicPresent", "Genic",
                          "proportionCodingPresent",  "numberCodingPresent", "ProteinCoding",
                          "proportion_intergenic_present", "intergenic_region_present", "Intergenic","ratioIntergenicCodingPresent")
  return(collectInfo)
}

### Plot the data
plotData <- function(dataFile){
  dataFile <- read.table(dataFile, header=TRUE, sep="\t")
  proportion_all_samples <- dataFile[c(3,6,9,14)]
  proportion_all_samples <- melt(proportion_all_samples)
  g1 <- ggplot(proportion_all_samples, aes(organism,value, fill = organism)) +
    geom_boxplot() + stat_smooth() +
    coord_cartesian(ylim = c(0, 100)) + facet_wrap(~variable) +
    geom_hline(yintercept=c(70,80), linetype="dashed", color = "black") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

  number_all_samples <- dataFile[c(4,7,10,14)]
  number_all_samples <- melt(number_all_samples)
  g2 <- ggplot(number_all_samples, aes(organism,value, fill = organism)) +
    geom_boxplot() + stat_smooth() +
    facet_wrap(~variable) + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

  pdf(file = file.path(output_folder, "Plot_after_calls.pdf"), width = 16, height = 10)
  grid.arrange(g1,g2,ncol = 1, nrow = 2)
  dev.off()
}

collectSamples <- data.frame()
for(species in unique(annotation$speciesId)){
  sum_by_species <- read.table(paste0(sum_species, "/sum_abundance_gene_level+fpkm+intergenic+classification_", species, ".tsv"), h=T, sep="\t")
  print(species)
  for(libraryId in annotation$libraryId[annotation$speciesId == species]){

        cat("Doing the library: ",libraryId, "\n")
        ## genic + intergenic
        kallisto_gene_counts <- read.table(file.path(cells_folder, libraryId, "abundance_gene_level+fpkm+intergenic.tsv"), h=T, sep="\t")
        ## genic
        gene_counts <- read.table(file.path(cells_folder, libraryId, "abundance_gene_level+new_tpm+new_fpkm.tsv"), h=T, sep = "\t")

        selected_coding <- kallisto_gene_counts$biotype %in% "protein_coding"
        max_intergenic <- collectmaxIntergenic(gaussian, sum_by_species)
        selected_intergenic <- (kallisto_gene_counts$type == "intergenic" & sum_by_species$tpm < max_intergenic)
        results <- calculate_r(kallisto_gene_counts, selected_coding, selected_intergenic, ratioValue = as.numeric(ratioValue))
        TPM_cutoff <- results[[1]]
        r_cutoff <- results[[2]]
        recalculateTPM <- calculateFinalTPM(kallisto_gene_counts, gene_counts, TPM_cutoff=TPM_cutoff)

        ## Setting the presence flags
        gene_counts$call <- ifelse(gene_counts$tpm >= recalculateTPM[[3]], "present", "-")
        kallisto_gene_counts$call <- ifelse(kallisto_gene_counts$tpm >= TPM_cutoff, "present", "-")

        ## Export cutoff infoormation file + new files with calls
        cutoff_info_file <- cutoff_info(libraryId, kallisto_gene_counts, max_intergenic, TPM_cutoff, a[[3]], r_cutoff)
        pathExport <- file.path(cells_folder, libraryId)
        write.table(gene_counts,file = file.path(pathExport, "abundance_gene_level+new_tpm+new_fpkm+calls.tsv"),quote=F, sep = "\t", col.names=T, row.names=F)
        write.table(kallisto_gene_counts,file = file.path(pathExport, "abundance_gene_level+fpkm+intergenic+calls.tsv"),quote=F, sep = "\t", col.names=T, row.names=F)
        write.table(t(t(cutoff_info_file)),file = file.path(pathExport, "cutoff_info_file.tsv"),quote=F, sep = "\t", col.names=F, row.names=T)

        ## add this to big data frame with all samples
        this_sample <- as.data.frame(t(cutoff_info_file), stringsAsFactors=F)
        this_sample[, "species"]  <- species
        this_sample[, "organism"] <- as.character(unique(annotation$organism[annotation$speciesId == species]))
        collectSamples <- rbind(collectSamples, this_sample)
  }
}
cat("libraryId\tcutoffTPM\tproportionGenicPresent\tnumberGenicPresent\tGenic\tproportionCodingPresent\tnumberCodingPresent\tProteinCoding\tproportion_intergenic_present\tintergenic_region_present\tIntergenic\tratioIntergenicCodingPresent\tspeciesID\torganism\n",
    file = file.path(output_folder, "All_samples.tsv"), sep = "\t")
write.table(collectSamples, file = file.path(output_folder, "All_samples.tsv"),col.names =F , row.names = F, append = T,quote = FALSE, sep = "\t")
plotData(file.path(output_folder, "All_samples.tsv"))
