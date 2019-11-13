## Julien Roux Mar 23, 2016
## This script loads expression at gene level for all libraries of each species, and sums the signal. This alows to better deconvolute the "real" intergenic regions, and "good" protein-coding genes to be used for presence / absence expression calls (see rna_seq_presence_absence.R)

## SFonsecaCosta Sep 17, 2019
## Re-calculate TPM and FPKM in sum_by_species!
## Instead of sum all est_counts, TPM and FPKM across all libraries from the same specie at gene level, here we sum the est_counts of all libraries at transcript level
## Followed be the recalculation of eff_length by weighted.mean, and then the TPM and FPKM are recalculated at transcript level 
## In the end the information is provided again at gene level (after sum the est_counts, TPM and FPKM from transcriptID that belong to the same geneID).


## Usage:
## R CMD BATCH --no-save --no-restore '--args rna_seq_sample_info="rna_seq_sample_info.txt" rna_seq_sample_excluded="rna_seq_sample_excluded.txt" kallisto_count_folder="all_results_bgee_v14" sum_by_species_folder="$(RNASEQ_VITALIT_SUM_RES)"' rna_seq_sum_by_species.R rna_seq_sum_by_species.Rout
## rna_seq_sample_info     - file with info on mapped libraries
## rna_seq_sample_excluded - file with excluded libraries
## kallisto_count_folder   - path to kallisto result folder
## sum_by_species_folder   - folder where to export the plots, summed data, and classification of coding / intergenic regions

## Session info
print(sessionInfo())

## libraries
library(dplyr)

## reading in arguments provided in command line
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed in command line
command_arg <- c("rna_seq_sample_info", "rna_seq_sample_excluded", "kallisto_count_folder", "sum_by_species_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read rna_seq_sample_info.txt file. If file not exists, script stops
if( file.exists(rna_seq_sample_info) ){
  sampleInfo <- read.table(rna_seq_sample_info, h=T, sep="\t", comment.char="")
  names(sampleInfo)[1] <- "libraryId"
} else {
  stop( paste("rna_seq_sample_info.txt file not found [", rna_seq_sample_info, "]\n"))
}
print(dim(sampleInfo))

## Remove excluded libraries from sample info file (mapping failed, or bad quality)
if( file.exists(rna_seq_sample_excluded) ){
  sampleExcluded <- read.table(rna_seq_sample_excluded, h=T, sep="\t", comment.char="")
  names(sampleExcluded)[1] <- "libraryId"
  sampleInfo <- sampleInfo[ !sampleInfo$libraryId %in% sampleExcluded$libraryId[sampleExcluded$excluded == TRUE], ]
} else {
  warning( paste("rna_seq_sample_excluded.txt file not found [", rna_seq_sample_excluded, "]\n"))
}
print(dim(sampleInfo))

###############################################################################
## Loop across species and sum the expression across libraries

## write header of file number_libraries.txt
cat("speciesId\tspeciesName\tnumberLibrariesUsed\tnumberLibraries\n", file = paste0(sum_by_species_folder, "/number_libraries.txt"), sep = "\t")
## write header of file gaussian_choice_by_species.txt
cat("speciesId\torganism\tnumberGaussiansCoding\tnumberGaussiansIntergenic\tselectedGaussianCoding\tselectionSideCoding\tselectedGaussianIntergenic\tselectionSideIntergenic\tcomment\tannotatorId\n", file = paste0(sum_by_species_folder, "/gaussian_choice_by_species_TO_FILL.txt"), sep = "\t")

for(species in unique(sampleInfo$speciesId)){
  cat(paste0("\nSumming data for ", as.character(unique(sampleInfo$organism[sampleInfo$speciesId == species])), " (species ID: ", species,")\n"))
  numLibs = 0
  
  raw_counts_info <- c()
  effec_length_info <- c()
  
  for(libraryId in sampleInfo$libraryId[sampleInfo$speciesId == species]){
    ## For each library use the file at transcriptID level
    file <- paste0(kallisto_count_folder, "/", libraryId, "/abundance+gene_id+fpkm+intergenic.tsv")
    
    if (file.exists(file)){
      cat("  Reading the ", libraryId, "\n")
      numLibs = numLibs + 1
      ## read transcript level data for each library
      kallisto_transcrip_counts <- read.table(file, h=T, sep="\t")
      raw_counts_info <- cbind(raw_counts_info, kallisto_transcrip_counts$est_counts)
      effec_length_info <- cbind(effec_length_info, kallisto_transcrip_counts$eff_length)
    } else {
      cat("  No data found from", libraryId, "\n")
    }
  }
  
  if (numLibs == 1){
    ## provide output directly from gene level
    file2 <- paste0(kallisto_count_folder,"/", libraryId, "/abundance_gene_level+fpkm+intergenic.tsv")
    kallisto_transcrip_counts <- read.table(file2, h=T, sep="\t")
    summed <- kallisto_transcrip_counts
  } else {
    summed <- c()
    ## sum the raw counts for all libraries
    summed$counts <- apply(raw_counts_info, 1, sum)
    ## calculate effect_length for each transcript by using weighted.mean...
    calculate_effect_length <- function(counts, effec_length){
      myWeightedMean <- c()
      ## transpose the matrix (this means each column is a transcript ID and each row is a library)
      rawcounts <- t(counts)
      raw_effeclength <- t(effec_length)
      
      for (i in 1:ncol(raw_effeclength)) {
        if (sum(rawcounts[,i]) == 0 ){
          ## provide the same weight for a transcript in case the est_count is always zero
          rawcounts[,i] <- ifelse(rawcounts[,i] == 0, 1)
          weightedMean <-  weighted.mean(x=raw_effeclength[,i], w=rawcounts[,i])
          myWeightedMean[i] <- weightedMean
        } else {
          weightedMean <-  weighted.mean(x=raw_effeclength[,i], w=rawcounts[,i])
          myWeightedMean[i] <- weightedMean
        }       
      }
      return(myWeightedMean)
    }
    ## each row is the weighted.mean of eff_length for each transcriptID
    summed$eff_length <- as.data.frame(calculate_effect_length(raw_counts_info, effec_length_info))
    
    ## re-calculate TPM and FPKM after collect the weighted.mean of eff_length for each transcriptID and after sum the est_count for each transcriptID
    estCount_to_tpm <- function(est_count, effec_length){
      rate <- log(est_count) - log(effec_length)
      denom <- log(sum(exp(rate)))
      exp(rate - denom + log(1e6))
    }
    summed$tpm <- estCount_to_tpm(summed$counts, summed$eff_length) 
    
    estCount_to_fpkm <- function(est_count, effec_length){
      N <- sum(est_count)
      exp( log(est_count) + log(1e9) - log(effec_length) - log(N) )
    }
    summed$fpkm <- estCount_to_fpkm(summed$counts, summed$eff_length) 
    
    summed <- data.frame(kallisto_transcrip_counts$target_id, kallisto_transcrip_counts$gene_id, kallisto_transcrip_counts$length, summed$eff_length, summed$counts, summed$tpm, summed$fpkm, kallisto_transcrip_counts$type, kallisto_transcrip_counts$biotype)
    colnames(summed) <- c("target_id", "gene_id", "length", "eff_length", "est_counts","tpm", "fpkm", "type", "biotype")
    
    ## select intergenic regions
    intergenic_regions <- summed[summed$type == "intergenic", c(1, 5, 6, 7, 8, 9)]
    names(intergenic_regions)[1] <- "gene_id"
    ## Sum est_counts, TPM and FPKM for genic regions
    sumGenic <- summed[summed$type == "genic", c(2, 5, 6, 7, 8, 9)]
    sumGenic <- aggregate(sumGenic[,2:4], list(sumGenic$gene_id), sum)
    names(sumGenic)[1] <- "gene_id"
    sumGenic$type <-  rep("genic", times=length(sumGenic$gene_id))
    select_genic <- dplyr::filter(summed, type == "genic")
    select_biotype_genic = distinct(select_genic, gene_id, .keep_all= TRUE)
    sumGenic$biotype <- select_biotype_genic$biotype
    ## Final Table with genic and intergenic information
    summed <- rbind(sumGenic, intergenic_regions)
 }
 
  ## Export number of libraries (and total number of libraries for this species) used in this species
  cat(c(species, as.character(unique(sampleInfo$organism[sampleInfo$speciesId == species])), numLibs, paste0(length(sampleInfo$libraryId[sampleInfo$speciesId == species]), "\n")), file = paste0(sum_by_species_folder, "/number_libraries.txt"), sep = "\t", append = TRUE)
  
  ## if no library ws found for this species
  if (numLibs == 0){
    cat("  No library found for this species, skipping it.")
    next
  }
  
  ## Else:
  cat("  Plotting density of aggregated data\n")
  ## Density plot of summed data
  pdf(file = paste0(sum_by_species_folder, "/distribution_TPM_genic_intergenic_sum_", species, ".pdf"), width = 6, height = 5)
  ## par(mar=c(5,6,1,1)) ## bottom, left, top and right margins
  ## density of log2(TPM) of summed data
  dens <- density(log2(na.omit(summed$tpm) + 10^-6))
  ## Subgroups densities. Visualization trick: we add an invisible set of points at x=-30, to make densities comparable
  ## genic regions
  dens_genic <- density(c(rep(-30, times=sum(summed$type != "genic")), log2(summed$tpm[summed$type == "genic"] + 10^-6)))
  ## protein-coding genes only (had to take care of NAs strange behavior)
  dens_coding <- density(c(rep(-30, times=sum(!summed$biotype %in% "protein_coding")), log2(summed$tpm[summed$biotype %in% "protein_coding"] + 10^-6)))
  ## intergenic
  dens_intergenic <- density(c(rep(-30, times=sum(summed$type != "intergenic")), log2(summed$tpm[summed$type == "intergenic"] + 10^-6)))
  ## Plot whole distribution
  plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-23, 21), lwd=2, main=paste0(as.character(unique(sampleInfo$organism[sampleInfo$speciesId == species])), " (", numLibs, " libraries)"), bty="n", axes=T, xlab="log2(TPM + 10^-6)")
  ## Add subgroups distributions (genic, intergenic, etc):
  ## genic
  lines(dens_genic, col="firebrick3", lwd=2)
  ## protein-coding genes
  lines(dens_coding, col="firebrick3", lwd=2, lty=2)
  ## intergenic
  lines(dens_intergenic, col="dodgerblue3", lwd=2)
  ## legend
  legend("topleft", c(paste0("all (", length(summed[,1]),")"), paste0("genic (", sum(summed$type == "genic"), ")"), paste0("coding (", sum(summed$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "firebrick3", "dodgerblue3"), lty=c(1, 1, 2, 1), bty="n")
  dev.off()
  
  ## Redo plot for log2(FPKMs) (probably not proportional anymore)
  pdf(file = paste0(sum_by_species_folder, "/distribution_FPKM_genic_intergenic_sum_", species, ".pdf"), width = 6, height = 5)
  ## par(mar=c(5,6,1,1)) ## bottom, left, top and right margins
  ## FPKM density of summed data
  dens <- density(log2(na.omit(summed$fpkm) + 10^-6))
  ## Subgroups densities. Visualization trick: we add an invisible set of points at x=-30, to make densities comparable
  ## genic regions
  dens_genic <- density(c(rep(-30, times=sum(summed$type != "genic")), log2(summed$fpkm[summed$type == "genic"] + 10^-6)))
  ## protein-coding genes only (had to take care of NAs strange behavior)
  dens_coding <- density(c(rep(-30, times=sum(!summed$biotype %in% "protein_coding")), log2(summed$fpkm[summed$biotype %in% "protein_coding"] + 10^-6)))
  ## intergenic
  dens_intergenic <- density(c(rep(-30, times=sum(summed$type != "intergenic")), log2(summed$fpkm[summed$type == "intergenic"] + 10^-6)))
  ## Plot whole distribution
  plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-23, 21), lwd=2, main=paste0(as.character(unique(sampleInfo$organism[sampleInfo$speciesId == species])), " (", numLibs, " libraries)"), bty="n", axes=T, xlab="log2(FPKM + 10^-6)")
  ## Add subgroups distributions (genic, intergenic, etc):
  ## genic
  lines(dens_genic, col="firebrick3", lwd=2)
  ## protein-coding genes
  lines(dens_coding, col="firebrick3", lwd=2, lty=2)
  ## intergenic
  lines(dens_intergenic, col="dodgerblue3", lwd=2)
  ## legend
  legend("topleft", c(paste0("all (", length(summed[,1]),")"), paste0("genic (", sum(summed$type == "genic"), ")"), paste0("coding (", sum(summed$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "firebrick3", "dodgerblue3"), lty=c(1, 1, 2, 1), bty="n")
  dev.off()
  
  ## Redo plot for log(read counts)
  pdf(file = paste0(sum_by_species_folder, "/distribution_counts_genic_intergenic_sum_", species, ".pdf"), width = 6, height = 5)
  ## par(mar=c(5,6,1,1)) ## bottom, left, top and right margins
  ## density of read counts on summed data
  dens <- density(log2(na.omit(summed$est_counts) + 1))
  ## Subgroups densities. Visualization trick: we add an invisible set of points at x=-10, to make densities comparable
  ## genic regions
  dens_genic <- density(c(rep(-10, times=sum(summed$type != "genic")), log2(summed$est_counts[summed$type == "genic"] + 1)))
  ## protein-coding genes only (had to take care of NAs strange behavior)
  dens_coding <- density(c(rep(-10, times=sum(!summed$biotype %in% "protein_coding")), log2(summed$est_counts[summed$biotype %in% "protein_coding"] + 1)))
  ## intergenic
  dens_intergenic <- density(c(rep(-10, times=sum(summed$type != "intergenic")), log2(summed$est_counts[summed$type == "intergenic"] + 1)))
  ## Plot whole distribution
  plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-1, 20), lwd=2, main=paste0(as.character(unique(sampleInfo$organism[sampleInfo$speciesId == species])), " (", numLibs, " libraries)"), bty="n", axes=T, xlab="log2(read counts + 1)")
  ## Add subgroups distributions (genic, intergenic, etc):
  ## genic
  lines(dens_genic, col="firebrick3", lwd=2)
  ## protein-coding genes
  lines(dens_coding, col="firebrick3", lwd=2, lty=2)
  ## intergenic
  lines(dens_intergenic, col="dodgerblue3", lwd=2)
  ## legend
  legend("topleft", c(paste0("all (", length(summed[,1]),")"), paste0("genic (", sum(summed$type == "genic"), ")"), paste0("coding (", sum(summed$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "firebrick3", "dodgerblue3"), lty=c(1, 1, 2, 1), bty="n")
  dev.off()
  
  
  ## Deconvolute TPM intergenic and genic distributions
  ## As in Hebenstreit 2011 Mol Syst Biol: use clustering approach
  ## Mclust: Normal Mixture Modelling for Model-Based Clustering, Classification, and Density Estimation
  ## We do not chose the number of gaussians, and let mclust choose
  cat("  Deconvoluting sub-distributions of genic and intergenic regions\n")
  library(mclust)
  
  ## Focus on regions with enough signal (remove TPM = 0 or very small)
  summed_filtered <- summed[summed$tpm > 10^-6, ]
  
  ## open PDF device
  pdf(file = paste0(sum_by_species_folder, "/distribution_TPM_genic_intergenic_sum_deconvolution_", species, ".pdf"), width = 6, height = 5)
  
  ## Coding regions
  mod1 = densityMclust(log2(summed_filtered$tpm[summed_filtered$biotype %in% "protein_coding"]))
  plot(mod1, what = "BIC")
  cat("    Protein-coding genes:\n")
  print(summary(mod1, parameters = TRUE))
  plot(mod1, what = "density", data = log2(summed_filtered$tpm[summed_filtered$biotype %in% "protein_coding"]), breaks = 100, xlab="log2(TPM) - protein-coding genes")
  
  ## Intergenic regions
  mod2 = densityMclust(log2(summed_filtered$tpm[summed_filtered$type == "intergenic"]))
  plot(mod2, what = "BIC")
  cat("    Intergenic regions:\n")
  print(summary(mod2, parameters = TRUE))
  plot(mod2, what = "density", data = log2(summed_filtered$tpm[summed_filtered$type == "intergenic"]), breaks = 100, xlab="log2(TPM) - intergenic")
  
  ## Plot the density of the original data, and the density of regions classified to different gaussians
  cat("  Plotting density of deconvoluted genic and intergenic regions\n")
  dens <- density(log2(summed_filtered$tpm))
  ## Plot whole distribution
  plot(dens, ylim=c(0, max(dens$y)*1.1), xlim=c(-20, 20), lwd=2, main=paste0(as.character(unique(sampleInfo$organism[sampleInfo$speciesId == species])), " (", numLibs, " libraries)"), bty="n", axes=T, xlab="log2(TPM)")
  
  ## protein-coding genes only (had to take care of NAs strange behavior)
  dens_coding <- density(log2(summed_filtered$tpm[summed_filtered$biotype %in% "protein_coding"]))
  ## Normalize density for number of observations
  dens_coding$y <- dens_coding$y * sum(summed_filtered$biotype %in% "protein_coding") / length(summed_filtered$tpm)
  lines(dens_coding, col="firebrick3", lwd=2, lty=2)
  ## Sub-distributions
  for (i in 1:mod1$G){
    ## if any point classified
    if (sum(mod1$classification == i) >= 2){
      dens_coding_sub <- density(log2(summed_filtered$tpm[summed_filtered$biotype %in% "protein_coding"][mod1$classification == i]))
      ## y-axis scaling
      dens_coding_sub$y <- dens_coding_sub$y * length(summed_filtered$tpm[summed_filtered$biotype %in% "protein_coding"][mod1$classification == i]) / length(summed_filtered$tpm)
      lines(dens_coding_sub, col=paste0("grey", trunc(100/(mod1$G+1))*i), lwd=2, lty=2)
      ## Print gaussian number on plot: at location of max value of gaussian (italics)
      text(dens_coding_sub$x[dens_coding_sub$y == max(dens_coding_sub$y)], 0.005, labels = i, col=paste0("grey", trunc(100/(mod1$G+1))*i), font=3)
    }
  }
  ## intergenic
  dens_intergenic <- density(log2(summed_filtered$tpm[summed_filtered$type == "intergenic"]))
  dens_intergenic$y <- dens_intergenic$y * sum(summed_filtered$type == "intergenic") / length(summed_filtered$tpm)
  lines(dens_intergenic, col="dodgerblue3", lwd=2)
  for (i in 1:mod2$G){
    ## if any point classified
    if (sum(mod2$classification == i) >= 2){
      dens_intergenic_sub <- density(log2(summed_filtered$tpm[summed_filtered$type == "intergenic"][mod2$classification == i]))
      ## y-axis scaling
      dens_intergenic_sub$y <- dens_intergenic_sub$y * length(summed_filtered$tpm[summed_filtered$type == "intergenic"][mod2$classification == i]) / length(summed_filtered$tpm)
      lines(dens_intergenic_sub, col=paste0("grey", trunc(100/(mod2$G+1))*i), lwd=2)
      ## Print gaussian number on plot: at location of max value of gaussian
      text(dens_intergenic_sub$x[dens_intergenic_sub$y == max(dens_intergenic_sub$y)], 0.005, labels = i, col=paste0("grey", trunc(100/(mod2$G+1))*i))
    }
  }
  ## legend
  legend("topleft", c(paste0("all (", length(summed_filtered[,1]),")"), paste0("coding (", sum(summed_filtered$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed_filtered$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "dodgerblue3"), lty=c(1, 2, 1), bty="n")
  dev.off()
  
  ## Export file with summed data and classification of intergenic and coding regions
  cat("  Exporting aggregated data and classification of coding and intergenic regions\n")
  ## Add new column to summed object
  summed$classification <- NA
  summed$classification[summed$tpm > 10^-6 & summed$biotype %in% "protein_coding"] <- paste("coding_", mod1$classification, sep="")
  summed$classification[summed$tpm > 10^-6 & summed$type == "intergenic"] <- paste("intergenic_", mod2$classification, sep="")
  write.table(summed, file = paste0(sum_by_species_folder, "/sum_abundance_gene_level+fpkm+intergenic+classification_", species, ".tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ## Export file with speciesId and number of coding and intergenic gaussians (to be filled manually with selected gaussians)
  cat(paste0(species, "\t", mod1$G, "\t", mod2$G, "\n"), file = paste0(sum_by_species_folder, "/gaussian_choice_by_species_TO_FILL.txt"), sep = "\t", append=T)
  
  rm(summed)
  rm(summed_filtered)
  rm(numLibs)
}

## TODO export gaussian parameters for each species
