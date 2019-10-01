## Julien Roux June 12, 2016
## This script computes cutoff for present/absent call and returns calls for all genes
## It loops through all libraries in rna_seq_sample_info.txt (removing those in rna_seq_sample_excluded.txt)
## Output is written in kallisto results folder

## Usage:
## R CMD BATCH --no-save --no-restore '--args rna_seq_sample_info="rna_seq_sample_info.txt" rna_seq_sample_excluded="rna_seq_sample_excluded.txt" kallisto_count_folder= "all_results_bgee_v15" sum_by_species_folder="$(RNASEQ_CLUSTER_SUM_RES)" gaussian_choice="$(RNASEQ_CLUSTER_GAUSSIAN_CHOICE)" out_folder="all_results_bgee_v15" desired_r_cutoff="r_cutoff_value" plot_only=FALSE' rna_seq_presence_absence.R rna_seq_presence_absence.Rout
## rna_seq_sample_info      - file with info on mapped libraries
## rna_seq_sample_excluded  - file with excluded libraries
## kallisto_count_folder    - path to kallisto result folder
## sum_by_species_folder    - path to folder where summed data where exported by rna_seq_sum_by_species.pl script.
## gaussian_choice          - path to file with manually chosen deconvoluted gaussians
## out_folder               - path to folder where to write results
## desired_r_cutoff         - desired cutoff value (proportion of intergenic, value between 0 and 1)
## plot_only                - specify if only plotting of final boxplots should be made (using .RDa file presence_absence_all_samples.RDa)

## Session info
print(sessionInfo())

## reading in arguments provided in command line
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed in command line
command_arg <- c("rna_seq_sample_info", "rna_seq_sample_excluded", "kallisto_count_folder", "sum_by_species_folder", "gaussian_choice", "out_folder", "desired_r_cutoff", "plot_only")
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

## File indicating gaussian_choice
if( file.exists(gaussian_choice) ){
  gaussian <- read.table(gaussian_choice, h=T, sep="\t", comment.char="")
} else {
  warning( paste("File with gaussian choices not found [", gaussian_choice, "]\n"))
}

###############################################################################
## Functions:
plot_distributions <- function(counts, selected_coding, selected_intergenic, cutoff){
  ## Plotting of the distribution of TPMs for coding and intergenic regions + cutoff
  ## Note: this code is largely similar to plotting section in rna_seq_analysis.R

  par(mar=c(5,6,1,1)) ## bottom, left, top and right margins
  dens <- density(log2(na.omit(counts$tpm) + 10^-6))

  ## protein-coding genes only (had to take care of NAs strange behavior)
  dens_coding <- density(log2(counts$tpm[selected_coding] + 10^-6))
  ## Normalize density for number of observations
  dens_coding$y <- dens_coding$y * sum(selected_coding) / length(counts$tpm)

  ## intergenic
  dens_intergenic <- density(log2(counts$tpm[selected_intergenic] + 10^-6))
  dens_intergenic$y <- dens_intergenic$y * sum(selected_intergenic) / length(counts$tpm)

  ## Plot whole distribution
  plot(dens, ylim=c(0, max(dens$y)*1.1), xlim=c(-23, 21), lwd=2, main=libraryId, bty="n", axes=F, xlab="")
  axis(2, las=1)
  ## Add 2 x-axes: TPMs and RPKMs
  ## See http://stackoverflow.com/questions/8443820/r-multiple-x-axis-with-annotations
  axis(1, at=seq(-30 , 30, by=10), line=0, mgp = c(3, 0.5, 0), cex.axis=0.8)
  mtext(expression(log[2]('TPM'+10^-6)), 1,  adj = 1, padj = 0, line=0.2, at=par("usr")[1], col="black", cex=0.8)
  ## To make FPKM scale, we need to know what log2(TPM + 10^-6) value corresponds to any log2(FPKM + 10^-6) value. We know that for any given gene, TPMg/FPKMg = coef
  ##    log2(FPKM + 10^-6) = x
  ## <=>              FPKM = exp(x*log(2)) - 10^-6
  ## <=>               TPM = coef * (exp(x*log(2)) - 10^-6)
  ## <=> log2(TPM + 10^-6) = log2( coef * (exp(x*log(2)) - 10^-6) + 10^-6)
  coef <- na.omit(counts$tpm / counts$fpkm)[1]
  ## We generate scale from -20 to 100 FPKMs
  axis(1, at=log2( coef * (exp(seq(-20 , 100, by=10)*log(2)) - 10^-6) + 10^-6), labels=seq(-20 , 100, by=10), line=2, mgp = c(3, 0.5, 0), cex.axis=0.8)
  mtext(expression(log[2]('FPKM'+10^-6)), 1,  adj = 1, padj = 0, line=2.2, at=par("usr")[1], col="black", cex=0.8)

  ## Plot the TPM cutoff
  ## abline(v=cutoff, col="gray", lty=1, lwd=2)
  arrows(log2(cutoff + 10e-6), par("usr")[3], log2(cutoff + 10e-6), par("usr")[4]/2, col="gray", lty=1, lwd=2, angle=160, length=0.1)

  ## Add subgroups distributions (coding, intergenic, etc):
  ## protein-coding genes
  lines(dens_coding, col="firebrick3", lwd=2, lty=2)
  ## intergenic
  lines(dens_intergenic, col="dodgerblue3", lwd=2, lty=2)

  ## legend
  legend("topright", c("all", "selected protein-coding genes", "selected intergenic regions"), lwd=2, col=c("black", "firebrick3", "dodgerblue3"), lty=c(1, 2, 2), bty="n")
  return();
}

cutoff_info <- function(library_id, counts, column, max_intergenic, TPM_cutoff, TPM_final_cutoff, FPKM_cutoff, FPKM_final_cutoff, r_cutoff){
  ## Calculate summary statistics to export in cutoff info file
  genic_present <- sum(counts[[column]][counts$type == "genic"] == "present")/sum(counts$type == "genic") * 100
  number_genic_present <- sum(counts[[column]][counts$type == "genic"] == "present")

  coding_present <- sum(counts[[column]][counts$biotype %in% "protein_coding"] == "present")/sum(counts$biotype %in% "protein_coding") * 100
  number_coding_present <- sum(counts[[column]][counts$biotype %in% "protein_coding"] == "present")

  intergenic_present <- sum(counts[[column]][counts$type == "intergenic"] == "present")/sum(counts$type == "intergenic") * 100
  number_intergenic_present <- sum(counts[[column]][counts$type == "intergenic"] == "present")

  ## Export cutoff_info_file
  to_export <- c(library_id,
                 max_intergenic,
                 TPM_cutoff, TPM_final_cutoff, FPKM_cutoff, FPKM_final_cutoff,
                 genic_present, number_genic_present, sum(counts$type == "genic"),
                 coding_present,  number_coding_present, sum(counts$biotype %in% "protein_coding"),
                 intergenic_present, number_intergenic_present, sum(counts$type == "intergenic"),
                 r_cutoff
                 )
  names(to_export) <- c("libraryId",
                        "max_intergenic",
                        "cutoffTPM", "cutoffGenicTPM", "cutoffFPKM", "cutoffGenicFPKM",
                        "proportionGenicPresent", "numberGenicPresent", "numberGenic",
                        "proportionCodingPresent", "numberPresentCoding", "numberCoding",
                        "proportionIntergenicPresent", "numberIntergenicPresent", "numberIntergenic",
                        "ratioIntergenicCodingPresent")
  return(to_export)
}

calculate_and_plot_r <- function(counts, selected_coding, selected_intergenic, desired_r_cutoff){
  ## r = (number of intergenic regions with TPM values higher than x * number of coding regions) /
  ##     (number of coding regions with TPM values higher than x * number of intergenic regions)
  ##   = 0.05
  ## What is value of x (cutoff)? calculate the distribution of r for a range of TPMs, then select the closest value to 0.05

  ## Counting how many intergenic regions have equal or higher value of TPM for every value of TPM
  ## For each gene's TPM (sorted), calculate r

  ##   r <- sapply(sort(unique(counts$tpm[selected_coding])), function(x){
  ##     return(
  ##            ( sum(counts$tpm[selected_intergenic] >= x) / sum(selected_intergenic) ) /
  ##            ( sum(counts$tpm[selected_coding] >= x) / sum(selected_coding) )
  ##            )
  ##   })
  ##   ## This is too long! Takes >30s for each library

  ## For each unique value of coding TPM, let's calculate the number of intergenic TPMs that are larger:
  ## ptm <- proc.time()
  summed_intergenic <- sapply(unique(sort(counts$tpm[selected_coding])), function(x){
    return( sum(counts$tpm[selected_intergenic] >= x) )
  })
  ## It is not necessary to do the same for coding regions: for a sorted vector the number of greater elements is equal to lenght(vector) - position + 1. Here, it is a bit trickier since we do not consider all coding TPM values, but the unique ones, so we use the rle function to know the lengths of the runs of similar values and sum them
  summed_coding <- c(0, cumsum(rle(sort(counts$tpm[selected_coding]))$lengths))
  summed_coding <- summed_coding[-(length(summed_coding))]
  summed_coding <- sum(selected_coding) - summed_coding

  ## Now we can calculate r
  r <- ( summed_intergenic / sum(selected_intergenic) ) /
       ( summed_coding / sum(selected_coding) )
  ## This is twice faster as code above!

  percent <- (1-desired_r_cutoff)*100

  ## Select the minimal value of TPM for which the ratio of genes and intergenic regions is equal to 0.05 or lower (first test if at least 1 TPM value has this property):
  if (sum(r < desired_r_cutoff) == 0){
    TPM_cutoff <- sort(unique(counts$tpm[selected_coding]))[which(r == min(r))[1]]
    r_cutoff <- min(r)
    cat(paste0("    There is no TPM cutoff for which " , percent,"%", " of the expressed genes would be coding. TPM cutoff is fixed at the first value with maximum coding/intergenic ratio. r=", r_cutoff, " at TPM=", TPM_cutoff,"\n"))
  } else {
    TPM_cutoff <- sort(unique(counts$tpm[selected_coding]))[which(r < desired_r_cutoff)[1]]
    r_cutoff <- desired_r_cutoff
    cat(paste0("    TPM cutoff for which " , percent,"%", " of the expressed genes are be coding found at TPM=", TPM_cutoff,"\n"))
  }

  ## Plot TPMs vs. ratio: should mostly go down
  plot(log2(sort(unique(counts$tpm[selected_coding]))+10e-6), r, pch=16, xlab="log2(TPM + 10^-6)", type="l")
  abline(h=desired_r_cutoff, lty=2, col="gray")
  if (r_cutoff > desired_r_cutoff){
    abline(h=desired_r_cutoff, lty=3, col="gray")
  }
  arrows(log2(TPM_cutoff + 10e-6), par("usr")[3], log2(TPM_cutoff + 10e-6), par("usr")[4]/2, col="gray", lty=1, lwd=2, angle=160, length=0.1)
  return(c(TPM_cutoff, r_cutoff))
}

###############################################################################
## data frame to store info for all samples
all_samples <- data.frame()

## Calculate presence / absence thresholds (only if plot_only=FALSE)
if ( !plot_only ){

  ## Loop through all species, with human last because they are so many samples...
  for(species in c(unique(sampleInfo$speciesId[sampleInfo$speciesId != 9606]), 9606)){
    ## Sum by species data
    file <- paste0(sum_by_species_folder, "/sum_abundance_gene_level+fpkm+intergenic+classification_", species, ".tsv")
    if (!(file.exists(file))){
      cat("Sum by species data file is missing for", species, "\n")
    } else {
      organism <- as.character(unique(sampleInfo$organism[sampleInfo$speciesId == species]))
      cat(paste0("\n\nTreating data for ", organism, " (species ID: ", species,")\n"))
      ## Read sum by species data
      sum_by_species <- read.table(file, h=T, sep="\t")

      for(libraryId in sampleInfo$libraryId[sampleInfo$speciesId == species]){
        ## Gene level file including intergenic regions
        file1 <- paste0(kallisto_count_folder, "/", libraryId, "/abundance_gene_level+fpkm+intergenic.tsv")
        ## Gene level file without intergenic regions
        file2 <- paste0(kallisto_count_folder, "/", libraryId, "/abundance_gene_level+new_tpm+new_fpkm.tsv")
        if (!(file.exists(file1) & file.exists(file2))){
          cat("  Missing data file for", libraryId, "\n")
        } else {
          cat("  Treating", libraryId, "\n")

          ## If out folder doesn't exist, create it
          if (!file.exists(file.path(out_folder, libraryId))){
            dir.create(file.path(out_folder, libraryId), mode = "0755")
          }

          ## read files
          kallisto_gene_counts <- read.table(file1, h=T, sep="\t")
          gene_counts <- read.table(file2, h=T, sep = "\t")

          #######################################################################################
          ## This section calculates the TPM cutoff used to decide if a gene is present or absent
          ## We use Marta's approach, using all coding regions and selected deconvoluted intergenic regions

          ## Plot the density of selected coding and intergenic regions, with cutoff
          pdf(file = paste0(out_folder, "/", libraryId, "/distribution_TPM_genic_intergenic+cutoff.pdf"), width = 6, height = 5)

          ## select all coding regions
          selected_coding <- kallisto_gene_counts$biotype %in% "protein_coding"

          ## For intergenic regions, select all regions that (in summed data) have TPM below threshold defined by selected gaussian (max TPM of regions from this gaussian, the gaussian is on the left of the threshold)
          if (gaussian$selectionSideIntergenic[gaussian$speciesId == species] == "Left"){
            max_intergenic <- max(sum_by_species$tpm[
                                                     sum_by_species$classification %in%
                                                     paste0("intergenic_", gaussian$selectedGaussianIntergenic[gaussian$speciesId == species])
                                                     ])
          } else if (gaussian$selectionSideIntergenic[gaussian$speciesId == species] == "Right") {
            ## In some case it is tricky to define the gaussians in the conventional way, so it is possible to select all regions that (in summed data) have TPM below threshold defined by selected gaussian, BUT taking min TPM of regions from this gaussian, the gaussian is on the right of the threshold!
            max_intergenic <- min(sum_by_species$tpm[
                                                     sum_by_species$classification %in%
                                                     paste0("intergenic_", gaussian$selectedGaussianIntergenic[gaussian$speciesId == species])
                                                     ])
          }
          ## TODO ideally this could be done only once per species

          selected_intergenic <- (kallisto_gene_counts$type == "intergenic"
                                  & sum_by_species$tpm <= max_intergenic)

          ## Print what the max level of summed TPM for intergenic regions kept is
          cat(paste0("    The max (summed by species) TPM level of selected intergenic is ", max_intergenic,"\n"))

          ## Function to calculate TPM cutoff at r of 5% or the minimal r across all TPM range
          results <- calculate_and_plot_r(kallisto_gene_counts, selected_coding, selected_intergenic, as.numeric(desired_r_cutoff))
          TPM_cutoff <- results[1]
          r_cutoff <- results[2]

          ## FPKM are proportional to TPM. The proportionality coefficient is taken from the first non-zero gene
          FPKM_cutoff <-  TPM_cutoff * na.omit(kallisto_gene_counts$fpkm / kallisto_gene_counts$tpm)[1]
          ## From FPKMs calculated over genic and intergenic regions to FPKMs calculated over genic regions only
          FPKM_final_cutoff <- FPKM_cutoff * sum(kallisto_gene_counts$est_counts)/sum(gene_counts$est_counts)
          ## Back to TPM
          TPM_final_cutoff <- FPKM_final_cutoff * na.omit(gene_counts$tpm / gene_counts$fpkm)[1]

          ## Plot distributions and cutoff
          plot_distributions(kallisto_gene_counts, selected_coding, selected_intergenic, TPM_cutoff)

          ## Setting the presence / absence flags
          gene_counts$call <- ifelse(gene_counts$tpm >= TPM_final_cutoff, "present", "absent")
          kallisto_gene_counts$call <- ifelse(kallisto_gene_counts$tpm >= TPM_cutoff, "present", "absent")
          ## Export cutoff info file
          cutoff_info_file <- cutoff_info(libraryId, kallisto_gene_counts, "call", max_intergenic, TPM_cutoff, TPM_final_cutoff, FPKM_cutoff, FPKM_final_cutoff, r_cutoff)

          ## close PDF device
          dev.off()

          ## Saving calls and stats to output folder:
          write.table(gene_counts,
                      file = paste0(out_folder, "/", libraryId, "/abundance_gene_level+new_tpm+new_fpkm+calls.tsv"),
                      ## This is the final file to be used to populate Bgee with TPM values and calls
                      quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
          write.table(kallisto_gene_counts,
                      file = paste0(out_folder, "/", libraryId, "/abundance_gene_level+fpkm+intergenic+calls.tsv"),
                      quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
          write.table(t(t(cutoff_info_file)),
                      file = paste0(out_folder, "/", libraryId, "/cutoff_info_file.tsv"),
                      quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)
                      ## t(t(cutoff_info_file)) is a solution to export a vector vertically

          ## add this to big data frame with all samples
          this_sample <- as.data.frame(t(cutoff_info_file), stringsAsFactors=F)
          this_sample[, "species"]  <- species
          this_sample[, "organism"] <- organism
          all_samples <- rbind(all_samples, this_sample)
        }
      }
    }
  }
  ## backup of data for plotting
  save(all_samples, file=paste0(out_folder, "/presence_absence_all_samples.RDa"))
  ## get library number here because after that "libraryId" field name becomes "#libraryId"
  libraryIdNbr <- length(unique(all_samples$libraryId))
  ## Add a # at beginning of header line
  names(all_samples)[1] <- paste0('#', names(all_samples)[1])
  ## export tsv file
  write.table(all_samples, file=paste0(out_folder, "/presence_absence_all_samples.txt"), sep="\t", quote=F, col.names=TRUE, row.names=FALSE)
  ## TODO give name of the output file as argument of this script

  cat(paste0("Done. ", libraryIdNbr, " libraries successfully treated in ", length(unique(all_samples$organism)), " species.\n"))
}

cat("Now plotting the presence / absence parameters and results\n")
if (plot_only){
  load(file=paste0(out_folder, "/presence_absence_all_samples.RDa"))
}

## PDF for all boxplot
pdf(file = paste0(out_folder, "/presence_absence_boxplots.pdf"), width = 12, height = 5)
## 14 boxplots to plot, 1 per column in all_samples
## Manually set the y-scale limits for most of the columns. NA are here for columsn thta need to be printed in log scale
ylimits <- list(NA, NA, NA, NA, c(0, 100), c(0, 70000), c(0, 70000), c(0, 100), c(0, 30000), c(0, 30000), c(0, 100), c(0, 30000), c(0, 30000), c(0, 0.2))

i = 1
for (column in c(3:16)){
  print(names(all_samples)[column])

  ## get the subset of data to use
  subsetNames <- names(all_samples)[column]
  subsetData <- as.numeric(as.character(all_samples[, column]))
  organism <- as.factor(as.character(all_samples$organism))

  ## reoder organism levels to get boxes ordered by median values in boxplot
  organism <- reorder(organism, subsetData, median)

  par(mar=c(12,4,2,1)) ## bottom, left, top and right margins
  if (is.na(ylimits[[i]][1])){
    boxplot(log2(subsetData+10^-6) ~ organism,
            pch=16, notch=T, outcex=0.5, boxwex=0.7,
            ylab=paste0("log2(", subsetNames, " + 10e-6)"),
            main="", las=2)
    text(x=1:length(unique(organism)),
         y=log2(max(subsetData)+10^-6),
         labels=table(organism),
         col="black", cex=0.8)
  } else {
    boxplot(subsetData ~ organism,
            pch=16, notch=T, outcex=0.5, boxwex=0.7,
            ylab=subsetNames,
            main="", las=2,
            ylim=ylimits[[i]])
    text(x=1:length(unique(organism)),
         y=ylimits[[i]][2],
         labels=table(organism),
         col="black", cex=0.8)
  }
  abline(v=1:length(unique(organism)), lty=2, col="lightgray")
  if (names(all_samples)[column] == "proportionCodingPresent"){
    abline(h=c(70, 80), lty=2, col="red")
  }
  i = i+1
}
## close PDF device
dev.off()
## TODO? Plotting could be made in a bit more elegant way, but this should work fine

## TODO? add the possibility to remove regions above / below a arbitrary TPM threshold in the summed data? This could be useful when the distributions deconvoluted clearly do not correspond to what seems to make sense

