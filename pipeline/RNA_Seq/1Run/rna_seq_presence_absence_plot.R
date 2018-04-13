## Julien Roux Jun 12, 2016
## This script plots summary statistics based on the present/absent calls for all genes
## IMPORTANT: this plotting script was merged to the end of rna_seq_presence_absence.pl, the current script is left for reference only and should be eventually removed

## Usage:
## R CMD BATCH --no-save --no-restore '--args rna_seq_sample_info="rna_seq_sample_info.txt" rna_seq_sample_excluded="rna_seq_sample_excluded.txt" presence_folder="presence_absence_bgee_v14"' rna_seq_presence_absence_plot.R rna_seq_presence_absence_plot.Rout
## rna_seq_sample_info      - file with info on mapped libraries
## rna_seq_sample_excluded  - file with excluded libraries
## presence_folder          - path to folder where presence_absence files are stored

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
command_arg <- c("rna_seq_sample_info", "rna_seq_sample_excluded", "presence_folder")
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
all_samples <- data.frame()

for(libraryId in sampleInfo$libraryId){
  ## Cutoff info file
  file <- paste0(presence_folder, "/", libraryId, "/cutoff_info_file.tsv")
  if (!(file.exists(file))){
    cat("  Missing data file for", libraryId, "\n")
  } else {
    cat("  Treating", libraryId, "\n")

    ## read file
    cutoff_info_file <- as.data.frame(t(read.table(file, h=F, sep="\t", row.names=1)))
    cutoff_info_file$species <- sampleInfo$speciesId[sampleInfo$libraryId == libraryId]
    cutoff_info_file$organism <- sampleInfo$organism[sampleInfo$libraryId == libraryId]

    ## add this to big data frame with all samples
    all_samples <- rbind(all_samples, cutoff_info_file)
  }
}
print(dim(all_samples))
print(length(unique(all_samples$libraryId)))
print(length(unique(all_samples$organism)))
save(all_samples, file=paste0(presence_folder, "/all_samples.RDa"))

## PDF for all boxplot
pdf(file = paste0(presence_folder, "/boxplots_4_methods_comparison.pdf"), width = 12, height = 5)
ylimits <- list(NA, NA, NA, NA, c(0, 100), c(0, 40000), c(0, 100), c(0, 20000), c(0, 100), c(0, 7000), c(0, 100))
i = 1
for (column in c(4:9, 11:12, 14:15)){
  print(names(all_samples)[column])
  method <- "method3"
  ## for (method in c("method1", "method2", "method3", "method4")){
    ## get the subset of data to use
    subsetNames <- names(all_samples)[column]
    subsetData <- as.numeric(as.character(all_samples[all_samples$method == method, column]))
    organism <- as.factor(as.character(all_samples$organism[all_samples$method == method]))
 
    ## reoder organism levels to get boxes ordered by median values in boxplot
    organism <- reorder(organism, subsetData, median)
  
    ## For each method, boxplot split by species of each column
    par(mar=c(12,4,2,1)) ## bottom, left, top and right margins
    if (is.na(ylimits[[i]][1])){
      boxplot(log2(subsetData+10^-6) ~ organism,
              pch=16, notch=T, outcex=0.5, boxwex=0.7,
              ylab=paste0("log2(", subsetNames, " + 10e-6)"),
              main=method, las=2)
      text(x=1:length(unique(organism)),
           y=log2(max(subsetData)+10^-6),
           labels=table(organism),
           col="black", cex=0.8)
    } else {
      boxplot(subsetData ~ organism,
              pch=16, notch=T, outcex=0.5, boxwex=0.7,
              ylab=subsetNames,
              main=method, las=2,
              ylim=ylimits[[i]])
      text(x=1:length(unique(organism)),
           y=ylimits[[i]][2],
           labels=table(organism),
           col="black", cex=0.8)
    }

    abline(v=1:length(unique(organism)), lty=2, col="lightgray")
    if (column == 11){
      abline(h=c(70, 80), lty=2, col="red")
    }
  ## }
  i = i+1
}   
## close PDF device        
dev.off()

## TODO - cross with annotation file to get tissue and redo same plot split by tissue?
##          look at testis / brain / whole body: more genes expressed?
##        - Automatically report those not included within 70-80%
##        - Look for samples with ratioIntergenicCodingPresent > 5%
##        - Check library with very overlapping distributions: SRX190962, Erinaceus europaeus
##        - Find library with completely non overlapping distribution. These should be fine as proportion of intergenic will go down quickly from 100% to 0%, and cutoff will be set here
##        - Split by tissue!

