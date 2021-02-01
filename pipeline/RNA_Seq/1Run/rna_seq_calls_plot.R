## Julien Wollbrett July 15, 2020
## This script plot presence/absence information for all RNA-Seq libraries grouped by species
## It loops through all libraries

## Usage:
## R CMD BATCH --no-save --no-restore '--args sample_info="/path/to/rna_seq_sample_info.tsv" calls_dir=$(RNASEQ_CLUSTER_BGEECALL_OUTPUT) rna_seq_calls_plot.R rna_seq_calls_plot.Rout
## bgeecall_sample_info       - path to file with info about libraries processed with BgeeCall
## calls_dir                  - path to folder where BgeeCall wrote the calls.


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
command_arg <- c("bgeecall_sample_info", "calls_dir")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

##load data
sample_info_data <- read.table(bgeecall_sample_info, sep = "\t", header = TRUE, comment.char = "")

## init variables
info_file <- "gene_cutoff_info_file.tsv"
all_samples <- data.frame(matrix(nrow = 0, ncol = 15))
samples_columns <- c("libraryId", "cutoffTPM", "proportionGenicPresent", "numberGenicPresent", "numberGenic", "proportionCodingPresent", "numberPresentCoding", "numberCoding", "proportionIntergenicPresent", "numberIntergenicPresent", "numberIntergenic", "pValueCutoff", "meanIntergenic", "sdIntergenic", "speciesId")
names(all_samples) <- samples_columns


## check that calls have been generated for all libraries
libraries_wo_calls <- 0
lib_dirs <- list.dirs(path = calls_dir, full.names = FALSE, recursive = FALSE)
message(nrow(sample_info_data), " libraries in the bgeecall info file")
for(line in seq(nrow(sample_info_data))) {
  library_id <- basename(as.character(sample_info_data$rnaseq_lib_path[line]))
  if(library_id %in% lib_dirs) {
    info_file_path <- file.path(calls_dir, library_id, info_file)
    if(file.exists(info_file_path)) {
      info <- read.table(info_file_path, row.names = 1)
      library_info <- data.frame(library_id, info["cutoffTPM",], info["proportionGenicPresent",],
        info["numberGenicPresent",], info["numberGenic",], info["proportionCodingPresent",], 
        info["numberPresentCoding",], info["numberCoding",], info["proportionIntergenicPresent",], 
        info["numberIntergenicPresent",], info["numberIntergenic",], info["pValueCutoff",], 
        info["meanIntergenic",], info["sdIntergenic",], sample_info_data[line,]$species)
      names(library_info) <- samples_columns
      all_samples <- rbind(all_samples, library_info)
    } else {
      warning(library_id, " : info file was not generated")
      libraries_wo_calls <- libraries_wo_calls + 1
    }
  } else {
    warning(library_id, " : library directory not created")
    libraries_wo_calls <- libraries_wo_calls + 1
  }
}

if (libraries_wo_calls > 0 ) {
  stop("Calls were not generated for ", libraries_wo_calls, " libraries.")
}
message(nrow(all_samples), "libraries found in the calls directory")
save(all_samples, file=file.path(calls_dir, "presence_absence_all_samples.RDa"))
write.table(all_samples, file = file.path(calls_dir, "presence_absence_all_samples.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

## PDF for all boxplot
pdf(file = paste0(calls_dir, "/presence_absence_boxplots.pdf"), width = 12, height = 5)
## 7 boxplots to plot, 1 per column in all_samples
## Manually set the y-scale limits for most of the columns. NA are here for columsn thta need to be printed in log scale
ylimits <- list(NA, c(0, 100), c(0, 70000), c(0, 70000), c(0, 100), c(0, 30000), c(0,30000),  c(0, 100), c(0, 30000), c(0,30000), c(0,0.3))

i = 1
for (column in c(2:12)){
  print(names(all_samples)[column])
  
  ## get the subset of data to use
  subsetNames <- names(all_samples)[column]
  subsetData <- as.numeric(as.character(all_samples[, column]))
  organism <- as.factor(as.character(all_samples$speciesId))
  
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
