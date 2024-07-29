## JW, Feb. 2024

## This script is used to generate a tsv file containing calls summary for all properly processed libraries

## Usage:
## R CMD BATCH --no-save --no-restore '--args metadata_file="metadata_file" calls_dir="calls_dir" summary_calls_file="summary_calls_file"' calls_summary.R calls_summary.Rout
## metadata_file            --> Path to the file containing metadata of all libraries for which calls have to be processed
## calls_dir                --> Path to the directory where calls have been generated
## summary_calls_file       --> Path to the tsv that will be created by this script and contains calls summary for all libraries/celltype

library(ggplot2)
library(gridExtra)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){
  stop("no arguments provided\n")
} else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("metadata_file", "calls_dir", "summary_calls_file")
for (c_arg in command_arg) {
  if (!exists(c_arg)) {
    stop(paste(c_arg,"command line argument not provided\n"))
  }
}
if (! file.exists(metadata_file)) {
  stop("metadata file ",metadata_file, " does not exist.")
}
metadata <- read.table(file = metadata_file, header = TRUE, sep = "\t")
summaryData <- c()
for (libraryId in sort(unique(metadata$library_id))) {
  libraryCallsFile <- file.path(calls_dir, libraryId, paste0(libraryId, "_stats.tsv"))
  if (! file.exists(libraryCallsFile)) {
    warning("an error happened in the generation of calls for library ", libraryId, ". The calls summary file for that library has not been generated.",
         " Please solve the issue or remove this library from the metadata file and rerun the rule.")
  } else {
    librarySummary <- read.table(file = libraryCallsFile, header = TRUE, sep = "\t")
    summaryData <- rbind(summaryData, librarySummary)
  }
}
write.table(x = summaryData, file = summary_calls_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
## final plot per species
pdf(file.path(calls_dir, paste0("All_libraries_stats_information_10X.pdf")), width=16, height=6)
g1 <- ggplot(summaryData, aes(x=organism, y=as.numeric(Proportion_genic_present))) + 
  geom_boxplot()+ylim(0,100)+xlab(" ")+ylab("% Genic Present")
g2 <- ggplot(summaryData, aes(x=organism, y=as.numeric(Proportion_coding_present))) + 
  geom_boxplot(notch=TRUE)+ylim(0,100)+xlab(" ")+ylab("% Protein Coding Present")
g3 <- ggplot(summaryData, aes(x=organism, y=as.numeric(Proportion_intergenic_present))) + 
  geom_boxplot(notch=TRUE)+ylim(0,100)+xlab(" ")+ylab("% Intergenic Present")
g4 <- ggplot(summaryData, aes(x=organism, y=as.numeric(meanRefIntergenic))) + 
  geom_boxplot(notch=TRUE)+ylim(0,1)+xlab(" ")+ylab("Mean Reference Intergenic")
g5 <- ggplot(summaryData, aes(x=organism, y=as.numeric(sdRefIntergenic))) + 
  geom_boxplot(notch=TRUE)+ylim(0,10)+xlab(" ")+ylab("SD Reference Intergenic")
g6 <- ggplot(summaryData, aes(x=organism, y=as.numeric(pValue_cutoff))) + 
  geom_boxplot(notch=TRUE)+ylim(0,1)+xlab(" ")+ylab("pValue_cutoff")
grid.arrange(g1,g2,g3,g4,g5,g6, nrow = 2, ncol = 3)
dev.off()
