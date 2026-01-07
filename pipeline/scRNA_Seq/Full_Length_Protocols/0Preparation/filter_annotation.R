## SFonsecaCosta 2019
## This script allow to verify the number of cells per cell-population regarding the: experimentID, species, cellTypeId, stageId, strain, uberonId and sex after the annotation process.
## Just cell-population that belongs to the same experimentID, species, cellTypeId, stageId, strain, uberonId and sex and have at least 50 cells will be keeped to continue in the pipeline.
## The output file generated (passScRNASeqLibrary.tsv) will be used to download the data and to continue in the pipeline.

## Julien W Apr. 2025
## This script now also filter libraries based on their protocolType and only keep full-length libraries.
## renamed the script to better describe its purpose

## parameters::
## scRNASeqLibrary --> Library annotation file from manual annotation
## scRNASeqExperiment --> Experiment annotation file from manual annotation
## filteredLibrariesFile --> output file with the libraries that passed the filtering
## filteredExperimentsFile --> output file with the experiments that passed the filtering
## cellsThreshold --> (optional) keep experiments with a number of libraries higher of equal to cellsThreshold. If not provided as argument then no filtering.

## Libraries
library(dplyr)
library(forcats)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if (length(cmd_args) == 0 ){ stop("no arguments provided\n")
} else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("scRNASeqLibrary", "scRNASeqExperiment", "filteredLibrariesFile", "filteredExperimentsFile")
for (c_arg in command_arg) {
  if (!exists(c_arg)) {
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read scRNASeqLibrary file, if file not exists, script stops
if (!file.exists(scRNASeqLibrary)) {
  stop( paste("scRNASeqLibrary file not found [", scRNASeqLibrary, "]"))
}
library_annotation <- read.table(scRNASeqLibrary, h=TRUE, sep="\t", comment.char="")
names(library_annotation)[1] <- "libraryId"
experiment_annotation <- read.table(scRNASeqExperiment, h=TRUE, sep="\t", comment.char="")
names(experiment_annotation)[1] <- "experimentId"
## replace all empty values per NA for sex and strain to avoid having NA and empty string
## considered as different condition
library_annotation$sex <- fct_na_value_to_level(library_annotation$sex, level = "NA")
library_annotation$strain <- fct_na_value_to_level(library_annotation$strain, level = "NA")

print(paste("Number of libraries before filtering: ", nrow(library_annotation)))

#create a new empty data.frame for libs that passed the threshold
passed <- data.frame(matrix(nrow = 0, ncol = ncol(library_annotation)))
colnames(passed) <- colnames(library_annotation)

filtered_lib_annotation <- library_annotation[library_annotation$protocolType == "Full-length",]
print(paste("Number of libraries with protocolType = Full-length: ", nrow(filtered_lib_annotation)))


# group by species, exp, condition and count number of libraries
group_by_libs_above_threshold <- filtered_lib_annotation %>% group_by(speciesId, experimentId, cellTypeId, stageId,
  strain, anatId, sex) %>% summarize(numberLibs = n(), .groups = "drop")
# filtering of cellpopulation is done only if the cellsThreshold argument were provided
if (exists("cellsThreshold")) {
  group_by_libs_above_threshold <- group_by_libs_above_threshold %>% filter(numberLibs >= as.integer(cellsThreshold))
} else {
  print("No filtering based on minimum number of cells per cell population")
}

for (rowNumber in seq(nrow(group_by_libs_above_threshold))) {
  subset <- filtered_lib_annotation[filtered_lib_annotation$speciesId %in% group_by_libs_above_threshold$speciesId[rowNumber]
                       & filtered_lib_annotation$experimentId %in% group_by_libs_above_threshold$experimentId[rowNumber]
                       & filtered_lib_annotation$cellTypeId %in% group_by_libs_above_threshold$cellTypeId[rowNumber]
                       & filtered_lib_annotation$stageId %in% group_by_libs_above_threshold$stageId[rowNumber]
                       & filtered_lib_annotation$strain %in% group_by_libs_above_threshold$strain[rowNumber]
                       & filtered_lib_annotation$anatId %in% group_by_libs_above_threshold$anatId[rowNumber]
                       & filtered_lib_annotation$sex %in% group_by_libs_above_threshold$sex[rowNumber],]
  passed <- rbind(passed, subset)
}
write.table(x = passed, file = filteredLibrariesFile, sep = "\t", col.names = TRUE, row.names = FALSE,
  quote = FALSE)
not_passed <- filtered_lib_annotation[! filtered_lib_annotation$libraryId %in% passed$libraryId,]
print(paste("Number of FL libraries that passed the threshold of minimum number of cells per celltype: ", nrow(passed)))
print(paste("Number of FL libraries that did not passed the threshold of minimum number of cells per celltype: ", nrow(not_passed)))
# now we need to filter the experiment annotation file
# we need to keep only the experiments that have at least one library that passed the filtering
experiment_annotation <- experiment_annotation[experiment_annotation$experimentId %in% passed$experimentId,]
# write the experiment annotation file
write.table(x = experiment_annotation, file = filteredExperimentsFile, sep = "\t", col.names = TRUE, row.names = FALSE,
  quote = FALSE)
print(paste("Number of experiments that passed the filtering: ", nrow(experiment_annotation)))
