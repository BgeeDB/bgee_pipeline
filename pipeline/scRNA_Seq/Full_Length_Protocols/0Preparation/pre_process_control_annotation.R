## SFonsecaCosta 2019
## This script allow to verify the number of cells per cell-population regarding the: experimentID, species, cellTypeId, stageId, strain, uberonId and sex after the annotation process.
## Just cell-population that belongs to the same experimentID, species, cellTypeId, stageId, strain, uberonId and sex and have at least 50 cells will be keeped to continue in the pipeline.
## The output file generated (passScRNASeqLibrary.tsv) will be used to download the data and to continue in the pipeline.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeqLibrary="scRNASeqLibrary.tsv" cellsThreshold="minimum number of cells" output_folder="output_folder"' pre_process_control_annotation.R pre_process_control_annotation.Rout
## scRNASeqLibrary --> File from manual annotation (scRNASeqFLLibrary.tsv)
## cellsThreshold --> (optional) keep experiments with a number of libraries higher of equal to cellsThreshold. If not provided as argument then no filtering.
## output_file_pass --> path to the output file containing libraries passing the cell threshold
## output_file_not_pass --> path to the output file containing libraries passing the cell threshold

## Libraries
library(dplyr)

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
command_arg <- c("scRNASeqLibrary", "output_file_pass", "output_file_not_pass")
for (c_arg in command_arg) {
  if (!exists(c_arg)) {
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read scRNASeqLibrary file, if file not exists, script stops
if (!file.exists(scRNASeqLibrary)) {
  stop( paste("scRNASeqLibrary file not found [", scRNASeqLibrary, "]"))
}
annotation <- read.table(scRNASeqLibrary, h=T, sep="\t", comment.char="")
names(annotation)[1] <- "libraryId"
## replace all empty values per NA for sex and strain to avoid having NA and empty string
## considered as different condition
annotation$sex <- forcats::fct_explicit_na(annotation$sex)
annotation$strain <- forcats::fct_explicit_na(annotation$strain)

#create a new empty data.frame for libs that passed the threshold
passed <- data.frame(matrix(nrow = 0, ncol = ncol(annotation)))
colnames(passed) <- colnames(annotation)

# group by species, exp, condition and count number of libraries
group_by_libs_above_threshold <- annotation %>% group_by(speciesId, experimentId, cellTypeId_abInitio, stageId,
  strain, anatId, sex) %>% summarize(numberLibs = n())
# filtering of cellpopulation is done only if the cellsThreshold argument were provided
if (exists("cellsThreshold")) {
  group_by_libs_above_threshold <- group_by_libs_above_threshold %>% filter(numberLibs >= as.integer(cellsThreshold))
} else {
  warning("No filtering based on minimum number of cells per cell population")
}

for (rowNumber in seq(nrow(group_by_libs_above_threshold))) {
  subset <- annotation[annotation$speciesId %in% group_by_libs_above_threshold$speciesId[rowNumber]
                       & annotation$experimentId %in% group_by_libs_above_threshold$experimentId[rowNumber]
                       & annotation$cellTypeId_abInitio %in% group_by_libs_above_threshold$cellTypeId_abInitio[rowNumber]
                       & annotation$stageId %in% group_by_libs_above_threshold$stageId[rowNumber]
                       & annotation$strain %in% group_by_libs_above_threshold$strain[rowNumber]
                       & annotation$anatId %in% group_by_libs_above_threshold$anatId[rowNumber]
                       & annotation$sex %in% group_by_libs_above_threshold$sex[rowNumber],]
  passed <- rbind(passed, subset)
}
write.table(x = passed, file = output_file_pass, sep = "\t", col.names = TRUE, row.names = FALSE,
  quote = FALSE)
not_passed <- annotation[! annotation$libraryId %in% passed$libraryId,]
write.table(x = not_passed, file = output_file_not_pass, sep = "\t", col.names = TRUE, row.names = FALSE,
            quote = FALSE)
