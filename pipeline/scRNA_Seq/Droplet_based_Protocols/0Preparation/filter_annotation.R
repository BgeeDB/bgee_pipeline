## Julien Wollbrett, Mar 15 2023

## This script takes as input raw annotation and produce annotation containing only
## library/experiment meant to be processed by the Bgee pipeline.

## TODO: As it is the same annotation file used for target based and full length single cell,
## this script will probably move to single_cell_utils.R and be used by both pipelines.

## July 2023: added possibility to filter based on species and library IDs. the filtering is on speciesIds AND libraryIds.
##            to to so use attributes speciesIds and libraryIds as a comma separated list without space (e.g speciesIds="7727,9606")

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("scRNASeqExperiment","scRNASeqTBLibrary", "acceptedProtocols", "strainMapping", "outputDir")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

############################ Functions #################################

readTsvFile<- function(fileToRead, header = TRUE, sep = "\t", quote = "",
    commentChar = "") {
  if (file.exists(fileToRead)) {
    data <- read.table(fileToRead, sep = sep, header = header, quote = quote,
      comment.char = commentChar)
  } else  {
    stop(fileToRead, " file does not exist.")
  }
  return(data)
}

#########################################################################

### First check that all files exist and read them
raw_lib_annot <- readTsvFile(scRNASeqTBLibrary, quote = "\"", commentChar="")
raw_exp_annot <- readTsvFile(scRNASeqExperiment, quote = "\"", commentChar="")
protocols <- readTsvFile(acceptedProtocols)
standardised_strains <- readTsvFile(strainMapping, quote = "\"", commentChar = "")
colnames(standardised_strains)[1] <- "sourceStrain"
target_based_protocols <- protocols[protocols$Library_construction == "3'end",]

### Then filter on protocol accepted in Bgee pipeline
filtered_exp_annot <- c()
for (rowNumber in seq(nrow(raw_exp_annot))) {
  protocols_used_in_exp <- unlist(strsplit(as.character(raw_exp_annot$protocol[rowNumber]), ", "))
  if(any(protocols_used_in_exp %in% target_based_protocols$Protocols)) {
    filtered_exp_annot <- rbind(filtered_exp_annot, raw_exp_annot[rowNumber,])
  }
}
filtered_lib_annot <- raw_lib_annot[raw_lib_annot$protocol %in% target_based_protocols$Protocols
  & raw_lib_annot$whiteList %in% target_based_protocols$target_whiteList,]

## Then update Strain to fit standardized strain names
for (i in seq(nrow(standardised_strains))) {
  filtered_lib_annot$strain[filtered_lib_annot$strain == standardised_strains[i,]$sourceStrain
    & filtered_lib_annot$speciesId == standardised_strains[i,]$speciesId] <- standardised_strains[i,]$targetStrain
}

# If speciesIds are provided then use them to filter libraries to process in the pipeline
if (exists('speciesIds')) {
  speciesIds <- as.list(strsplit(x = speciesIds, split = ",")[[1]])
  filtered_lib_annot <- filtered_lib_annot[filtered_lib_annot$speciesId %in% speciesIds,]
}
# If libraryIds are provided then use them to filter libraries to process in the pipeline
if (exists('libraryIds')) {
  libraryIds <- as.list(strsplit(x = libraryIds, split = ",")[[1]])
  filtered_lib_annot <- filtered_lib_annot[filtered_lib_annot$libraryId %in% libraryIds,]
}
### Finally write filtered files in the output dir
colnames(filtered_exp_annot)[1] <- "experimentId"
colnames(filtered_lib_annot)[1] <- "libraryId"

write.table(x = filtered_lib_annot, file.path(outputDir, basename(scRNASeqTBLibrary)), col.names = TRUE,
  row.names = FALSE, sep = "\t", quote = FALSE)
write.table(x = filtered_exp_annot, file.path(outputDir, basename(scRNASeqExperiment)), col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)
