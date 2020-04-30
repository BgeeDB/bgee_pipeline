## SFonsecaCosta, Mars 2020

## This script is used to add the aleatory cut-off (used as q-values) and to make the adj_calls in the MAS5 output files.

## Usage:
## R CMD BATCH --no-save --no-restore '--args affy_anno="affy_annotation_file" mas5_path="mas5_path_out" cut_off="cut_off_value"' adj_calls_mas5.R adj_calls_mas5.Rout
## affy_anno --> annotation file from bgee
## mas5_path --> path to the output files from MAS5 (files that contain: probeId, expression/calls or calls/expression). Note the script deal with different order of the columns.
## cut_off --> correspond to the q-value cutoff that will be applied to call present or absent probes in the adj_calls column.

## libraries used
library(data.table)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed in command line
command_arg <- c("affy_anno","mas5_path", "cut_off")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## reading annotation file
if( file.exists(affy_anno) ){
  annotation <- read.table(affy_anno, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "chipId"
  annotation <- dplyr::filter(annotation, normalizationTypeId == "1")
} else {
  stop( paste("Affymetrix annotation file not found [", affy_anno, "]\n"))
}

#########################################################################
## function to add qValue and than do adj_call
mas5 <- function(chipID, cut_off){
  
  chip <- fread(chipID)
  ## get character column with P,M,A info to attribute q-value
  columnCalls <- chip[, sapply(chip[,2:3], class) == 'character']
  columnCalls <- names(columnCalls[columnCalls == TRUE])
  
  ## target column with call info and based on that add qValue and than adjusted call
  if (columnCalls == "V2"){
    colnames(chip) <- c("probeId", "call", "expression")
    ## add aleatory q-value
    chip$qValue <- ifelse(chip$call == "A" , 0.5, ifelse(chip$call == "P" , 0.01, 0.05 ))
    ## adj_calls
    chip$adjusted_call <- ifelse(chip$qValue <= cut_off , "P", "A" )
  } else {
    colnames(chip) <- c("probeId", "expression", "call")
    ## add aleatory q-value
    chip$qValue <- ifelse(chip$call == "A" , 0.5, ifelse(chip$call == "P" , 0.01, 0.05 ))
    ## adj_calls
    chip$adjusted_call <- ifelse(chip$qValue <= cut_off , "P", "A" )
  }
  return(chip)
}

## loop through all chips in each experimentId
for (experimentId in unique(annotation$experimentId)) {
  
  filePath <- file.path(mas5_path, experimentId)
  
  if (file.exists(filePath) == "TRUE"){
    AllFiles <- list.files(filePath, full.names = TRUE, recursive = TRUE)
    
    for (chip in AllFiles) {
      nameChipFile <- basename(chip)
      adjCalls <- mas5(chipID = chip, cut_off = cut_off)
      ## re-write the file with adj_call info
      write.table(adjCalls, file = file.path(filePath, nameChipFile), sep = "\t", row.names = FALSE)
    }
    cat("Adjusted calls added to all chipId files in", experimentId, "experiment.")
    
  } else {
    cat("Folder not exit for this experiment,", experimentId, "\n")
  }
}

