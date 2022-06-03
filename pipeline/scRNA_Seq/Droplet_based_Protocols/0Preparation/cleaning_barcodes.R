## SFonsecaCosta, June 2022

## This script is used to remove doublets and multiplets from the annotation file (information collected from the authors)
## Note: if a barcode is detected >=2 in same experiment/library this barcode is removed from the barcode annotation file and will not be used during the analyis

## Usage:
## R CMD BATCH --no-save --no-restore '--args barcodesFolder="path_barcodes_folder" output="output_folder"' cleaning_barcodes.R cleaning_barcodes.Rout
## barcodesFolder --> Folder where are the barcodes files annotated by Bgee for each experimentID
## output --> Path where should be saved the output results after removing the duplicate barcodes per experimentID/LibraryID

## libraries used
library(data.table)
library(dplyr)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("barcodesFolder", "output")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}
########################################################################################################################################################
getBarcodes_files <-  list.files(barcodesFolder, pattern = "scRNASeq_barcode_", full.names = TRUE)

for (barcodeFile in getBarcodes_files) {
  
  readBarcodeFile <- fread(barcodeFile, header = TRUE, sep="\t")
  fileName <- basename(barcodeFile)
  
  saveUniqueBarcodes <- c()
  saveDuplicatedBarcodes <- c()
  for (libraryID in unique(readBarcodeFile$library)) {
    
    selectBarcodes <- dplyr::filter(readBarcodeFile, readBarcodeFile$library == libraryID)
    
    getUniqueBarcodes <- selectBarcodes %>% 
      group_by(barcode) %>% 
      filter(n()==1)
    
    getBarcodesDuplicates <- selectBarcodes %>% 
      group_by(barcode) %>% 
      filter(n() != 1)
  
    saveUniqueBarcodes <- rbind(saveUniqueBarcodes, getUniqueBarcodes)
    saveDuplicatedBarcodes <- rbind(saveDuplicatedBarcodes, getBarcodesDuplicates)
  }
  write.table(saveUniqueBarcodes, file = file.path(output, fileName), sep = "\t", row.names = FALSE, quote = FALSE)
  if (nrow(saveDuplicatedBarcodes) != 0){
    message("Exist barcode duplicates in 1 or more libraryID of the correspondent experimentID: ", fileName)
    write.table(saveDuplicatedBarcodes, file = file.path(output, paste0("DuplicateBarcodes_", fileName)), sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    message("Not exist barcode duplicates in any libraryID of the correspondent experimentID: ", fileName)
  }
}
