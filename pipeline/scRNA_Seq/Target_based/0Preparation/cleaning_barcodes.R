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
if(!dir.exists(output)) {
  dir.create(output)
}

barcodes_files_path <-  list.files(barcodesFolder, pattern = "scRNASeq_barcode_", full.names = TRUE)

for (barcode_file_path in barcodes_files_path) {
  
  barcodes <- fread(barcode_file_path, header = TRUE, sep="\t")
  barcode_file_name <- basename(barcode_file_path)
  
  all_unique_barcodes <- c()
  all_duplicated_barcodes <- c()
  for (library_id in unique(barcodes$library)) {
    
    select_barcodes <- dplyr::filter(barcodes, barcodes$library == library_id)
    
    unique_barcodes <- select_barcodes %>% 
      group_by(barcode) %>% 
      filter(n()==1)
    
    barcodes_duplicates <- select_barcodes %>% 
      group_by(barcode) %>% 
      filter(n() != 1)
  
    all_unique_barcodes <- rbind(all_unique_barcodes, unique_barcodes)
    saveDuplicatedBarcodes <- rbind(all_duplicated_barcodes, barcodes_duplicates)
  }
  write.table(all_unique_barcodes, file = file.path(output, barcode_file_name), sep = "\t", row.names = FALSE, quote = FALSE)
  if (isTRUE(all_duplicated_barcodes) && nrow(all_duplicated_barcodes) != 0){
    message("Exist barcode duplicates in 1 or more libraryID of the correspondent experimentID: ", barcode_file_name)
    write.table(all_duplicated_barcodes, file = file.path(output, paste0("DuplicateBarcodes_", barcode_file_name)), sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    message("Not exist barcode duplicates in any libraryID of the correspondent experimentID: ", barcode_file_name)
  }
}
