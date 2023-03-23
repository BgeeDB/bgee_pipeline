## Julien Wollbrett, Mar 23 2023

## This script takes as input metadata allowing to download all FCA data and generate
## a metadata file similar to the one generated from Bgee library annotation file

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("FCAMetadata","outputMetadataFile")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

data_FCA <- read.table(FCAMetadata, sep = "\t", header = TRUE)
header <- c("sample_accession", "experiment_id", "library_id", "run_accession", "read_count", "tax_id", "scientific_name", "instrument_model", "library_layout", "fastq_ftp", "submitted_ftp", "source")
data_FCA$experimentId <- "ERP129698"
data_FCA$tax_id <- "7227"
data_FCA$scientific_name <- "Drosophila melanogaster"
data_FCA$fastq <- paste(data_FCA[,56], data_FCA[,58], sep = ";")
data_FCA$source <- "FCA"
reorder_FCA <-data_FCA[,c(3, 63, 2, 53, 1, 64, 65, 1, 22, 66, 66, 67)]
colnames(reorder_FCA) <- header
reorder_FCA$read_count <- ""
reorder_FCA$instrument_model <- ""
write.table(x = reorder_FCA, outputMetadataFile, row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)
