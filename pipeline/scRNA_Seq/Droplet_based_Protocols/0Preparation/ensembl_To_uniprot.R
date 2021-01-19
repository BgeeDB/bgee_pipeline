## SFonsecaCosta, Jan 2021

## This script is used to:
## Retrieve ensembl ID from UNIPROT ID for each species

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNA_Seq_info_target.txt" output="output"' ensembl_To_uniprot.R ensembl_To_uniprot.Rout
## scRNASeq_Info --> File that results from annotation and metadata (libraries downloaded and with extra information as readlength or SRR)
## output --> Folder where we should save the results

## Libraries
library("biomaRt")

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}
## checking if all necessary arguments were passed....
command_arg <- c("scRNASeq_Info","output")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}
## Read scRNASeq_Info file. If file not exists, script stops
if( file.exists(scRNASeq_Info) ){
  scRNASeqAnnotation <- fread(scRNASeq_Info, h=T, sep="\t")
} else {
  stop( paste("The scRNASeq_Info file not found [", scRNASeq_Info, "]\n"))
}

for (species in unique(scRNASeqAnnotation$scientific_name)) {
  
  ## for the moment we just have 3 species with data for droplet-based pipeline
  if (species == "Homo sapiens"){
    speciesInfo <- "hsapiens"
  } else if (species == "Mus musculus"){
    speciesInfo <- "mmusculus"
  } else if (species == "Heterocephalus glaber") {
    ## use female version
    speciesInfo <- "hgfemale"
  } else {
    message("Species still not introduzed!")
  }
  
  mart <- useMart('ENSEMBL_MART_ENSEMBL')
  mart <- useDataset(paste0(speciesInfo,"_gene_ensembl"), mart)
  
  annotLookup <- getBM(mart = mart, attributes = c('ensembl_gene_id','uniprot_gn_id'),uniqueRows=TRUE)
  colnames(annotLookup)[2] <- c("markerGene_ID_UniProt_Ensembl")
  
  write.table(annotLookup, file = paste0(output, "/ensembl_Uniprot_", gsub(" ","_",species),".tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
}
