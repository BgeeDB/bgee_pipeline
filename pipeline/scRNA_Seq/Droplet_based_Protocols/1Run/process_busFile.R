## SFonsecaCosta, 2 Oct 2019

## This script is used to process the bus files to an equivalent-class-UMI and to a gene-UMI
## As suggested by bustools this step is divide in 3 sub-steps:

## 1) Correct the barcodes using bustools correct: fix the barcodes that are within one hamming distance of the barcodes in the whitelist using whitelist.txt
## 2) Sort the busfile using bustools sort: organize the busfile by barcode, UMI, set and multiplicity.
## 3) Count records in the BUS with bustools count: generate the UMI count matrix using transcripts_to_genes.txt.

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNA_Seq_info_TargetBased.txt" folder_data="folder_data" folderSupport="folderSupport" whiteList_Path="whiteList_Path"' process_busFile.R process_busFile.Rout
## scRNASeq_Info --> File that results from annotation and metadata (libraries downloaded and with extra information as SRR) 
## folder_data --> Folder where are all the libraries in fastq format
## folderSupport --> Folder where is placed the informative files as: transcriptomes index + gtf_all + transcript_to_gene 
## whiteList_Path --> Folder where is located the barcode_whitelist files

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

## checking if all necessary arguments were passed.
command_arg <- c("scRNASeq_Info", "folder_data", "folderSupport", "whiteList_Path")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read scRNA-Seq info file. If file not exists, script stops
if( file.exists(scRNASeq_Info) ){
  scRNAInfo <- fread(scRNASeq_Info)
} else {
  stop( paste("The annotation file not found [", scRNASeq_Info, "]\n"))
}
###############################################################################################
## for each library 
for (library in unique(scRNAInfo$libraryId)) {
  
  ## verify if library exist
  libraryInfo <- file.exists(paste0(folder_data, "/", library, "/"))
  
  if (libraryInfo == TRUE){
    
    ## verify if busOutput exist for the library
    pathBusOut <-  paste0(folder_data, library, "/busOutput")
    if (!dir.exists(pathBusOut)){
      cat("The Kallisto bus folder doesn't exist for the library: ", library, "\n")
    } else {
      cat("Start correction, sort and counts for the library: ", library, "\n")
      
      ## Note: the whiteList we use in this pipeline are the files provided by 10X platform (add to source files)
      collectWhitelist <- as.character(scRNAInfo$whiteList[scRNAInfo$libraryId == library]) 
      selectedWhitheList <- paste0("10X_",collectWhitelist)
      cat("whitelist:  ", selectedWhitheList, "\n")
      
      ## step 1 --> correct the barcodes
      cat("Correct barcodes......", "\n")
      system(sprintf('%s -w %s -o %s %s', paste0(bustools, " ", "correct"), paste0(whiteList_Path, "barcode_whitelist_", selectedWhitheList,".txt"), paste0(pathBusOut, "/output.correct.bus"), paste0(pathBusOut, "/output.bus")))
      
      ## step 2 --> sort the bus file
      cat("Sort bus file......", "\n")
      system(sprintf('%s -t 4 -o %s %s', paste0(bustools, " ", "sort"), paste0(pathBusOut, "/output.correct.sort.bus"), paste0(pathBusOut, "/output.correct.bus")))
      
      ## Creat folders to export the information per TCC and gene matrix (counts)
      tcc_counts <- paste0(pathBusOut, "/tcc_counts")
      if (!dir.create(tcc_counts)){
        dir.create(tcc_counts)
      } else {
        print("File already exist.....")
      }
      
      gene_counts <- paste0(pathBusOut, "/gene_counts")
      if (!dir.create(gene_counts)){
        dir.create(gene_counts)
      } else {
        print("File already exist.....")
      }
      
      collectSpecies <- as.character(scRNAInfo$scientific_name[scRNAInfo$libraryId == library]) 
      collectSpecies <- gsub(" ", "_", collectSpecies)
      
      ## step 3 --> count with bustools count
      ## TCC level
      cat("TCC level......", "\n")
      system(sprintf('%s -o %s -g %s -e %s -t %s %s', paste0(bustools, " ", "count"), paste0(tcc_counts,"/tcc"), paste0(folderSupport, "/transcript_to_gene_with_intergenic_", collectSpecies, ".tsv"),paste0(pathBusOut, "/matrix.ec"), paste0(pathBusOut, "/transcripts.txt"), paste0(pathBusOut, "/output.correct.sort.bus")))
      ## GENE level
      cat("Gene level......", "\n")
      system(sprintf('%s -o %s -g %s -e %s -t %s %s %s', paste0(bustools, " ", "count"), paste0(gene_counts,"/gene"), paste0(folderSupport, "/transcript_to_gene_with_intergenic_", collectSpecies, ".tsv"), paste0(pathBusOut, "/matrix.ec"), paste0(pathBusOut, "/transcripts.txt"), paste0("--genecounts"), paste0(pathBusOut, "/output.correct.sort.bus")))
      
    }
  } else {
    cat("Library ", library ," not present in the folder to process.", "\n")
  }
}
