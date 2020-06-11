## SFC, 20 Sep 2019
## Modified May 2020

## This script is used to run Kallisto bus in all libraries from all species and experiments

## Usage:
## R CMD BATCH --no-save --no-restore '--args metadata_file="metadata_info_10X.txt" annotation_file="scRNASeqLibrary.tsv" folder_data="folder_data" folderSupport="folderSupport" output="output"' Kallisto_bus.R Kallisto_bus.Rout
## metadata_file --> Metadata file exported after run the script retrieve_metadata.R
## annotation_file --> Bgee annotation file with information about all libraries 
## folder_data --> Folder where all the libraries are located
## folderSupport --> Folder where is placed the informative files as transcriptomes index + gtf_all 
## output --> Path where should be saved the .Rout file

## libraries used
library(stringr)
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
command_arg <- c("metadata_file","annotation_file", "folder_data", "folderSupport", "output")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read metadata file. If file not exists, script stops
if( file.exists(metadata_file) ){
  metadata <- read.table(metadata_file, h=T, sep="\t", comment.char="")
} else {
  stop( paste("The metadata file not found [", metadata_file, "]\n"))
}

## Read annotation file. If file not exists, script stops
if( file.exists(annotation_file) ){
  annotation <- read.table(annotation_file, h=T, sep="\t", comment.char="")
  colnames(annotation)[1] <- "libraryId"
} else {
  stop( paste("The annotation file not found [", annotation_file, "]\n"))
}
##########################################################################################################################################################
## first merge information from the metadata (like: reads information and SRR) and annotation file (from Bgee)
metadataCollect <- metadata[c(1,2,3,5,7)]
colnames(metadataCollect)[1] <- "libraryId"
targetBased <- data.frame(dplyr::filter(annotation, protocol == "10X Genomics" & protocolType == "3'end"))

## generate the informative file of the target based protocols with all the information!
scRNASeqInfo <- merge(targetBased, metadataCollect, by="libraryId", incomparables = NaN)
write.table(scRNASeqInfo, file = paste0(output, "/scRNA_Seq_info_TargetBased.txt"), col.names = TRUE, row.names = FALSE , quote = FALSE, sep = "\t")

for (species in unique(scRNASeqInfo$scientific_name)) {
  cat("Species:", species, "\n")
  
  ## collect species info
  specieID <- gsub(" ", "_", species)
  specieID <- list.files(folderSupport, pattern = paste0("^", specieID, ".*", "transcriptome.idx", "$"))
  print(specieID)
  
  ## collect all libraries from the species
  collectLibrary <- scRNASeqInfo$libraryId[scRNASeqInfo$scientific_name == species]
  
  for (i in collectLibrary) {
    cat("Treating library: ", i, "\n")
    
    ## verify if library exist
    libraryInfo <- file.exists(paste0(folder_data, "/", i, "/"))
    
    if (libraryInfo == TRUE){
      
      ## collect info about whitelist
      whiteLInfo <- as.character(scRNASeqInfo$whiteList[scRNASeqInfo$libraryId==i]) 
      cat("The whitelist used is:", whiteLInfo, "\n")
      
      ## verify if exist FASTQ (with or without recursive folder) or SRR folder with fastq.gz files for the library
      detectFastqPath_1 <- file.exists(paste0(folder_data, "/", i, "/", "FASTQ"))
      detectFastqPath_2 <- list.dirs(paste0(folder_data, "/", i, "/", "FASTQ"), recursive=TRUE)[-1]
      detectFastqPath_2 <- rlang::is_empty(detectFastqPath_2)
      
      
      if (detectFastqPath_1 == TRUE & detectFastqPath_2 == TRUE){
        
        print("This library come from HCA repository.")
        ## select all fastq.gz files
        detectFastqPath <- list.dirs(paste0(folder_data, "/", i, "/", "FASTQ"), recursive=TRUE)
        detectFastqFiles <- list.files(path=detectFastqPath, pattern = "\\.fastq.gz$")
        cat("FASTQ Files detected: ", detectFastqFiles, "\n") 
        
        ## Note: libraries from HCA have just one R1 and one R2
        ReadBarcodes <- list.files(path=detectFastqPath, pattern = "*_R1_001.fastq.gz$")
        ReadSeq <- list.files(detectFastqPath,"*_R2_001.fastq.gz")
        
        ## collect all fastq.gz files
        filesKallisto <- rbind(paste0(detectFastqPath,"/",ReadBarcodes), paste0(detectFastqPath,"/",ReadSeq))
        filesKallisto <- toString(filesKallisto)
        filesKallisto <- gsub(",", " " ,filesKallisto)
        cat("Files to pass to Kallisto: ", "\n", filesKallisto, "\n")
        
        ## create directory for bus_output for each library
        busOutput <- paste0(folder_data, "/", i, "/busOutput")
        if (!dir.create(busOutput)){
          dir.create(busOutput)
        } else {
          print("File already exist.....")
        }
        
        ## RUN Kallisto bus
        system(sprintf('%s -i %s -o %s -x %s -t 4 %s', paste0("kallisto bus"), paste0(folderSupport, "/", specieID), paste0(busOutput), paste0("10x",whiteLInfo), paste0(filesKallisto)))
        
      } else if (detectFastqPath_1 == TRUE & detectFastqPath_2 == FALSE){
        
        print("This library come from EBI repository. Downloaded bam files that were converted to fastq.gz.")
        ## verify if exist multiplefiles in the folder (multiple fastq.gz represent lanes)
        detectFastqPath <- list.dirs(paste0(folder_data, "/", i, "/", "FASTQ"), recursive=TRUE)[-1]
        detectFastqFiles <- list.files(path=detectFastqPath, pattern = "\\.fastq.gz$")
        cat("FASTQ Files detected: ", detectFastqFiles, "\n") 
        
        ReadBarcodes <- str_subset(detectFastqFiles,"bamtofastq_S1_L0\\d+_R1_\\d+")
        cat("Length of barcode files detected: ", length(ReadBarcodes), "\n")
        ReadSeq <- str_subset(detectFastqFiles,"bamtofastq_S1_L0\\d+_R2_\\d+")
        cat("Length sequence files detected: ", length(ReadSeq), "\n")
        
        ## collect all fastq.gz files
        filesKallisto <- rbind(paste0(detectFastqPath,"/",ReadBarcodes), paste0(detectFastqPath,"/",ReadSeq))
        filesKallisto <- toString(filesKallisto)
        filesKallisto <- gsub(",", " " ,filesKallisto)
        cat("Files to pass to Kallisto: ", "\n", filesKallisto, "\n")
        
        ## create directory for bus_output for each library
        busOutput <- paste0(folder_data, "/", i, "/busOutput")
        if (!dir.create(busOutput)){
          dir.create(busOutput)
        } else {
          print("File already exist.....")
        }
        
        ## RUN Kallisto bus
        system(sprintf('%s -i %s -o %s -x %s -t 4 %s', paste0("kallisto bus"), paste0(folderSupport, "/", specieID), paste0(busOutput), paste0("10x",whiteLInfo), paste0(filesKallisto)))
        
      } else {
        
        print("This library come from SRA repository. fastq.gz files are saved directlly in a SRR folder")
        detectFastqPath <- list.dirs(paste0(folder_data, "/", i), recursive=TRUE)[-1]
        detectFastqFiles <- list.files(path=detectFastqPath, pattern = "\\.fastq.gz$")
        cat("FASTQ Files detected: ", detectFastqFiles, "\n") 
        
        ReadBarcodes <- list.files(path=detectFastqPath, pattern = "*_R1.fastq.gz$")
        ReadSeq <- list.files(detectFastqPath,"*_R2.fastq.gz")
        
        ## collect all fastq.gz files
        filesKallisto <- rbind(paste0(detectFastqPath,"/",ReadBarcodes), paste0(detectFastqPath,"/",ReadSeq))
        filesKallisto <- toString(filesKallisto)
        filesKallisto <- gsub(",", " " ,filesKallisto)
        cat("Files to pass to Kallisto: ", "\n", filesKallisto, "\n")
        
        ## create directory for bus_output for each library
        busOutput <- paste0(folder_data, "/", i, "/busOutput")
        if (!dir.create(busOutput)){
          dir.create(busOutput)
        } else {
          print("File already exist.....")
        }
        
        ## RUN Kallisto bus
        system(sprintf('%s -i %s -o %s -x %s -t 4 %s', paste0("kallisto bus"), paste0(folderSupport, "/", specieID), paste0(busOutput), paste0("10x",whiteLInfo), paste0(filesKallisto)))
        
      }
    } else {
      cat("Library not present in the folder to run Kallisto bus!", "\n")
    }
  }
}
