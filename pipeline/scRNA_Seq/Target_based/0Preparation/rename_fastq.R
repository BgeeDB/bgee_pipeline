# Julien Wollbrett
# created: 27/02/2024

## This script allows to rename fastq files. Different logic of renaming can be run depending on
## how the fastq file was generated
## 1. generated from bamtofastq
##   When bam files are downloaded from SRA we use the tool bamtofastq to generate fastq files.
##   The name of those files does not follow our standard naming so we have to rename them.
##   More than just the name, it hapens that one file type (I1, R1, R2) exists for several lanes (L001, L002, ...)
##   for the same run. We merge those lanes in one file.
## 2. downloaded using fastq-dump
##   In that case file names are RUNID_1.fastq, RUNID_2.fastq and sometimes also RUNID_3.fastq
##   In the high majority of cases _1, _2 and _3 in fastq file name can be used to properly detect I1, R1 and R2 files
##   BUT there are some exceptions. We tried to use read length with the assumption that I1 < R1 < R2 to properly rename
##   fastq files but as described in https://www.biostars.org/p/9588563/ it created a bug. In order to solve that issue
##   we implemented a new approach following the logic of the fasterqParseR package (https://github.com/Nusob888/fasterqParseR)
##   We decided to reimplement that logic to allow to rename already gzipped fastq files and to be able to update the code in
##   case a new protocol requires to update the logic of this renaming.

library(readr)

## reading arguments
cmd_args = commandArgs(TRUE)
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed.
command_arg <- c("fastqPath", "runId", "renaming", "whitelistFolder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

##################### FUNCTIONS ######################

## rename fastq files downloaded from fastq-dump
## this function is based on the assignSRAReads function from the package fasterqParseR
## https://github.com/Nusob888/fasterqParseR
renameFastqFilesFromFastqDump <- function(fastqPath = NULL, whitelistFolder = NULL, runId = NULL) {
  
  renamingFileInfo <- file.path(fastqPath, "renaming.info")
  
  #read whitelists
  ##TODO do not hardcode name of whitelists
  whitelists <- list(
    readr::read_csv(file = file.path(whitelistFolder, "barcode_whitelist_10X_Genomics_V2.txt.zip"), col_names = FALSE),
    readr::read_csv(file = file.path(whitelistFolder, "barcode_whitelist_10X_Genomics_V3.txt.zip"), col_names = FALSE))
  names(whitelists) <- c("10xV2", "10xV3")
  
  summary <- c()
  fastqFiles <- list.files(fastqPath, pattern="fastq.gz$|fastq$" ,full.names = TRUE)
  for (fastqFile in fastqFiles) {
    # read fastq files with cat or zcat depending on file extension
    catCmd <- ifelse(grepl(".gz", fastqFile), "zcat", "cat")
    # calculate mean read length to distinguish I1 to R1 and R2
    meanReadLength <- mean(as.numeric(system(paste0(catCmd, " ", fastqFile, " 2>/dev/null | head -1000 ",
                                                    "| awk '{if(NR%4==2) print length($1)}'"), intern = T)))
    if(meanReadLength %in% seq(5,10)){
      assignedRead <- "I1"
      chemistry <- NA
    } else if(meanReadLength > 10){
      sequences <- system(paste0(catCmd, " ", fastqFile, " 2>/dev/null | head -40000 | awk '{if(NR%4==2) print /^@/ ? $1 : ",
                                 " substr($0,1,16)}'"), intern = TRUE)
      
      whitelistCounts <- cbind(lapply(whitelists, function(x){
        sum(sequences %in% x$X1)
      }))
      
      if (sum(unlist(whitelistCounts)) > 1000) {
        assignedRead <- "R1"
        chemistry <- names(whitelistCounts[which.max(whitelistCounts),])
      } else{
        assignedRead <- "R2"
        chemistry <- NA
      }
    }
    newFastqName <- file.path(fastqPath, paste0(runId, "_", assignedRead,
                                                ifelse(catCmd == "cat", ".fastq", ".fastq.gz")))
    assigned <- data.frame(fastqFile, meanReadLength, assignedRead, chemistry, newFastqName)
    summary <- rbind(summary, assigned)
  }
  if (nrow(summary[summary$assignedRead == "R2",]) != 1 | nrow(summary[summary$assignedRead == "R1",]) != 1) {
    stop("can not rename fastq files as more or less than one R1 or R2 files have been detected for the run ", runId)
  }
  
  
  if (! identical(summary$fastqFile, summary$newFastqName)) {
    for (i in seq_len(nrow(summary))) {
      message("mv ",as.character(summary$fastqFile[i]), " ", paste0(as.character(summary$newFastqName[i]),"_new"))
      file.rename(from = as.character(summary$fastqFile[i]), to = paste0(as.character(summary$newFastqName[i]),"_new"))
    }
    for (i in seq_len(nrow(summary))) {
      message("mv ",paste0(summary$fastqFile[i],"_new"), " ", summary$newFastqName[i])
      file.rename(paste0(as.character(summary$newFastqName[i]),"_new"), as.character(summary$newFastqName[i]))
    }
    write.table(x = summary, file = renamingFileInfo, sep = "\t", quote = FALSE,
      col.names = TRUE, row.names = FALSE)
      message("properly renamed files for the run ", runId)
  } else {
      message("fastq files already named properly")
  }
}

renameFastqFilesFromBamToFastq <- function(fastqFileType = NULL, fastqPath=NULL, fileIsMandatory = NULL, runId = NULL){
  patternRegex <- paste0("_", fastqFileType, "_|_", fastqFileType, ".f")
  filelist <- list.files(input_dir, pattern = patterRegex, full.names = TRUE)
  renamedFile = paste0(fastqPath,"/", runId, "_", fastqFileType, ".fastq.gz")
  if (fileIsMandatory && length(filelist) == 0) {
    stop("No ",fastqFileType, " fastq file found for that run")
  }
  # rename R1 file in only one file
  if (length(filelist) == 1) {
    file.rename(from = filelist[1], to = renamedFile)
    # merge files into a new one called with the proper name if several lanes
  } else if (length(filelist) > 1) {
    mergedCommand = "cat "
    removeCmds = c()
    for (file in filelist) {
      mergedCommand <- paste0(mergedCommand, file, " ")
      removeCmds <- append(removeCmds, paste0("rm ", file))
    }
    mergedCommand = paste0(mergedCommand,  "> ", renamedFile)
    system(mergedCommand)
    for (removeCmd in removeCmds) {
      system(removeCmd)
    }
  }
}

########################## MAIN ###########################

if (renaming != "bamtofastq" && renaming != "fastqdump") {
  stop("renaming approach called ", renaming, " is not recognized.")
} else if (renaming == "bamtofastq") {
  i1 <- renameFastqFilesFromBamToFastq("I1", fastqPath, 0)
  r1 <- renameFastqFilesFromBamToFastq("R1", fastqPath, 1)
  r2 <- renameFastqFilesFromBamToFastq("R2", fastqPath, 1)
  create.table(x = rbind(i1, r1, r2), file = renamingFileInfo, sep ="\t", quote = FALSE,
               col.names = TRUE, row.names = FALSE)
} else if (renaming == "fastqdump") {
  renameFastqFilesFromFastqDump(fastqPath, whitelistFolder, runId)
}

