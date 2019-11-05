## SFonsecaCosta Nov 4, 2019

## This script is used to do the quality control of RNA-Seq libraries.
## The QC criteria are based on: 1) % reads mapped to transcriptome 2) coverage 3) p_pseudoaligned and 4) % protein coding genes detected (as exemple, for polyA)
## In the end, 2 files are exported with information about the libraries that pass or not the QC and a summary stats: rna_seq_sample_info_QC.txt and SummaryInformation_QC.txt

## Usage:
## R CMD BATCH --no-save --no-restore '--args rna_seq_sample_info="rna_seq_sample_info.txt" kallisto_count_folder="kallisto_count_folder" argumentsList="argumentsList" output="output""' rna_seq_QR.R rna_seq_QR.Rout
## rna_seq_sample_info --> file with info on mapped libraries
## kallisto_count_folder --> path to kallisto result folder
## argumentsList --> file with cutoff information
## output --> path where the files of genome_stats + QC and .Rout will be saved

## Libraries used
library(RCurl)
library(R.utils)
library(rjson)

######################################## READING ARGUMENTS #########################################################
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("rna_seq_sample_info", "kallisto_count_folder", "argumentsList", "output")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read annotation file. If file not exists, script stops
if( file.exists(rna_seq_sample_info) ){
  annotation <- read.table(rna_seq_sample_info, header = TRUE , sep="\t", comment.char = "", quote = "")
} else {
  stop( paste("The rna_seq_sample_info file not found [", rna_seq_sample_info, "]\n"))
}
## Read argumentsList file. If file not exists, script stops
if( file.exists(argumentsList) ){
  arguments <- read.table(argumentsList, header = TRUE , sep="\t", comment.char = "", quote = "")
} else {
  stop( paste("The rna_seq_sample_info file not found [", argumentsList, "]\n"))
}

######################################## FUNCTIONS  #################################################################
## Function to collect genome_stats from ensembl
collectStats <- function(speciesID){
  setwd(output)
  metazoaDro <- grepl("drosophila_",speciesID)
  metazoaCel <- grepl("caenorhabditis_",speciesID)
  
  if (metazoaDro == "TRUE" || metazoaCel == "TRUE"){
    url1 <- "ftp://ftp.ensemblgenomes.org/pub/metazoa/current/mysql/"
    cat("Download from metazoa ensembl!", "\n")
    result <- getURL(url1,verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE)
    } else {
    cat("Download from ensembl!", "\n")
    url2 <- "ftp://ftp.ensembl.org/pub/current_mysql/"
    result <- getURL(url2,verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE)
    }
  
  result <- gsub("\n", " ", result, fixed=TRUE)
  collectInfo <- as.vector(unlist(strsplit(result," ")),mode="list") 
  collectInfo <- do.call(rbind.data.frame, collectInfo)
  colnames(collectInfo) <- "species"
  result2 <- grep(paste0(speciesID,"_core_"),collectInfo$species)
  speciesBgeeNme <- collectInfo[paste0(result2),]
  genomeStats.gz <- system(sprintf('%s %s', paste0("wget"), paste0(url2, speciesBgeeNme, "/", "genome_statistics.txt.gz")))
  gzFiles <- list.files(path = output, pattern = "*.gz$", full.names = TRUE)
  unzipFile <- gunzip(gzFiles, remove=TRUE)
  system(sprintf('%s %s %s', paste0("mv"), paste0(output, "genome_statistics.txt"), paste0(output, "genome_statistics_", speciesID, ".txt")))
}

## Function to collect information for each library using the fastp.json output file
collectInformationFASTP <- function(annotation, kallisto_count_folder, library){
  
  libraryInfo <- annotation$libraryType[annotation$X.libraryId == library]
  fastpInfo <- paste0(kallisto_count_folder, library)
  
  if(libraryInfo == "SINGLE"){
    fastpFile <- list.files(path=fastpInfo, pattern = "\\.fastp.json.xz$")
    fastpFile <- paste0(fastpInfo, "/", fastpFile)
  } else {
    fastpFile <- list.files(path=fastpInfo, pattern = "\\_1.fastp.json.xz$")
    fastpFile <- paste0(fastpInfo, "/", fastpFile)
  }
  ## If lanes exist, the reads would be summed and calculated the average for reads length 
  readsMapInfo <- 0
  readLength <- c()
  cgInfo <- c()
  for (files in fastpFile) {
    fastpName <- sub("\\..*", "", files)
    unZipFiles <- system(sprintf('%s < %s > %s', paste0("unxz"), paste0(files), paste0(fastpName, ".fastp.json")))
    ## unzip files
    filesInfo <- paste0(fastpName, ".fastp.json")
    for (i in filesInfo) {
      fastpjson <- fromJSON(file = i)
      readsInfo <- fastpjson$summary$before_filtering$total_reads
      readsLengthInfo <- fastpjson$summary$before_filtering$read1_mean_length
      CG_content <- fastpjson$summary$before_filtering$gc_content
      readsMapInfo <- readsMapInfo+readsInfo
      readLength <- rbind(readLength,readsLengthInfo)
      cgInfo <- rbind(cgInfo,CG_content)
    }
  }
  readLength <- mean(readLength)
  cgInfo <- mean(cgInfo)
  
  ## Export info about library
  libraryInfo <- data.frame(readsMapInfo, readLength, cgInfo)
  return(libraryInfo)
}  

## Function to do the QC per library
qc <- function(fastpLibrary, kallistoInfo, abundaceFile, genomeSize, library, arguments){
  
  libraryInfo <- annotation[annotation$X.libraryId == library,]
  protocol <- annotation$RNASeqProtocol[annotation$X.libraryId == library]
  libraryType <- annotation$libraryType[annotation$X.libraryId == library]
  readsMap <- fastpLibrary$readsMapInfo
  ## Coverage = read Length * Number of reads / genome size
  if (libraryType == "SINGLE"){
    coverage <- fastpLibrary$readLength * readsMap / genomeSize
  } else {
    coverage <- (2*fastpLibrary$readLength) * readsMap / genomeSize 
  }
  CG <- fastpLibrary$cgInfo
  p_alignmentKallisto <- kallistoInfo$p_pseudoaligned
  
  ## calculate propotion (for protein coding, lncRNA or miRNA)
  genicRegion <- nrow(dplyr::filter(abundaceFile, type == "genic"))
  proteinCoding <- dplyr::filter(abundaceFile, type == "genic" & biotype == "protein_coding")
  miRNA <- dplyr::filter(abundaceFile, type == "genic" & biotype == "miRNA")
  lncRNA <- dplyr::filter(abundaceFile, type == "genic" & biotype == "lncRNA")
  
  if (protocol == "polyA"){
    cat("Library ", library, "is PolyA", "\n")
    ## calculate proportion of protein_coding with TPM value higher then zero
    sizeProteinCoding <- nrow(proteinCoding)
    proteinCodingFilter <- nrow(proteinCoding[proteinCoding$tpm > 0, ])
    detected <- (proteinCodingFilter/sizeProteinCoding)*100
    ## collec all info for the library
    collectInfo <- data.frame(readsMap, coverage, CG, p_alignmentKallisto, detected)
    argumentsUsed <- dplyr::filter(arguments, library == "polyA")
    
    } else if (protocol == "miRNA"){
      cat("Library is: miRNA", "\n")
      ## calculate proportion of miRNA with TPM value higher then zero
      sizemiRNA <- nrow(miRNA)
      miRNAFilter <- nrow(miRNA[miRNA$tpm > 0, ])
      detected <- (miRNAFilter/sizemiRNA)*100
      ## collec all info for the library
      collectInfo <- data.frame(readsMap, coverage, CG, p_alignmentKallisto, detected)
      argumentsUsed <- dplyr::filter(arguments, library == "miRNA")
    
      } else if (protocol == "lncRNA"){
        cat("Library is: lncRNA", "\n")
      ## calculate proportion of lncRNA with TPM value higher then zero
      sizelncRNA <- nrow(lncRNA)
      lncRNAFilter <- nrow(lncRNA[lncRNA$tpm > 0, ])
      detected <- (lncRNAFilter/sizemlncRNA)*100
      ## collec all info for the library
      collectInfo <- data.frame(readsMap, coverage, CG, p_alignmentKallisto, detected)
      argumentsUsed <- dplyr::filter(arguments, library == "lncRNA")
      
      } else {
    cat("Protocol not recognized!", "\n")
  }
  collectInfo$InfoQC <- ""
  collectInfo <- cbind(libraryInfo, collectInfo)
  
  ## filter library based on the arguments cutoff
  if (collectInfo$readsMap < argumentsUsed$reads && collectInfo$coverage < argumentsUsed$coverage){
    collectInfo$InfoQC <- "1"
    collectInfo$description <- "reads and coverage are low than the cutoff"
  } else if (collectInfo$p_alignmentKallisto < argumentsUsed$p_alignment){
    collectInfo$InfoQC <- "2"
    collectInfo$description <- "p_alignment is lower than the cutoff"
  } else if (collectInfo$detected < argumentsUsed$RNA){
    collectInfo$InfoQC <- "3"
    collectInfo$description <- "proportion of tpm higher then zero is lower than the cutoff"
  } else {
    collectInfo$InfoQC <- "0"
    collectInfo$description <- "pass QC"
  }
  return(collectInfo)
}

######################################### CREATING OUTPUT FILES  ####################################################
## Export new rna_seq_sample_info file with information about quality control (pass or not pass library)
rna_seq_sample_info_QC <- paste0(output, "rna_seq_sample_info_QC.txt")
if (!file.exists(rna_seq_sample_info_QC)){
  file.create(rna_seq_sample_info_QC)
  cat("libraryId\texperimentId\tspeciesId\torganism\tgenomeFilePath\tdatabase\tplatform\tlibraryType\tlibraryInfo\treadLength\trunIds\tRNASeqProtocol\tInfoQC\n",file = paste0(output, "/","rna_seq_sample_info_QC.txt"), sep = "\t")	
} else {
  print("File already exist.....", "\n")
}
## Export detailed file with summary stats
SummaryInformation_QC <- paste0(output, "SummaryInformation_QC.txt")
if (!file.exists(SummaryInformation_QC)){
  file.create(SummaryInformation_QC)
  cat("libraryId\texperimentId\tspeciesId\torganism\tlibraryType\tlibraryInfo\treadLength\tRNASeqProtocol\treadsLib\tcoverage\tCG\tp_alignmentKallisto\tdetectedRNA\tInfoQC\tdescription\n",file = paste0(output, "/","SummaryInformation_QC.txt"), sep = "\t")	
} else {
  print("File already exist.....", "\n")
}

######################################### APPLY FOR EACH LIBRARY ######################################################
## Download for all species present in the annotation file (rna_seq_sample_info)!
for (species in unique(annotation$organism)) {
  cat("Collect genome stats from ensembl!", "\n")
  speciesInfo <- gsub(" ","_",species)
  speciesID <- tolower(speciesInfo)
  collectStats(speciesID = speciesID)
}

for (speciesID in unique(annotation$organism)) {
  speciesInfo <- gsub(" ","_",speciesID)
  speciesInfo <- tolower(speciesInfo)
 
  ## read download file from ensembl to collect genome size
  genomeSize <- read.table(paste0(output, "genome_statistics_", speciesInfo, ".txt"), header=FALSE, sep="\t")
  genomeSize <- dplyr::filter(genomeSize, V2 == "ref_length")
  genomeSize <- genomeSize$V3
  
  ## collect libraries that belong to same species!
  for (libraryID in annotation$X.libraryId[annotation$organism == speciesID]) {
    
    fastpInfo <- paste0(kallisto_count_folder, libraryID)
    ## collect stats from fastp.json file
    fastpLibrary <- collectInformationFASTP(annotation = annotation, kallisto_count_folder = kallisto_count_folder, library = libraryID)
    kallistoInfo <- fromJSON(file = paste0(kallisto_count_folder, libraryID, "/" , "run_info.json"))
    abundaceFilePath <- paste0(kallisto_count_folder, libraryID, "/" , "abundance_gene_level+fpkm+intergenic.tsv")
    abundaceFile <- read.table(abundaceFilePath, header=TRUE, sep="\t")
    ## collect stats + QC information for the library
    finalInfo <- qc(fastpLibrary = fastpLibrary, kallistoInfo = kallistoInfo, abundaceFile = abundaceFile, genomeSize = genomeSize, library = libraryID, arguments = arguments)
    write.table(finalInfo[c(1:12,18)], file = paste0(output, "rna_seq_sample_info_QC.txt"), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(finalInfo[c(1:4,8:10,12:19)], file = paste0(output, "SummaryInformation_QC.txt"), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
    
  }
}
