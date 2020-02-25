## SFonsecaCosta Nov 4, 2019

## This script is used to do the quality control of RNA-Seq libraries.
## The QC criteria are based on: 1) reads depth  2) coverage 3) % p_pseudoaligned and 4) % protein coding genes detected (as exemple, for polyA)
## In the end, 2 files are exported with information about the libraries that pass or not the QC and a summary stats: rna_seq_sample_info_QC.txt and SummaryInformation_QC.txt

## Usage:
## R CMD BATCH --no-save --no-restore '--args rna_seq_sample_info="rna_seq_sample_info.txt" kallisto_count_folder="kallisto_count_folder" argumentsList="argumentsList" output="output" ensembl_version=98 metazoa_version=45' rna_seq_QR.R rna_seq_QR.Rout
## rna_seq_sample_info --> file with info on mapped libraries
## kallisto_count_folder --> path to kallisto result folder
## argumentsList --> file with cutoff information
## output --> path where the files of genome_stats + QC and .Rout will be saved
## ensembl_version --> version of ensembl.org
## metazoa_version --> version of ensembl metazoa

## Note: To run this script you should run Kallisto version > 0.44 or higher

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
command_arg <- c("rna_seq_sample_info", "kallisto_count_folder", "argumentsList", "output", "ensembl_version", "metazoa_version")
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
collectStats <- function(speciesID, database){

  if (database == "EnsemblMetazoa"){
    url <- paste0("ftp://ftp.ensemblgenomes.org/pub/metazoa/release-",metazoa_version,"/mysql/")
    message("Download from metazoa ensembl!", "\n")
  } else if (database == "Ensembl"){
    message("Download from ensembl!", "\n")
    url <- paste0("ftp://ftp.ensembl.org/pub/release-",ensembl_version,"/mysql/")
  } else {
    exit("unknown database to collect genome statistics ",database)
  }
  message(url)
  result <- getURL(url,verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE)
  collectInfo <- as.vector(unlist(strsplit(result,"\n", fixed = TRUE)),mode="list")
  collectInfo <- do.call(rbind.data.frame, collectInfo)
  colnames(collectInfo) <- "species"
  result2 <- grep(paste0(speciesID,"_core_"),collectInfo$species)
  speciesBgeeNme <- collectInfo[result2,]
  genomeStatsFileName <- "genome_statistics.txt.gz"
  download.file(url = file.path(url, speciesBgeeNme, "genome_statistics.txt.gz"),
                destfile = file.path(output, genomeStatsFileName),
                method = "wget", quiet = FALSE)
  gzFiles <- list.files(path = output, pattern = "*.gz$", full.names = TRUE)
  unzipFile <- gunzip(gzFiles, remove=TRUE)
  renamed_file <- file.rename(file.path(output, "genome_statistics.txt"),
                              file.path(output, paste0("genome_statistics_", speciesID, ".txt")))
}

## Function to collect information for each library using the fastp.json output file
collectInformationFASTP <- function(annotation, kallisto_count_folder, library){

  libraryInfo <- annotation$libraryType[annotation$X.libraryId == library]
  fastpInfo <- paste0(kallisto_count_folder, library)

  fastpFile <- list.files(path=fastpInfo, pattern = "\\.fastp.json.xz$")
  fastpFile <- file.path(fastpInfo, fastpFile)
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
    message("Library ", library, " is PolyA", "\n")
    ## calculate proportion of protein_coding with TPM value higher then zero
    sizeProteinCoding <- nrow(proteinCoding)
    proteinCodingFilter <- nrow(proteinCoding[proteinCoding$tpm > 0, ])
    detected <- (proteinCodingFilter/sizeProteinCoding)*100
    ## collec all info for the library
    collectInfo <- data.frame(readsMap, coverage, CG, p_alignmentKallisto, detected)
    argumentsUsed <- dplyr::filter(arguments, library == "polyA")

    } else if (protocol == "miRNA"){
      message("Library ", library, " is miRNA", "\n")
      ## calculate proportion of miRNA with TPM value higher then zero
      sizemiRNA <- nrow(miRNA)
      miRNAFilter <- nrow(miRNA[miRNA$tpm > 0, ])
      detected <- (miRNAFilter/sizemiRNA)*100
      ## collec all info for the library
      collectInfo <- data.frame(readsMap, coverage, CG, p_alignmentKallisto, detected)
      argumentsUsed <- dplyr::filter(arguments, library == "miRNA")

      } else if (protocol == "lncRNA"){
      message("Library ", library, " is lncRNA", "\n")
      ## calculate proportion of lncRNA with TPM value higher then zero
      sizelncRNA <- nrow(lncRNA)
      lncRNAFilter <- nrow(lncRNA[lncRNA$tpm > 0, ])
      detected <- (lncRNAFilter/sizemlncRNA)*100
      ## collec all info for the library
      collectInfo <- data.frame(readsMap, coverage, CG, p_alignmentKallisto, detected)
      argumentsUsed <- dplyr::filter(arguments, library == "lncRNA")

      } else {
    message("Protocol not recognized!", "\n")
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
rna_seq_sample_info_QC <- file.path(output, "rna_seq_sample_info_QC.txt")
if (!file.exists(rna_seq_sample_info_QC)){
  file.create(rna_seq_sample_info_QC)
  cat("libraryId\texperimentId\tspeciesId\torganism\tgenomeFilePath\tdatabase\tplatform\tlibraryType\tlibraryInfo\treadLength\trunIds\tRNASeqProtocol\tInfoQC\n",file = rna_seq_sample_info_QC, sep = "\t")
} else {
  message("File ", rna_seq_sample_info_QC, " already exists.....\n")
}
## Export detailed file with summary stats
summaryInformation_QC <- file.path(output, "SummaryInformation_QC.txt")
if (!file.exists(summaryInformation_QC)){
  file.create(summaryInformation_QC)
  cat("libraryId\texperimentId\tspeciesId\torganism\tlibraryType\tlibraryInfo\treadLength\tRNASeqProtocol\treadsLib\tcoverage\tCG\tp_alignmentKallisto\tdetectedRNA\tInfoQC\tdescription\n",file = summaryInformation_QC, sep = "\t")
} else {
  message("File already exist.....", "\n")
}

######################################### APPLY FOR EACH LIBRARY ######################################################
## Download for all species present in the annotation file (rna_seq_sample_info)!
for (species in unique(annotation$organism)) {
  message("Collect genome stats from ensembl!")
  database <- as.character(unique(annotation$database[annotation$organism == species]))
  lengthName <- sapply(strsplit(species, " "), length)
  if (lengthName == 2){
    speciesInfo <- gsub(" ","_",species)
    speciesInfo <- tolower(speciesInfo)
  } else if (lengthName > 2) {
    speciesInfo <- strsplit(species, "\\s+")[[1]]
    speciesInfo <- paste0(speciesInfo[1], "_", speciesInfo[3])
    speciesInfo <- tolower(speciesInfo)
  }
  collectStats(speciesID = speciesInfo, database = database)

  message("Read genome stats")
  ## read download file from ensembl to collect genome size
  genomeSize <- read.table(file.path(output, paste0("genome_statistics_", speciesInfo, ".txt")), header=FALSE, sep="\t")
  genomeSize <- dplyr::filter(genomeSize, V2 == "ref_length")
  genomeSize <- genomeSize$V3

  ## collect libraries that belong to same species!
  for (libraryID in annotation$X.libraryId[annotation$organism == species]) {
      message("Treating ", libraryID)
    ## verify if library exist
    if (!dir.exists(file.path(kallisto_count_folder, libraryID))){
      message("The folder for this library is not present!")
      ## check if all files necessary are present
    } else if (dir.exists(file.path(kallisto_count_folder, libraryID)) &&  !file.exists(file.path(kallisto_count_folder, libraryID, "abundance_gene_level+fpkm+intergenic.tsv")) || !file.exists(file.path(kallisto_count_folder, libraryID, "run_info.json")) || length(list.files(path = file.path(kallisto_count_folder, libraryID), pattern = "\\.fastp.json.xz$")) == 0){
      message("For this library: ", libraryID, " missing files to run the QC!")
    } else {
      ## collect stats from fastp.json file
      fastpLibrary <- collectInformationFASTP(annotation = annotation, kallisto_count_folder = kallisto_count_folder, library = libraryID)
      kallistoInfo <- fromJSON(file = file.path(kallisto_count_folder, libraryID, "run_info.json"))
      abundaceFilePath <- file.path(kallisto_count_folder, libraryID, "abundance_gene_level+fpkm+intergenic.tsv")
      abundaceFile <- read.table(abundaceFilePath, header=TRUE, sep="\t")
      ## collect stats + QC information for the library
      finalInfo <- qc(fastpLibrary = fastpLibrary, kallistoInfo = kallistoInfo, abundaceFile = abundaceFile, genomeSize = genomeSize, library = libraryID, arguments = arguments)
      write.table(finalInfo[c(1:12,18)], file = paste0(output, "rna_seq_sample_info_QC.txt"), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(finalInfo[c(1:4,8:10,12:19)], file = paste0(output, "SummaryInformation_QC.txt"), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
    }
  }
}

