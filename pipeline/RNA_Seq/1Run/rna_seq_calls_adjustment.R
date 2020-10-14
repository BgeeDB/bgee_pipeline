## SFonsecaCosta Mars 2020
## This script is used to make the calls adjustment for all libraries of RNA-Seq pipeline.

## Usage:
## R CMD BATCH --no-save --no-restore '--args rna_seq_sample_info="rna_seq_sample_info.txt" RNASeq_folder_calls="RNASeq_folder_calls" sum_by_species_folder="sum_by_species_folder" gaussianFile="gaussian.txt" cutoff="cutoff"' rna_seq_calls_adjustment.R rna_seq_calls_adjustment.Rout
## rna_seq_sample_info --> file with info on mapped libraries
## RNASeq_folder_calls --> path to folder where presence/absence files are stored per library
## sum_by_species_folder --> Folder were are the files sum_by_species
## gaussianFile --> file with selected gaussian
## cutoff --> qValue cut-off used to call present/absent genes

## Libraries used
library(data.table)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}
## checking if all necessary arguments were passed....
command_arg <- c("rna_seq_sample_info", "RNASeq_folder_calls", "sum_by_species_folder", "gaussianFile" ,"cutoff")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}
## reading the input file
if( file.exists(rna_seq_sample_info) ){
  annotation <- read.table(rna_seq_sample_info, h=T, sep="\t", comment.char="")
  colnames(annotation)[1] <- "libraryId"
} else {
  stop( paste("rna_seq_sample_info file not found [", rna_seq_sample_info, "]\n"))
}
## reading the input file
if( file.exists(rna_seq_sample_info) ){
  gaussian <- read.table(gaussianFile, h=T, sep="\t", comment.char="")
} else {
  stop( paste("rna_seq_sample_info file not found [", rna_seq_sample_info, "]\n"))
}

##################################################################################
## function to retrieve max value of the intergenic 
max_intergenic <- function(sum_by_species, gaussian, species){
  ## select max TPM value of intergenic for the species having in consideration the gaussian selection
  if (gaussian$selectionSideIntergenic[gaussian$speciesId == species] == "Left"){
    max_intergenic <- max(sum_by_species$tpm[sum_by_species$classification %in%
                                               paste0("intergenic_", gaussian$selectedGaussianIntergenic[gaussian$speciesId == species])])
    print(paste0("Left side....", max_intergenic))
  } else if (gaussian$selectionSideIntergenic[gaussian$speciesId == species] == "Right") {
    max_intergenic <- min(sum_by_species$tpm[sum_by_species$classification %in%
                                               paste0("intergenic_", gaussian$selectedGaussianIntergenic[gaussian$speciesId == species])])
    print(paste0("Right side....", max_intergenic))
  }
  return(max_intergenic)
}

## function to calculate qValue and than make the adjust call
libraryAdjCall <- function(libraryFile, cutoff, sum_by_species, max_intergenic){
  
  cutoff <- as.numeric(cutoff)
  
  readCallFile <- fread(libraryFile)
  ## density of genic and intergenic region
  selected_genic <- readCallFile$type %in% "genic"
  ## select just reference intergenic
  selected_intergenic <- readCallFile$type %in% "intergenic" & sum_by_species$tpm <= max_intergenic
  
  #dens <- density(log2(na.omit(readCallFile$tpm) + 10^-6))
  dens_genic <- density(log2(readCallFile$tpm[selected_genic] + 10^-6))
  #dens_genic$y <- dens_genic$y * sum(selected_genic) / length(readCallFile$tpm)
  dens_intergenic <- density(log2(readCallFile$tpm[selected_intergenic] + 10^-6))
  #dens_intergenic$y <- dens_intergenic$y * sum(selected_intergenic) / length(readCallFile$tpm)
  
  ## linear interpolation for each type
  genicRegion <- approxfun(dens_genic$x, dens_genic$y)
  intergenicRegion <- approxfun(dens_intergenic$x, dens_intergenic$y)
  ## numerical integration
  numInt_geniRegion <- integrate(genicRegion, min(dens_genic$x), max(dens_genic$x), subdivisions=1000)$value
  numInt_intergenicRegion <- integrate(intergenicRegion, min(dens_intergenic$x), max(dens_intergenic$x), subdivisions=1000)$value
  
  ## for each TPM value collect the genic and intergenic linear interpolation
  interpolationInfo <- c()
  for (i in 1:nrow(readCallFile)) {
    log2TpmValue <- log2(readCallFile$tpm[i])
    genicY <- genicRegion(log2TpmValue)
    intergenicY <- intergenicRegion(log2TpmValue)
    ## creat a info table
    geneInfo <- c(log2TpmValue, genicY,intergenicY)
    interpolationInfo <- rbind(interpolationInfo, geneInfo)
  }
  interpolationInfo <- data.frame(interpolationInfo)
  colnames(interpolationInfo) <- c("log2TpmValue", "genicY","intergenicY")
  interpolationInfo$geneId <- readCallFile$gene_id
  
  ## calculate the qValue for each geneId in the library
  qValueInfo <- c() 
  
  for (i in 1:nrow(interpolationInfo)) {
    
    log2TpmValue <- interpolationInfo$log2TpmValue[i]
    genicY_info <- interpolationInfo$genicY[i]
    intergenicY_info <- interpolationInfo$intergenicY[i]
    
    ## attribute qValue = 1 to inf log2TPM values and to values with Na to genicY_info and intergenicY_info (peak on left side of the plot)
    if ( log2TpmValue == "-Inf" | log2TpmValue != "-Inf" & is.na(genicY_info) == TRUE & is.na(intergenicY_info) == TRUE){
      qValue <- "1"
      qValueInfo <- rbind(qValueInfo,qValue)
      
    } else if (log2TpmValue != "-Inf" & (is.na(intergenicY_info) == TRUE & is.na(genicY_info) == FALSE) | (is.na(intergenicY_info) == FALSE & is.na(genicY_info) == TRUE)){
      
      ## retrieve log2TPM value where was possible to calculate the linear interpolation for genic and intergenic
      maxNumIntegration <- dplyr::filter(interpolationInfo, genicY != "NaN" & intergenicY != "NaN" )
      
      ## if log2TPM value is a negative value (means peak more at left side and we have value to genicY_info and Na to intergenicY_info)
      ## we attribute qValue based on the minimum where was possible to calculate the linear interpolation for both of the types
      if(log2TpmValue < 0){
        
        log2TpmValue <- min(maxNumIntegration$log2TpmValue)
        ## integrate for genic and intergenic region
        unscaled_genic <- integrate(genicRegion, log2TpmValue, max(dens_genic$x), subdivisions=1000, stop.on.error = FALSE)$value
        scaled_genic <- unscaled_genic / numInt_geniRegion
        unscaled_intergenic <- integrate(intergenicRegion, log2TpmValue, max(dens_intergenic$x), subdivisions=1000, stop.on.error = FALSE)$value
        scaled_intergenic <- unscaled_intergenic / numInt_intergenicRegion
        ## calculate qValue for target gene
        qValue <- scaled_intergenic / (scaled_intergenic + scaled_genic)
        qValueInfo <- rbind(qValueInfo,qValue) 
        
      } else {
        
        ## if log2TPM value is positive and the genicY_info is numeric and intergenicY_info is NA
        ## this represent the values in the tail on the right side of the plot.
        ## we attribute qValue for this cases based on the last log2TPM value where was possible to calculate the linear interpolation for both of the types
        log2TpmValue <- max(maxNumIntegration$log2TpmValue)
        ## integrate for genic and intergenic region
        unscaled_genic <- integrate(genicRegion, log2TpmValue, max(dens_genic$x), subdivisions=1000, stop.on.error = FALSE)$value
        scaled_genic <- unscaled_genic / numInt_geniRegion
        unscaled_intergenic <- integrate(intergenicRegion, log2TpmValue, max(dens_intergenic$x), stop.on.error = FALSE, subdivisions=1000)$value
        scaled_intergenic <- unscaled_intergenic / numInt_intergenicRegion
        ## calculate qValue for target gene
        qValue <- scaled_intergenic / (scaled_intergenic + scaled_genic)
        qValueInfo <- rbind(qValueInfo,qValue)
      }
    } else {
      
      ## integrate for genic and intergenic region for a particular log2TPM and than calculate qValue
      ## note in some really particular cases the argument stop.on.error is used in the function integrate().
      unscaled_genic <- integrate(genicRegion, log2TpmValue, max(dens_genic$x), subdivisions=1000, stop.on.error = FALSE)$value
      scaled_genic <- unscaled_genic / numInt_geniRegion
      unscaled_intergenic <- integrate(intergenicRegion, log2TpmValue, max(dens_intergenic$x), subdivisions=1000, stop.on.error = FALSE)$value
      scaled_intergenic <- unscaled_intergenic / numInt_intergenicRegion
      ## calculate q-Value
      qValue <- scaled_intergenic / (scaled_intergenic + scaled_genic)
      qValueInfo <- rbind(qValueInfo,qValue)
    }
  }
  ## add qValue and adj_call columns to the final table
  readCallFile$qValue <- as.numeric(qValueInfo[,1])
  readCallFile$adjusted_call <- ifelse(readCallFile$qValue <= cutoff , "present", "absent")
  return(readCallFile)
}

## function to retrieve information about all libraries before and after the adj of the calls
collectInfoAdjCalls <- function(adjCalls, annotation, libraryId){
  
  genic <- nrow(dplyr::filter(adjCalls, type == "genic" ))
  proteinCoding <- nrow(dplyr::filter(adjCalls, biotype == "protein_coding" ))
  intergenic <- nrow(dplyr::filter(adjCalls, type == "intergenic" ))
  
  ## % genic before and after correction calls
  genicPresent_call <- nrow(dplyr::filter(adjCalls, type == "genic" & call == "present"))
  genicInfo_call <-  (genicPresent_call/genic)*100 
  genicPresent_adjCall <- nrow(dplyr::filter(adjCalls, type == "genic" & adjusted_call == "present"))
  genicInfo_adjCall <-  (genicPresent_adjCall/genic)*100 
  ## % Protein coding before and after correction calls
  proteinCodingPresent_call <- nrow(dplyr::filter(adjCalls, biotype == "protein_coding" & call == "present"))
  proteinCodingInfo_call <- (proteinCodingPresent_call/proteinCoding)*100 
  proteinCodingPresent_adjCall <- nrow(dplyr::filter(adjCalls, biotype == "protein_coding" & adjusted_call == "present"))
  proteinCodingInfo_adjCall <- (proteinCodingPresent_adjCall/proteinCoding)*100 
  ## % Intergenic before and after correction calls
  intergenicPresent_call <- nrow(dplyr::filter(adjCalls, type == "intergenic" & call == "present"))
  intergnicInfo_call <- (intergenicPresent_call/intergenic)*100 
  intergenicPresent_adjcall <- nrow(dplyr::filter(adjCalls, type == "intergenic" & adjusted_call == "present"))
  intergnicInfo_adjcall <- (intergenicPresent_adjcall/intergenic)*100
  
  tpmCutoff <-  min(adjCalls$tpm[adjCalls$adjusted_call == "present"])
  speciesId <- annotation$speciesId[annotation$libraryId == libraryId]
  organism <- as.character(unique(annotation$organism[annotation$libraryId == libraryId]))
  infoLibrary <- c(libraryId, genic, proteinCoding, intergenic, genicInfo_call, genicInfo_adjCall, proteinCodingInfo_call, proteinCodingInfo_adjCall,
                   intergnicInfo_call, intergnicInfo_adjcall, tpmCutoff, speciesId, organism)
}

infoFileLibraries <- file.path(RNASeq_folder_calls, "Info_File_for_all_libraries_after_calls_adjustment.tsv")
if (!file.exists(infoFileLibraries)){
  file.create(infoFileLibraries)
  cat("libraryId\tNumber_genic\tNumber_proteinCoding\tNumber_intergenic\tPercentage_genicInfo_call\tPercentage_genicInfo_adjCall\tPercentage_proteinCodingInfo_call\tPercentage_proteinCodingInfo_adjCall\tPercentage_intergnicInfo_call\tPercentage_intergnicInfo_adjcall\ttpmCutoff_after_adj\tspeciesId\torganism\n", file = file.path(RNASeq_folder_calls, "Info_File_for_all_libraries_after_calls_adjustment.tsv"), sep = "\t")
} else {
  print("File already exist.....")
}

## apply libraryAdjCall() to all libraries and retrieve info from all libraries
for (library in unique(annotation$libraryId)) {
  
  pathLibrary <- file.path(RNASeq_folder_calls, library) 
  libraryCalls <- list.files(pathLibrary,  pattern="abundance_gene_level\\+fpkm\\+intergenic\\+calls.tsv$", full.names=T, recursive = TRUE)
  
  if(file.exists(pathLibrary) & length(libraryCalls) != 0){
    nameFile <- basename(libraryCalls)
    ## retrieve just reference intergenic per species
    species <- annotation$speciesId[annotation$libraryId == library]
    sum_by_species <- fread(paste0(sum_by_species_folder, "/sum_abundance_gene_level+fpkm+intergenic+classification_",species,".tsv"))
    max_intergenic_Value <- max_intergenic(sum_by_species = sum_by_species, gaussian = gaussian, species = species)
    
    ## adj the calls
    callsAdjustment <- libraryAdjCall(libraryFile = libraryCalls, cutoff = cutoff, sum_by_species = sum_by_species, max_intergenic = max_intergenic_Value)
    ## retrieve info
    infoLibrary <- collectInfoAdjCalls(adjCalls = callsAdjustment, annotation = annotation, libraryId = library)
    ## re-write the calls file with extra information: qValue + adj_call
    write.table(callsAdjustment, file = paste0(pathLibrary, "/", nameFile), row.names = F, col.names = T, quote = FALSE, sep = "\t")
    ## write info file
    write.table(t(infoLibrary), file = infoFileLibraries, row.names = F, col.names = F, quote = FALSE, append = T, sep = "\t")
    cat("DONE library ", library, "\n")
  } else if (file.exists(pathLibrary) & length(libraryCalls) == 0){
    cat("File with calls not found for the library ", library , "\n")
  } else{
    cat("Folder not exist for this library", library , "\n")
  }
}
