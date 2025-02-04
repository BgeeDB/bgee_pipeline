## SFonsecaCosta, May 2019

## This script is used to do the Sum of raw counts for scRNA-Seq data.
## Mainly this sum all est_counts for each transcript ID per cell population that belongs to the same: experiment, species, cellTypeId, stageId, strain, uberonId, sex, followed by the calculation of the effective length by using the weighted.mean

## Usage:
## R CMD BATCH --no-save --no-restore '--args NEW_scRNASeq_sample_info="NEW_scRNASeq_sample_info.tsv" cells_folder="cells_folder" output_folder="output_folder"' Sum_RawCounts_cellPopulation.R Sum_RawCounts_cellPopulation.Rout
## NEW_scRNASeq_sample_info --> info of samples per library/experiment/species.
## cells_folder --> where we is located all the libraries/cells after Kallisto (treated data)
## output_folder --> folder where we should save each new abundance.tsv file that belong to each different cell type and experimemnt.

## Libraries used
library(plyr)
library(dplyr)
library(reshape2)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed in command line
command_arg <- c("NEW_scRNASeq_sample_info", "cells_folder", "output_folder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read NEW_scRNASeq_sample_info file. If file not exists, script stops
if( file.exists(NEW_scRNASeq_sample_info) ){
  annotation <- read.table(NEW_scRNASeq_sample_info, h=T, sep="\t", comment.char="")
  names(annotation)[1] <- "libraryId"
  ## remove special characters in strain
  annotation$strain <- gsub('\\/', '_', annotation$strain)
  annotation$strain <- gsub('\\s+', '_', annotation$strain)
} else {
  stop( paste("NEW_scRNASeq_info_file file not found [", NEW_scRNASeq_sample_info, "]\n"))
}
##########################################################################################################################################
## create the output file -> sum_scRNASeq_sample_info.txt
sum_infoFile <- file.path(output_folder, "sum_scRNASeq_sample_info.txt")
if (!file.exists(sum_infoFile)){
  file.create(sum_infoFile)
} else {
  print("File already exist.....")
}
cat("libraryId\texperimentId\tcellTypeName\tcellTypeId\tspeciesId\tplataform\twhiteList\tprotocol\tprotocolType\tlibraryType\tinfoOrgan\tstageId\tuberonId\tsex\tstrain\treadLength\torganism\n",
    file = file.path(output_folder, "sum_scRNASeq_sample_info.txt"), sep = "\t")

collectInfoSum <- function(scrna_seq_sample_infoFile){
  libraryId <- paste0(experiment, "_", cellId, "_", stageId, "_", strain, "_", uberonId, "_", sex)
  experimentId <- paste0(experiment, "_", cellId, "_", stageId, "_", strain, "_", uberonId, "_", sex)
  cellTypeName <- as.character(unique(annotation$cellTypeName[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex]))
  cellTypeId <- paste0(cellId)
  speciesId <- species
  plataform <- as.character(unique(annotation$platform[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex]))
  whiteList <- as.character(unique(annotation$whiteList[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex]))
  protocol <- as.character(unique(annotation$protocol[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex]))
  protocolType <- as.character(unique(annotation$protocolType[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex]))
  libraryType <- as.character(unique(annotation$libraryType[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex]))[1]
  infoOrgan <- as.character(unique(annotation$infoOrgan[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex]))
  stageId <- paste0(stageId)
  uberonId <- paste0(uberonId)
  sex <- paste0(sex)
  strain <- paste0(strain)
  readLength <- "NA"
  organism <- as.character(unique(annotation$organism[annotation$speciesId == species]))

  ## Export information of the sum
  collectInfo <- c(libraryId, experimentId, cellTypeName, cellTypeId, speciesId, plataform,
                   whiteList, protocol, protocolType, libraryType,  infoOrgan, stageId,
                   uberonId, sex, strain, readLength, organism)
  names(collectInfo) <- c("libraryId", "experimentId", "cellTypeName", "cellTypeId", "speciesId", "plataform",
                          "whiteList", "protocol", "protocolType", "libraryType", "infoOrgan", "stageId",
                          "uberonId", "sex", "strain", "readLength", "organism")
  return(collectInfo)
}

## Loop through the cell-type/Experiments to sum all the Raw est_counts per cell type ....
for (species in unique(annotation$speciesId)) {
  cat("Species:", species, "\n")
  for (experiment in unique(annotation$experimentId[annotation$speciesId == species])){
    cat("Name of experiments that belongs to the species:", experiment, "\n")
    for (cellId in unique(annotation$cellTypeId[annotation$experimentId == experiment])){
      cat("CellId info:", cellId, "\n")
      for (stageId in unique(annotation$stageId[annotation$cellTypeId == cellId])){
        cat("StageId info:", stageId, "\n")
        for (strain in unique(annotation$strain[annotation$stageId == stageId])){
          cat("Strain info:", strain, "\n")
          for (uberonId in unique(annotation$uberonId[annotation$strain == strain])){
            cat("UberonId info:", uberonId, "\n")
            for (sex in unique(annotation$sex[annotation$uberonId == uberonId])){
              cat("sex info:", sex, "\n")


              infoLib <- annotation$libraryId[annotation$speciesId == species & annotation$experimentId == experiment & annotation$cellTypeId == cellId & annotation$stageId == stageId & annotation$strain == strain & annotation$uberonId == uberonId & annotation$sex == sex]
              ## print information
              print(length(infoLib))
              if (length(infoLib) != 0){

                file <- file.path(cells_folder, infoLib)
                AllFiles <- list.files(file, pattern="abundance.tsv", full.names=T, recursive = TRUE)
                print(AllFiles)

                AllFiles <- lapply(AllFiles, read.delim)
                DATA <- do.call("cbind", AllFiles)
                ## select target info and length
                select_info <- DATA[,1:2]

                ## select all columns with -> est_count and obtain the sum
                select_est_count = DATA[, grepl("^est_count", names( DATA))]
                select_est_count_sum <- apply(select_est_count, 1, sum)
                ## select all columns with -> eff_length
                select_effec_length = DATA[, grepl("^eff_length", names( DATA))]

                ## calculate effect_lengh by using weighted.mean for each transcrip...
                calculate_effect_length <- function(counts, effec_length){
                  myWeightedMean <- c()
                  ## transpose the matrix (this means each column is a transcript ID and each row a cell)
                  rawcounts <- t(counts)
                  raw_effeclength <- t(effec_length)

                  for (i in 1:ncol(raw_effeclength)) {
                    if (sum(rawcounts[,i]) == 0 ){
                      ## provide the same weight for all transcripts in case the est_count is always zero
                      rawcounts[,i] <- ifelse(rawcounts[,i] == 0, 1)
                      weightedMean <-  weighted.mean(x=raw_effeclength[,i], w=rawcounts[,i])
                      myWeightedMean[i] <- weightedMean
                    } else {
                      weightedMean <-  weighted.mean(x=raw_effeclength[,i], w=rawcounts[,i])
                      myWeightedMean[i] <- weightedMean
                    }
                  }
                  return(myWeightedMean)
                }

                ## export result where each row is the weighted.mean of eff_length of each transcript
                select_effec_length_Wmean <- as.data.frame(calculate_effect_length(select_est_count, select_effec_length))
                colnames(select_effec_length_Wmean) <- "select_effec_length_Wmean"
                ## re-build the information by using weighted.mean...
                final_Final_weightMean <- data.frame(select_info, select_effec_length_Wmean, select_est_count_sum)

                ## calculate TPM after collect the weighted.mean of the eff_length and sum of est_count
                estCount_to_tpm <- function(est_count, effec_length){
                  rate <- log(est_count) - log(effec_length)
                  denom <- log(sum(exp(rate)))
                  exp(rate - denom + log(1e6))
                }
                tpmValues_wMean <- estCount_to_tpm(final_Final_weightMean$select_est_count_sum, final_Final_weightMean$select_effec_length_Wmean)

                ## export the final file by using weighted.mean of effective_Lenght....
                AllCells_sum_wMean <- data.frame(final_Final_weightMean, tpmValues_wMean)
                colnames(AllCells_sum_wMean) <- c("target_id","length", "eff_length","est_counts", "tpm")
                ## create a sub-directory for each experimentID
                output_dir <- paste0(output_folder, "/", experiment, "_", cellId, "_", stageId, "_", strain, "_", uberonId, "_", sex)

                if (!dir.exists(output_dir)){
                  dir.create((output_dir))
                } else {
                  print("Directory already exists.....")
                }

                write.table(AllCells_sum_wMean, file = file.path(output_dir, "abundance.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
                collectInfo <- collectInfoSum(scrna_seq_sample_infoFile=annotation)
                write.table(t(collectInfo), file = sum_infoFile,col.names =F , row.names = F, append = T,quote = FALSE, sep = "\t")

              } else {
                cat("Combination with zero cells!")
              }
            }
          }
        }
      }
    }
  }
}

