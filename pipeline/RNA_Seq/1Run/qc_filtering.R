library(dplyr)

#reports from kallisto for all libraries to insert
reportsFile = "../../generated_files/RNA_Seq/reports_info_all_samples.txt"
reportsKallistoLibs <- read.table(file = reportsFile, header = T, sep = "\t")
#sample info file used by the pipeline
sampleFile = "../../generated_files/RNA_Seq/rna_seq_sample_info.txt"
sampleInfo <- read.table(file = sampleFile, header = T, sep = "\t", quote = "", comment.char = "!")
#file summarizing present/absent calls proportions
presenceAbsenceFile = "../../generated_files/RNA_Seq/presence_absence_all_samples.txt"
presenceAbsence <- read.table(file = presenceAbsenceFile, header = T, sep = "\t")
#annotation file coming from Bgee curators
libAnnot <- read.table(file = "../../generated_files/RNA_Seq/RNASeqLibrary_full.tsv", header = TRUE, sep = "\t", quote = "\"")
#update column name because R does not properly manage the # as 1st character of the header
colnames(sampleInfo)[1] <- "libraryId"
colnames(libAnnot)[1] <- "libraryId"

#security check that the calls have been generated for all libraries kallisto was run on
if(nrow(presenceAbsence) != nrow(reportsKallistoLibs)) {
  stop("kallisto report file and calls summary file should contain the same number of libraries")
}

#remove from kallisto reports all libraries with less than than 100'000 reads or less than 5% of reads mapped to the transcriptome
filteredReportsKallistoLibs <- reportsKallistoLibs %>% filter(reads > 100000 & prop_aligned > 5)
message("filtering on number of reads and proportion of reads mapped removed ",
  nrow(reportsKallistoLibs)-nrow(filteredReportsKallistoLibs), " libraries.")

#now remove all libraries with less than 13% of protein coding gene present from the file summarizing calls proportion
#to do so we first consider only libraries that target protein coding genes
nonProtCodingLibAnnot <- libAnnot %>% filter(RNASelection %in% c("lncRNA", "miRNA", "circRNA"))
nonProtCodingPresenceAbsence <- presenceAbsence %>% filter(libraryId %in% nonProtCodingLibAnnot$libraryId)
filteredPresenceAbsence <- presenceAbsence %>% filter(!libraryId %in% nonProtCodingPresenceAbsence$libraryId & proportionCodingPresent > 13)

# now apply both filtering to all files
finalReportsKallistoLibs <- filteredReportsKallistoLibs %>% filter(libraryId %in% nonProtCodingLibAnnot$libraryId |
  libraryId %in% filteredPresenceAbsence$libraryId)
message("filtering on proportion of protein coding gene present removed ",
        nrow(filteredReportsKallistoLibs)-nrow(finalReportsKallistoLibs), " libraries.")
finalPresenceAbsence <- presenceAbsence %>% filter(libraryId %in% finalReportsKallistoLibs$libraryId)

#summary of libraries removed by bith filtering
removedLibraries <-  reportsKallistoLibs %>% filter(!libraryId %in% finalReportsKallistoLibs$libraryId)
removedLibrariesWithSpeciesAndExp <- merge(x = removedLibraries, y = sampleInfo[,c("libraryId", "experimentId", "speciesId")], by = "libraryId", all.x = T)
message("Summary of libraries removed applying both filtering.")
removedLibrariesWithSpeciesAndExp %>% group_by(speciesId) %>% summarize (count = n()) %>% arrange(count)

#remove qc filtered libraries from the sample_info file
finalSampleInfo <- sampleInfo %>% filter(!libraryId %in% removedLibraries$libraryId)

#update read length of sample_info file for library that have no value for that field
absent <- 0
for(i in seq_len(nrow(finalSampleInfo))) {
  if (is.na(finalSampleInfo$readLength[i]) || finalSampleInfo$readLength[i] == 0) {
    print(finalSampleInfo$libraryId[i])
    if (finalSampleInfo$libraryId[i] %in% finalReportsKallistoLibs$libraryId) {
      readLength <- as.numeric(finalReportsKallistoLibs$max_read_size[finalReportsKallistoLibs$libraryId == finalSampleInfo$libraryId[i]])
      finalSampleInfo$readLength[i] <- if_else(finalSampleInfo$libraryType[i] == "SINGLE", readLength, readLength*2)
      absent <- absent + 1
      finalSampleInfo$readLength[i] <- readLength
    }
  }
}
message("read length has been corrected for ", absent, " libraries")

#save final version of the files
write.table(x = finalReportsKallistoLibs, file = reportsFile, quote = FALSE, sep = "\t", col.names = TRUE,
            row.names = FALSE)
write.table(x = finalPresenceAbsence, file = presenceAbsenceFile, quote = FALSE, sep = "\t", col.names = TRUE,
            row.names = FALSE)
write.table(x = finalSampleInfo, file = sampleFile, quote = FALSE, sep = "\t", col.names = TRUE,
            row.names = FALSE)
# save libraries removed by the filtering
write.table(x = removedLibrariesWithSpeciesAndExp[,c("libraryId", "experimentId", "speciesId")],
  file = "../../generated_files/RNA_Seq/libraries_qc_filtered.tsv", quote = FALSE, sep = "\t", col.names = TRUE,
  row.names = FALSE)

#remove libraries that are not protein coding from the libraries to use to generate the plot
presenceAbsencePlot <- finalPresenceAbsence %>% filter(!libraryId %in% nonProtCodingLibAnnot$libraryId)

all_samples <- presenceAbsencePlot
pdf(file = "~/Downloads/presence_absence_boxplots.pdf", width = 12, height = 5)
## 7 boxplots to plot, 1 per column in all_samples
## Manually set the y-scale limits for most of the columns. NA are here for columsn thta need to be printed in log scale
ylimits <- list(NA, c(0, 100), c(0, 70000), c(0, 70000), c(0, 100), c(0, 30000), c(0,30000),  c(0, 100), c(0, 30000), c(0,30000), c(0,0.3))

i = 1 
for (column in c(2:12)){
  print(names(all_samples)[column])
  
  ## get the subset of data to use
  subsetNames <- names(all_samples)[column]
  subsetData <- as.numeric(as.character(all_samples[, column]))
  organism <- as.factor(as.character(all_samples$speciesId))
  
  ## reoder organism levels to get boxes ordered by median values in boxplot
  organism <- reorder(organism, subsetData, median)
  
  par(mar=c(12,4,2,1)) ## bottom, left, top and right margins
  if (is.na(ylimits[[i]][1])){
    boxplot(log2(subsetData+10^-6) ~ organism,
            pch=16, notch=T, outcex=0.5, boxwex=0.7,
            ylab=paste0("log2(", subsetNames, " + 10e-6)"),
            main="", las=2)
    text(x=1:length(unique(organism)),
         y=log2(max(subsetData)+10^-6),
         labels=table(organism),
         col="black", cex=0.8)
  } else {
    boxplot(subsetData ~ organism,
            pch=16, notch=T, outcex=0.5, boxwex=0.7,
            ylab=subsetNames,
            main="", las=2,
            ylim=ylimits[[i]])
    text(x=1:length(unique(organism)),
         y=ylimits[[i]][2],
         labels=table(organism),
         col="black", cex=0.8)
  }
  abline(v=1:length(unique(organism)), lty=2, col="lightgray")
  if (names(all_samples)[column] == "proportionCodingPresent"){
    abline(h=c(70, 80), lty=2, col="red")
  }
  i = i+1
}
## close PDF device
dev.off()
