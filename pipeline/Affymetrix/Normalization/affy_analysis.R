## Julien Roux 15/05/08  ##Marta Rosikiewicz 18/06/12

## options(echo=FALSE)

## First read in the arguments listed at the command line
args=(commandArgs(TRUE))

## args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  q()
} else {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

## source("http://www.bioconductor.org/biocLite.R")
## biocLite()
library("affy")
library("genefilter")
library("gcrma")
library("data.table")
library("fdrtool")

## adding additional path to variable specifying where R is looking for packages
## this code allows to install packages without root access on the server
Rlib_path<-".."
Rlib_folder<-"Rlib_folder"
Rlib_folder_path<-file.path(Rlib_path,Rlib_folder)
if(!file.exists(Rlib_folder_path)){
    if(!dir.create(Rlib_folder_path)){stop("can't create Rlib_folder_path")}
}

.libPaths(c(Rlib_folder_path,.libPaths()))



## filenames=dir(celpath)

cat("Reading files:\n",paste(filenames,collapse="\n"),"\n")
##read one .CEL file and compute affinities
data <- ReadAffy(filenames=filenames, celfile.path=celpath)

###################
## normalization ##
###################

## compute expression measure by gcRMA
## data.gcrma<-gcrma(data, type="affinities")

## test if the affinity info is not already computed
affinities <- dir(affin)

if (sum(unlist(strsplit(affinities, ".RData")) == array) == 1){
  load(paste(affin,array,".RData",sep=""))
  print("Affinities loaded")
} else {
  ai <- compute.affinities(cdfName(data))
  ## ai <- compute.affinities.local(data, Array=NULL)
  save(ai,file=paste(affin,array,".RData",sep=""))
  print("Affinities calculated and saved")
}

## From: Zhijin Wu
## The "fullmodel" option uses both MM intensities and PM probe sequences;
## The  "affinities" option uses MM probes only as negative control probes to
## estimate the relationship of NSB with sequence. IT does not use the MMs as
## paried measurement of background for each PM probe. The idea here is that
## a set of "non-PM" probes are needed but not one for each PM.
## The "MM" option uses only the MM measurements for background adjustment.

## gcrma normalization
## in case there is only one celfile the quantile normalization is turned off
if(length(filenames)==1){
    ##gcrma background correction
    ##no quantile nomrlization in case of one cel file
    ##GSB.adjust=FALSE: no comparisons between chips in case of only  cel file
    data_bg_cor<-bg.adjust.gcrma(data,affinity.info=ai,type="affinities",GSB.adjust=FALSE)
    data.gcrma<-rma(data_bg_cor, subset=NULL, verbose=TRUE, destructive=TRUE, normalize=FALSE,
    background=FALSE)
}else {
    data.gcrma <- gcrma(data, affinity.info=ai, type="affinities")
}


################################
## detection calls (Schuster) ##
################################

## find probesets that are very likely to have absent target transcripts
## probesets that are very likely to have an absent target transcript
## are defined as having MAS5 present/absent P values for all replicates > 0.50
raw.calls <- mas5calls(data)

f1 <- kOverA(length(sampleNames(data)), 0.50)
ff1 <- filterfun(f1)
ab.50 <- genefilter(assayData(raw.calls)[["se.exprs"]], ff1)
gn <- rownames(assayData(raw.calls)[["se.exprs"]])
ab.50 <- gn[ab.50]

## GC-RMA PM probes using affinities specific to each replicate sample
## requires "ai" object previously computed
## ai <- compute.affinities.local(data, Array=NULL)
## pipepline adjusted for experiments with single celfile
if(length(filenames)==1){
    ##GSB.adjust=FALSE: no comparisons between chips in case of only  cel file
    data.PM <- bg.adjust.gcrma(data, affinity.info=ai, type="fullmodel", GSB.adjust=FALSE,fast=FALSE)
}else{
    data.PM <- bg.adjust.gcrma(data, affinity.info=ai, type="fullmodel", fast=FALSE)
}

## replace MM probe values with a PM threshold value
for(i in 1:length(sampleNames(data))){
mm(data.PM)[,i] <- mean(pm(data.PM, ab.50)[,i], trim=0.02)}

## calculate present/absent calls
## alpha values set the "present" and "marginal" P value cutoffs
## alpha1 is set to 0.03 and the expected false discover rate at this P value is roughly 5% to 8%
## alpha1 is set to 0.12 and the expected false discover rate at this P value is roughly 10% to 13%

thres.calls <- mas5calls(data.PM, tau=0.015, alpha1 = 0.024, alpha2 = 0.111 )

################################################################################################
## Output table with expression values + pValues + calls + qValues + adjustment of the calls  ##
################################################################################################
qValue_cutoff <- "0.01"

for (i in 1:length(filenames)){
  # probeset names are taken from the row.names of exprs(data.gcrma)[,i] & exprs(thres.calls)[,i]
  normalizedTable <- cbind(exprs(data.gcrma)[,i],assayData(thres.calls)[["se.exprs"]][,i],exprs(thres.calls)[,i])
  normalizedTable <- data.frame(as.numeric(normalizedTable[,1]), as.numeric(normalizedTable[,2]),normalizedTable[,3])
  colnames(normalizedTable) <- c("expression", "pValue", "call")
  normalizedTable <- setDT(normalizedTable, keep.rownames = TRUE)[]
  colnames(normalizedTable)[1] <- "probeId"
  
  ## select probes with value of expression == minimum value of expression AND probes with value of expression != minimum value of expression
  minExpression <- min(normalizedTable$expression)
  probeMinExp <- dplyr::filter(normalizedTable, expression == minExpression)
  probesets <- dplyr::filter(normalizedTable, expression != minExpression)

  ## calculate q-value just using the p-values from probes that don't have a value of expression equal to the minimum of expression
  qValue <- fdrtool(probesets$pValue, statistic="pvalue", plot = FALSE)
  probesets <-cbind(probesets, qValue$qval)
  colnames(probesets)[5] <- "qValue"
  
  ## attribute qValue to probes with minimum value of expression
  probeMinExp$qValue <- 1
  
  ## final table with all information
  finalTable <- rbind(probesets, probeMinExp)
  finalTable$adjusted_call <- ifelse(finalTable$qValue <= qValue_cutoff, "P", "A")
  write.table(finalTable, file=paste(processed,filenames[i],".out",sep=""), sep="\t",quote=F,col.names=T,row.names=F)
}

