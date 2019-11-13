# Marta Rosikiewicz
# created: 30/10/2012
# modified: 21/10/2013: removed intronic regions from analysis; size of intergenic regions limited to min 2000 max 20000
# modified 13/05/2014: the script was corrected to work with new version of ENSEMBL gtf format

# Julien Roux
# modified 10/12/2015: updated and simplified script for better readability; checked compatibility with kallisto-based pipeline; removed length calculation and export; added export of table to map transcript IDs to gene IDs; modified to read directly compressed GTF file
# modified 22/02/2016: clarified and simplified even more the script (removed variables used only once); added export of gene biotypes table; fixed small bug in coordinates of intergenic region: because of using round(), there were too short by 1bp. Now using ceiling() and floor()

# Julien Wollbrett
# modified 07/01/2019: remove potential scientific notation for start or stop position of intergenic regions

## From GTF file downloaded from Ensembl, this script prepares a new GTF file (gtf_all) with:
## - exonic regions from all transcripts of each gene
## - intergenic regions. The 500nt flanking genes are excluded (the minimal distance from start or stop of intergenic region to boundary of the nearest gene: 500nt). The intergenic region length is limited to min 2000 max 20000. Intergenic regions with lower size are discarded. Larger regions are limited to +/- 10000 around the center of the intergenic region.
## A conversion table (transcript ID to gene ID) is exported: extension .gene2transcript
## A gene biotype table is exported: extension .gene2biotype

## Invoking:
# R CMD BATCH --no-save --no-restore '--args gene_gtf_path="gene_gtf_path" output_gtf_path="output_gtf_path"' prepare_gtf.R Rout_path

## Example:
# R CMD BATCH --no-save --no-restore '--args gene_gtf_path="~/Desktop/RNAseq/pipeline_Ensembl_73/Mus_musculus.GRCm38.73.gtf.gz" output_gtf_path="~/Desktop/RNAseq/pipeline_Ensembl_73/gtf_folder/Mus_musculus.GRCm38.73"' prepare_GTF.R

## Arguments to provide:
# "gene_gtf_path" - full path to input gene gtf file
# "output_gtf_path" - full path to output folder + base name for output files
# Rout_path - path to  output .Rout file (optional)

##################################################################################
## adding additional path to variable specifying where R is looking for packages
## this code allows to install packages without root access on the server
Rlib_path <- ".."
Rlib_folder <- "Rlib_folder"
Rlib_folder_path <- file.path(Rlib_path,Rlib_folder)
if(!file.exists(Rlib_folder_path)){
  dir.create(Rlib_folder_path)
}
.libPaths(c(Rlib_folder_path,.libPaths()))

library("GenomicFeatures")
library("chipseq")
## TODO Are these packages really needed? We seem to use functions from packages that depend on them, not from them

## Session info
print(sessionInfo())

## reading in command line arguments
cmd_args = commandArgs(TRUE);
print ("command arguments\n")
print (cmd_args)

if(length(cmd_args)==0){
  stop("no arguments provided")
} else {
  for(i in 1:length(cmd_args)){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessarily arguments were provided properly
if (!exists("gene_gtf_path")){ stop("gene_gtf_path not defined") }
if (!exists("output_gtf_path")){ stop("output_gtf_path not defined") }

########################################################################
## Function for obtaining the part of the annotation field from gtf file
## Input: whole field, already splitted, for example: gene_id "FBgn0264003"; gene_name "mir-5613"; gene_source "FlyBase"; gene_biotype "pre_miRNA"
## Output: could be for example gene_id: FBgn0264003
get_annot_value <- function(split_annotation, field_name){
  ## find the right field
  field_all <- split_annotation[grep(field_name, split_annotation, fixed=T)];

  ## split the field
  field_value <- strsplit(field_all, ' ', fixed=T)[[1]][2];

  ## remove the last ';' if necessairy
  field_value <- sub(';', '', field_value,fixed=T)

  return(field_value)
}

###################################################################
## reading in gtf file (gzipped, no need to uncompress)
cat("Reading GTF file...\n")
gene_gtf <- as.matrix(read.table(file=gzfile(gene_gtf_path, "r"), sep="\t", strip.white=TRUE, as.is=TRUE, colClasses="character", comment.char='#'))
## Header lines starting with # are not read

## selecting exon lines
gene_gtf_exon <- gene_gtf[gene_gtf[,3]=="exon",]
## selecting gene lines
gene_gtf <- gene_gtf[gene_gtf[,3]=="gene",]

##  selecting only genes from assembled chromosome sequence (potentially useful code to select only fully asembled genome sequence)
##  chromosomes_all <- c(1:100, "X", "Y", "MT")
##  gene_gtf_exon <- gene_gtf_exon[gene_gtf_exon[,1] %in% chromosomes_all,]

## getting gene start, stop, chromosome, strand, biotype for each gene (from the exon GTF)
## GTF format reports 1-based coordinates (start and end)
cat("Extracting gene informations...\n")
## splitting the annotation field using single space "; " as a pattern
split_annotation_list <- strsplit(gene_gtf_exon[,9], "; ",  fixed =T)

## getting the vector of the gene IDs (1 for every exon)
gene_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_id') })

## getting the vector of the transcript IDs (1 for every exon)
transcript_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'transcript_id') })

## getting the table with mappings between transcript and gene IDs (for export)
gene_transcript_ids <- unique(cbind(gene_ids, transcript_ids), MARGIN=1)

## getting the vector of gene_biotypes: splitting gene gtf
split_annotation_list <- strsplit(gene_gtf[,9], "; ",  fixed =T)
gene_biotypes <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_biotype') })
names(gene_biotypes) <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_id') })

## getting the chromosome, start and the end of the gene
## For start, take the minimum exon start position
## For end, take the maximum exon end position
gene_start <- sapply(split(as.numeric(gene_gtf_exon[,4]), gene_ids), function(x){ sort(as.numeric(x))[1] })
gene_stop <- sapply(split(as.numeric(gene_gtf_exon[,5]), gene_ids), function(x){ rev(sort(as.numeric(x)))[1] })
gene_chr <- sapply(split(gene_gtf_exon[,1], gene_ids), function(x){ x[1] })

## chromosome/contig names from given gtf files
chromosomes <- unique(gene_gtf_exon[,1])
## removing patch contigs before selecting intergenic regions
chromosomes <- chromosomes[grep('PATCH', chromosomes, invert=TRUE, ignore.case=TRUE)]

###################################################################
## Select the set of intergenic regions
cat("Selecting set of intergenic regions...\n")

## This object will include the coordinates of the selected intergenic regions
intergenic_regions <- matrix(ncol=3, nrow=0)
colnames(intergenic_regions) <- c("chr", "start", "end")

# To avoid scientific notation we change the value of option "scipen" to 999. At the end of the script we will change to its initial value
scipen_initial_value <- getOption("scipen")
options(scipen = 999) 

for(chr in chromosomes){
    cat(chr, " ")
    ## keeping genes from selected chromosome
    ## skip chromosome/contig if 1 or less gene
    if(( sum(gene_chr==as.character(chr)) <= 1 )){ next }

    ## constructing coverage map of the chromosomes using start and stop coordinates of the genes and selecting regions with 0 coverage (intergenic)
    gene_IR <- IRanges(as.numeric(gene_start[gene_chr==as.character(chr)]), as.numeric(gene_stop[gene_chr==as.character(chr)]))
    inter_IR <- slice(coverage(gene_IR), lower=0, upper=0, rangesOnly=TRUE)
    inter_gene_data <- as.data.frame(cbind(inter_IR@start, end(inter_IR), inter_IR@width))
    colnames(inter_gene_data) <- c("start", "end", "width")

    ## if selected chromosome/contig has only no intergenic regions with length min 2000nt then skip this contig
    if( sum(inter_gene_data[,3] > 2000) == 0 ){ next }

    ## restrict to regions larger than 2000nt
    inter_gene_data <- inter_gene_data[inter_gene_data[,3] > 2000,]

    ## finding the center of the intergenic region
    inter_gene_data$center <- round((inter_gene_data[,1]+inter_gene_data[,2])/2)

    ## getting the size of "usable" intergenic regions around the center point
    ## width - 500nt around each flank
    inter_gene_data$size <- inter_gene_data[,3] - 1000
    inter_gene_data$size[inter_gene_data$size > 20000] <- 20000 ## limit the size of usable regions to max 20000

    ## storing information (chr, start, stop, size) for selected intergenic regions on this chromosome
    intergenic_regions <- rbind(intergenic_regions, cbind(chr, inter_gene_data$center - ceiling(inter_gene_data$size/2) + 1 , inter_gene_data$center + floor(inter_gene_data$size/2)))
}
options(scipen = scipen_initial_value)

## preparing intergenic gtf data
cat("\nPreparing intergenic GTF data...\n")
intergenic_regions_gtf <- matrix(ncol=9, nrow=nrow(intergenic_regions))
intergenic_regions_gtf[,1] <- intergenic_regions[,1]
intergenic_regions_gtf[,2] <- "intergenic"
intergenic_regions_gtf[,3] <- "exon"
intergenic_regions_gtf[,c(4,5)] <- intergenic_regions[,c(2,3)]
intergenic_regions_gtf[,c(6,8)] <- "."
## Strand is chosen randomly
set.seed(12) ## setting seed for random number generator
intergenic_regions_gtf[,7] <- sample(c("+", "-"), nrow(intergenic_regions), replace =TRUE)

## intergenic_id - chr "_" start "_" stop
intergenic_id <- apply(intergenic_regions[,1:3], 1, function(x){ paste(x, collapse="_") })
gene_id <- paste("gene_id ", paste(intergenic_id, ";", sep=""), sep="")
transcript_id <- paste("transcript_id ", paste(intergenic_id, ";", sep=""), sep="")
intergenic_regions_gtf[,9] <- apply(cbind(gene_id, transcript_id), 1, function(x){ paste(x, collapse=" ") })

####################################
## Output:
## GTF file with both genic exons and intergenic regions
write.table(rbind(gene_gtf_exon, intergenic_regions_gtf), file=paste(output_gtf_path, ".gtf_all", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Table with mappings between transcript and gene IDs
write.table(gene_transcript_ids, file=paste(output_gtf_path, ".gene2transcript", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Table with biotype for all gene IDs
write.table(gene_biotypes, file=paste(output_gtf_path, ".gene2biotype", sep=""), sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
