# Marta Rosikiewicz
# created: 30/10/2012
# modified: 21/10/2013: removed intronic regions from analysis; size of intergenic regions limited to min 2000 max 20000
# modified 13/05/2014: the script was corrected to work with new version of ENSEMBL gtf format

# Julien Roux
# modified 10/12/2015: updated and simplified script for better readability; checked compatibility with kallisto-based pipeline; removed length calculation and export; added export of table to map transcript IDs to gene IDs; modified to read directly compressed GTF file
# modified 22/02/2016: clarified and simplified even more the script (removed variables used only once); added export of gene biotypes table; fixed small bug in coordinates of intergenic region: because of using round(), there were too short by 1bp. Now using ceiling() and floor()

# Julien Wollbrett
# modified 07/01/2019: Correct a bug of scientific notation for start or stop of intergenic region. Scientific notation is disabled.
# modified 11/12/2018: Remove N from intergenic regions. Need to import fasta file to access to the sequences.

## From GTF file downloaded from Ensembl, this script prepares a new GTF file (gtf_all) with:
## - exonic regions from all transcripts of each gene
## - intergenic regions. The 500nt flanking genes are excluded (the minimal distance from start or stop of intergenic region to boundary of the nearest gene: 500nt). The intergenic region length is limited to min 2000 max 20000. Intergenic regions with lower size are discarded. Larger regions are limited to +/- 10000 around the center of the intergenic region.
## - remove Ns from intergenic sequences. Genome is used to retrieve sequence of previously detected intergenic regions. If the sequence contains more than `block_size` consucutive Ns the sequence is splitted in 2. At the end all intergenic sequences without block of Ns but with a higher proportion of N than `N_proportion` OR shorter than 1000bp are discarded
## - A summary of N removal is provided: extension .Nremoval

## Invoking:
# R CMD BATCH --no-save --no-restore --slave '--args gene_gtf_path = "gene_gtf_path" genome_fasta_path = "genome_fasta_path" N_block_size = 31 N_proportion = 0.05 output_gtf_path = "output_gtf_path"' prepare_gtf.R Rout_path

## Example:
# R CMD BATCH --no-save --no-restore --slave'--args gene_gtf_path = "~/Desktop/RNAseq/pipeline_Ensembl_73/Mus_musculus.GRCm38.73.gtf.gz" genome_fasta_path = "~/Desktop/RNAseq/pipeline_Ensembl_73/Mus_musculus.GRCm38.73.genome.fa" N_block_size = 31 N_proportion = 0.05 output_gtf_path = "~/Desktop/RNAseq/pipeline_Ensembl_73/gtf_folder/Mus_musculus.GRCm38.73"' prepare_GTF.R

## Arguments to provide:
# "gene_gtf_path" - full path to input gene gtf file
# "genome_fasta_path" - full path to input genome fasta file
# "N_block_size" - number of successive N from which it is considered as a block of N (and removed)
# "N_proportion" - higher proportion of N than this threshold results in removing the sequence
# "output_gtf_path" - full path to output folder + base name for output files
# Rout_path - path to  output .Rout file (optional)

##################################################################################
## adding additional path to variable specifying where R is looking for packages
## this code allows to install packages without root access on the server
Rlib_path <- ".."
Rlib_folder <- "Rlib_folder"
Rlib_folder_path <- file.path(Rlib_path,Rlib_folder)
if (!file.exists(Rlib_folder_path)) {
  dir.create(Rlib_folder_path)
}
.libPaths(c(Rlib_folder_path,.libPaths()))

library("GenomicFeatures")
library("Biostrings")

## Session info
print(sessionInfo())

## reading in command line arguments
cmd_args = commandArgs(TRUE);
print("command arguments\n")
print(cmd_args)

if (length(cmd_args) ==  0) {
  stop("no arguments provided")
} else {
  for (i in 1:length(cmd_args)) {
    eval(parse(text = cmd_args[i]))
  }
}

## checking if all necessarily arguments were provided properly
if (!exists("gene_gtf_path")) { stop("gene_gtf_path not defined") }
if (!exists("genome_fasta_path")) { stop("genome_fasta_path not defined") }
if (!exists("output_gtf_path")) { stop("output_gtf_path not defined") }
if (!exists("N_block_size")) { stop("Threshold of consecutive N in a sequence not defined") }
if (!exists("N_proportion")) { stop("Threshold of N proportion in a sequence not defined") }

######################## Functions #############################

## Function for obtaining the part of the annotation field from gtf file
## Input: whole field, already splitted, for example: gene_id "FBgn0264003"; gene_name "mir-5613"; gene_source "FlyBase"; gene_biotype "pre_miRNA"
## Output: could be for example gene_id: FBgn0264003
get_annot_value <- function(split_annotation, field_name){
  ## find the right field
  field_all <- split_annotation[grep(field_name, split_annotation, fixed = T)];

  ## split the field
  field_value <- strsplit(field_all, ' ', fixed = T)[[1]][2];

  ## remove the last ';' if necessairy
  field_value <- sub(';', '', field_value,fixed = T)

  return(field_value)
}

## Function used to remove Ns from intergenic regions.
## Historicaly intergenic regions were detected only using the gtf file.
## These regions should contain ATGC bp but sometimes it happens that they contain a lot of N (sequencing errors).
## These N are mainly present in big blocks (more than `max_block_size bp) and are pseudorandomly transformed to A, T, G or C during the index generation step of kallisto.
## The purpose of this function is to remove N using some rules:
## - Remove all regions only composed of N if the size of this region is bigger than a defined value (we choose a size of 31 because it correspond to standard kmer size in kallisto).
##   Removing these region can result to more intergenic regions with a shortest length. These "splited" regions are kept if there size is bigger than `min_intergenic_length
## - Remove all intergenic sequences containing more than a certain proportion of N.
## Attributs:
## - chr_number: the id of the chromosome/contig
## - chr_intergenic_regions: intergenic regions defined using gtf file
## - max_block_size : threshold on maximum size of a block of N. If a sequence contains a block of N bigger or equals to this threshold we remove the block and split the
##	 sequence in 2. It will potentially create 2 intergenic sequences (if the block is not at the beginning nor the end of the sequence).
## - min_intergenic_length : minimum length of an intergenic region we want to keep. We remove all intergenic sequences smaller than this minimmum length.
## - max_proportion_N : maximum proportion of N in a sequence. If the proportion of N in a sequence is bigger than this threshold we do not keep the sequence.
## TODO: we should maybe also take into account the composition (GC content, ...?) or size distribution of genic regions to create/filter these intergenic regions.

remove_Ns_from_intergenic <- function (chr_number, chr_intergenic_regions, max_block_size = 31, max_proportion_N = 0.05, min_intergenic_length = 1000) {
  ## Intergenic regions after removing unwanted N
  intergenic_regions_without_N <- matrix(ncol = 4, nrow = 0)
  colnames(intergenic_regions_without_N) <- c("chr", "start", "end", "sequence")

  # For each intergenic sequence
  for (line in 1 : nrow(chr_intergenic_regions)) {

    # absolut start and end positions in the chromosome
    chr_start <- as.numeric(chr_intergenic_regions[line,"start"])
    chr_end <- as.numeric(chr_intergenic_regions[line,"end"])

    # start position in the intergenic sequence. Allows to retrieve the sequence
    intergenic_start <- 1

    # retrieve intergenic sequence from fasta file using subseq function from BioStrings library
    intergenic_sequence <- chr_intergenic_regions[line,"sequence"]

    # if less N than the minimum size of a block we add the intergenic sequence to corrected intergenic regions
    if (max_block_size >=  count_number_of_occurences("N", intergenic_sequence)) {
      intergenic_regions_without_N <- rbind(intergenic_regions_without_N, chr_intergenic_regions[line,])

      # if no block of N and proportion of N lower than threshold then we add the sequence to corrected intergenic regions
    } else if (count_number_of_occurences(strrep("N", max_block_size), intergenic_sequence) ==  0
               & proportion_of_N(intergenic_sequence) <=  max_proportion_N) {
      intergenic_regions_without_N <- rbind(intergenic_regions_without_N, chr_intergenic_regions[line,])

    # Now need to parse the sequence in order to find potential block of N
    } else {

      # for each bp of the sequence
      seq_position <- 0
      while (seq_position + 1 <=  nchar(intergenic_sequence)) {
        seq_position <- seq_position + 1

        # start a new block of N if the current bp is a N
        if (substr(intergenic_sequence, seq_position, seq_position) ==  "N") {
          block_size = 1

          # continue to increase size of the block of N until we are at a real bp position or at the end of the sequence
          while (seq_position + 1 <=  nchar(intergenic_sequence) && substr(intergenic_sequence,seq_position + 1,seq_position + 1) ==  "N") {
            block_size <- block_size + 1
            seq_position <- seq_position + 1
          }

          # check the size of the block of N to know if it should be removed
          if (block_size >=  max_block_size) {

            #absolute position in the chromosome
            current_chr_stop <- as.numeric(chr_intergenic_regions[line,"start"]) + seq_position - (block_size + 1)
            # position in the intergenic sequence (the first one in this addition correspond to the start position of the sequence)
            current_intergenic_stop <- seq_position - block_size
            current_size <- (current_chr_stop - chr_start) + 1
            # check size and proportion of N of the subsequence
            if (current_size >=  min_intergenic_length && proportion_of_N(substr(intergenic_sequence, intergenic_start, current_intergenic_stop)) <=  max_proportion_N) {
              intergenic_regions_without_N <- rbind(intergenic_regions_without_N,
                                                    cbind(chr_number, chr_start, current_chr_stop, substr(intergenic_sequence, intergenic_start, current_intergenic_stop))[1,])
            }
            intergenic_start <- 1 + seq_position
            chr_start <- as.numeric(chr_intergenic_regions[line,"start"]) + seq_position
          }
        }
      }
      # test if it remains one intergenic region to add
      last_portion_size <- (seq_position - intergenic_start) + 1
      # check size and proportion of N of the subsequence
      if (last_portion_size >=  min_intergenic_length && proportion_of_N(substr(intergenic_sequence, intergenic_start, seq_position)) <=  max_proportion_N) {
        intergenic_regions_without_N <- rbind(intergenic_regions_without_N,
                                              cbind(chr_number, chr_start, chr_end, substr(intergenic_sequence, intergenic_start, seq_position))[1,])
      }
    }

  }
  return(intergenic_regions_without_N)
}

# Function used to calculate the proportion of N in a sequence
proportion_of_N <- function (sequence) {
  return (countPattern("N", sequence) / nchar(sequence))
}

# Function allowing to count number of times one string is present in an other string (e.g count_number_of_occurences("N", "ATGNNCTN") ==  3)
# - pattern : pattern that could be present in the string
# - string : string in which the pattern could be present
count_number_of_occurences <- function (pattern, string) {
  return (lengths(regmatches(string, gregexpr(pattern, string))))
}

###################################################################

## reading in fasta file (gzipped, no need to uncompress)
cat("Reading FASTA file...\n")
ref_fasta_genome <- readDNAStringSet(genome_fasta_path)

# header of chromosomes/contigs has to be changed in order to map name of chr/contig present in gtf file.
# Only the name of chs/contig has to be kept
# regex allowing to retrieve name of the chs/contig (keep all characters before the first space)
chr_contig_regex <- "^([^ ]+).+"
# change header of each sequence
names(ref_fasta_genome) <- sub(chr_contig_regex, "\\1", names(ref_fasta_genome))

## reading in gtf file (gzipped, no need to uncompress)
cat("Reading GTF file...\n")
gene_gtf <- as.matrix(read.table(file = gzfile(gene_gtf_path, "r"), sep = "\t", strip.white = TRUE, as.is = TRUE, colClasses = "character", comment.char = '#'))
## Header lines starting with # are not read

## selecting exon lines
gene_gtf_exon <- gene_gtf[gene_gtf[,3] == "exon",]
## selecting gene lines
gene_gtf <- gene_gtf[gene_gtf[,3] == "gene",]

##  selecting only genes from assembled chromosome sequence (potentially useful code to select only fully asembled genome sequence)
##  chromosomes_all <- c(1:100, "X", "Y", "MT")
##  gene_gtf_exon <- gene_gtf_exon[gene_gtf_exon[,1] %in% chromosomes_all,]

## getting gene start, stop, chromosome, strand, biotype for each gene (from the exon GTF)
## GTF format reports 1-based coordinates (start and end)
cat("Extracting gene informations...\n")
## splitting the annotation field using single space "; " as a pattern
split_annotation_list <- strsplit(gene_gtf_exon[,9], "; ",  fixed  = T)

## getting the vector of the gene IDs (1 for every exon)
gene_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_id') })

## getting the vector of the transcript IDs (1 for every exon)
transcript_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'transcript_id') })

## getting the table with mappings between gene IDs and transcript (for tximport)
tx2gene_ids <- unique(cbind(transcript_ids, gene_ids), MARGIN = 1)

## getting the vector of gene_biotypes: splitting gene gtf
split_annotation_list <- strsplit(gene_gtf[,9], "; ",  fixed  = T)
gene_biotypes <- cbind(
    sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_id') }),
    sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_biotype') }),
    "genic")

## getting the chromosome, start and the end of the gene
## For start, take the minimum exon start position
## For end, take the maximum exon end position
gene_start <- sapply(split(as.numeric(gene_gtf_exon[,4]), gene_ids), function(x){ sort(as.numeric(x))[1] })
gene_stop <- sapply(split(as.numeric(gene_gtf_exon[,5]), gene_ids), function(x){ rev(sort(as.numeric(x)))[1] })
gene_chr <- sapply(split(gene_gtf_exon[,1], gene_ids), function(x){ x[1] })

## chromosome/contig names from given gtf files
chromosomes <- unique(gene_gtf_exon[,1])
## removing patch contigs before selecting intergenic regions
chromosomes <- chromosomes[grep('PATCH', chromosomes, invert = TRUE, ignore.case = TRUE)]

###################################################################
## Select the set of intergenic regions
cat("Selecting set of intergenic regions...\n")

## This object will include the coordinates of the selected intergenic regions
final_intergenic_regions <- matrix(ncol = 3, nrow = 0)
colnames(final_intergenic_regions) <- c("chr", "start", "end")

## Matrix summarizing impact of N removal
summary_N_removal <- matrix(ncol = 4, nrow = 2, data = 0)
colnames(summary_N_removal) <- c("intergenic_regions", "total_bp", "total_N", "proportion_N")
rownames(summary_N_removal) <- c("before", "after")

# To avoid scientific notation we change the value of option "scipen" to 999. At the end of the script we will change to its initial value
scipen_initial_value <- getOption("scipen")
options(scipen = 999)

for(chr in chromosomes){

  ## Variable used to store intergenic regions of this chr with all N
  chr_intergenic_regions <- matrix(ncol = 4, nrow = 0)
  colnames(chr_intergenic_regions) <- c("chr", "start", "end", "sequence")

  cat(paste0("start generation of intergenic regions for chromosome ", chr, "\n"))
  ## keeping genes from selected chromosome
  ## skip chromosome/contig if 1 or less gene
  if(( sum(gene_chr == as.character(chr)) <=  1 )){ next }

  ## constructing coverage map of the chromosomes using start and stop coordinates of the genes and selecting regions with 0 coverage (intergenic)
  gene_IR <- IRanges(as.numeric(gene_start[gene_chr == as.character(chr)]), as.numeric(gene_stop[gene_chr == as.character(chr)]))
  inter_IR <- slice(coverage(gene_IR), lower = 0, upper = 0, rangesOnly = TRUE)
  inter_gene_data <- as.data.frame(cbind(inter_IR@start, end(inter_IR), inter_IR@width))
  colnames(inter_gene_data) <- c("start", "end", "width")

  ## if selected chromosome/contig has only no intergenic regions with length min 2000nt then skip this contig
  if( sum(inter_gene_data[,3] > 2000) ==  0 ){ next }

  ## restrict to regions larger than 2000nt
  inter_gene_data <- inter_gene_data[inter_gene_data[,3] > 2000,]

  ## finding the center of the intergenic region
  inter_gene_data$center <- round((inter_gene_data[,1]+inter_gene_data[,2])/2)

  ## getting the size of "usable" intergenic regions around the center point
  ## width - 500nt around each flank
  inter_gene_data$size <- inter_gene_data[,3] - 1000
  inter_gene_data$size[inter_gene_data$size > 20000] <- 20000 ## limit the size of usable regions to max 20000

  ## get the sequence of the chromosome
  chr_sequence <- as.character(ref_fasta_genome[[chr]])

  ## get corrected start and stop position for intergenic regions size <=  20000
  inter_gene_data$corrected_start <- inter_gene_data$center - ceiling(inter_gene_data$size/2) + 1
  inter_gene_data$corrected_end <- inter_gene_data$center + floor(inter_gene_data$size/2)

  ## storing information (chr, start, stop) for selected intergenic regions on this chromosome
  chr_intergenic_regions <- rbind(chr_intergenic_regions, cbind(chr,  inter_gene_data$corrected_start, inter_gene_data$corrected_end,
                                                                apply(inter_gene_data, 1, function(x) substr(chr_sequence, x["corrected_start"], x["corrected_end"]))))
  # Keep information of number of N, bp and number of intergenic regions before removing blocks of N
  summary_N_removal["before","total_N"] <- summary_N_removal["before","total_N"] + sum(apply(chr_intergenic_regions, 1, function(x) count_number_of_occurences("N", x["sequence"])))
  summary_N_removal["before","intergenic_regions"] <- summary_N_removal["before","intergenic_regions"] + nrow(chr_intergenic_regions)
  summary_N_removal["before","total_bp"] <- summary_N_removal["before","total_bp"] + sum(as.numeric(inter_gene_data$size))


  # Remove blocks of N and intergenic regions with big proportion of N
  chr_intergenic_regions_after_N_removal <- remove_Ns_from_intergenic(chr, chr_intergenic_regions, as.numeric(N_block_size), as.numeric(N_proportion))

  # Keep information of number of N, bp and number of intergenic regions before removing blocks of N
  summary_N_removal["after","total_N"] <- summary_N_removal["after","total_N"] + sum(apply(chr_intergenic_regions_after_N_removal, 1, function(x) count_number_of_occurences("N", x["sequence"])))
  summary_N_removal["after","intergenic_regions"] <- summary_N_removal["after","intergenic_regions"] + nrow(chr_intergenic_regions_after_N_removal)
  summary_N_removal["after","total_bp"] <- summary_N_removal["after","total_bp"] + sum(as.numeric(chr_intergenic_regions_after_N_removal[,"end"]) - as.numeric(chr_intergenic_regions_after_N_removal[,"start"]) + 1)
  final_intergenic_regions <- rbind(final_intergenic_regions, chr_intergenic_regions_after_N_removal[,c(1:3)])

}
options(scipen = scipen_initial_value)

## preparing intergenic gtf data
cat("\nPreparing intergenic GTF data...\n")
intergenic_regions_gtf <- matrix(ncol = 9, nrow = nrow(final_intergenic_regions))
intergenic_regions_gtf[,1] <- final_intergenic_regions[,1]
intergenic_regions_gtf[,2] <- "intergenic"
intergenic_regions_gtf[,3] <- "exon"
intergenic_regions_gtf[,c(4,5)] <- final_intergenic_regions[,c(2,3)]
intergenic_regions_gtf[,c(6,8)] <- "."
## Strand is chosen randomly
set.seed(12) ## setting seed for random number generator
intergenic_regions_gtf[,7] <- sample(c("+", "-"), nrow(final_intergenic_regions), replace  = TRUE)

## intergenic_id - chr "_" start "_" stop
intergenic_id <- apply(final_intergenic_regions[,1:3], 1, function(x){ paste(x, collapse = "_") })
gene_id <- paste("gene_id ", paste(intergenic_id, ";", sep = ""), sep = "")
transcript_id <- paste("transcript_id ", paste(intergenic_id, ";", sep = ""), sep = "")
intergenic_regions_gtf[,9] <- apply(cbind(gene_id, transcript_id), 1, function(x){ paste(x, collapse = " ") })

## Caluclate proportion of N before and after N removal
summary_N_removal[,"proportion_N"] <- summary_N_removal[,"total_N"] / summary_N_removal[,"total_bp"]

####################################
output_file_path <- file.path(output_gtf_path, basename(substring(gene_gtf_path, 1, (nchar(gene_gtf_path) -7))))

## Output:
## GTF file with both genic exons and intergenic regions
message("write file : ", paste(output_file_path, ".gtf_all", sep = ""))
write.table(x = rbind(gene_gtf_exon, intergenic_regions_gtf),
            file = paste(output_file_path, ".gtf_all", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

## GTF file with both genic exons and intergenic regions
message("write file : ", paste(output_file_path, ".gtf_transcriptome", sep = ""))
write.table(x = gene_gtf_exon,
            file = paste(output_file_path, ".gtf_transcriptome", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

## Table summaryzing the step of block of Ns removal
message("write file : ", paste(output_file_path, ".Nremoval", sep = ""))
write.table(x = summary_N_removal,
            file = paste(output_file_path, ".Nremoval", sep = ""),
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE)

## Table of mapping between transcript_id and gene_id
intergenic_tx2gene_ids <- cbind(intergenic_id, intergenic_id)
all_tx2gene_ids <- rbind(tx2gene_ids, intergenic_tx2gene_ids)
message("write file : ", paste(output_file_path, ".tx2gene", sep = ""))
write.table(x = all_tx2gene_ids,
            file = paste(output_file_path, ".tx2gene", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = c("TXNAME", "GENEID"),
            quote = FALSE)

## Table of mapping between gene_id and both biotype and type (genic or intergenic)
intergenic_gene_biotypes <- cbind(
    intergenic_id,
    NA,
    "intergenic")
gene_biotypes <- rbind(gene_biotypes, intergenic_gene_biotypes)
message("write file : ", paste(output_file_path, ".gene2biotype", sep = ""))
write.table(x = gene_biotypes,
            file = paste(output_file_path, ".gene2biotype", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = c("id", "biotype", "type"),
            quote = FALSE)
