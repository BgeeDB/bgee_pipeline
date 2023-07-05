## Julien Wollbrett, May 10 2023

## highly inspired from script process_busFile.R written by S. Fonseca but slightly improved to be
## able to parallelize bustools step per library

## This script is used to process the bus files to an equivalent-class-UMI and to a gene-UMI
## As suggested by bustools this step is divide in 3 sub-steps:

## 1) Correct the barcodes using bustools correct: fix the barcodes that are within one hamming distance of the barcodes in the whitelist using whitelist.txt
## 2) Sort the busfile using bustools sort: organize the busfile by barcode, UMI, set and multiplicity.
## 3) Count records in the BUS with bustools count: generate the UMI count matrix using transcripts_to_genes.txt.
## 4) Normalized CPM expression level for all genes without intergenic regions

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNA_Seq_info_TargetBased.txt" kallisto_bus_results="kallisto_bus_results" folder_gtf="folder_gtf" whiteList_Path="whiteList_Path"' process_busFile.R process_busFile.Rout
## libraryId            --> the ID of the library for which bustools has to be run
## scRNASeq_Info        --> File that results from annotation and metadata (libraries downloaded and with extra information as SRR)
## kallisto_bus_results --> Folder where are all kallisto bus results are saved
## folder_gtf           --> Folder containing informative files as: transcriptomes index + gtf_all + transcript_to_gene
## whiteList_Path       --> Folder containing the barcode_whitelist files

## libraries used
library(data.table)
library(BUSpaRse)
library(Matrix)
library(Seurat)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed.
command_arg <- c("libraryId", "scRNASeq_Info", "kallisto_bus_results", "folder_gtf", "whiteList_Path")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read scRNA-Seq info file. If file not exists, script stops
if( file.exists(scRNASeq_Info) ){
  scRNAInfo <- fread(scRNASeq_Info)
} else {
  stop( paste("The annotation file not found [", scRNASeq_Info, "]\n"))
}

generateGenesCpm <- function(tx2geneFile, cpm_dir, gene_counts, bustoolsGeneMatrixName) {
  # create cpm directory
  if (!dir.exists(cpm_dir)) {
    dir.create(cpm_dir)
  }
  #load transcript to gene mapping file
  tx2gene <- read.table(tx2geneFile, sep = "\t", header = FALSE);
  # load counts from bustools using a function from BUSpaRse
  counts <- read_count_output(gene_counts, bustoolsGeneMatrixName, FALSE)
  # intergenic have same name in  both transcript and gene columns.
  # we use this caracteristic to keep only counts from genic region
  tx2geneWithoutIntergenic <- tx2gene[!tx2gene[,1] %in% tx2gene[,2],]
  genesWithoutIntergenic <- unique(tx2geneWithoutIntergenic[,2])
  subset_counts <- counts[counts@Dimnames[[1]][counts@Dimnames[[1]] %in% genesWithoutIntergenic],]

  # now calculate cpm using a function from Seurat
  subset_cpm <- NormalizeData(subset_counts, normalization.method="RC", scale.factor=1e06)

  # Write output files
  cpmMatrixFile <- file.path(cpm_dir, "cpm_counts.mtx")
  Matrix::writeMM(obj = subset_cpm, file =cpmMatrixFile)
  write.table(x = subset_cpm@Dimnames[[1]], file = file.path(cpm_dir, "cpm_counts.genes.txt"),
    sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(x = subset_cpm@Dimnames[[2]], file = file.path(cpm_dir, "cpm_counts.barcodes.txt"),
    sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  # the writeMM function does not save the same .mtx than bustools. Columns 1 and 2
  # are switched. It means that when reading the file created by writeMM with the BUSpaRse
  # read function, rows and columns are exchanged. In order to keep the same structure of
  # matrix, we manually exchange column 1 and 2 from the .mtx file
  #XXX Did not find more elegant solution...
  # transform the matrix to be as the one from bustools
  matrixTsv <- read.table(cpmMatrixFile, sep =" ", header = FALSE, comment.char = "%")
  transformedMatrixTsv <- matrixTsv[,c(2,1,3)]
  write.table(x = transformedMatrixTsv, file = cpmMatrixFile, sep = " ", col.names = FALSE,
    quote = FALSE, row.names = FALSE)
}
###############################################################################################

if (dir.exists(file.path(kallisto_bus_results, libraryId))){

  ## verify if busOutput exist for the library
  pathBusOut <-  file.path(kallisto_bus_results, libraryId)
  if (!file.exists(file.path(pathBusOut, "run_info.json"))) {
    warning("kallisto was not run for the library : ", libraryId)
  } else {
    ## the last step of this script for each library is to generate the gene cpm count matrix.
    ## The presence of the directory containing this count matrix is used as a marker to check
    ## if the library has already been processed
    cpm_dir <- file.path(pathBusOut, "cpm_counts")
    if (dir.exists(cpm_dir)) {
      message("library ", libraryId, " already processed. bustools is not run again")
    } else {  
      message("Start correction, sort and counts for the library: ", libraryId)

      ## Note: the whiteList we use in this pipeline are the files provided by 10X platform (add to source files)
      collectWhitelist <- unique(as.character(scRNAInfo$whiteList[scRNAInfo$libraryId == libraryId]))
      selectedWhitheList <- paste0("10X_",collectWhitelist)
      message("whitelist:  ", selectedWhitheList)

      ## step 1 --> correct the barcodes
      message("Correct barcodes")
      system(sprintf('%s -w %s -o %s %s', "bustools correct", paste0(whiteList_Path, "barcode_whitelist_", selectedWhitheList,".txt"), paste0(pathBusOut, "/output.correct.bus"), paste0(pathBusOut, "/output.bus")))

      ## step 2 --> sort the bus file
      message("Sort bus file")
      system(sprintf('%s -t 4 -o %s %s', "bustools sort", paste0(pathBusOut, "/output.correct.sort.bus"), paste0(pathBusOut, "/output.correct.bus")))

      ## Create folders to export the information per TCC and gene matrix (counts)
      tcc_counts <- file.path(pathBusOut, "tcc_counts")
      if (!dir.exists(tcc_counts)){
        dir.create(tcc_counts)
      } else {
        print("File already exist.....")
      }

      gene_counts <- file.path(pathBusOut, "gene_counts")
      if (!dir.exists(gene_counts)){
        dir.create(gene_counts)
      } else {
        print("File already exist.....")
      }

      collectSpecies <- unique(as.character(scRNAInfo$scientific_name[scRNAInfo$libraryId == libraryId]))
      collectSpecies <- gsub(" ", "_", collectSpecies)

      tx2gene_file = file.path(folder_gtf, paste0(collectSpecies, "_transcript_to_gene_with_intergenic.tsv"))
      print(tx2gene_file)
      if (!file.exists(tx2gene_file)) {
        if(file.exists(paste0(tx2gene_file, ".xz"))) {
          system(sprintf("unxz %s", paste0(tx2gene_file, ".xz")))
        } else {
          stop("transcript to gene file [", tx2gene_file, "] not found.")
        }
      }
      ## step 3 --> count with bustools count
      ## TCC level
      message("TCC level")
      ## the --em parameter estimates gene abundances using an EM algorithm for reads that pseudoalign to multiple genes.
      system(sprintf('bustools count --em -o %s -g %s -e %s -t %s %s', file.path(tcc_counts, "tcc"),
        tx2gene_file, file.path(pathBusOut, "matrix.ec"), file.path(pathBusOut, "transcripts.txt"),
        file.path(pathBusOut, "output.correct.sort.bus")))
      ## GENE level
      message("Gene level")
      bustoolsGeneMatrix <- "gene"
      ## the --em parameter estimates gene abundances using an EM algorithm for reads that pseudoalign to multiple genes.
      system(sprintf('bustools count --em -o %s -g %s -e %s -t %s --genecounts %s',
        file.path(gene_counts, bustoolsGeneMatrix), tx2gene_file, file.path(pathBusOut, "/matrix.ec"),
        file.path(pathBusOut, "transcripts.txt"), file.path(pathBusOut, "output.correct.sort.bus")))
      ## CPM level
      message("CPM for genes without intergenic")
      generateGenesCpm(tx2gene_file, cpm_dir, gene_counts, bustoolsGeneMatrix)
    }
  }
} else {
  message("Library ", library ," not present in the folder to process.")
}

