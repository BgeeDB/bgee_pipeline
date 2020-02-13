## Julien Roux Feb 2, 2016
## Partially inspired from previous pipeline scripts written by Marta Rosikiewicz
## This script analyzes counts and TPMs returned by Kallisto at the transcript level and computes counts, RKPMs and TPMs at gene level

## Usage:
## R CMD BATCH --no-save --no-restore '--args  kallisto_count_folder= "kallisto_count_folder" gene2transcript_file="gene2transcript_file" gene2biotype_file="gene2biotype_file" library_id="library_id" rna_seq_analysis.R library_id.Rout
## kallisto_count_folder - path to kallisto result file
## gene2transcript_file  - path to file with genes to transcripts mapping
## gene2biotype_file     - path to file with biotype for all genes
## library_id            - Id of library to treat

## Session info
print(sessionInfo())
library(tximport)

#################### PRELIMINARY STEPS ####################

## reading in arguments provided in command line
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed in command line
command_arg <- c("kallisto_count_folder", "tx2gene_file", "gene2biotype_file", "library_id")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## reading kallisto's output. If file not exists, script stops
kallisto_count_tsv_file <- paste0(kallisto_count_folder, "/abundance.tsv")
if( !file.exists(kallisto_count_file) ){
	kallisto_count <- read.table(kallisto_count_file, h=T, sep="\t")
} else {
  stop( paste("Kallisto results file not found [", kallisto_count_file, "]\n"))
}
## Basic sanity check:
if ( sum(is.na(kallisto_count$tpm)) == length(kallisto_count$tpm)){
  stop( paste("Kallisto results include only NAs, please check for a problem (k-mer size too small?) [", kallisto_count_file, "]\n"))
} else if ( sum(is.na(kallisto_count$tpm)) > 0.2 * length(kallisto_count$tpm)){
  warning( paste("Kallisto results include >20% NAs, please check for a problem [", kallisto_count_file, "]\n"))
}

## If transcript to gene file does not exists, script stops
if( file.exists(tx2gene_file) ){
  tx2gene <- read.table(tx2gene_file, header = TRUE, sep="\t")
} else {
  stop( paste("transcript to gene file not found [", tx2gene_file, "]\n"))
}

## reading gene biotype file. If file not exists, script stops
if( file.exists(gene2biotype_file) ){
	gene2biotype <- read.table(gene2biotype_file, h=T, sep="\t")
  gene2biotype <- gene2biotype[order(gene2biotype$id), ] ## order by gene Id
} else {
  stop( paste("Gene biotype file not found [", gene2biotype_file, "]\n"))
}

#################### FUNCTIONS ####################

## calculate FPKMs, using functions from
## https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expressio
tpmToFpkm <- function(tpm, counts, effLen){
  exp(log(tpm) + log(sum(counts/effLen)) - log(sum(counts)) + log(1e3))
}

###############################################################################

tximportObject <- tximport(kallisto_count_tsv_file, type = "kallisto", tx2gene = tx2gene)
tx_df <- as.data.frame(tximportObject)
tx_df$id <- rownames(tx_df)
tx_df$countsFromAbundance <- NULL
kallisto_gene_count <- merge(tx_df, gene2biotype, by = "id", all = FALSE)
kallisto_gene_count$fpkm <- tpmToFpkm(abundance$abundance, counts = abundance$counts, effLen = abundance$length)
kallisto_gene_count <- kallisto_gene_count[order(kallisto_gene_count[,6], kallisto_gene_count[,1]),]
# order first by type (genic first) and then by gene id
# the pipeline did not use tximport to summarize abundance at gene level. It uses an homemade script.
# We modify the output of tximport in order to create the same file than the old homemade script
kallisto_gene_count <- kallisto_gene_count[,c(1, 3, 2, 7, 6, 5)]
names(kallisto_gene_count) <- c("gene_id", "est_counts", "tpm", "fpkm", "type", "biotype")

write.table(kallisto_gene_count, file = paste0(kallisto_count_folder, "/abundance_gene_level+fpkm+intergenic.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


###############################################################################
## Plotting of the distribution of TPMs for genic and intergenic regions
pdf(file = paste0(kallisto_count_folder, "/distribution_TPM_genic_intergenic.pdf"), width = 6, height = 5)
par(mar=c(5,6,1,1)) ## bottom, left, top and right margins

dens <- density(log2(na.omit(kallisto_gene_count$tpm) + 10^-6))
## Subgroups densities. Visualization trick: we add an invisible set of points at x=-30, to make densities comparable
## genic regions
dens_genic <- density(c(rep(-30, times=sum(kallisto_gene_count$type != "genic")), log2(kallisto_gene_count$tpm[kallisto_gene_count$type == "genic"] + 10^-6)))
## protein-coding genes only (had to take care of NAs strange behavior)
dens_coding <- density(c(rep(-30, times=sum(!kallisto_gene_count$biotype %in% "protein_coding")), log2(kallisto_gene_count$tpm[kallisto_gene_count$biotype %in% "protein_coding"] + 10^-6)))
## intergenic
dens_intergenic <- density(c(rep(-30, times=sum(kallisto_gene_count$type != "intergenic")), log2(kallisto_gene_count$tpm[kallisto_gene_count$type == "intergenic"] + 10^-6)))

## Plot whole distribution
plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-23, 21), lwd=2, main="", bty="n", axes=F, xlab="")
axis(2, las=1)
## Add 2 x-axes: TPMs and RPKMs
## See http://stackoverflow.com/questions/8443820/r-multiple-x-axis-with-annotations
axis(1, at=seq(-30 , 30, by=10), line=0, mgp = c(3, 0.5, 0), cex.axis=0.8)
mtext(expression(log[2]('TPM'+10^-6)), 1,  adj = 1, padj = 0, line=0.2, at=par("usr")[1], col="black", cex=0.8)
## To make FPKM scale, we need to know what log2(TPM + 10^-6) value corresponds to any log2(FPKM + 10^-6) value. We know that for any given gene, TPMg/FPKMg = coef
##    log2(FPKM + 10^-6) = x
## <=>              FPKM = exp(x*log(2)) - 10^-6
## <=>               TPM = coef * (exp(x*log(2)) - 10^-6)
## <=> log2(TPM + 10^-6) = log2( coef * (exp(x*log(2)) - 10^-6) + 10^-6)
coef <- na.omit(kallisto_gene_count$tpm / kallisto_gene_count$fpkm)[1]
## We generate scale from -20 to 100 FPKMs
axis(1, at=log2( coef * (exp(seq(-20 , 100, by=10)*log(2)) - 10^-6) + 10^-6), labels=seq(-20 , 100, by=10), line=2, mgp = c(3, 0.5, 0), cex.axis=0.8)
mtext(expression(log[2]('FPKM'+10^-6)), 1,  adj = 1, padj = 0, line=2.2, at=par("usr")[1], col="black", cex=0.8)

## Add subgroups distributions (genic, intergenic, etc):
## genic
lines(dens_genic, col="firebrick3", lwd=2)
## protein-coding genes
lines(dens_coding, col="firebrick3", lwd=2, lty=2)
## intergenic
lines(dens_intergenic, col="dodgerblue3", lwd=2)

## legend
legend("topright", c("all", "genic", "protein-coding genes", "intergenic"), lwd=2, col=c("black", "firebrick3", "firebrick3", "dodgerblue3"), lty=c(1, 1, 2, 1), bty="n")
dev.off()
