source('prepare_GTF.R')

context('testing preparation of GTF file')
## intergenic_sequence <- "TIGCCCCCETTINNNNTIIJPOFPITINNNTIENRTINNNTIETI"
## contains 4 sequences without block of N bigger than 4bps. One of these sequences contains more than 5% of N and is discarded.

test_that("Well remove Ns", {
    intergenic_regions_without_N <- matrix(ncol=4, nrow=0)
    colnames(intergenic_regions_without_N) <- c("chr", "start", "end", "sequence")
    chr_intergenic_regions <- matrix(ncol=4, nrow=0)
    colnames(chr_intergenic_regions) <- c("chr", "start", "end", "sequence")
    chr_intergenic_regions <- rbind(chr_intergenic_regions, cbind("4", "1", "45", "TIGCCCCCETTINNNNTIIJPOFPITINNNTIERRTINNNTIETI"))
    intergenic_regions_without_N <- remove_Ns ("4", chr_intergenic_regions, max_block_size = 2, min_intergenic_length = 5, max_proportion_N = 0.05)
	  #test type
    expect_type(intergenic_regions_without_N, "character")
	  #test values
    test_intergenic_regions_without_N <- matrix(ncol=4, nrow=3,data = "")
    colnames(test_intergenic_regions_without_N) <- c("chr", "start", "end", "sequence")
    test_intergenic_regions_without_N[1,] <-cbind("4","1","12","TIGCCCCCETTI")[1,]
    test_intergenic_regions_without_N[2,] <-cbind("4","17","27","TIIJPOFPITI")[1,]
    test_intergenic_regions_without_N[3,] <-cbind("4","41","45","TIETI")[1,]
	  expect_equal(intergenic_regions_without_N, test_intergenic_regions_without_N)
})
   
