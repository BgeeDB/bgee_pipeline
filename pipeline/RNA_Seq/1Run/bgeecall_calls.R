## Julien Wollbrett Jun 12, 2020
## This script generate all present/absent gene expression calls using BgeeCall

## Usage as for Bgee15:
## R CMD BATCH --no-save --no-restore '--args bgeecall_input_file=$(RNASEQ_BGEECALL_FILE) account="account" time="2:00:00" partition="partition" working_path="/path/to/R/session/directory/"' bgeecall_calls.R bgeecall_calls.Rout
## bgeecall_input_file      - file with info on RNA-Seq libraries
## account                  - name of the account on the cluster
## time                     - maximum time to run each job
## partition                - partition on which jobs have to be run
## working_path             - path to directory where species specific files will be created by BgeeCall

## Session info
print(sessionInfo())

library(BgeeCall)

## reading in arguments provided in command line
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed in command line
command_arg <- c("bgeecall_input_file", "account", "time", "partition")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}
#detect number of rows in input file in order to select nodes number
libraries <- nrow(read.table(bgeecall_input_file, sep="\t", quote = "", header = TRUE))
#maximum size for an Array in slurm Jura cluster is 1000.
if(libraries>1000) {
  libraries <- 1000
}
#specific to UNIL clusters: load R module
#TODO: Maybe load kallisto from UNIL modules
modules <- c("module add Bioinformatics/Software/vital-it;", "module add R/3.6.1;", "module add UHTS/Analysis/kallisto/0.46.0;")
#generate BgeeCall objects
kallistoMetadata <- new("KallistoMetadata", download_kallisto=FALSE)
userMetadata <- new("UserMetadata", working_path = working_path)
# use local version of intergenic sequences
bgeeMetadata <- new("BgeeMetadata", intergenic_release="custom")
# slurm options for index generation
slurm_options_index <- list(account = account, time = time, partition = partition, mem = "60G")
# generate calls
generate_slurm_calls(userFile=bgeecall_input_file, slurm_options = slurm_options_index, modules = modules, kallistoMetadata = kallistoMetadata, bgeeMetadata = bgeeMetadata, userMetadata = userMetadata, nodes=libraries)
