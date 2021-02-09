## Julien Wollbrett Jun 12, 2020
## This script generate all required kallisto indexes from the BgeeCall input tsv file.

## Usage as for Bgee15:
## R CMD BATCH --no-save --no-restore '--args bgeecall_input_file=$(RNASEQ_BGEECALL_FILE) account="account" time="2:00:00" partition="partition" working_path="/path/to/R/session/directory/"' bgeecall_index.R bgeecall_index.Rout
## bgeecall_input_file      - file with info on RNA-Seq libraries
## account                  - name of the account on the cluster
## time                     - maximum time to run each job
## partition                - partition on which jobs have to be run
## working_path             - path to directory where species specific files will be created by BgeeCall

library(BgeeCall)

## Session info
print(sessionInfo())

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
#specific to UNIL clusters: load R module
modules <- c("module add Bioinformatics/Software/vital-it;", "module add R/3.6.1;", "module add UHTS/Analysis/kallisto/0.46.0;")
#generate BgeeCall objects to download kallisto and define path where data at species level will be stored
kallistoMetadata <- new("KallistoMetadata")
userMetadata <- new("UserMetadata", working_path = working_path)
bgeeMetadata <- new("BgeeMetadata", intergenic_release="custom")
# slurm options for index generation. 30G is enough for the majority of species. However as for Bgee15 it had to be increased to 90G for few species
slurm_options_index <- list(account = account, time = time, partition = partition, mem = "60G")
# generate indexes
generate_slurm_indexes(userFile=bgeecall_input_file, slurm_options = slurm_options_index, modules = modules, kallistoMetadata = kallistoMetadata, bgeeMetadata=bgeeMetadata, userMetadata = userMetadata, nodes=50)
