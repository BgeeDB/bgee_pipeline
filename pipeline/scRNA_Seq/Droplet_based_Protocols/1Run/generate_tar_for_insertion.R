## Script used to create a tarball used to move all results from the sensitive cluster
## to bgee servers in order to insert data on the DB

## example how to use it :

md_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("pathToDirToTar", "pathtoTarFile")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

pathToDirToTar <- "/home/julien/Documents/temp/testTargetBase/"

allFiles <- list.files(path = pathToDirToTar, full.names = TRUE, recursive = TRUE, all.files = TRUE)
## remove big files not mandatory for data insertion
filesToCompress <- allFiles[grep(pattern = "output.bus$", x = allFiles, invert = TRUE)]
filesToCompress <- filesToCompress[grep(pattern = "output.correct.*$", x = filesToCompress, invert = TRUE)]

relativePathsFilesToCompress <- gsub(pattern = pathToDirToTar, replacement = "./", x = filesToCompress)

tar(tarfile = file.path(pathtoTarFile, "all_results_except_bus.tar.gz"), files = filesToCompress,
    compression = "gzip", compression_level = 9)
