startDir = "/scratch/jwollbre/SRP125768/"
allDir <- list.dirs(path = startDir, full.names = TRUE, recursive = TRUE)

for (dir in allDir) {
  message("dir : ", dir)
  files <- list.files(path = dir, full.names = FALSE, , pattern = "_[1-2].fastq.gz", recursive = FALSE)
  message("files before check length : ",files)
  if (length(files) == 2) {
    for(file in files) {
      message("file before suffix checking : ", file)
      suffix <- sub(".*(_[0-9].fastq.gz)", "\\1", file)
      if (suffix == "_1.fastq.gz") {
        new_file_name <- gsub(pattern = "_1.fastq.gz", replacement = "_R1.fastq.gz", x = file)
        message("old file name: ", file.path(dir,file), " new file name: ", file.path(dir, new_file_name))
        file.rename(from = file.path(dir,file), to = file.path(dir, new_file_name))
      }else if (suffix == "_2.fastq.gz") {
        new_file_name <- gsub(pattern = "_2.fastq.gz", replacement = "_R2.fastq.gz", x = file)
        message("old file name: ", file.path(dir,file), " new file name: ", file.path(dir, new_file_name))
        file.rename(from = file.path(dir,file), to = file.path(dir, new_file_name))
      }
    }
  }
}
