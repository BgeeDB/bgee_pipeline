## SFonsecaCosta, June 2022

## This script is used to :
##    - remove doublets and multiplets from the annotation file (information collected from the authors)
##    - remove barcodes with no cell-type ID
##    - if a cluster ID but no free text cell-type annotation is provided by the author, then use the author cluster ID as free text annotation
##    - create internal cluster IDs used to generate calls. The only goal of these IDs is to avoid using free text author annotation to name calls files
## Note: if a barcode is detected >=2 in same experiment/library this barcode is removed from the barcode annotation file and will not be used during the analyis

## Usage:
## R CMD BATCH --no-save --no-restore '--args barcodesFolder="path_barcodes_folder" output="output_folder"' cleaning_barcodes.R cleaning_barcodes.Rout
## barcodesFolder --> Folder where are the barcodes files annotated by Bgee for each experimentID
## output --> Path where should be saved the output results after removing the duplicate barcodes per experimentID/LibraryID

## libraries used
library(dplyr)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("barcodesFolder", "output")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}
########################################################################################################################################################
if(!dir.exists(output)) {
  dir.create(output)
}

barcodes_files_path <-  list.files(barcodesFolder, pattern = "scRNASeq_barcode.*.tsv", full.names = TRUE)
print(paste("Found", length(barcodes_files_path), "barcode files in", barcodesFolder))
for (barcode_file_path in barcodes_files_path) {
  
  barcodes <- tryCatch(
    read.table(barcode_file_path, header = TRUE, sep = "\t"),
    error = function(e) {
      stop(paste("Error reading file:", barcode_file_path, "\n", e$message))
    }
  )
  #first we remove all rows with cellTypeId empty or NA as annotation of a barcode with a term coming from an ontology
  # is mandatory to insert data in Bgee
  barcodes <- barcodes %>% filter(!is.na(cellTypeId) & cellTypeId != "")
  barcode_file_name <- basename(barcode_file_path)
  unique_barcodes_list <- list()
  duplicated_barcodes_list <- list()
  # Initialize statistics tracking
  stats_same_annotation <- list()
  stats_different_annotation <- list()
  
  for (library_id in unique(barcodes$library)) {
    cluster_id <- 1
    select_barcodes <- dplyr::filter(barcodes, barcodes$library == library_id)
    
    # Count duplicates with same annotation (cellTypeId and cell_type)
    duplicates_same_annotation <- select_barcodes %>%
      group_by(barcode) %>%
      filter(n() > 1 & n_distinct(cellTypeId) == 1 & n_distinct(cell_type) == 1) %>%
      ungroup()
    
    num_barcodes_same_annotation <- n_distinct(duplicates_same_annotation$barcode)
    stats_same_annotation[[library_id]] <- num_barcodes_same_annotation
    
    # Get the extra occurrences of duplicates with same annotation (all but the first one)
    # These will be moved to the Duplicated file
    duplicates_same_annotation_extra <- select_barcodes %>%
      group_by(barcode) %>%
      filter(n() > 1 & n_distinct(cellTypeId) == 1 & n_distinct(cell_type) == 1) %>%
      slice(-1) %>%  # Remove first occurrence, keep all others
      ungroup()
    
    #select only unique barcodes/library. This check should be done at annotation level before running the pipeline
    unique_barcodes <- select_barcodes %>%
      group_by(barcode) %>% 
      # keep barcodes that are unique or that are duplicated but with the same cellTypeId and cell_type.
      # Keep only one occurance of deplicated barcodes with same cellTypeId
      filter(n_distinct(cellTypeId) == 1 & n_distinct(cell_type) == 1) %>%
      slice(1)

    uniq_celltype_freetext <- unique(unique_barcodes$cell_type)
    # As the cell_type column is used to detect clusters of barcodes, we can not have cell_type that is NA
    contains_celltype_freeText <- TRUE
    if (is.null(uniq_celltype_freetext) || all(is.na(uniq_celltype_freetext)) || (length(uniq_celltype_freetext) == 1 && uniq_celltype_freetext == "")) {
      unique_barcodes$cell_type <- ""
      message("free text cell-type information is always NA or empty. Will check if cluster info are available for file ", barcode_file_name, " and library ", library_id)
      contains_celltype_freeText <- FALSE
    }
    contains_cluster <- TRUE
    uniq_cluster_id <- unique(unique_barcodes$cluster)
    if (is.null(uniq_cluster_id) || all(is.na(uniq_cluster_id)) || (length(uniq_cluster_id) == 1 && uniq_cluster_id == "")) {
      unique_barcodes$cluster <- ""
      contains_cluster <- FALSE
    }
    if (contains_cluster && !contains_celltype_freeText) {
      # Check if each cluster ID maps to exactly one cellTypeId (1-to-1 relationship)
      num_unique_clusters <- length(uniq_cluster_id)
      num_cluster_celltype_combinations <- nrow(unique(unique_barcodes[, c("cluster", "cellTypeId")]))
      
      if (num_unique_clusters != num_cluster_celltype_combinations) {
        message("cluster IDs are not unique for a combination of library and celltypeId for file ", barcode_file_name,
          ". Removes cluster IDs from that file in order not to store it as author celltype info in the Bgee database")
        unique_barcodes$cluster <- ""
        contains_cluster <- FALSE
      } else {
        message("No free text celltype annotation but an author cluster ID exists. Will use this cluster IDs as free text annotation. File : ", barcode_file_name, " and library ", library_id)
        unique_barcodes$cell_type <- paste0("cluster ", unique_barcodes$cluster)
        contains_celltype_freeText <- TRUE
      }
    }

    # we now generate an internal cluster ID. The only use of that ID is to be able to name count matrices files
    # generated by the pipeline without using the free text celltype annotation from the authors.
    # This cluster ID is not inserted in the database and is only internal to the pipeline.
    unique_barcodes$internal_cluster_id <- -1
    # if free text annotation are provided by the author, or cluster IDs have been used to generate free text annotation,
    # then they are used to generate the internal cluster IDs
    if (contains_celltype_freeText) {
      for (celltype_author_annotation in unique(unique_barcodes$cell_type)) {
        unique_barcodes$internal_cluster_id[unique_barcodes$cell_type == celltype_author_annotation] <- cluster_id
        cluster_id <- cluster_id + 1
      }
    # if not free text annotation or unique cluster IDs are provided by the author then cellTypeIds are used to generate
    # internal cluster IDs
    } else {
      for (celltype_id in unique(unique_barcodes$cellTypeId)) {
        message("cell type ID for library ", library_id, " is ", celltype_id)
        unique_barcodes$internal_cluster_id[unique_barcodes$cellTypeId == celltype_id] <- cluster_id
        cluster_id <- cluster_id + 1
      }
    }

    # Get duplicates with different annotation - ALL occurrences go to duplicated file
    barcodes_duplicates_different <- select_barcodes %>% 
      group_by(barcode) %>% 
      filter(n() > 1 & (n_distinct(cellTypeId) > 1 | n_distinct(cell_type) > 1))
    
    # Combine both types of duplicates for the Duplicated file
    barcodes_duplicates <- bind_rows(duplicates_same_annotation_extra, barcodes_duplicates_different)
    
    # Count duplicates with different annotation
    num_barcodes_different_annotation <- n_distinct(barcodes_duplicates_different$barcode)
    stats_different_annotation[[library_id]] <- num_barcodes_different_annotation
    
    # Report statistics for this library
    if (num_barcodes_same_annotation > 0) {
      message("Library ", library_id, ": ", num_barcodes_same_annotation, 
              " barcodes were present more than once with same cellTypeId and cell_type (keeping one occurrence, others in Duplicated file)")
    }
    if (num_barcodes_different_annotation > 0) {
      message("Library ", library_id, ": ", num_barcodes_different_annotation, 
              " barcodes were present more than once with different annotations (all moved to Duplicated file)")
    }

    ## Append results to lists
    unique_barcodes_list[[library_id]] <- unique_barcodes
    duplicated_barcodes_list[[library_id]] <- barcodes_duplicates

  }

    ## Combine results from all libraries using bind_rows
  all_unique_barcodes <- bind_rows(unique_barcodes_list)
  all_duplicated_barcodes <- bind_rows(duplicated_barcodes_list)
  
  # Summary statistics for the entire file
  total_same_annotation <- sum(unlist(stats_same_annotation))
  total_different_annotation <- sum(unlist(stats_different_annotation))
  if (total_same_annotation > 0 || total_different_annotation > 0) {
    message("====== SUMMARY for ", barcode_file_name, " ======")
    message("Total barcodes with duplicates (same annotation): ", total_same_annotation)
    message("Total barcodes with duplicates (different annotation): ", total_different_annotation)
    message("================================================")
  }
  
  tryCatch(
    write.table(all_unique_barcodes, file = file.path(output, barcode_file_name), sep = "\t", row.names = FALSE, quote = TRUE),
    error = function(e) {
      stop(paste("Error writing file:", barcode_file_name, "\n", e$message))
    }
  )

  if (nrow(all_duplicated_barcodes) > 0){
    write.table(all_duplicated_barcodes, file = file.path(output, paste0("DuplicateBarcodes_", barcode_file_name)), sep = "\t", row.names = FALSE, quote = FALSE)
  }
}
