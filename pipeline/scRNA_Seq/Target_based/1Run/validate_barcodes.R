## Script to validate barcodes in Raw_Counts files against annotation files
## Checks if barcodes in Raw_Counts_* files match expected barcodes from scRNASeq_barcode_EXP_ID.tsv

## Usage:
## R CMD BATCH --no-save --no-restore '--args scRNASeq_Info="scRNA_Seq_info_target.txt" infoFolder="folder_with_barcodes_files" outputFolder="output_folder" ncores=4' validate_barcodes.R validate_barcodes.Rout

library(dplyr)
library(data.table)
library(foreach)
library(doParallel)

## reading arguments
cmd_args = commandArgs(TRUE)
print(cmd_args)
if(length(cmd_args) == 0) {
  stop("no arguments provided\n")
} else {
  for(i in 1:length(cmd_args)) {
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed
command_arg <- c("scRNASeq_Info", "infoFolder", "outputFolder")
for(c_arg in command_arg) {
  if(!exists(c_arg)) {
    stop(paste(c_arg, "command line argument not provided\n"))
  }
}

## Set number of cores (default to all available - 1)
if(!exists("ncores")) {
  ncores <- max(1, parallel::detectCores() - 1)
} else {
  ncores <- as.integer(ncores)
}
message("Using ", ncores, " cores for parallel processing")

## Read scRNASeq_Info file
if(file.exists(scRNASeq_Info)) {
  scRNASeqAnnotation <- fread(scRNASeq_Info, header = TRUE, sep = "\t", quote = "\"", data.table = FALSE)
} else {
  stop(paste("The scRNASeq_Info file not found [", scRNASeq_Info, "]\n"))
}

## Get unique experiments and libraries
experiments <- unique(scRNASeqAnnotation$experiment_id)
libraries <- unique(scRNASeqAnnotation$library_id)

message("Found ", length(libraries), " libraries across ", length(experiments), " experiments")

## Setup parallel backend
cl <- makeCluster(ncores)
registerDoParallel(cl)

## Results storage
libraries_with_issues <- list()

## Process experiments in parallel
all_issues <- foreach(experimentID = experiments, 
                      .combine = list,
                      .multicombine = TRUE,
                      .packages = c("dplyr", "data.table"),
                      .errorhandling = "pass",
                      .verbose = FALSE) %dopar% {

  # Read barcode annotation file once per experiment
  barcodeFile <- file.path(infoFolder, paste0("scRNASeq_barcode_", experimentID, ".tsv"))
  
  # Get all libraries for this experiment (needed for diagnostics)
  experiment_libraries <- unique(scRNASeqAnnotation$library_id[
    scRNASeqAnnotation$experiment_id == experimentID
  ])
  
  # Initialize diagnostics structure
  experiment_diagnostics <- list(
    experiment_id = experimentID,
    library_count = length(experiment_libraries),
    missing_folders = character(0),
    missing_barcodes = character(0),
    missing_raw_files = character(0),
    barcode_file_not_found = FALSE
  )
  
  # fread can read .gz files directly - check both versions
  if(!file.exists(barcodeFile)) {
    barcodeFileGz <- paste0(barcodeFile, ".gz")
    if(file.exists(barcodeFileGz)) {
      barcodeFile <- barcodeFileGz
    } else {
      # Mark that barcode file wasn't found and return diagnostics
      experiment_diagnostics$barcode_file_not_found <- TRUE
      return(list(issues = list(), diagnostics = experiment_diagnostics))
    }
  }

  barcodeExperiment <- fread(barcodeFile, header = TRUE, sep = "\t", quote = "\"", data.table = FALSE)

  # Storage for this experiment's issues
  experiment_issues <- list()

  # Loop through each library of this experiment
  for(libraryId in experiment_libraries) {

    library_folder <- file.path(outputFolder, libraryId)

    # Check if library folder exists
    if(!dir.exists(library_folder)) {
      experiment_diagnostics$missing_folders <- c(experiment_diagnostics$missing_folders, library_folder)
      next  # Skip to next library
    }

    # Filter barcodes for this specific library
    barcodeLibrary <- dplyr::filter(barcodeExperiment, library == libraryId)

    if(nrow(barcodeLibrary) == 0) {
      experiment_diagnostics$missing_barcodes <- c(experiment_diagnostics$missing_barcodes, libraryId)
      next  # Skip to next library
    }

    # Find all Raw_Counts files for this library
    raw_count_files <- list.files(library_folder, pattern = "^Raw_Counts_.*\\.tsv$", full.names = TRUE)

    if(length(raw_count_files) == 0) {
      experiment_diagnostics$missing_raw_files <- c(experiment_diagnostics$missing_raw_files, libraryId)
      next  # Skip to next library
    }

    # Storage for this library's issues
    library_issues <- list()

    # Check each Raw_Counts file
    for(raw_file in raw_count_files) {

      # Extract internal cluster ID from filename
      filename <- basename(raw_file)
      cluster_id <- gsub("Raw_Counts_(.*)\\.tsv", "\\1", filename)

      # Fast column name extraction - read only header (0 rows)
      barcode_cols <- colnames(fread(raw_file, nrows = 0, sep = "\t"))
      
      # Get barcode columns (exclude metadata columns)
      metadata_cols <- c("gene_id", "biotype", "type", "internalClusterId", "cellTypeId")
      barcode_cols <- setdiff(barcode_cols, metadata_cols)

      # Get expected barcodes for this cluster
      expected_barcodes_cluster <- barcodeLibrary$barcode[
        barcodeLibrary$internal_cluster_id == cluster_id &
        barcodeLibrary$cellTypeName != "Unassigned"
      ]

      if(length(expected_barcodes_cluster) == 0) {
        next
      }

      # Compare barcodes (use unique values since setdiff automatically deduplicates)
      found_in_file <- barcode_cols
      total_for_cluster <- expected_barcodes_cluster
      expected_for_cluster <- unique(expected_barcodes_cluster)

      missing_barcodes <- setdiff(expected_for_cluster, found_in_file)
      extra_barcodes <- setdiff(found_in_file, expected_for_cluster)
      matched_barcodes <- intersect(expected_for_cluster, found_in_file)

      # Check for mismatches and store detailed information
      if(length(missing_barcodes) > 0 || length(extra_barcodes) > 0) {
        issue_info <- list(
          file = filename,
          cluster_id = cluster_id,
          missing_count = length(missing_barcodes),
          extra_count = length(extra_barcodes),
          matched_count = length(matched_barcodes),
          missing_barcodes = missing_barcodes,
          extra_barcodes = extra_barcodes,
          total_count = length(total_for_cluster),
          expected_uniq_count = length(expected_for_cluster),
          found_count = length(found_in_file)
        )
        library_issues[[filename]] <- issue_info
      }
    }

    # Add library to list if it has issues
    if(length(library_issues) > 0) {
      experiment_issues[[libraryId]] <- library_issues
    }
  }
  
  # Return experiment issues and diagnostics (will be combined by foreach)
  list(issues = experiment_issues, diagnostics = experiment_diagnostics)
}

## Stop parallel cluster
stopCluster(cl)

## Debug: Check structure of results
message("DEBUG: Length of all_issues: ", length(all_issues))
message("DEBUG: Class of all_issues: ", class(all_issues))
if(length(all_issues) > 0) {
  message("DEBUG: Structure of first element: ", paste(names(all_issues[[1]]), collapse = ", "))
}

## Separate issues and diagnostics from parallel results
all_diagnostics <- list()
all_experiment_issues <- list()
for(i in seq_along(all_issues)) {
  if(!is.null(all_issues[[i]]$diagnostics)) {
    all_diagnostics[[length(all_diagnostics) + 1]] <- all_issues[[i]]$diagnostics
  }
  if(!is.null(all_issues[[i]]$issues) && length(all_issues[[i]]$issues) > 0) {
    all_experiment_issues <- c(all_experiment_issues, all_issues[[i]]$issues)
  }
}

message("DEBUG: Collected ", length(all_diagnostics), " diagnostics")
message("DEBUG: Collected ", length(all_experiment_issues), " library issues")

## Flatten the nested list structure
libraries_with_issues <- all_experiment_issues

## Report processing summary and diagnostics
message("\n========================================")
message("Processing completed for ", length(experiments), " experiments")
message("Total libraries in metadata: ", length(libraries))
message("========================================")

# Report all experiments (brief summary)
for(diag in all_diagnostics) {
  has_issues <- diag$barcode_file_not_found || 
                length(diag$missing_folders) > 0 || 
                length(diag$missing_barcodes) > 0 || 
                length(diag$missing_raw_files) > 0
  
  # Check if this experiment has barcode validation issues
  exp_libs <- names(libraries_with_issues)
  experiment_id <- diag$experiment_id
  libs_with_issues_in_exp <- grep(paste0("^", experiment_id), exp_libs, value = TRUE)
  
  has_validation_issues <- length(libs_with_issues_in_exp) > 0
  
  status_icon <- if(has_issues || has_validation_issues) "⚠" else "✓"
  message(status_icon, " Experiment ", diag$experiment_id, " (", diag$library_count, " libraries)")
}

message("\n========================================")
message("DIAGNOSTIC DETAILS")
message("========================================")

# Report detailed diagnostics only for experiments with issues
experiments_with_diagnostics <- 0
for(diag in all_diagnostics) {
  has_issues <- diag$barcode_file_not_found || 
                length(diag$missing_folders) > 0 || 
                length(diag$missing_barcodes) > 0 || 
                length(diag$missing_raw_files) > 0
  
  # Check if this experiment has barcode validation issues
  exp_libs <- names(libraries_with_issues)
  experiment_id <- diag$experiment_id
  libs_with_validation_issues <- character(0)
  for(lib in exp_libs) {
    # Extract experiment ID from library ID (assuming format like SRX123456)
    if(grepl(paste0("^", substr(experiment_id, 1, 3)), lib)) {
      libs_with_validation_issues <- c(libs_with_validation_issues, lib)
    }
  }
  
  has_validation_issues <- length(libs_with_validation_issues) > 0
  
  if(has_issues || has_validation_issues) {
    experiments_with_diagnostics <- experiments_with_diagnostics + 1
    message("\nExperiment ", diag$experiment_id, ":")
    
    if(diag$barcode_file_not_found) {
      message("  ⚠ Barcode file not found")
    }
    
    if(length(diag$missing_folders) > 0) {
      message("  ⚠ ", length(diag$missing_folders), " library folders not found")
      for(folder in head(diag$missing_folders, 3)) {
        message("    - ", folder)
      }
      if(length(diag$missing_folders) > 3) {
        message("    ... and ", length(diag$missing_folders) - 3, " more")
      }
    }
    
    if(length(diag$missing_barcodes) > 0) {
      message("  ⚠ ", length(diag$missing_barcodes), " libraries with no barcodes: ", 
              paste(head(diag$missing_barcodes, 5), collapse = ", "),
              if(length(diag$missing_barcodes) > 5) paste0(" ... +", length(diag$missing_barcodes) - 5) else "")
    }
    
    if(length(diag$missing_raw_files) > 0) {
      message("  ⚠ ", length(diag$missing_raw_files), " libraries with no raw count files: ", 
              paste(head(diag$missing_raw_files, 5), collapse = ", "),
              if(length(diag$missing_raw_files) > 5) paste0(" ... +", length(diag$missing_raw_files) - 5) else "")
    }
    
    if(has_validation_issues) {
      message("  ⚠ ", length(libs_with_validation_issues), " libraries with barcode validation issues: ",
              paste(head(libs_with_validation_issues, 5), collapse = ", "),
              if(length(libs_with_validation_issues) > 5) paste0(" ... +", length(libs_with_validation_issues) - 5) else "")
    }
  }
}

if(experiments_with_diagnostics == 0) {
  message("\n✓ All experiments processed successfully without diagnostic issues")
} else {
  message("\n", experiments_with_diagnostics, " experiments had diagnostic issues (shown above)")
}
message("========================================\n")

## Output results
if(length(libraries_with_issues) > 0) {
  cat("\n========================================\n")
  cat("BARCODE VALIDATION ISSUES FOUND\n")
  cat("========================================\n\n")

  for(lib in names(libraries_with_issues)) {
    cat("Library:", lib, "\n")
    cat(paste0(rep("-", 60), collapse=""), "\n")

    lib_issues <- libraries_with_issues[[lib]]

    for(file_info in lib_issues) {
      cat("  File:", file_info$file, "\n")
      cat("  Cluster ID:", file_info$cluster_id, "\n")
      cat("  Expected barcodes:", file_info$total_count, "\n")
      cat("  Expected uniq barcodes:", file_info$expected_uniq_count, "\n")
      cat("  Found barcodes:", file_info$found_count, "\n")
      cat("  Matched barcodes:", file_info$matched_count, "\n")

      if(file_info$missing_count > 0) {
        cat("  ⚠ MISSING BARCODES:", file_info$missing_count, "\n")
        if(file_info$missing_count <= 10) {
          cat("    ", paste(file_info$missing_barcodes, collapse = ", "), "\n")
        } else {
          cat("    First 10:", paste(head(file_info$missing_barcodes, 10), collapse = ", "), "\n")
          cat("    ... and", file_info$missing_count - 10, "more\n")
        }
      }

      if(file_info$extra_count > 0) {
        cat("  ⚠ UNEXPECTED BARCODES:", file_info$extra_count, "\n")
        if(file_info$extra_count <= 10) {
          cat("    ", paste(file_info$extra_barcodes, collapse = ", "), "\n")
        } else {
          cat("    First 10:", paste(head(file_info$extra_barcodes, 10), collapse = ", "), "\n")
          cat("    ... and", file_info$extra_count - 10, "more\n")
        }
      }

      cat("\n")
    }
  }

  cat("========================================\n")
  cat("SUMMARY: Total libraries with barcode issues:", length(libraries_with_issues), "\n")
  cat("========================================\n")

  # Also output the warning messages for backwards compatibility
  for(lib in names(libraries_with_issues)) {
    warning("Library with barcode issues: ", lib, call. = FALSE)
  }
  warning("Total libraries with barcode issues: ", length(libraries_with_issues), call. = FALSE)

} else {
  message("\n✓ No libraries with barcode issues found. All barcodes match!")
}
