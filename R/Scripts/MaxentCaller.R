# =============================================================================
# MaxentCaller.R - Wrapper for executing MaxEnt via java -jar
# =============================================================================

MaxentCaller <- new.env()

MaxentCaller$maxent_jar <- "C:/maxent/maxent.jar"
MaxentCaller$java_memory <- "2048"

MaxentCaller$setMaxentJar <- function(path) {
  MaxentCaller$maxent_jar <- path
}

MaxentCaller$runMaxent <- function(species_path, args = character(0), output_dir = NULL) {
  if (!file.exists(MaxentCaller$maxent_jar)) {
    stop("MaxEnt jar not found: ", MaxentCaller$maxent_jar)
  }
  
  if (is.null(output_dir)) {
    output_dir <- file.path(PathManager$getModelsPath(),
                            tools::file_path_sans_ext(basename(species_path)))
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cmd_args <- c(
    paste0("-mx", MaxentCaller$java_memory, "m"),
    "-jar", MaxentCaller$maxent_jar,
    paste0("samplesfile=", species_path),
    paste0("environmentallayers=", PathManager$getLayersPath()),
    paste0("outputdirectory=", output_dir),
    "redoifexists=true", "autorun=true",
    "visible=false", "warnings=false"
  )
  cmd_args <- c(cmd_args, args)
  
  cat("Running MaxEnt for:", basename(species_path), "\n")
  cat("Output:", output_dir, "\n")
  
  # Write the exact command to a file so we can inspect / replay
  cmd_log <- file.path(output_dir, "maxent_command.txt")
  writeLines(c("java", cmd_args), cmd_log)
  
  result <- tryCatch({
    system2("java", args = cmd_args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    cat("  MaxEnt execution error:", e$message, "\n")
    return(NULL)
  })
  
  # Save full stdout/stderr
  if (!is.null(result)) {
    log_file <- file.path(output_dir, "maxent_stdout.log")
    writeLines(as.character(result), log_file)
    
    err_lines <- grep("error|exception|failed|cannot|missing",
                      result, ignore.case = TRUE, value = TRUE)
    if (length(err_lines) > 0) {
      cat("  !! MaxEnt reported problems:\n")
      for (ln in head(err_lines, 10)) cat("     ", ln, "\n")
    }
  }
  
  # Verify an .asc actually appeared
  ascs <- list.files(output_dir, pattern = "\\.asc$", full.names = FALSE)
  ascs <- ascs[!grepl("_clamping|_novel", ascs)]
  if (length(ascs) == 0) {
    cat("  !! No prediction .asc found in", output_dir, "\n")
    cat("     See:", cmd_log, "\n")
    cat("     See:", file.path(output_dir, "maxent_stdout.log"), "\n")
  }
  
  return(output_dir)
}


MaxentCaller$buildLayerArgs <- function(layer_names, layers_path = NULL) {
  if (is.null(layers_path)) layers_path <- PathManager$getLayersPath()
  all_layers <- list.files(layers_path, pattern = "\\.asc$", full.names = FALSE)
  all_layer_names <- tools::file_path_sans_ext(all_layers)
  
  args <- character(0)
  for (layer in all_layer_names) {
    if (!(layer %in% layer_names)) {
      args <- c(args, paste0("togglelayertype=", layer))
    }
  }
  return(args)
}

MaxentCaller$getFolderNames <- function(directory) {
  CSVPermutationAnalyzer$getFolderNames(directory)
}

MaxentCaller$getFileNames <- function(directory) {
  CSVPermutationAnalyzer$getFileNames(directory)
}

