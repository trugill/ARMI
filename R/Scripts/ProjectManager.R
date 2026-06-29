# =============================================================================
# ProjectManager.R - Orchestrates the 13-step workflow
# =============================================================================

ProjectManager <- new.env()
ProjectManager$bestVariables <- NULL
ProjectManager$bestMultiplier <- 1.0

ProjectManager$runMethodology <- function(config, progress_callback = NULL) {
  
  update_progress <- function(step, name) {
    if (!is.null(progress_callback)) progress_callback(step, name)
    cat("\n[STEP", step, "]", name, "\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
  }
  
  WorkflowConfiguration_validate(config)
  
  # Steps 1-2
  if (config$step1$enabled) {
    update_progress(1, "Loading location data")
    species_file <- config$step1$occurrence_file
    if (config$step1$use_raw_gbif) {
      species_file <- CSVSpeciesExtractor$processCSVFile(species_file, FALSE)
    }
    PathManager$initialize(species_file)
    
    if (config$step1$use_iucn_shp && nchar(config$step1$shp_path) > 0) {
      update_progress(2, "Clipping points to shapefile")
      RFunctions$clipToShp(species_file, config$step1$shp_path)
    }
  }
  
  # Step 3
  if (config$step3$enabled) {
    update_progress(3, "Trimming duplicate occurrence points")
    asc_files <- list.files(PathManager$getLayersPath(), pattern = "\\.asc$", full.names = TRUE)
    if (length(asc_files) == 0) {
      cat("  No .asc files yet - will run after step 4\n")
    } else {
      RFunctions$trimDupes()
      PathManager$trimDupePath()
    }
  }
  
  # Steps 4-5
  if (config$step4$enabled) {
    update_progress(4, "Defining spatial extent and clipping rasters")
    if (config$step4$clip_to_extent) {
      pts <- read.csv(PathManager$getSpeciesPath())
      extent_array <- c(min(pts$longitude) - 1, max(pts$longitude) + 1,
                        min(pts$latitude) - 1, max(pts$latitude) + 1)
    } else {
      extent_array <- c(config$step4$min_lon, config$step4$max_lon,
                        config$step4$min_lat, config$step4$max_lat)
    }
    cat("  Extent:", paste(round(extent_array, 2), collapse = ", "), "\n")
    RFunctions$clipRasters(extent_array = extent_array)
    
    if (config$step3$enabled) {
      species_name <- tools::file_path_sans_ext(basename(PathManager$getSpeciesPath()))
      if (!grepl("trimmed", species_name)) {
        RFunctions$trimDupes()
        PathManager$trimDupePath()
      }
    }
  }
  
  # Step 6
  global_model_path <- NULL
  if (config$step6$enabled) {
    update_progress(6, "Running global Maxent model")
    global_output <- file.path(PathManager$getModelsPath(), "Global")
    if (!dir.exists(global_output)) dir.create(global_output, recursive = TRUE)
    MaxentCaller$runMaxent(PathManager$getSpeciesPath(), 
                            args = config$step6$extra_args,
                            output_dir = global_output)
    global_model_path <- file.path(global_output, "maxentResults.csv")
  }
  
  # Step 7
  if (config$step7$enabled) {
    update_progress(7, "Identifying top variables")
    if (is.null(global_model_path)) {
      global_model_path <- file.path(PathManager$getModelsPath(), "Global", "maxentResults.csv")
    }
    top_vars <- CSVPermutationAnalyzer$getTopPermutationImportance(
      global_model_path, config$step7$threshold, config$step7$top_count)
    cat("  Top variables:", paste(top_vars, collapse = ", "), "\n")
  }
  
  # Step 8
  if (config$step8$enabled) {
    update_progress(8, "Calculating variable correlation")
    if (is.null(global_model_path)) {
      global_model_path <- file.path(PathManager$getModelsPath(), "Global", "maxentResults.csv")
    }
    RFunctions$correlation(global_model_path, config$step7$threshold, config$step7$top_count)
  }
  
  # Step 9
  if (config$step9$enabled) {
    update_progress(9, "Generating valid variable combinations")
    RFunctions$combos(config$step9$correlation_threshold)
  }
  
  # Step 10
  if (config$step10$enabled) {
    update_progress(10, "Running Maxent for all variable combinations")
    combo_folders <- CSVPermutationAnalyzer$getFolderNames(PathManager$getModelsPath())
    combo_folders <- combo_folders[!combo_folders %in% c("Global", "Best_10x", "BetaOptim")]
    cat("  Running", length(combo_folders), "permutations\n")
    CSVPermutationAnalyzer$refreshCSV(PathManager$getModelSelectionPath())
    
    for (i in seq_along(combo_folders)) {
      folder <- combo_folders[i]
      cat("\n  [", i, "/", length(combo_folders), "] ", folder, "\n", sep = "")
      
      if (length(config$step10$required_tifs) > 0) {
        layers_in_combo <- strsplit(folder, "_")[[1]]
        if (!all(config$step10$required_tifs %in% layers_in_combo)) {
          cat("    Skipping (required TIFs missing)\n")
          next
        }
      }
      
      output_dir <- file.path(PathManager$getModelsPath(), folder)
      layer_names <- strsplit(folder, "_")[[1]]
      maxent_args <- MaxentCaller$buildLayerArgs(layer_names)
      MaxentCaller$runMaxent(PathManager$getSpeciesPath(), args = maxent_args, output_dir = output_dir)
      
      species_name <- PathManager$getSpeciesName()
      asc_file <- file.path(output_dir, paste0(species_name, ".asc"))
      lambdas_file <- file.path(output_dir, paste0(species_name, ".lambdas"))
      CSVPermutationAnalyzer$addToFirstEmptyRow(PathManager$getModelSelectionPath(),
        c(PathManager$getSpeciesPath(), asc_file, lambdas_file))
    }
  }
  
  # Step 11
  if (config$step11$enabled) {
    update_progress(11, "Identifying best model by information criterion")
    results_path <- file.path(PathManager$getCSVPath(), "modelSelectionResults.csv")
    RFunctions$processModelSelection(PathManager$getModelSelectionPath(), results_path)
    
    best <- CSVPermutationAnalyzer$findLowestScoreVariables(results_path, config$step11$selection_criterion)
    if (!is.null(best)) {
      ProjectManager$bestVariables <- strsplit(best$variables, "_")[[1]]
      cat("  Best variables:", paste(ProjectManager$bestVariables, collapse = ", "), "\n")
      cat("  Best", config$step11$selection_criterion, "score:", best$score, "\n")
    }
  }
  
  # Step 12
  if (config$step12$enabled) {
    update_progress(12, "Optimizing regularization (beta multiplier)")
    if (is.null(ProjectManager$bestVariables)) stop("No best variables. Run Step 11 first.")
    
    beta_values <- seq(1.0, config$step12$max_beta, by = config$step12$beta_increment)
    beta_dir <- file.path(PathManager$getModelsPath(), "BetaOptim")
    if (!dir.exists(beta_dir)) dir.create(beta_dir, recursive = TRUE)
    beta_csv <- file.path(PathManager$getCSVPath(), "betaSelection.csv")
    CSVPermutationAnalyzer$refreshCSV(beta_csv)
    
    for (beta in beta_values) {
      cat("\n  Testing beta =", beta, "\n")
      output_dir <- file.path(beta_dir, paste0("beta_", beta))
      layer_args <- MaxentCaller$buildLayerArgs(ProjectManager$bestVariables)
      beta_args <- c(layer_args, paste0("betamultiplier=", beta))
      MaxentCaller$runMaxent(PathManager$getSpeciesPath(), args = beta_args, output_dir = output_dir)
      
      species_name <- PathManager$getSpeciesName()
      asc_file <- file.path(output_dir, paste0(species_name, ".asc"))
      lambdas_file <- file.path(output_dir, paste0(species_name, ".lambdas"))
      CSVPermutationAnalyzer$addToFirstEmptyRow(beta_csv,
        c(PathManager$getSpeciesPath(), asc_file, lambdas_file))
    }
    
    beta_results_path <- file.path(PathManager$getCSVPath(), "betaResults.csv")
    RFunctions$processModelSelection(beta_csv, beta_results_path)
    best_beta <- CSVPermutationAnalyzer$findLowestScoreVariables(beta_results_path, config$step12$selection_criterion)
    
    if (!is.null(best_beta)) {
      beta_folder <- basename(dirname(best_beta$path))
      ProjectManager$bestMultiplier <- as.numeric(gsub("beta_", "", beta_folder))
      cat("  Best beta multiplier:", ProjectManager$bestMultiplier, "\n")
    }
  }
  
  # Step 13
  if (config$step13$enabled) {
    update_progress(13, "Running final optimized model")
    if (is.null(ProjectManager$bestVariables)) stop("No best variables. Run Step 11 first.")
    
    final_dir <- file.path(PathManager$getModelsPath(), "Best_10x")
    if (!dir.exists(final_dir)) dir.create(final_dir, recursive = TRUE)
    
    layer_args <- MaxentCaller$buildLayerArgs(ProjectManager$bestVariables)
    final_args <- c(layer_args,
                    paste0("betamultiplier=", ProjectManager$bestMultiplier),
                    paste0("replicates=", config$step13$replicates),
                    "replicatetype=crossvalidate")
    
    for (proj in config$step13$projection_layers) {
      if (dir.exists(proj)) {
        final_args <- c(final_args, paste0("projectionlayers=", proj))
      }
    }
    
    MaxentCaller$runMaxent(PathManager$getSpeciesPath(), args = final_args, output_dir = final_dir)
    cat("\n=== FINAL MODEL COMPLETE ===\n")
    cat("Output:", final_dir, "\n")
  }
  
  if (!is.null(progress_callback)) progress_callback(14, "Complete!")
}

