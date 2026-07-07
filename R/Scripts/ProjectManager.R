# =============================================================================
# ProjectManager.R - Orchestrates the 13-step workflow
# =============================================================================
ProjectManager <- new.env()
ProjectManager$bestVariables <- NULL
ProjectManager$bestMultiplier <- 1.0

# =============================================================================
# WorkflowConfiguration - constructor for the master config object
# =============================================================================
WorkflowConfiguration <- function() {
  list(
    step0  = NULL,
    step1  = NULL,
    step3  = NULL,
    step4  = NULL,
    step6  = NULL,
    step7  = NULL,
    step8  = NULL,
    step9  = NULL,
    step10 = NULL,
    step11 = NULL,
    step12 = NULL,
    step13 = NULL
  )
}

# =============================================================================
# WorkflowConfiguration_validate - sanity check each step config
# =============================================================================
WorkflowConfiguration_validate <- function(config) {
  if (is.null(config) || !is.list(config)) {
    stop("WorkflowConfiguration_validate: config must be a list")
  }
  
  step_names <- c("step0", "step1", "step3", "step4", "step6",
                  "step7", "step8", "step9", "step10", "step11",
                  "step12", "step13")
  
  for (sn in step_names) {
    step <- config[[sn]]
    if (is.null(step)) next  # missing step is OK -- treated as disabled
    
    if (exists("StepConfig") && is.function(StepConfig$validate)) {
      tryCatch(
        StepConfig$validate(step),
        error = function(e) {
          stop(sprintf("Validation failed for %s: %s", sn, e$message))
        }
      )
    }
  }
  
  return(TRUE)
}

ProjectManager$runMethodology <- function(config, progress_callback = NULL) {
  
  update_progress <- function(step, name) {
    if (!is.null(progress_callback)) progress_callback(step, name)
    cat("\n[STEP", step, "]", name, "\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
  }
  
  WorkflowConfiguration_validate(config)
  
  # Uniform helper: treats missing/null step configs as disabled
  step_on <- function(step) {
    !is.null(step) && isTRUE(step$enabled)
  }
  
  # =========================================================================
  # Step 0 - Data download (environmental layers + GBIF)
  # =========================================================================
  if (step_on(config$step0)) {
    update_progress(0, "Downloading environmental data")
    
    if (isTRUE(config$step0$download_worldclim)) {
      cat("  Downloading WorldClim...\n")
      tryCatch(
        DataDownloader$downloadWorldClim(config$step0$target_res_arcmin),
        error = function(e) cat("  WorldClim skipped:", e$message, "\n")
      )
    }
    if (isTRUE(config$step0$download_spectre)) {
      cat("  Downloading SPECTRE...\n")
      tryCatch(
        DataDownloader$downloadSPECTRE(config$step0$target_res_arcmin),
        error = function(e) cat("  SPECTRE skipped:", e$message, "\n")
      )
    }
    if (isTRUE(config$step0$download_landcover)) {
      cat("  Downloading Land Cover...\n")
      tryCatch(
        DataDownloader$downloadLandcover(config$step0$target_res_arcmin),
        error = function(e) cat("  LandCover skipped:", e$message, "\n")
      )
    }
    if (isTRUE(config$step0$download_footprint)) {
      cat("  Downloading Human Footprint...\n")
      tryCatch(
        DataDownloader$downloadFootprint(config$step0$target_res_arcmin),
        error = function(e) cat("  Footprint skipped:", e$message, "\n")
      )
    }
    if (isTRUE(config$step0$use_custom_layers)) {
      cat("  Processing custom layers (EnvironmentalLayers/custom)...\n")
      tryCatch(
        DataDownloader$processCustomLayers(config$step0$target_res_arcmin),
        error = function(e) cat("  Custom layers skipped:", e$message, "\n")
      )
    }
    
    cat("  Converting resampled TIF layers to ASCII (TIF/ -> ASC/)...\n")
    tryCatch(
      DataDownloader$convertTifToAsc(),
      error = function(e) cat("  TIF -> ASC conversion skipped:", e$message, "\n")
    )
    
    # GBIF download + save as tab-delimited (raw GBIF format)
    if (isTRUE(config$step0$download_gbif) &&
        !is.null(config$step1) && nchar(config$step1$occurrence_file) > 0) {
      
      species_tag  <- config$step1$occurrence_file
      species_name <- gsub("_", " ", species_tag)
      
      species_dir <- file.path(MAIN_PATH, "output", species_tag, "CSV")
      csv_path    <- file.path(species_dir, paste0(species_tag, ".csv"))
      
      skip_existing <- isTRUE(ConfigManager$current$gbif$skip_existing)
      
      if (skip_existing && file.exists(csv_path) && file.size(csv_path) > 100) {
        cat("  GBIF: existing CSV found for", species_name,
            "- skipping re-download.\n")
        cat("    ", csv_path, "\n")
        
        if (isTRUE(config$step0$extract_iucn_range)) {
          shp_source_dir <- file.path(MAIN_PATH, "SHP", "IUCN")
          out_shp_dir    <- file.path(MAIN_PATH, "output", species_tag, "SHP")
          range_shp      <- file.path(out_shp_dir,
                                      paste0(species_tag, "_range.shp"))
          
          if (file.exists(range_shp)) {
            cat("  IUCN range shapefile already exists - skipping extraction.\n")
          } else if (dir.exists(shp_source_dir)) {
            cat("  Extracting IUCN range polygon (no cached GBIF taxonomy",
                "- using name as-is)...\n")
            SHPExtractor$extractSpeciesRange(
              species_tag    = species_tag,
              sci_name       = species_name,
              group          = NA_character_,
              shp_source_dir = shp_source_dir,
              out_shp_dir    = out_shp_dir
            )
          } else {
            cat("  IUCN range extraction requested but source dir not found:\n")
            cat("    ", shp_source_dir, "\n")
          }
        }
        
      } else {
        if (skip_existing && file.exists(csv_path)) {
          cat("  GBIF: existing CSV for", species_name,
              "is too small (", file.size(csv_path),
              " bytes) - re-downloading.\n")
        }
        
        cat("  Downloading GBIF data for:", species_name, "\n")
        
        gbif_result <- tryCatch(
          DataDownloader$downloadGBIF(
            species_tag = species_name,
            max_records = ConfigManager$current$gbif$max_records,
            page_size   = ConfigManager$current$gbif$page_size
          ),
          error = function(e) {
            cat("  GBIF download failed:", e$message, "\n")
            NULL
          }
        )
        
        if (!is.null(gbif_result) &&
            !is.null(gbif_result$data) &&
            nrow(gbif_result$data) > 0) {
          
          if (!dir.exists(species_dir)) dir.create(species_dir, recursive = TRUE)
          
          df <- gbif_result$data
          
          if (!"species" %in% names(df)) {
            if ("scientificName" %in% names(df)) df$species <- df$scientificName
            else df$species <- gbif_result$taxon$matched_name %||% species_name
          }
          if (!"decimalLatitude" %in% names(df) && "latitude" %in% names(df)) {
            df$decimalLatitude <- df$latitude
          }
          if (!"decimalLongitude" %in% names(df) && "longitude" %in% names(df)) {
            df$decimalLongitude <- df$longitude
          }
          
          tryCatch({
            write.table(df, csv_path,
                        sep = "\t", row.names = FALSE, quote = FALSE, na = "")
            cat("  Saved", nrow(df),
                "records (tab-delimited) to:\n    ", csv_path, "\n")
          }, error = function(e) {
            cat("  Failed to save CSV:", e$message, "\n")
          })
          
          if (isTRUE(config$step0$extract_iucn_range)) {
            shp_source_dir <- file.path(MAIN_PATH, "SHP", "IUCN")
            out_shp_dir    <- file.path(MAIN_PATH, "output", species_tag, "SHP")
            
            if (dir.exists(shp_source_dir)) {
              cat("  Extracting IUCN range polygon...\n")
              
              group <- gbif_result$taxon$group %||% NA_character_
              match_name <- gbif_result$taxon$canonical_name %||%
                gbif_result$taxon$matched_name  %||% species_name
              
              SHPExtractor$extractSpeciesRange(
                species_tag    = species_tag,
                sci_name       = match_name,
                group          = group,
                shp_source_dir = shp_source_dir,
                out_shp_dir    = out_shp_dir
              )
            } else {
              cat("  IUCN range extraction requested but source dir not found:\n")
              cat("    ", shp_source_dir, "\n")
            }
          }
          
        } else {
          cat("  GBIF returned no records for", species_name, "\n")
        }
      }
    }
  }  # end Step 0
  
  # =========================================================================
  # Steps 1-2 - Load location data, optional SHP clip
  # =========================================================================
  if (step_on(config$step1)) {
    update_progress(1, "Loading location data")
    
    species_tag <- config$step1$occurrence_file
    
    species_file <- file.path(MAIN_PATH, "output", species_tag, "CSV",
                              paste0(species_tag, ".csv"))
    
    if (!file.exists(species_file)) {
      stop("Occurrence CSV not found at: ", species_file,
           "\n  (Did Step 0 GBIF download succeed? Or place a manual CSV there.)")
    }
    
    cat("  Using CSV:", species_file, "\n")
    
    if (config$step1$use_raw_gbif) {
      processed <- CSVSpeciesExtractor$processCSVFile(species_file, FALSE)
      if (!is.null(processed) && is.character(processed) && file.exists(processed)) {
        species_file <- processed
        cat("  Processed species CSV:", species_file, "\n")
      }
    }
    
    PathManager$initialize(species_file)
    
    # ---- Wipe Models/ folder ONLY if not reusing anything ----
    reusing_anything <-
      isTRUE(ConfigManager$current$advanced$skip_global_if_exists) ||
      isTRUE(ConfigManager$current$advanced$reuse_clipped_layers) ||
      isTRUE(ConfigManager$current$advanced$reuse_clipped_mask)
    
    models_dir <- PathManager$getModelsPath()
    if (dir.exists(models_dir) && !reusing_anything) {
      cat("  Wiping existing Models/ folder for clean run:\n    ", models_dir, "\n")
      unlink(models_dir, recursive = TRUE, force = TRUE)
    } else if (dir.exists(models_dir)) {
      cat("  Preserving existing Models/ folder (reuse flags set)\n")
    }
    dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (!reusing_anything) {
      for (stale in c(PathManager$getCorrelationPath(),
                      PathManager$getModelSelectionPath(),
                      file.path(PathManager$getCSVPath(), "modelSelectionResults.csv"),
                      file.path(PathManager$getCSVPath(), "betaSelection.csv"),
                      file.path(PathManager$getCSVPath(), "betaResults.csv"))) {
        if (file.exists(stale)) {
          cat("  Removing stale:", basename(stale), "\n")
          file.remove(stale)
        }
      }
    }
    
    if (config$step1$use_iucn_shp && nchar(config$step1$shp_path) > 0) {
      update_progress(2, "Clipping points to shapefile")
      RFunctions$clipToShp(species_file, config$step1$shp_path)
    }
  }
  
  # =========================================================================
  # Step 3 - Trim duplicate occurrence points
  # =========================================================================
  if (step_on(config$step3)) {
    update_progress(3, "Trimming duplicate occurrence points")
    asc_files <- list.files(ASC_DIR, pattern = "\\.asc$", full.names = TRUE)
    if (length(asc_files) == 0) {
      cat("  No resampled .asc layers in", ASC_DIR, "yet - will retry after Step 4\n")
    } else {
      use_range_clip <- isTRUE(config$step1$use_iucn_shp) ||
        isTRUE(ConfigManager$current$species$use_iucn_range)
      
      if (use_range_clip) {
        species_tag <- basename(PathManager$getProjectPath())
        range_shp <- file.path(PathManager$getShapefilePath(),
                               paste0(species_tag, "_range.shp"))
        
        if (file.exists(range_shp)) {
          cat("  Range polygon found - clipping .trimmed.csv to polygon + cell dedup\n")
          RFunctions$trimDupesWithRange(range_shp_path = range_shp)
        } else {
          cat("  'Clip occurrences outside of range polygon' is enabled but no ",
              "range shapefile exists at:\n    ", range_shp, "\n",
              "  Falling back to cell dedup only.\n", sep = "")
          RFunctions$trimDupes()
        }
      } else {
        RFunctions$trimDupes()
      }
      PathManager$trimDupePath()
    }
  }
  
  # =========================================================================
  # Steps 4-5 - Clip rasters to study extent (with reuse support)
  # =========================================================================
  if (step_on(config$step4)) {
    update_progress(4, "Defining spatial extent and clipping rasters")
    
    layers_dir <- PathManager$getLayersPath()
    mask_dir   <- PathManager$getClippedMaskPath()
    
    reuse_layers <- isTRUE(config$step4$reuse_existing_layers) &&
      dir.exists(layers_dir) &&
      length(list.files(layers_dir, pattern = "\\.asc$")) > 0
    
    reuse_mask <- isTRUE(config$step4$reuse_existing_mask) &&
      dir.exists(mask_dir) &&
      length(list.files(mask_dir, pattern = "\\.asc$")) > 0
    
    if (reuse_layers) {
      n_existing <- length(list.files(layers_dir, pattern = "\\.asc$"))
      cat("  Reusing existing ClippedLayers/ (", n_existing,
          " .asc files) - skipping training-extent clip\n", sep = "")
    }
    if (reuse_mask) {
      n_existing <- length(list.files(mask_dir, pattern = "\\.asc$"))
      cat("  Reusing existing ClippedMask/ (", n_existing,
          " .asc files) - skipping mask clip\n", sep = "")
    }
    
    need_training_clip <- !reuse_layers
    need_mask_clip     <- isTRUE(config$step4$use_mask) && !reuse_mask
    
    if (need_training_clip || need_mask_clip) {
      
      # Determine strategy
      strategy <- if (isTRUE(config$step4$clip_to_extent)) "auto_extent"
      else if (!is.null(config$step4$max_lat) && !is.na(config$step4$max_lat)) "manual_bbox"
      else "iucn_range"
      
      # Read species points (needed for auto_extent & iucn_range fallback)
      pts <- read.csv(PathManager$getSpeciesPath())
      lon_col <- if ("longitude" %in% names(pts)) "longitude"
      else if ("decimalLongitude" %in% names(pts)) "decimalLongitude" else NA
      lat_col <- if ("latitude" %in% names(pts)) "latitude"
      else if ("decimalLatitude" %in% names(pts)) "decimalLatitude" else NA
      
      point_extent <- function() {
        if (is.na(lon_col) || is.na(lat_col)) {
          stop("Could not find longitude/latitude columns in ",
               PathManager$getSpeciesPath())
        }
        c(min(pts[[lon_col]], na.rm = TRUE) - 1,
          max(pts[[lon_col]], na.rm = TRUE) + 1,
          min(pts[[lat_col]], na.rm = TRUE) - 1,
          max(pts[[lat_col]], na.rm = TRUE) + 1)
      }
      
      extent_array <- if (strategy == "manual_bbox") {
        c(config$step4$min_lon, config$step4$max_lon,
          config$step4$min_lat, config$step4$max_lat)
      } else if (strategy == "iucn_range") {
        species_tag <- basename(dirname(PathManager$getCSVPath()))
        range_shp <- file.path(MAIN_PATH, "output", species_tag, "SHP",
                               paste0(species_tag, "_range.shp"))
        
        if (file.exists(range_shp)) {
          shp_extent <- SHPExtractor$getExtentFromShapefile(range_shp, buffer_deg = 1.0)
          if (!is.null(shp_extent)) {
            cat("  Using extent from IUCN range shapefile:\n    ", range_shp, "\n")
            shp_extent
          } else {
            cat("  IUCN shapefile exists but extent read failed -- falling back to point extent.\n")
            point_extent()
          }
        } else {
          cat("  No IUCN range shapefile at:\n    ", range_shp, "\n")
          cat("  Falling back to point extent.\n")
          point_extent()
        }
      } else {
        point_extent()
      }
      
      cat("  Strategy:", strategy, "\n")
      cat("  Extent:", paste(round(extent_array, 2), collapse = ", "), "\n")
      
      input_folder <- ASC_DIR
      active_dirs <- ConfigManager$current$data_sources$active_layer_dirs
      if (!is.null(active_dirs) && length(active_dirs) > 0) {
        cat("  active_layer_dirs configured:", paste(active_dirs, collapse = ", "),
            "(only applies to nested raw folders; ASC/ is flat)\n")
      }
      cat("  Using all layers in:", input_folder, "\n")
      
      if (!dir.exists(input_folder) || length(list.files(input_folder, pattern = "\\.asc$")) == 0) {
        stop("No resampled ASCII layers found in: ", input_folder,
             "\n  Run Step 0 first (download/resample layers -> TIF/ -> ASC/),",
             "\n  or place a manual .tif under EnvironmentalLayers/custom and re-run Step 0.")
      }
      
      if (need_training_clip) {
        RFunctions$clipRasters(
          input_folder = input_folder,
          extent_array = extent_array,
          output_folder = layers_dir
        )
      }
      
      if (need_mask_clip) {
        mask_extent <- c(config$step4$mask_min_lon, config$step4$mask_max_lon,
                         config$step4$mask_min_lat, config$step4$mask_max_lat)
        if (any(is.na(mask_extent))) {
          cat("  use_mask is TRUE but mask_min/max_lat/lon are incomplete - skipping mask clip\n")
        } else {
          cat("  Clipping layers to mask extent (for projection) ->", mask_dir, "\n")
          RFunctions$clipRasters(
            input_folder = input_folder,
            extent_array = mask_extent,
            output_folder = mask_dir
          )
        }
      }
    }  # end (need_training_clip || need_mask_clip)
    
    # Post-step-4 duplicate trimming (only if Step 3 also enabled and not already trimmed)
    if (step_on(config$step3)) {
      species_name <- tools::file_path_sans_ext(basename(PathManager$getSpeciesPath()))
      if (!grepl("trimmed", species_name)) {
        use_range_clip <- isTRUE(config$step1$use_iucn_shp) ||
          isTRUE(ConfigManager$current$species$use_iucn_range)
        
        if (use_range_clip) {
          species_tag <- basename(PathManager$getProjectPath())
          range_shp <- file.path(PathManager$getShapefilePath(),
                                 paste0(species_tag, "_range.shp"))
          if (file.exists(range_shp)) {
            RFunctions$trimDupesWithRange(range_shp_path = range_shp)
          } else {
            RFunctions$trimDupes()
          }
        } else {
          RFunctions$trimDupes()
        }
        PathManager$trimDupePath()
      }
    }
  }  # end Step 4
  
  # =========================================================================
  # Step 6 - Global Maxent model (with skip-if-exists support)
  # =========================================================================
  global_model_path <- file.path(PathManager$getModelsPath(), "Global", "maxentResults.csv")
  
  if (step_on(config$step6)) {
    global_output <- file.path(PathManager$getModelsPath(), "Global")
    
    # Only skip if the checkbox is checked AND the file exists
    skip_global <- isTRUE(config$step6$skip_if_exists) &&
      file.exists(global_model_path) &&
      file.size(global_model_path) > 100
    
    if (skip_global) {
      update_progress(6, "Reusing existing global Maxent model")
      cat("  Found existing global model:\n    ", global_model_path, "\n")
      cat("  Skipping rerun (skip_global_if_exists = TRUE)\n")
    } else {
      update_progress(6, "Running global Maxent model")
      if (!dir.exists(global_output)) dir.create(global_output, recursive = TRUE)
      MaxentCaller$runMaxent(PathManager$getSpeciesPath(),
                             args = config$step6$extra_args,
                             output_dir = global_output)
    }
  }
  
  # =========================================================================
  # Step 7 - Identify top variables
  # =========================================================================
  if (step_on(config$step7)) {
    update_progress(7, "Identifying top variables")
    top_vars <- CSVPermutationAnalyzer$getTopPermutationImportance(
      global_model_path, config$step7$threshold, config$step7$top_count)
    cat("  Top variables:", paste(top_vars, collapse = ", "), "\n")
  }
  
  # =========================================================================
  # Step 8 - Calculate correlation
  # =========================================================================
  if (step_on(config$step8)) {
    update_progress(8, "Calculating variable correlation")
    thr <- if (step_on(config$step7)) config$step7$threshold else 0.05
    cnt <- if (step_on(config$step7)) config$step7$top_count else 10
    RFunctions$correlation(global_model_path, thr, cnt)
  }
  
  # =========================================================================
  # Step 9 - Generate variable combinations
  # =========================================================================
  if (step_on(config$step9)) {
    update_progress(9, "Generating valid variable combinations")
    RFunctions$combos(config$step9$correlation_threshold)
  }
  
  # =========================================================================
  # Step 10 - Run all permutations
  # =========================================================================
  if (step_on(config$step10)) {
    update_progress(10, "Running Maxent for all variable combinations")
    combo_folders <- CSVPermutationAnalyzer$getFolderNames(PathManager$getModelsPath())
    combo_folders <- combo_folders[!combo_folders %in% c("Global", "Best_10x", "BetaOptim")]
    cat("  Running", length(combo_folders), "permutations\n")
    CSVPermutationAnalyzer$refreshCSV(PathManager$getModelSelectionPath())
    
    full_layer_dir <- PathManager$getLayersPath()
    name_map <- RFunctions$buildLayerNameMap(full_layer_dir)
    
    for (i in seq_along(combo_folders)) {
      folder <- combo_folders[i]
      cat("\n  [", i, "/", length(combo_folders), "] ", folder, "\n", sep = "")
      
      combo_dir_path <- file.path(PathManager$getModelsPath(), folder)
      layer_names <- RFunctions$getComboVariables(combo_dir_path, folder)
      cat("    Raw combo variables:", paste(layer_names, collapse = ", "), "\n")
      
      if (length(config$step10$required_tifs) > 0) {
        if (!all(config$step10$required_tifs %in% layer_names)) {
          cat("    Skipping (required TIFs missing)\n")
          next
        }
      }
      
      wanted_ids <- character(0)
      unresolved <- character(0)
      for (ln in layer_names) {
        src <- RFunctions$resolveLayerFile(ln, name_map, full_layer_dir)
        if (is.null(src)) {
          unresolved <- c(unresolved, ln)
        } else {
          wanted_ids <- c(wanted_ids,
                          tools::file_path_sans_ext(basename(src)))
        }
      }
      
      if (length(unresolved) > 0) {
        cat("    WARNING: unresolved layers:", paste(unresolved, collapse = ", "), "\n")
        cat("    Skipping this combo\n")
        next
      }
      
      output_dir <- file.path(PathManager$getModelsPath(), folder)
      
      all_asc <- list.files(full_layer_dir, pattern = "\\.asc$",
                            full.names = FALSE)
      all_layer_ids <- tools::file_path_sans_ext(all_asc)
      unwanted_ids  <- setdiff(all_layer_ids, wanted_ids)
      
      toggle_args <- as.vector(rbind("-N", unwanted_ids))
      
      maxent_args <- c(
        paste0("environmentallayers=", full_layer_dir),
        "prefixes=false",
        toggle_args
      )
      
      cat("    Wanted:  ", paste(wanted_ids, collapse = ", "), "\n")
      cat("    Disabled:", paste(unwanted_ids, collapse = ", "), "\n")
      
      MaxentCaller$runMaxent(PathManager$getSpeciesPath(),
                             args = maxent_args,
                             output_dir = output_dir)
      
      asc_files <- list.files(output_dir, pattern = "\\.asc$", full.names = TRUE)
      asc_files <- asc_files[!grepl("_clamping|_novel", basename(asc_files))]
      
      if (length(asc_files) == 0) {
        cat("    WARNING: no .asc output produced in", output_dir, "\n")
        next
      }
      asc_file <- asc_files[1]
      
      lambdas_files <- list.files(output_dir, pattern = "\\.lambdas$", full.names = TRUE)
      if (length(lambdas_files) == 0) {
        cat("    WARNING: no .lambdas output produced in", output_dir, "\n")
        next
      }
      lambdas_file <- lambdas_files[1]
      
      CSVPermutationAnalyzer$addToFirstEmptyRow(
        PathManager$getModelSelectionPath(),
        c(PathManager$getSpeciesPath(), asc_file, lambdas_file)
      )
    }
  }
  
  # =========================================================================
  # Step 11 - Identify best model
  # =========================================================================
  if (step_on(config$step11)) {
    update_progress(11, "Identifying best model by information criterion")
    results_path <- file.path(PathManager$getCSVPath(), "modelSelectionResults.csv")
    RFunctions$processModelSelection(PathManager$getModelSelectionPath(), results_path)
    
    best <- CSVPermutationAnalyzer$findLowestScoreVariables(results_path, config$step11$selection_criterion)
    if (!is.null(best)) {
      best_combo_dir <- dirname(best$path)
      ProjectManager$bestVariables <- RFunctions$getComboVariables(best_combo_dir, best$variables)
      cat("  Best variables:", paste(ProjectManager$bestVariables, collapse = ", "), "\n")
      cat("  Best", config$step11$selection_criterion, "score:", best$score, "\n")
    }
  }
  
  # =========================================================================
  # Step 12 - Beta optimization
  # =========================================================================
  if (step_on(config$step12)) {
    update_progress(12, "Optimizing regularization (beta multiplier)")
    if (is.null(ProjectManager$bestVariables)) stop("No best variables. Run Step 11 first.")
    
    beta_values <- seq(1.0, config$step12$max_beta, by = config$step12$beta_increment)
    beta_dir <- file.path(PathManager$getModelsPath(), "BetaOptim")
    if (!dir.exists(beta_dir)) dir.create(beta_dir, recursive = TRUE)
    beta_csv <- file.path(PathManager$getCSVPath(), "betaSelection.csv")
    CSVPermutationAnalyzer$refreshCSV(beta_csv)
    
    full_layer_dir <- PathManager$getLayersPath()
    name_map <- RFunctions$buildLayerNameMap(full_layer_dir)
    
    wanted_ids <- character(0)
    unresolved <- character(0)
    for (ln in ProjectManager$bestVariables) {
      src <- RFunctions$resolveLayerFile(ln, name_map, full_layer_dir)
      if (is.null(src)) {
        unresolved <- c(unresolved, ln)
      } else {
        wanted_ids <- c(wanted_ids, tools::file_path_sans_ext(basename(src)))
      }
    }
    if (length(unresolved) > 0) {
      stop("Step 12: could not resolve variables: ", paste(unresolved, collapse = ", "))
    }
    
    all_asc <- list.files(full_layer_dir, pattern = "\\.asc$", full.names = FALSE)
    all_layer_ids <- tools::file_path_sans_ext(all_asc)
    unwanted_ids  <- setdiff(all_layer_ids, wanted_ids)
    
    cat("  Wanted:  ", paste(wanted_ids, collapse = ", "), "\n")
    cat("  Disabled:", paste(unwanted_ids, collapse = ", "), "\n")
    
    toggle_args <- as.vector(rbind("-N", unwanted_ids))
    
    base_args <- c(
      paste0("environmentallayers=", full_layer_dir),
      "prefixes=false",
      toggle_args
    )
    
    for (beta in beta_values) {
      cat("\n  Testing beta =", beta, "\n")
      output_dir <- file.path(beta_dir, paste0("beta_", beta))
      beta_args <- c(base_args, paste0("betamultiplier=", beta))
      
      MaxentCaller$runMaxent(PathManager$getSpeciesPath(),
                             args = beta_args,
                             output_dir = output_dir)
      
      asc_files <- list.files(output_dir, pattern = "\\.asc$", full.names = TRUE)
      asc_files <- asc_files[!grepl("_clamping|_novel", basename(asc_files))]
      if (length(asc_files) == 0) {
        cat("    WARNING: no .asc output produced in", output_dir, "\n")
        next
      }
      asc_file <- asc_files[1]
      
      lambdas_files <- list.files(output_dir, pattern = "\\.lambdas$", full.names = TRUE)
      if (length(lambdas_files) == 0) {
        cat("    WARNING: no .lambdas output produced in", output_dir, "\n")
        next
      }
      lambdas_file <- lambdas_files[1]
      
      CSVPermutationAnalyzer$addToFirstEmptyRow(beta_csv,
                                                c(PathManager$getSpeciesPath(),
                                                  asc_file, lambdas_file))
    }
    
    beta_results_path <- file.path(PathManager$getCSVPath(), "betaResults.csv")
    RFunctions$processModelSelection(beta_csv, beta_results_path)
    best_beta <- CSVPermutationAnalyzer$findLowestScoreVariables(
      beta_results_path, config$step12$selection_criterion)
    
    if (!is.null(best_beta)) {
      beta_folder <- basename(dirname(best_beta$path))
      ProjectManager$bestMultiplier <- as.numeric(gsub("beta_", "", beta_folder))
      cat("  Best beta multiplier:", ProjectManager$bestMultiplier, "\n")
    }
  }
  
  # =========================================================================
  # Step 13 - Final model
  # =========================================================================
  if (step_on(config$step13)) {
    update_progress(13, "Running final optimized model")
    if (is.null(ProjectManager$bestVariables)) stop("No best variables. Run Step 11 first.")
    
    final_dir <- file.path(PathManager$getModelsPath(), "Best_10x")
    if (!dir.exists(final_dir)) dir.create(final_dir, recursive = TRUE)
    
    full_layer_dir <- PathManager$getLayersPath()
    name_map <- RFunctions$buildLayerNameMap(full_layer_dir)
    
    wanted_ids <- character(0)
    unresolved <- character(0)
    for (ln in ProjectManager$bestVariables) {
      src <- RFunctions$resolveLayerFile(ln, name_map, full_layer_dir)
      if (is.null(src)) {
        unresolved <- c(unresolved, ln)
      } else {
        wanted_ids <- c(wanted_ids, tools::file_path_sans_ext(basename(src)))
      }
    }
    if (length(unresolved) > 0) {
      stop("Step 13: could not resolve variables: ", paste(unresolved, collapse = ", "))
    }
    
    all_asc <- list.files(full_layer_dir, pattern = "\\.asc$", full.names = FALSE)
    all_layer_ids <- tools::file_path_sans_ext(all_asc)
    unwanted_ids  <- setdiff(all_layer_ids, wanted_ids)
    
    cat("  Wanted:  ", paste(wanted_ids, collapse = ", "), "\n")
    cat("  Disabled:", paste(unwanted_ids, collapse = ", "), "\n")
    
    toggle_args <- as.vector(rbind("-N", unwanted_ids))
    
    final_args <- c(
      paste0("environmentallayers=", full_layer_dir),
      "prefixes=false",
      toggle_args,
      paste0("betamultiplier=", ProjectManager$bestMultiplier),
      paste0("replicates=", config$step13$replicates),
      "replicatetype=crossvalidate"
    )
    
    proj_layers <- config$step13$projection_layers
    if (length(proj_layers) == 0 &&
        dir.exists(PathManager$getClippedMaskPath()) &&
        length(list.files(PathManager$getClippedMaskPath(), pattern = "\\.asc$")) > 0) {
      cat("  No projection_layers configured - defaulting to ClippedMask/:",
          PathManager$getClippedMaskPath(), "\n")
      proj_layers <- PathManager$getClippedMaskPath()
    }
    
    for (proj in proj_layers) {
      if (dir.exists(proj)) {
        final_args <- c(final_args, paste0("projectionlayers=", proj))
      }
    }
    
    MaxentCaller$runMaxent(PathManager$getSpeciesPath(),
                           args = final_args,
                           output_dir = final_dir)
    
    cat("\n=== FINAL MODEL COMPLETE ===\n")
    cat("Output:", final_dir, "\n")
  }
  
  if (!is.null(progress_callback)) progress_callback(14, "Complete!")
}