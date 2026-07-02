# =============================================================================
# ConfigManager.R - JSON-based configuration management
# Looks in docs/config/ (new location) with fallback to config/ (legacy)
# =============================================================================

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
library(jsonlite)

ConfigManager <- new.env()

# Default fallback path
ConfigManager$DEFAULT_MAIN_PATH <- "D:/GithubProjects/ARMIGithub"

# ---- Path detection ----
ConfigManager$detectMainPath <- function() {
  candidates <- c(getwd(), ConfigManager$DEFAULT_MAIN_PATH)
  
  # Walk up from CWD
  current <- getwd()
  for (i in 1:5) {
    candidates <- c(candidates, current)
    current <- dirname(current)
  }
  
  # Look for config in either location
  for (path in unique(candidates)) {
    # Prefer new location: docs/config/
    if (file.exists(file.path(path, "docs", "config", "user_config.json")) ||
        file.exists(file.path(path, "docs", "config", "default_config.json"))) {
      return(path)
    }
  }
  # Fall back to default
  return(ConfigManager$DEFAULT_MAIN_PATH)
}

# ---- Load configuration ----
ConfigManager$load <- function(config_path = NULL) {
  main_path <- ConfigManager$detectMainPath()
  
  if (is.null(config_path)) {
    # NEW LOCATION: docs/config/
    new_user    <- file.path(main_path, "docs", "config", "user_config.json")
    new_default <- file.path(main_path, "docs", "config", "default_config.json")
    
    # LEGACY LOCATION: config/
    old_user    <- file.path(main_path, "config", "user_config.json")
    old_default <- file.path(main_path, "config", "default_config.json")
    
    # Try in priority order, skipping empty files
    candidates <- c(new_user, new_default, old_user, old_default)
    config_path <- NULL
    
    for (cand in candidates) {
      if (file.exists(cand) && file.size(cand) > 50) {
        config_path <- cand
        break
      }
    }
    
    if (is.null(config_path)) {
      stop("No valid config file found. Searched:\n  ",
           paste(candidates, collapse = "\n  "),
           "\n\nAll were either missing or empty (<50 bytes).")
    }
  }
  
  # Validate file size
  if (file.size(config_path) < 50) {
    stop("Config file appears empty or corrupted (size = ",
         file.size(config_path), " bytes): ", config_path)
  }
  
  # Parse JSON with a helpful error message
  config <- tryCatch(
    jsonlite::fromJSON(config_path, simplifyVector = TRUE,
                       simplifyDataFrame = FALSE,
                       simplifyMatrix = FALSE),
    error = function(e) {
      stop("Failed to parse JSON config: ", e$message,
           "\nFile: ", config_path)
    }
  )
  
  # Inject computed main_path
  config$paths$main_path_resolved <- main_path
  
  ConfigManager$current <- config
  ConfigManager$current_path <- config_path
  
  cat("Loaded config:", config_path, "\n")
  return(config)
}

# ---- Save configuration ----
ConfigManager$save <- function(config = NULL, path = NULL) {
  if (is.null(config)) config <- ConfigManager$current
  if (is.null(config)) stop("No config to save")
  
  if (is.null(path)) path <- ConfigManager$current_path
  if (is.null(path)) {
    main_path <- ConfigManager$detectMainPath()
    config_dir <- file.path(main_path, "docs", "config")
    if (!dir.exists(config_dir)) {
      dir.create(config_dir, recursive = TRUE)
    }
    path <- file.path(config_dir, "user_config.json")
  }
  
  # Strip computed fields
  to_save <- config
  to_save$paths$main_path_resolved <- NULL
  
  json_str <- jsonlite::toJSON(to_save, pretty = TRUE, auto_unbox = TRUE,
                               null = "null", na = "null")
  writeLines(json_str, path)
  cat("Config saved to:", path, "\n")
}

# ---- Reset to defaults ----
ConfigManager$resetToDefaults <- function() {
  main_path <- ConfigManager$detectMainPath()
  default_path <- file.path(main_path, "docs", "config", "default_config.json")
  user_path    <- file.path(main_path, "docs", "config", "user_config.json")
  
  if (!file.exists(default_path)) {
    stop("default_config.json not found at: ", default_path)
  }
  
  file.copy(default_path, user_path, overwrite = TRUE)
  cat("Reset to defaults\n")
  ConfigManager$load(user_path)
}

# ---- Get/Set helpers ----
ConfigManager$get <- function(path) {
  keys <- strsplit(path, "\\.")[[1]]
  val <- ConfigManager$current
  for (k in keys) {
    val <- val[[k]]
    if (is.null(val)) return(NULL)
  }
  return(val)
}

ConfigManager$set <- function(path, value) {
  keys <- strsplit(path, "\\.")[[1]]
  cfg <- ConfigManager$current
  
  expr <- "cfg"
  for (k in keys) {
    expr <- paste0(expr, "[['", k, "']]")
  }
  eval(parse(text = paste0(expr, " <- value")))
  
  ConfigManager$current <- cfg
}

# ---- Load species list ----
ConfigManager$loadSpeciesList <- function(list_name = "default") {
  main_path <- ConfigManager$detectMainPath()
  
  candidates <- c(
    file.path(main_path, "docs", "config", "species_lists", paste0(list_name, ".json")),
    file.path(main_path, "docs", "config", "species_lists", "default.json")
  )
  
  for (loc in candidates) {
    if (file.exists(loc)) {
      data <- jsonlite::fromJSON(loc)
      if (is.list(data) && !is.null(data$lists)) {
        if (!is.null(data$lists[[list_name]])) {
          return(data$lists[[list_name]])
        }
        return(data$lists[[1]])
      }
      return(data)
    }
  }
  
  stop("Species list not found: ", list_name)
}

# ---- Apply config to global environment (for compatibility) ----
ConfigManager$applyToGlobals <- function() {
  cfg <- ConfigManager$current
  if (is.null(cfg)) stop("Config not loaded - call ConfigManager$load() first")
  
  # Paths
  MAIN_PATH  <<- cfg$paths$main_path_resolved
  MAXENT_JAR <<- file.path(MAIN_PATH, cfg$paths$maxent_jar)
  
  # Standard derived paths
  SCRIPTS_DIR        <<- file.path(MAIN_PATH, "R", "Scripts")
  ENV_LAYERS_DIR     <<- file.path(MAIN_PATH, "EnvironmentalLayers")
  TIF_DIR            <<- file.path(MAIN_PATH, "TIF")
  ASC_DIR            <<- file.path(MAIN_PATH, "ASC")
  OUTPUT_DIR         <<- file.path(MAIN_PATH, "output")
  BIN_DIR            <<- file.path(MAIN_PATH, "bin")
  LOGS_DIR           <<- file.path(MAIN_PATH, "logs")
  
  # Layers downloaded by R go under EnvironmentalLayers/R/<source>;
  # user-supplied layers go under EnvironmentalLayers/custom/
  ENV_LAYERS_R_DIR      <<- file.path(ENV_LAYERS_DIR, "R")
  ENV_LAYERS_CUSTOM_DIR <<- file.path(ENV_LAYERS_DIR, "custom")
  
  WORLDCLIM_RAW_DIR  <<- file.path(ENV_LAYERS_R_DIR, "WorldClim")
  SPECTRE_RAW_DIR    <<- file.path(ENV_LAYERS_R_DIR, "SPECTRE")
  LANDCOVER_RAW_DIR  <<- file.path(ENV_LAYERS_R_DIR, "Landcover")
  FOOTPRINT_RAW_DIR  <<- file.path(ENV_LAYERS_R_DIR, "Footprint")
  CUSTOM_RAW_DIR     <<- ENV_LAYERS_CUSTOM_DIR
  SHP_SOURCE_DIR     <<- file.path(MAIN_PATH, "SHP")   # Top-level SHP/ folder
  MODELS_DIR         <<- OUTPUT_DIR
  
  # Species
  SPECIES_LIST       <<- cfg$species$species_list
  USE_RAW_GBIF       <<- cfg$species$use_raw_gbif
  TRIM_OCCURRENCES   <<- cfg$species$trim_occurrences
  
  # Data sources
  DOWNLOAD_WORLDCLIM <<- cfg$data_sources$download_worldclim
  DOWNLOAD_SPECTRE   <<- cfg$data_sources$download_spectre
  DOWNLOAD_LANDCOVER <<- cfg$data_sources$download_landcover
  DOWNLOAD_FOOTPRINT <<- cfg$data_sources$download_footprint
  USE_CUSTOM_LAYERS  <<- cfg$data_sources$use_custom_layers
  DOWNLOAD_GBIF      <<- cfg$data_sources$download_gbif
  EXTRACT_IUCN_RANGE <<- cfg$data_sources$extract_iucn_range
  
  # GBIF
  GBIF_MAX_RECORDS    <<- cfg$gbif$max_records
  GBIF_PAGE_SIZE      <<- cfg$gbif$page_size
  GBIF_MIN_CONFIDENCE <<- cfg$gbif$min_confidence
  
  # Resolution
  RESOLUTION_MODE          <<- cfg$resolution$mode
  TARGET_RES_ARCMIN        <<- cfg$resolution$target_arcmin
  USER_TEMPLATE_RASTER     <<- cfg$resolution$user_template_raster
  USER_TEMPLATE_USE_EXTENT <<- cfg$resolution$user_template_use_extent
  
  # Extent
  CLIP_TO_EXTENT <<- cfg$extent$clip_to_extent
  MIN_LON <<- cfg$extent$min_lon
  MAX_LON <<- cfg$extent$max_lon
  MIN_LAT <<- cfg$extent$min_lat
  MAX_LAT <<- cfg$extent$max_lat
  
  # Mask
  USE_MASK     <<- cfg$mask$use_mask
  MASK_MIN_LON <<- cfg$mask$min_lon
  MASK_MAX_LON <<- cfg$mask$max_lon
  MASK_MIN_LAT <<- cfg$mask$min_lat
  MASK_MAX_LAT <<- cfg$mask$max_lat
  
  # Model parameters
  PERM_IMPORTANCE_THRESHOLD <<- cfg$model_parameters$perm_importance_threshold
  TOP_VARIABLE_COUNT        <<- cfg$model_parameters$top_variable_count
  CORRELATION_THRESHOLD     <<- cfg$model_parameters$correlation_threshold
  SELECTION_CRITERION       <<- cfg$model_parameters$selection_criterion
  MAX_COMBO_SIZE            <<- cfg$model_parameters$max_combo_size
  
  # Beta
  BETA_MIN  <<- cfg$beta_optimization$min
  BETA_MAX  <<- cfg$beta_optimization$max
  BETA_STEP <<- cfg$beta_optimization$step
  
  # Final model
  REPLICATES        <<- cfg$final_model$replicates
  REPLICATE_TYPE    <<- cfg$final_model$replicate_type
  PROJECTION_LAYERS <<- cfg$final_model$projection_layers
  
  # Cleanup
  USE_SHP             <<- cfg$cleanup$use_shp
  KEEP_CLIPPED_LAYERS <<- cfg$cleanup$keep_clipped_layers
  KEEP_CLIPPED_MASK   <<- cfg$cleanup$keep_clipped_mask
  
  # Advanced
  JAVA_MEMORY_MB       <<- cfg$advanced$java_memory_mb
  SKIP_EXISTING_MODELS <<- cfg$advanced$skip_existing_models
  REQUIRED_TIFS        <<- cfg$advanced$required_tifs
  
  invisible(NULL)
}

# ---- Per-species path builder ----
get_species_paths <- function(species_tag) {
  run_root <- file.path(OUTPUT_DIR, species_tag)
  list(
    species_tag         = species_tag,
    run_root            = run_root,
    csv_dir             = file.path(run_root, "CSV"),
    models_dir          = file.path(run_root, "Models"),
    shp_dir             = file.path(run_root, "SHP"),
    txt_dir             = file.path(run_root, "TXT"),
    clipped_layers_dir  = file.path(run_root, "ClippedLayers"),
    clipped_mask_dir    = file.path(run_root, "ClippedMask"),
    species_csv         = file.path(run_root, "CSV", paste0(species_tag, ".csv")),
    species_trimmed_csv = file.path(run_root, "CSV", paste0(species_tag, ".trimmed.csv")),
    species_clipped_csv = file.path(run_root, "CSV", paste0(species_tag, ".trimmed.clipped.csv")),
    correlation_csv     = file.path(run_root, "CSV", "correlation.csv"),
    model_selection_csv = file.path(run_root, "CSV", "modelSelection.csv"),
    range_shp           = file.path(run_root, "SHP", paste0(species_tag, "_range.shp"))
  )
}