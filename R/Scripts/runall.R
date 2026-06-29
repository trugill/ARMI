# =============================================================================
# MAXENT-R SYSTEM - ALL-IN-ONE INSTALLER & RUNNER
# =============================================================================
# This script contains ALL R files for the Maxent-R System.
# 
# USAGE:
#   1. Save this file as "MaxentRSystem_AllInOne.R"
#   2. Place it in your desired project directory
#   3. Run: source("MaxentRSystem_AllInOne.R")
#   
# On first run, it will automatically split itself into individual .R files
# in the same directory, then source them and launch the GUI.
# On subsequent runs, it will detect existing files and just launch the GUI.
#
# To force re-extraction, delete the generated .R files or set:
#   FORCE_REEXTRACT <- TRUE before sourcing
# =============================================================================

# ---- Configuration ----
if (!exists("FORCE_REEXTRACT")) FORCE_REEXTRACT <- FALSE

# Determine the directory this script is in
get_script_dir <- function() {
  # Try to get from sys.frames first (when source()'d)
  frames <- sys.frames()
  for (frame in rev(frames)) {
    if (!is.null(frame$ofile)) {
      return(normalizePath(dirname(frame$ofile), mustWork = FALSE))
    }
  }
  # Fall back to commandArgs for Rscript usage
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    path <- sub("^--file=", "", file_arg[1])
    return(normalizePath(dirname(path), mustWork = FALSE))
  }
  # Final fallback: current working directory
  return(normalizePath(getwd(), mustWork = FALSE))
}

SCRIPT_DIR <- get_script_dir()
cat("Working directory:", SCRIPT_DIR, "\n")

# List of files that will be generated
GENERATED_FILES <- c(
  "PathManager.R",
  "CSVPermutationAnalyzer.R",
  "CSVSpeciesExtractor.R",
  "RFunctions.R",
  "MaxentCaller.R",
  "StepConfig.R",
  "ProjectManager.R",
  "MaxentRSystemGUI.R",
  "main.R"
)

# =============================================================================
# FILE CONTENTS (as character vectors)
# =============================================================================

FILE_CONTENTS <- list()

# -----------------------------------------------------------------------------
FILE_CONTENTS[["PathManager.R"]] <- '# =============================================================================
# PathManager.R - Centralized path management
# =============================================================================

PathManager <- new.env()

PathManager$initialize <- function(species_file_path) {
  species_file_path <- normalizePath(species_file_path, mustWork = FALSE)
  
  PathManager$speciesPath <- species_file_path
  PathManager$speciesName <- tools::file_path_sans_ext(basename(species_file_path))
  PathManager$projectPath <- dirname(dirname(species_file_path))
  
  PathManager$csvPath <- file.path(PathManager$projectPath, "CSV")
  PathManager$layersPath <- file.path(PathManager$projectPath, "ClippedLayers")
  PathManager$tifPath <- file.path(PathManager$projectPath, "TIF")
  PathManager$modelsPath <- file.path(PathManager$projectPath, "Models", 
                                      PathManager$speciesName)
  PathManager$correlationPath <- file.path(PathManager$csvPath, "correlation.csv")
  PathManager$modelSelectionPath <- file.path(PathManager$csvPath, "modelSelection.csv")
  PathManager$shapefilePath <- file.path(PathManager$projectPath, "Shapefiles")
  PathManager$maskPath <- file.path(PathManager$projectPath, "Mask")
  PathManager$speciesDir <- file.path(PathManager$projectPath, "Species")
  
  dirs_to_create <- c(PathManager$csvPath, PathManager$layersPath, 
                      PathManager$modelsPath, PathManager$speciesDir)
  for (d in dirs_to_create) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  PathManager$savePaths()
}

PathManager$savePaths <- function() {
  paths_file <- file.path(PathManager$projectPath, "paths.txt")
  paths_content <- c(
    paste0("speciesPath=", PathManager$speciesPath),
    paste0("speciesName=", PathManager$speciesName),
    paste0("projectPath=", PathManager$projectPath),
    paste0("csvPath=", PathManager$csvPath),
    paste0("layersPath=", PathManager$layersPath),
    paste0("tifPath=", PathManager$tifPath),
    paste0("modelsPath=", PathManager$modelsPath),
    paste0("correlationPath=", PathManager$correlationPath),
    paste0("modelSelectionPath=", PathManager$modelSelectionPath),
    paste0("shapefilePath=", PathManager$shapefilePath),
    paste0("maskPath=", PathManager$maskPath)
  )
  writeLines(paths_content, paths_file)
}

PathManager$trimDupePath <- function() {
  base_name <- tools::file_path_sans_ext(basename(PathManager$speciesPath))
  new_name <- paste0(base_name, ".trimmed.csv")
  PathManager$speciesPath <- file.path(dirname(PathManager$speciesPath), new_name)
  PathManager$speciesName <- tools::file_path_sans_ext(new_name)
  PathManager$savePaths()
  return(PathManager$speciesPath)
}

PathManager$getSpeciesPath <- function() PathManager$speciesPath
PathManager$getSpeciesName <- function() PathManager$speciesName
PathManager$getLayersPath <- function() PathManager$layersPath
PathManager$getTifPath <- function() PathManager$tifPath
PathManager$getModelsPath <- function() PathManager$modelsPath
PathManager$getCSVPath <- function() PathManager$csvPath
PathManager$getCorrelationPath <- function() PathManager$correlationPath
PathManager$getModelSelectionPath <- function() PathManager$modelSelectionPath
PathManager$getShapefilePath <- function() PathManager$shapefilePath
PathManager$getMaskPath <- function() PathManager$maskPath
PathManager$getProjectPath <- function() PathManager$projectPath
'

# -----------------------------------------------------------------------------
FILE_CONTENTS[["CSVPermutationAnalyzer.R"]] <- '# =============================================================================
# CSVPermutationAnalyzer.R - CSV analysis utilities
# =============================================================================

CSVPermutationAnalyzer <- new.env()

CSVPermutationAnalyzer$getTopPermutationImportance <- function(csv_file_path, 
                                                                threshold, 
                                                                top_count) {
  if (!file.exists(csv_file_path)) stop("CSV file not found: ", csv_file_path)
  
  data <- read.csv(csv_file_path, stringsAsFactors = FALSE)
  perm_cols <- grep("permutation.importance", names(data), ignore.case = TRUE, value = TRUE)
  
  if (length(perm_cols) == 0) {
    warning("No permutation importance columns found")
    return(character(0))
  }
  
  values <- as.numeric(data[1, perm_cols])
  names(values) <- gsub("\\\\.permutation\\\\.importance", "", perm_cols, ignore.case = TRUE)
  
  threshold_pct <- threshold * 100
  filtered <- values[!is.na(values) & values >= threshold_pct]
  sorted <- sort(filtered, decreasing = TRUE)
  top_n <- head(sorted, top_count)
  return(names(top_n))
}

CSVPermutationAnalyzer$findLowestScoreVariables <- function(csv_file_path, selection_type) {
  if (!file.exists(csv_file_path)) stop("CSV file not found: ", csv_file_path)
  
  data <- read.csv(csv_file_path, stringsAsFactors = FALSE)
  
  score_col <- switch(toupper(selection_type),
                      "AIC" = "aic_score",
                      "AICC" = "aicc_score",
                      "BIC" = "bic_score",
                      "aicc_score")
  
  if (!(score_col %in% names(data))) stop("Score column not found: ", score_col)
  
  scores <- as.numeric(data[[score_col]])
  valid_idx <- which(!is.na(scores))
  
  if (length(valid_idx) == 0) return(NULL)
  
  min_idx <- valid_idx[which.min(scores[valid_idx])]
  path <- data$ascii_file_path[min_idx]
  folder_name <- basename(dirname(path))
  
  return(list(variables = folder_name, score = scores[min_idx], path = path))
}

CSVPermutationAnalyzer$addToFirstEmptyRow <- function(file_path, data) {
  if (!file.exists(file_path)) {
    df <- as.data.frame(t(data), stringsAsFactors = FALSE)
    write.csv(df, file_path, row.names = FALSE)
  } else {
    line <- paste(data, collapse = ",")
    write(line, file = file_path, append = TRUE)
  }
}

CSVPermutationAnalyzer$clearDirectory <- function(directory_path) {
  if (!dir.exists(directory_path)) return(TRUE)
  files <- list.files(directory_path, full.names = TRUE, recursive = FALSE)
  for (f in files) {
    if (dir.exists(f)) unlink(f, recursive = TRUE, force = TRUE)
    else file.remove(f)
  }
  return(TRUE)
}

CSVPermutationAnalyzer$refreshCSV <- function(path) {
  if (file.exists(path)) file.remove(path)
  file.create(path)
}

CSVPermutationAnalyzer$deleteDirectory <- function(path) {
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
}

CSVPermutationAnalyzer$createDirectory <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

CSVPermutationAnalyzer$getFolderNames <- function(directory) {
  if (!dir.exists(directory)) return(character(0))
  list.dirs(directory, full.names = FALSE, recursive = FALSE)
}

CSVPermutationAnalyzer$getFileNames <- function(directory) {
  if (!dir.exists(directory)) return(character(0))
  files <- list.files(directory, full.names = FALSE, recursive = FALSE)
  files[!file.info(file.path(directory, files))$isdir]
}
'

# -----------------------------------------------------------------------------
FILE_CONTENTS[["CSVSpeciesExtractor.R"]] <- '# =============================================================================
# CSVSpeciesExtractor.R - Extract species from GBIF tab-delimited files
# =============================================================================

CSVSpeciesExtractor <- new.env()

CSVSpeciesExtractor$processCSVFile <- function(input_file_path, is_temp = FALSE) {
  if (!file.exists(input_file_path)) stop("Input file does not exist: ", input_file_path)
  
  cat("Processing GBIF file:", input_file_path, "\\n")
  
  data <- tryCatch({
    read.delim(input_file_path, sep = "\\t", stringsAsFactors = FALSE,
               quote = "", fill = TRUE, header = TRUE)
  }, error = function(e) stop("Error reading file: ", e$message))
  
  col_names <- tolower(names(data))
  species_col <- which(grepl("species", col_names))[1]
  lat_col <- which(grepl("decimallatitude|latitude", col_names))[1]
  lon_col <- which(grepl("decimallongitude|longitude", col_names))[1]
  
  if (is.na(species_col) || is.na(lat_col) || is.na(lon_col)) {
    stop("Required columns not found. Need: species, latitude, longitude")
  }
  
  output_dir <- if (is_temp) {
    file.path(dirname(input_file_path), "temp_species")
  } else {
    file.path(dirname(input_file_path), "Species")
  }
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  species_names <- unique(data[[species_col]])
  species_names <- species_names[!is.na(species_names) & nchar(species_names) > 0]
  
  cat("Found", length(species_names), "unique species\\n")
  
  last_output <- NULL
  for (sp in species_names) {
    clean_name <- gsub("[^a-zA-Z0-9]", "_", sp)
    sp_rows <- data[data[[species_col]] == sp & 
                    !is.na(data[[lat_col]]) & 
                    !is.na(data[[lon_col]]), ]
    if (nrow(sp_rows) == 0) next
    
    out_df <- data.frame(
      species = clean_name,
      latitude = as.numeric(sp_rows[[lat_col]]),
      longitude = as.numeric(sp_rows[[lon_col]]),
      stringsAsFactors = FALSE
    )
    out_df <- out_df[complete.cases(out_df), ]
    if (nrow(out_df) == 0) next
    
    output_file <- file.path(output_dir, paste0(clean_name, ".csv"))
    write.csv(out_df, output_file, row.names = FALSE)
    last_output <- output_file
    cat("  Saved:", clean_name, "(", nrow(out_df), "records)\\n")
  }
  
  return(last_output)
}
'

# -----------------------------------------------------------------------------
FILE_CONTENTS[["RFunctions.R"]] <- '# =============================================================================
# RFunctions.R - Native R implementations (replaces ScriptEditor + RIntegration)
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(raster)
  library(sp)
  library(sf)
  library(dplyr)
  library(tools)
})

has_enmtools <- requireNamespace("ENMTools", quietly = TRUE)

RFunctions <- new.env()

# ---- TrimDupes ----
RFunctions$trimDupes <- function(points_file_path = NULL, 
                                  mask_file_path = NULL, 
                                  output_file_path = NULL) {
  if (is.null(points_file_path)) points_file_path <- PathManager$getSpeciesPath()
  
  if (is.null(mask_file_path)) {
    asc_files <- list.files(PathManager$getLayersPath(), pattern = "\\\\.asc$", 
                            full.names = TRUE)
    if (length(asc_files) == 0) stop("No .asc files found for mask")
    mask_file_path <- asc_files[1]
  }
  
  if (is.null(output_file_path)) {
    base_name <- tools::file_path_sans_ext(basename(points_file_path))
    output_file_path <- file.path(dirname(points_file_path),
                                  paste0(base_name, ".trimmed.csv"))
  }
  
  cat("Trimming duplicates from:", points_file_path, "\\n")
  
  original_points <- read.csv(points_file_path, stringsAsFactors = FALSE)
  species_name_from_file <- tools::file_path_sans_ext(basename(points_file_path))
  original_col_names <- colnames(original_points)
  cat("Original number of points:", nrow(original_points), "\\n")
  
  original_points_vect <- terra::vect(original_points, 
                                       geom = c(original_col_names[3], 
                                                original_col_names[2]))
  mask_raster <- terra::rast(mask_file_path)
  
  if (has_enmtools) {
    trimmed_points <- ENMTools::trimdupes.by.raster(points = original_points_vect, 
                                                     mask = mask_raster)
  } else {
    cells <- terra::cells(mask_raster, original_points_vect)[, "cell"]
    keep_idx <- !duplicated(cells) & !is.na(cells)
    trimmed_points <- original_points_vect[keep_idx, ]
  }
  
  coords_from_trimmed <- terra::crds(trimmed_points)
  final_output_df <- data.frame(
    species = species_name_from_file,
    latitude = coords_from_trimmed[, 2],
    longitude = coords_from_trimmed[, 1]
  )
  
  cat("Trimmed number of points:", nrow(final_output_df), "\\n")
  write.csv(final_output_df, file = output_file_path, row.names = FALSE)
  cat("Saved to:", output_file_path, "\\n")
  return(output_file_path)
}

# ---- Correlation ----
RFunctions$correlation <- function(maxent_result_path, threshold, count, 
                                    species_path = NULL) {
  cat("Calculating correlation between top variables...\\n")
  
  top_vars <- CSVPermutationAnalyzer$getTopPermutationImportance(
    maxent_result_path, threshold, count)
  
  if (length(top_vars) == 0) stop("No top variables identified")
  cat("Top variables:", paste(top_vars, collapse = ", "), "\\n")
  
  layers_path <- PathManager$getLayersPath()
  list_file <- file.path(layers_path, paste0(top_vars, ".asc"))
  list_file <- list_file[file.exists(list_file)]
  
  if (length(list_file) == 0) stop("No matching .asc files found")
  
  ac <- raster::stack(list_file)
  ac.brick <- raster::brick(ac)
  ac.brick[ac.brick < -100] <- NA
  rc2.corr <- raster::layerStats(ac.brick, "pearson", na.rm = TRUE)
  
  save_file <- PathManager$getCorrelationPath()
  write.csv(rc2.corr, save_file)
  cat("Correlation matrix saved to:", save_file, "\\n")
  return(save_file)
}

# ---- Combos ----
RFunctions$combos <- function(threshold = 0.8, csv_path = NULL, model_path = NULL) {
  if (is.null(csv_path)) csv_path <- PathManager$getCorrelationPath()
  if (is.null(model_path)) model_path <- PathManager$getModelsPath()
  
  cat("Generating variable combinations from:", csv_path, "\\n")
  
  corr <- read.csv(csv_path, header = TRUE)
  corr <- corr[, -c(1, ncol(corr))]
  colnames(corr) <- gsub("^pearson.correlation.coefficient.", "", colnames(corr))
  
  if (!dir.exists(model_path)) dir.create(model_path, recursive = TRUE)
  
  total_created <- 0
  ncols <- ncol(corr)
  max_combo_size <- min(7, ncols)
  
  for (size in 2:max_combo_size) {
    if (size > ncols) break
    combos <- combn(ncols, size, simplify = FALSE)
    
    for (combo in combos) {
      is_valid <- TRUE
      for (i in 1:(length(combo) - 1)) {
        for (j in (i + 1):length(combo)) {
          if (abs(corr[combo[i], combo[j]]) >= threshold) {
            is_valid <- FALSE
            break
          }
        }
        if (!is_valid) break
      }
      
      if (is_valid) {
        combo_name <- paste(colnames(corr)[combo], collapse = "_")
        combo_dir <- file.path(model_path, combo_name)
        if (!dir.exists(combo_dir)) dir.create(combo_dir, recursive = TRUE)
        combo_csv <- file.path(model_path, paste0("combinations_", size, ".csv"))
        write(combo_name, file = combo_csv, append = TRUE)
        total_created <- total_created + 1
      }
    }
  }
  
  cat("Generated", total_created, "variable combinations\\n")
  return(total_created)
}

# ---- Clip Rasters ----
RFunctions$clipRasters <- function(input_folder = NULL, extent_array, 
                                    output_folder = NULL) {
  if (is.null(input_folder)) input_folder <- PathManager$getTifPath()
  if (is.null(output_folder)) output_folder <- PathManager$getLayersPath()
  
  if (!dir.exists(input_folder)) stop("Input folder does not exist: ", input_folder)
  if (length(extent_array) != 4) stop("Extent array must have 4 values")
  
  left <- extent_array[1]; right <- extent_array[2]
  south <- extent_array[3]; north <- extent_array[4]
  
  if (left >= right || south >= north) stop("Invalid extent")
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  extensions <- c("\\\\.tif$", "\\\\.tiff$", "\\\\.img$", "\\\\.grd$", "\\\\.nc$")
  raster_files <- c()
  for (ext in extensions) {
    raster_files <- c(raster_files, 
                      list.files(input_folder, pattern = ext, 
                                 full.names = TRUE, ignore.case = TRUE))
  }
  
  if (length(raster_files) == 0) stop("No raster files found")
  
  cat("Found", length(raster_files), "raster files\\n")
  clip_extent <- terra::ext(left, right, south, north)
  
  successful <- 0
  for (i in seq_along(raster_files)) {
    input_file <- raster_files[i]
    file_name <- tools::file_path_sans_ext(basename(input_file))
    output_file <- file.path(output_folder, paste0(file_name, ".asc"))
    cat("  Processing", i, "of", length(raster_files), ":", basename(input_file), "\\n")
    
    tryCatch({
      r <- terra::rast(input_file)
      if (is.na(terra::crs(r)) || terra::crs(r) == "") {
        terra::crs(r) <- "EPSG:4326"
      } else if (!grepl("4326|WGS.*84", terra::crs(r))) {
        r <- terra::project(r, "EPSG:4326")
      }
      clipped <- terra::crop(r, clip_extent)
      if (terra::ncell(clipped) == 0) return()
      terra::NAflag(clipped) <- -9999
      terra::writeRaster(clipped, output_file, filetype = "AAIGrid",
                         overwrite = TRUE, NAflag = -9999)
      successful <<- successful + 1
    }, error = function(e) {
      cat("    Error:", e$message, "\\n")
    })
  }
  
  cat("Successfully clipped", successful, "of", length(raster_files), "rasters\\n")
  return(successful)
}

# ---- Clip to Shapefile ----
RFunctions$clipToShp <- function(points_path = NULL, shape_path, output_path = NULL) {
  if (is.null(points_path)) points_path <- PathManager$getSpeciesPath()
  if (is.null(output_path)) {
    base_name <- tools::file_path_sans_ext(basename(points_path))
    output_path <- file.path(dirname(points_path), paste0(base_name, ".clipped.csv"))
  }
  
  if (!file.exists(points_path)) stop("Points file not found: ", points_path)
  if (!file.exists(shape_path)) stop("Shapefile not found: ", shape_path)
  
  cat("Reading points from:", points_path, "\\n")
  points_df <- read.csv(points_path, stringsAsFactors = FALSE)
  
  required_cols <- c("species", "latitude", "longitude")
  if (!all(required_cols %in% names(points_df))) {
    stop("CSV must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  points_df <- points_df[complete.cases(points_df[c("latitude", "longitude")]), ]
  points_sf <- sf::st_as_sf(points_df, coords = c("longitude", "latitude"), crs = 4326)
  
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(TRUE))
  
  polygons_sf <- sf::st_read(shape_path, quiet = TRUE)
  invalid <- !sf::st_is_valid(polygons_sf)
  if (any(invalid)) {
    polygons_sf <- sf::st_make_valid(polygons_sf)
    polygons_sf <- polygons_sf[sf::st_is_valid(polygons_sf), ]
  }
  
  if (sf::st_crs(points_sf) != sf::st_crs(polygons_sf)) {
    points_sf <- sf::st_transform(points_sf, sf::st_crs(polygons_sf))
  }
  
  clipped_points <- sf::st_filter(points_sf, polygons_sf)
  if (sf::st_crs(clipped_points)$input != "EPSG:4326") {
    clipped_points <- sf::st_transform(clipped_points, 4326)
  }
  
  coords <- sf::st_coordinates(clipped_points)
  result_df <- sf::st_drop_geometry(clipped_points)
  result_df$longitude <- coords[, 1]
  result_df$latitude <- coords[, 2]
  result_df <- result_df[, c("species", "latitude", "longitude")]
  
  write.csv(result_df, output_path, row.names = FALSE)
  cat("Clipped", nrow(result_df), "points (saved to:", output_path, ")\\n")
  return(output_path)
}

# ---- Model Selection (AIC/AICc/BIC) ----
RFunctions$read_lambdas <- function(lambdas_path) {
  if (!file.exists(lambdas_path)) stop("Lambdas file not found: ", lambdas_path)
  raw_data <- read.csv(lambdas_path, header = FALSE, stringsAsFactors = FALSE,
                       comment.char = "", quote = "", strip.white = TRUE)
  raw_data$V2_numeric <- suppressWarnings(as.numeric(raw_data$V2))
  metadata_keywords <- c("linearPredictorNormalizer", "densityNormalizer",
                          "numBackgroundPoints", "entropy")
  is_metadata <- raw_data$V1 %in% metadata_keywords
  filtered <- raw_data[!is_metadata & !is.na(raw_data$V2_numeric), ]
  active <- filtered[filtered$V2_numeric != 0.0, ]
  return(active)
}

RFunctions$calculate_loglikelihood <- function(datapoints_path, raster_path) {
  tryCatch({
    datapoints <- read.csv(datapoints_path, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(datapoints) < 3) stop("Need 3+ columns")
    colnames(datapoints)[1:3] <- c("species", "latitude", "longitude")
    datapoints$latitude <- as.numeric(datapoints$latitude)
    datapoints$longitude <- as.numeric(datapoints$longitude)
    datapoints <- datapoints[!is.na(datapoints$latitude) & !is.na(datapoints$longitude), ]
    
    prediction_raster <- raster::raster(raster_path)
    raster_values <- raster::getValues(prediction_raster)
    valid_vals <- raster_values[!is.na(raster_values) & raster_values != -9999]
    probsum <- sum(valid_vals, na.rm = TRUE)
    
    if (probsum <= 0) {
      return(list(loglikelihood = NA, valid_points = 0,
                  total_points = nrow(datapoints), sample_size = 0,
                  probsum = probsum, error = "Probsum non-positive"))
    }
    
    points_sp <- sp::SpatialPoints(
      coords = datapoints[, c("longitude", "latitude")],
      proj4string = sp::CRS(sp::proj4string(prediction_raster))
    )
    extracted <- raster::extract(prediction_raster, points_sp)
    valid_idx <- which(!is.na(extracted) & extracted > 0)
    valid_vals <- extracted[valid_idx]
    
    if (length(valid_vals) == 0) {
      return(list(loglikelihood = 0, valid_points = 0,
                  total_points = nrow(datapoints), sample_size = 0, probsum = probsum))
    }
    
    probabilities <- valid_vals / probsum
    loglikelihood <- sum(log(probabilities))
    
    return(list(loglikelihood = loglikelihood,
                valid_points = length(valid_vals),
                total_points = nrow(datapoints),
                sample_size = length(valid_vals),
                probsum = probsum))
  }, error = function(e) {
    return(list(loglikelihood = NA, valid_points = NA, total_points = NA,
                sample_size = NA, probsum = NA, error = e$message))
  })
}

RFunctions$count_parameters <- function(lambdas_path) {
  tryCatch({
    lambdas <- RFunctions$read_lambdas(lambdas_path)
    return(nrow(lambdas))
  }, error = function(e) return(NA))
}

RFunctions$calculate_information_criteria <- function(ll, k, n) {
  if (is.na(ll) || is.na(k) || is.na(n) || n <= 0) {
    return(list(aic = NA, aicc = NA, bic = NA))
  }
  aic <- -2 * ll + 2 * k
  aicc <- if (n <= k + 1) NA else aic + (2 * k * (k + 1)) / (n - k - 1)
  bic <- -2 * ll + k * log(n)
  return(list(aic = aic, aicc = aicc, bic = bic))
}

RFunctions$processModelSelection <- function(csv_path, output_path) {
  input_data <- read.csv(csv_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(input_data)[1:3] <- c("datapoints_path", "raster_path", "lambdas_path")
  cat("Processing", nrow(input_data), "models for selection\\n")
  
  results <- data.frame(
    points_path = input_data$datapoints_path,
    ascii_file_path = input_data$raster_path,
    loglikelihood = NA_real_, parameter_count = NA_real_,
    sample_size = NA_integer_, aic_score = NA_real_,
    aicc_score = NA_real_, bic_score = NA_real_,
    probsum = NA_real_, valid_points = NA_integer_,
    total_points = NA_integer_, error_message = "",
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(input_data)) {
    cat("  Model", i, "of", nrow(input_data), "\\n")
    res <- RFunctions$calculate_loglikelihood(input_data$datapoints_path[i],
                                               input_data$raster_path[i])
    results$loglikelihood[i] <- res$loglikelihood
    results$sample_size[i] <- ifelse(is.null(res$sample_size), NA, res$sample_size)
    results$valid_points[i] <- ifelse(is.null(res$valid_points), NA, res$valid_points)
    results$total_points[i] <- ifelse(is.null(res$total_points), NA, res$total_points)
    results$probsum[i] <- ifelse(is.null(res$probsum), NA, res$probsum)
    results$error_message[i] <- ifelse(is.null(res$error), "", res$error)
    
    k <- RFunctions$count_parameters(input_data$lambdas_path[i])
    results$parameter_count[i] <- k
    ic <- RFunctions$calculate_information_criteria(res$loglikelihood, k, res$sample_size)
    results$aic_score[i] <- ic$aic
    results$aicc_score[i] <- ic$aicc
    results$bic_score[i] <- ic$bic
  }
  
  write.csv(results, output_path, row.names = FALSE)
  cat("Model selection saved to:", output_path, "\\n")
  return(results)
}
'

# -----------------------------------------------------------------------------
FILE_CONTENTS[["MaxentCaller.R"]] <- '# =============================================================================
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
  
  cat("Running MaxEnt for:", basename(species_path), "\\n")
  cat("Output:", output_dir, "\\n")
  
  result <- tryCatch({
    system2("java", args = cmd_args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    cat("MaxEnt execution error:", e$message, "\\n")
    return(NULL)
  })
  return(output_dir)
}

MaxentCaller$buildLayerArgs <- function(layer_names, layers_path = NULL) {
  if (is.null(layers_path)) layers_path <- PathManager$getLayersPath()
  all_layers <- list.files(layers_path, pattern = "\\\\.asc$", full.names = FALSE)
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
'

# -----------------------------------------------------------------------------
FILE_CONTENTS[["StepConfig.R"]] <- '# =============================================================================
# StepConfig.R - Workflow step configurations
# =============================================================================

StepConfig <- new.env()

StepConfig$Step1Config <- function(enabled = TRUE, occurrence_file = "", 
                                     use_raw_gbif = FALSE, use_iucn_shp = FALSE,
                                     shp_path = "") {
  list(enabled = enabled, occurrence_file = occurrence_file,
       use_raw_gbif = use_raw_gbif, use_iucn_shp = use_iucn_shp,
       shp_path = shp_path, type = "Step1Config")
}

StepConfig$Step3Config <- function(enabled = TRUE) {
  list(enabled = enabled, type = "Step3Config")
}

StepConfig$Step4Config <- function(enabled = TRUE, clip_to_extent = TRUE,
                                     max_lat = NA, max_lon = NA, 
                                     min_lat = NA, min_lon = NA,
                                     use_mask = FALSE,
                                     mask_max_lat = NA, mask_max_lon = NA,
                                     mask_min_lat = NA, mask_min_lon = NA) {
  list(enabled = enabled, clip_to_extent = clip_to_extent,
       max_lat = max_lat, max_lon = max_lon, min_lat = min_lat, min_lon = min_lon,
       use_mask = use_mask,
       mask_max_lat = mask_max_lat, mask_max_lon = mask_max_lon,
       mask_min_lat = mask_min_lat, mask_min_lon = mask_min_lon,
       type = "Step4Config")
}

StepConfig$Step6Config <- function(enabled = TRUE, extra_args = character(0)) {
  list(enabled = enabled, extra_args = extra_args, type = "Step6Config")
}

StepConfig$Step7Config <- function(enabled = TRUE, threshold = 0.05, top_count = 10) {
  list(enabled = enabled, threshold = threshold, top_count = top_count, type = "Step7Config")
}

StepConfig$Step8Config <- function(enabled = TRUE) {
  list(enabled = enabled, type = "Step8Config")
}

StepConfig$Step9Config <- function(enabled = TRUE, correlation_threshold = 0.8) {
  list(enabled = enabled, correlation_threshold = correlation_threshold, type = "Step9Config")
}

StepConfig$Step10Config <- function(enabled = TRUE, required_tifs = character(0)) {
  list(enabled = enabled, required_tifs = required_tifs, type = "Step10Config")
}

StepConfig$Step11Config <- function(enabled = TRUE, selection_criterion = "AICc") {
  list(enabled = enabled, selection_criterion = selection_criterion, type = "Step11Config")
}

StepConfig$Step12Config <- function(enabled = TRUE, max_beta = 5.0, 
                                      beta_increment = 0.5, selection_criterion = "AICc") {
  list(enabled = enabled, max_beta = max_beta, beta_increment = beta_increment,
       selection_criterion = selection_criterion, type = "Step12Config")
}

StepConfig$Step13Config <- function(enabled = TRUE, replicates = 10,
                                      projection_layers = character(0)) {
  list(enabled = enabled, replicates = replicates,
       projection_layers = projection_layers, type = "Step13Config")
}

StepConfig$validate <- function(config) {
  switch(config$type,
    "Step7Config" = {
      if (config$threshold < 0 || config$threshold > 1) stop("Step 7: threshold must be 0-1")
      if (config$top_count < 1) stop("Step 7: top_count must be >= 1")
    },
    "Step9Config" = {
      if (config$correlation_threshold < 0 || config$correlation_threshold > 1)
        stop("Step 9: correlation threshold must be 0-1")
    },
    "Step11Config" = {
      if (!(config$selection_criterion %in% c("AIC", "AICc", "BIC")))
        stop("Step 11: selection criterion must be AIC, AICc, or BIC")
    },
    "Step12Config" = {
      if (config$max_beta < 1) stop("Step 12: max_beta must be >= 1")
      if (config$beta_increment <= 0) stop("Step 12: beta_increment must be > 0")
    },
    "Step13Config" = {
      if (config$replicates < 1) stop("Step 13: replicates must be >= 1")
    }
  )
  return(TRUE)
}

WorkflowConfiguration <- function() {
  list(
    step1 = StepConfig$Step1Config(),
    step3 = StepConfig$Step3Config(),
    step4 = StepConfig$Step4Config(),
    step6 = StepConfig$Step6Config(),
    step7 = StepConfig$Step7Config(),
    step8 = StepConfig$Step8Config(),
    step9 = StepConfig$Step9Config(),
    step10 = StepConfig$Step10Config(),
    step11 = StepConfig$Step11Config(),
    step12 = StepConfig$Step12Config(),
    step13 = StepConfig$Step13Config()
  )
}

WorkflowConfiguration_validate <- function(wf) {
  for (step_name in names(wf)) {
    StepConfig$validate(wf[[step_name]])
  }
  return(TRUE)
}
'

# -----------------------------------------------------------------------------
FILE_CONTENTS[["ProjectManager.R"]] <- '# =============================================================================
# ProjectManager.R - Orchestrates the 13-step workflow
# =============================================================================

ProjectManager <- new.env()
ProjectManager$bestVariables <- NULL
ProjectManager$bestMultiplier <- 1.0

ProjectManager$runMethodology <- function(config, progress_callback = NULL) {
  
  update_progress <- function(step, name) {
    if (!is.null(progress_callback)) progress_callback(step, name)
    cat("\\n[STEP", step, "]", name, "\\n")
    cat(paste(rep("=", 60), collapse = ""), "\\n")
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
    asc_files <- list.files(PathManager$getLayersPath(), pattern = "\\\\.asc$", full.names = TRUE)
    if (length(asc_files) == 0) {
      cat("  No .asc files yet - will run after step 4\\n")
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
    cat("  Extent:", paste(round(extent_array, 2), collapse = ", "), "\\n")
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
    cat("  Top variables:", paste(top_vars, collapse = ", "), "\\n")
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
    cat("  Running", length(combo_folders), "permutations\\n")
    CSVPermutationAnalyzer$refreshCSV(PathManager$getModelSelectionPath())
    
    for (i in seq_along(combo_folders)) {
      folder <- combo_folders[i]
      cat("\\n  [", i, "/", length(combo_folders), "] ", folder, "\\n", sep = "")
      
      if (length(config$step10$required_tifs) > 0) {
        layers_in_combo <- strsplit(folder, "_")[[1]]
        if (!all(config$step10$required_tifs %in% layers_in_combo)) {
          cat("    Skipping (required TIFs missing)\\n")
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
      cat("  Best variables:", paste(ProjectManager$bestVariables, collapse = ", "), "\\n")
      cat("  Best", config$step11$selection_criterion, "score:", best$score, "\\n")
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
      cat("\\n  Testing beta =", beta, "\\n")
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
      cat("  Best beta multiplier:", ProjectManager$bestMultiplier, "\\n")
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
    cat("\\n=== FINAL MODEL COMPLETE ===\\n")
    cat("Output:", final_dir, "\\n")
  }
  
  if (!is.null(progress_callback)) progress_callback(14, "Complete!")
}
'

# -----------------------------------------------------------------------------
FILE_CONTENTS[["MaxentRSystemGUI.R"]] <- '# =============================================================================
# MaxentRSystemGUI.R - Main GUI using tcltk
# =============================================================================

library(tcltk)
library(tcltk2)

create_maxent_gui <- function() {
  
  main_window <- tktoplevel()
  tkwm.title(main_window, "R Automated Maxent-R System")
  tkwm.geometry(main_window, "750x700")
  
  title_frame <- tkframe(main_window)
  title_label <- tklabel(title_frame, text = "R Automated Maxent-R System",
                         font = tkfont.create(family = "Arial", size = 16, weight = "bold"))
  tkpack(title_label, pady = 10)
  tkpack(title_frame, side = "top", fill = "x")
  
  notebook <- tkwidget(main_window, "ttk::notebook")
  
  # TAB 1: Input Data
  tab1 <- tkframe(notebook)
  tkadd(notebook, tab1, text = "  Input Data  ")
  
  occ_frame <- ttklabelframe(tab1, text = "Occurrence Data", padding = 15)
  occ_path_frame <- tkframe(occ_frame)
  occ_label <- tklabel(occ_path_frame, text = "Occurrence File:")
  tkpack(occ_label, side = "left", padx = 5)
  occ_entry <- tkentry(occ_path_frame, width = 40)
  tkpack(occ_entry, side = "left", padx = 5, fill = "x", expand = TRUE)
  browse_btn <- tkbutton(occ_path_frame, text = "Browse...", command = function() {
    fp <- tclvalue(tkgetOpenFile())
    if (nchar(fp) > 0) {
      tkdelete(occ_entry, 0, "end")
      tkinsert(occ_entry, 0, fp)
    }
  })
  tkpack(browse_btn, side = "left", padx = 5)
  tkpack(occ_path_frame, fill = "x", pady = 5)
  
  use_raw_gbif <- tclVar("0")
  raw_check <- tkcheckbutton(occ_frame, text = "Use raw GBIF data", variable = use_raw_gbif)
  tkpack(raw_check, anchor = "w", pady = 5)
  
  use_iucn <- tclVar("0")
  iucn_check <- tkcheckbutton(occ_frame, text = "Use IUCN data (SHP)", variable = use_iucn)
  tkpack(iucn_check, anchor = "w", pady = 5)
  
  shp_frame <- tkframe(occ_frame)
  shp_label <- tklabel(shp_frame, text = "SHP path:")
  tkpack(shp_label, side = "left", padx = 5)
  shp_entry <- tkentry(shp_frame, width = 35)
  tkpack(shp_entry, side = "left", padx = 5, fill = "x", expand = TRUE)
  shp_browse <- tkbutton(shp_frame, text = "Browse...", command = function() {
    fp <- tclvalue(tkgetOpenFile(filetypes = "{{Shapefile} {.shp}}"))
    if (nchar(fp) > 0) {
      tkdelete(shp_entry, 0, "end")
      tkinsert(shp_entry, 0, fp)
    }
  })
  tkpack(shp_browse, side = "left", padx = 5)
  tkpack(shp_frame, fill = "x", pady = 5)
  
  tkpack(occ_frame, fill = "both", expand = TRUE, padx = 20, pady = 20)
  
  # TAB 2: Geographic Extent
  tab2 <- tkframe(notebook)
  tkadd(notebook, tab2, text = "  Geographic Extent  ")
  
  clip_frame <- ttklabelframe(tab2, text = "Clipping Strategy", padding = 15)
  clip_to_extent <- tclVar("1")
  
  coord_frame <- tkframe(clip_frame)
  max_row <- tkframe(coord_frame)
  tkpack(tklabel(max_row, text = "Max Lat:"), side = "left", padx = 5)
  max_lat_entry <- tkentry(max_row, width = 10, state = "disabled")
  tkpack(max_lat_entry, side = "left", padx = 5)
  tkpack(tklabel(max_row, text = "Max Lon:"), side = "left", padx = 15)
  max_lon_entry <- tkentry(max_row, width = 10, state = "disabled")
  tkpack(max_lon_entry, side = "left", padx = 5)
  tkpack(max_row, pady = 5)
  
  min_row <- tkframe(coord_frame)
  tkpack(tklabel(min_row, text = "Min Lat:"), side = "left", padx = 5)
  min_lat_entry <- tkentry(min_row, width = 10, state = "disabled")
  tkpack(min_lat_entry, side = "left", padx = 5)
  tkpack(tklabel(min_row, text = "Min Lon:"), side = "left", padx = 15)
  min_lon_entry <- tkentry(min_row, width = 10, state = "disabled")
  tkpack(min_lon_entry, side = "left", padx = 5)
  tkpack(min_row, pady = 5)
  
  clip_check <- tkcheckbutton(clip_frame, text = "Clip to extent of points",
                               variable = clip_to_extent,
                               command = function() {
                                 state <- if(tclvalue(clip_to_extent) == "0") "normal" else "disabled"
                                 tkconfigure(max_lat_entry, state = state)
                                 tkconfigure(max_lon_entry, state = state)
                                 tkconfigure(min_lat_entry, state = state)
                                 tkconfigure(min_lon_entry, state = state)
                               })
  tkpack(clip_check, anchor = "w", pady = 5)
  tkpack(coord_frame, fill = "x", pady = 10)
  tkpack(clip_frame, fill = "x", padx = 20, pady = 10)
  
  # TAB 3: Model Parameters
  tab3 <- tkframe(notebook)
  tkadd(notebook, tab3, text = "  Model Parameters  ")
  
  var_frame <- ttklabelframe(tab3, text = "Variable Selection", padding = 15)
  
  vi_frame <- tkframe(var_frame)
  tkpack(tklabel(vi_frame, text = "Variable Importance Threshold:"), side = "left", padx = 5)
  var_imp_val <- tclVar("0.05")
  vi_slider <- tkscale(vi_frame, from = 0, to = 1, orient = "horizontal",
                       length = 250, variable = var_imp_val, resolution = 0.01, showvalue = FALSE)
  tkpack(vi_slider, side = "left", padx = 5)
  vi_label <- tklabel(vi_frame, textvariable = var_imp_val, font = tkfont.create(weight = "bold"))
  tkpack(vi_label, side = "left", padx = 5)
  tkpack(vi_frame, fill = "x", pady = 5)
  
  nv_frame <- tkframe(var_frame)
  tkpack(tklabel(nv_frame, text = "Number of Important Variables:"), side = "left", padx = 5)
  num_vars_entry <- tkentry(nv_frame, width = 10)
  tkinsert(num_vars_entry, 0, "10")
  tkpack(num_vars_entry, side = "left", padx = 5)
  tkpack(nv_frame, fill = "x", pady = 5)
  
  rt_frame <- tkframe(var_frame)
  tkpack(tklabel(rt_frame, text = "Required TIFs (comma-sep, optional):"), side = "left", padx = 5)
  req_tif_entry <- tkentry(rt_frame, width = 30)
  tkpack(req_tif_entry, side = "left", padx = 5)
  tkpack(rt_frame, fill = "x", pady = 5)
  
  tkpack(var_frame, fill = "x", padx = 20, pady = 10)
  
  ms_frame <- ttklabelframe(tab3, text = "Model Selection", padding = 15)
  
  crit_frame <- tkframe(ms_frame)
  tkpack(tklabel(crit_frame, text = "Selection Criterion:"), side = "left", padx = 5)
  selection_criterion <- tclVar("AICc")
  tkpack(tkradiobutton(crit_frame, text = "AIC", variable = selection_criterion, value = "AIC"), side = "left", padx = 5)
  tkpack(tkradiobutton(crit_frame, text = "AICc", variable = selection_criterion, value = "AICc"), side = "left", padx = 5)
  tkpack(tkradiobutton(crit_frame, text = "BIC", variable = selection_criterion, value = "BIC"), side = "left", padx = 5)
  tkpack(crit_frame, fill = "x", pady = 5)
  
  corr_frame <- tkframe(ms_frame)
  tkpack(tklabel(corr_frame, text = "Correlation Threshold (r):"), side = "left", padx = 5)
  corr_val <- tclVar("0.80")
  corr_slider <- tkscale(corr_frame, from = 0, to = 1, orient = "horizontal",
                         length = 250, variable = corr_val, resolution = 0.01, showvalue = FALSE)
  tkpack(corr_slider, side = "left", padx = 5)
  corr_label <- tklabel(corr_frame, textvariable = corr_val, font = tkfont.create(weight = "bold"))
  tkpack(corr_label, side = "left", padx = 5)
  tkpack(corr_frame, fill = "x", pady = 5)
  
  beta_frame <- tkframe(ms_frame)
  tkpack(tklabel(beta_frame, text = "Max Beta:"), side = "left", padx = 5)
  max_beta_entry <- tkentry(beta_frame, width = 8)
  tkinsert(max_beta_entry, 0, "5.0")
  tkpack(max_beta_entry, side = "left", padx = 5)
  tkpack(tklabel(beta_frame, text = "Beta Increment:"), side = "left", padx = 15)
  beta_inc_entry <- tkentry(beta_frame, width = 8)
  tkinsert(beta_inc_entry, 0, "0.5")
  tkpack(beta_inc_entry, side = "left", padx = 5)
  tkpack(beta_frame, fill = "x", pady = 5)
  
  rep_frame <- tkframe(ms_frame)
  tkpack(tklabel(rep_frame, text = "Replicates (final model):"), side = "left", padx = 5)
  reps_entry <- tkentry(rep_frame, width = 10)
  tkinsert(reps_entry, 0, "10")
  tkpack(reps_entry, side = "left", padx = 5)
  tkpack(rep_frame, fill = "x", pady = 5)
  
  tkpack(ms_frame, fill = "x", padx = 20, pady = 10)
  
  # TAB 4: Step Selection
  tab4 <- tkframe(notebook)
  tkadd(notebook, tab4, text = "  Step Selection  ")
  
  info_lab <- tklabel(tab4, 
    text = "Select which steps to execute. Note: Some steps depend on previous steps.",
    justify = "left", foreground = "blue")
  tkpack(info_lab, anchor = "w", padx = 20, pady = 10)
  
  step_names <- c(
    "Step 1-2: Get location data and clip to extent",
    "Step 3: Remove duplicates",
    "Step 4-5: Define spatial extent / clip rasters",
    "Step 6: Run global model",
    "Step 7: Identify top variables",
    "Step 8: Calculate correlation",
    "Step 9: Generate permutations",
    "Step 10: Run all permutations",
    "Step 11: Identify top model",
    "Step 12: Optimize regularization",
    "Step 13: Run final model"
  )
  step_keys <- c("step1","step3","step4","step6","step7","step8",
                 "step9","step10","step11","step12","step13")
  step_vars <- list()
  
  btn_frame <- tkframe(tab4)
  tkpack(tkbutton(btn_frame, text = "Select All", command = function() {
    for (v in step_vars) tclvalue(v) <- "1"
  }), side = "left", padx = 5)
  tkpack(tkbutton(btn_frame, text = "Deselect All", command = function() {
    for (v in step_vars) tclvalue(v) <- "0"
  }), side = "left", padx = 5)
  tkpack(btn_frame, anchor = "w", padx = 20, pady = 10)
  
  for (i in seq_along(step_names)) {
    step_vars[[step_keys[i]]] <- tclVar("1")
    chk <- tkcheckbutton(tab4, text = step_names[i], variable = step_vars[[step_keys[i]]])
    tkpack(chk, anchor = "w", padx = 30, pady = 2)
  }
  
  tkpack(notebook, fill = "both", expand = TRUE, padx = 10, pady = 10)
  
  # Progress
  prog_frame <- tkframe(main_window)
  prog_label <- tklabel(prog_frame, text = "Ready", font = tkfont.create(size = 9))
  tkpack(prog_label, anchor = "w", padx = 20)
  progress_bar <- tkwidget(prog_frame, "ttk::progressbar", length = 700, mode = "determinate", maximum = 14)
  tkpack(progress_bar, padx = 20, pady = 5)
  tkpack(prog_frame, side = "top", fill = "x")
  
  update_progress <- function(step, name) {
    tkconfigure(progress_bar, value = step)
    tkconfigure(prog_label, text = sprintf("Step %d: %s", step, name))
    tcl("update")
  }
  
  # Buttons
  button_panel <- tkframe(main_window)
  
  process_btn <- tkbutton(button_panel, text = "Process", 
                          font = tkfont.create(size = 11, weight = "bold"),
                          command = function() {
    occ_file <- tclvalue(tkget(occ_entry))
    if (nchar(occ_file) == 0) {
      tkmessageBox(title = "Error", message = "Please specify the occurrence file.", icon = "error")
      return()
    }
    
    config <- WorkflowConfiguration()
    config$step1 <- StepConfig$Step1Config(
      enabled = tclvalue(step_vars$step1) == "1",
      occurrence_file = occ_file,
      use_raw_gbif = tclvalue(use_raw_gbif) == "1",
      use_iucn_shp = tclvalue(use_iucn) == "1",
      shp_path = tclvalue(tkget(shp_entry))
    )
    config$step3 <- StepConfig$Step3Config(enabled = tclvalue(step_vars$step3) == "1")
    
    clip_ext <- tclvalue(clip_to_extent) == "1"
    config$step4 <- StepConfig$Step4Config(
      enabled = tclvalue(step_vars$step4) == "1",
      clip_to_extent = clip_ext,
      max_lat = if (clip_ext) NA else as.numeric(tclvalue(tkget(max_lat_entry))),
      max_lon = if (clip_ext) NA else as.numeric(tclvalue(tkget(max_lon_entry))),
      min_lat = if (clip_ext) NA else as.numeric(tclvalue(tkget(min_lat_entry))),
      min_lon = if (clip_ext) NA else as.numeric(tclvalue(tkget(min_lon_entry)))
    )
    config$step6 <- StepConfig$Step6Config(enabled = tclvalue(step_vars$step6) == "1")
    config$step7 <- StepConfig$Step7Config(
      enabled = tclvalue(step_vars$step7) == "1",
      threshold = as.numeric(tclvalue(var_imp_val)),
      top_count = as.integer(tclvalue(tkget(num_vars_entry)))
    )
    config$step8 <- StepConfig$Step8Config(enabled = tclvalue(step_vars$step8) == "1")
    config$step9 <- StepConfig$Step9Config(
      enabled = tclvalue(step_vars$step9) == "1",
      correlation_threshold = as.numeric(tclvalue(corr_val))
    )
    
    req_tifs_str <- tclvalue(tkget(req_tif_entry))
    req_tifs <- if (nchar(req_tifs_str) > 0) trimws(strsplit(req_tifs_str, ",")[[1]]) else character(0)
    
    config$step10 <- StepConfig$Step10Config(
      enabled = tclvalue(step_vars$step10) == "1", required_tifs = req_tifs)
    config$step11 <- StepConfig$Step11Config(
      enabled = tclvalue(step_vars$step11) == "1",
      selection_criterion = tclvalue(selection_criterion))
    config$step12 <- StepConfig$Step12Config(
      enabled = tclvalue(step_vars$step12) == "1",
      max_beta = as.numeric(tclvalue(tkget(max_beta_entry))),
      beta_increment = as.numeric(tclvalue(tkget(beta_inc_entry))),
      selection_criterion = tclvalue(selection_criterion))
    config$step13 <- StepConfig$Step13Config(
      enabled = tclvalue(step_vars$step13) == "1",
      replicates = as.integer(tclvalue(tkget(reps_entry))))
    
    answer <- tkmessageBox(title = "Confirm",
                           message = "Begin processing with these settings?",
                           icon = "question", type = "okcancel")
    if (tclvalue(answer) != "ok") return()
    
    result <- tryCatch({
      ProjectManager$runMethodology(config, update_progress)
      "success"
    }, error = function(e) paste("Error:", e$message))
    
    if (result == "success") {
      tkmessageBox(title = "Complete", message = "Workflow completed successfully!", icon = "info")
    } else {
      tkmessageBox(title = "Error", message = result, icon = "error")
    }
  })
  tkpack(process_btn, side = "left", padx = 5)
  
  help_btn <- tkbutton(button_panel, text = "Help", command = function() {
    help_text <- paste(
      "R Automated Maxent-R System Help\\n\\n",
      "INPUT DATA TAB:\\n",
      "  - Specify your species occurrence CSV file\\n",
      "  - Check Use raw GBIF data if it is untrimmed GBIF data\\n",
      "  - Optionally clip to IUCN SHP boundary\\n\\n",
      "GEOGRAPHIC EXTENT TAB:\\n",
      "  - Auto-clip to point extent, or specify bounds manually\\n\\n",
      "MODEL PARAMETERS TAB:\\n",
      "  - Variable importance threshold (0-1)\\n",
      "  - Number of top variables to consider\\n",
      "  - Correlation threshold for combinations\\n",
      "  - Selection criterion: AIC, AICc, or BIC\\n",
      "  - Beta multiplier optimization range\\n\\n",
      "STEP SELECTION TAB:\\n",
      "  - Toggle individual workflow steps on/off\\n",
      sep = "")
    
    hw <- tktoplevel()
    tkwm.title(hw, "Help")
    ht <- tktext(hw, width = 70, height = 25, wrap = "word")
    tkinsert(ht, "end", help_text)
    tkconfigure(ht, state = "disabled")
    tkpack(ht, padx = 10, pady = 10)
  })
  tkpack(help_btn, side = "left", padx = 5)
  
  exit_btn <- tkbutton(button_panel, text = "Exit", command = function() tkdestroy(main_window))
  tkpack(exit_btn, side = "left", padx = 5)
  
  tkpack(button_panel, side = "bottom", pady = 15)
  tkfocus(main_window)
  invisible(main_window)
}
'

# -----------------------------------------------------------------------------
FILE_CONTENTS[["main.R"]] <- '# =============================================================================
# main.R - Application entry point
# =============================================================================

# Source all modules in correct order
script_dir <- tryCatch({
  if (!is.null(sys.frames()[[1]]$ofile)) {
    dirname(normalizePath(sys.frames()[[1]]$ofile))
  } else getwd()
}, error = function(e) getwd())

source(file.path(script_dir, "PathManager.R"))
source(file.path(script_dir, "CSVPermutationAnalyzer.R"))
source(file.path(script_dir, "CSVSpeciesExtractor.R"))
source(file.path(script_dir, "RFunctions.R"))
source(file.path(script_dir, "MaxentCaller.R"))
source(file.path(script_dir, "StepConfig.R"))
source(file.path(script_dir, "ProjectManager.R"))
source(file.path(script_dir, "MaxentRSystemGUI.R"))

# Configure MaxEnt jar path - EDIT THIS FOR YOUR SYSTEM
MAXENT_JAR_PATH <- "C:/maxent/maxent.jar"
if (file.exists(MAXENT_JAR_PATH)) {
  MaxentCaller$setMaxentJar(MAXENT_JAR_PATH)
} else {
  cat("WARNING: MaxEnt jar not found at:", MAXENT_JAR_PATH, "\\n")
  cat("Set the correct path using: MaxentCaller$setMaxentJar(\\"path/to/maxent.jar\\")\\n")
}

cat("\\n========================================\\n")
cat("R Automated Maxent-R System\\n")
cat("========================================\\n\\n")
cat("Launching GUI...\\n")

create_maxent_gui()
'

# =============================================================================
# EXTRACTION LOGIC
# =============================================================================

extract_files <- function() {
  cat("\n=== Extracting R files ===\n")
  for (filename in names(FILE_CONTENTS)) {
    filepath <- file.path(SCRIPT_DIR, filename)
    writeLines(FILE_CONTENTS[[filename]], filepath)
    cat("  Created:", filename, "\n")
  }
  cat("Extraction complete!\n\n")
}

# Check if files already exist
files_exist <- all(file.exists(file.path(SCRIPT_DIR, GENERATED_FILES)))

if (!files_exist || FORCE_REEXTRACT) {
  if (FORCE_REEXTRACT) {
    cat("FORCE_REEXTRACT is TRUE - regenerating all files\n")
  } else {
    cat("First run detected - extracting component files\n")
  }
  extract_files()
} else {
  cat("All component files found - skipping extraction\n")
  cat("(Set FORCE_REEXTRACT <- TRUE before sourcing to regenerate)\n\n")
}

# =============================================================================
# LAUNCH APPLICATION
# =============================================================================

cat("=== Loading modules ===\n")

# Source each file in the correct dependency order
source(file.path(SCRIPT_DIR, "PathManager.R"))
source(file.path(SCRIPT_DIR, "CSVPermutationAnalyzer.R"))
source(file.path(SCRIPT_DIR, "CSVSpeciesExtractor.R"))
source(file.path(SCRIPT_DIR, "RFunctions.R"))
source(file.path(SCRIPT_DIR, "MaxentCaller.R"))
source(file.path(SCRIPT_DIR, "StepConfig.R"))
source(file.path(SCRIPT_DIR, "ProjectManager.R"))
source(file.path(SCRIPT_DIR, "MaxentRSystemGUI.R"))

# Configure MaxEnt jar - edit this for your system
MAXENT_JAR_PATH <- "C:/maxent/maxent.jar"
if (file.exists(MAXENT_JAR_PATH)) {
  MaxentCaller$setMaxentJar(MAXENT_JAR_PATH)
  cat("MaxEnt jar configured:", MAXENT_JAR_PATH, "\n")
} else {
  cat("\nWARNING: MaxEnt jar not found at:", MAXENT_JAR_PATH, "\n")
  cat("Update the MAXENT_JAR_PATH variable or call:\n")
  cat("  MaxentCaller$setMaxentJar(\"your/path/to/maxent.jar\")\n\n")
}

cat("\n========================================\n")
cat(" R Automated Maxent-R System\n")
cat("========================================\n\n")
cat("Launching GUI...\n\n")

create_maxent_gui()