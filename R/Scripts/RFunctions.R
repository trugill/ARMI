# =============================================================================
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
    asc_files <- list.files(PathManager$getLayersPath(), pattern = "\\.asc$", 
                            full.names = TRUE)
    if (length(asc_files) == 0) stop("No .asc files found for mask")
    mask_file_path <- asc_files[1]
  }
  
  if (is.null(output_file_path)) {
    base_name <- tools::file_path_sans_ext(basename(points_file_path))
    output_file_path <- file.path(dirname(points_file_path),
                                  paste0(base_name, ".trimmed.csv"))
  }
  
  cat("Trimming duplicates from:", points_file_path, "\n")
  
  original_points <- read.csv(points_file_path, stringsAsFactors = FALSE)
  species_name_from_file <- tools::file_path_sans_ext(basename(points_file_path))
  original_col_names <- colnames(original_points)
  cat("Original number of points:", nrow(original_points), "\n")
  
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
  
  cat("Trimmed number of points:", nrow(final_output_df), "\n")
  write.csv(final_output_df, file = output_file_path, row.names = FALSE)
  cat("Saved to:", output_file_path, "\n")
  return(output_file_path)
}

# ---- Correlation ----
RFunctions$correlation <- function(maxent_result_path, threshold, count, 
                                   species_path = NULL) {
  cat("Calculating correlation between top variables...\n")
  
  top_vars <- CSVPermutationAnalyzer$getTopPermutationImportance(
    maxent_result_path, threshold, count)
  
  if (length(top_vars) == 0) stop("No top variables identified")
  cat("Top variables:", paste(top_vars, collapse = ", "), "\n")
  
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
  cat("Correlation matrix saved to:", save_file, "\n")
  return(save_file)
}

# ---- Combos ----
RFunctions$combos <- function(threshold = 0.8, csv_path = NULL, model_path = NULL) {
  if (is.null(csv_path)) csv_path <- PathManager$getCorrelationPath()
  if (is.null(model_path)) model_path <- PathManager$getModelsPath()
  
  cat("Generating variable combinations from:", csv_path, "\n")
  
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
  
  cat("Generated", total_created, "variable combinations\n")
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
  
  extensions <- c("\\.tif$", "\\.tiff$", "\\.img$", "\\.grd$", "\\.nc$", "\\.asc$")
  raster_files <- c()
  for (ext in extensions) {
    raster_files <- c(raster_files, list.files(input_folder,
                                               pattern = ext, full.names = TRUE, ignore.case = TRUE, recursive = TRUE))
  }
  
  # Honor active_layer_dirs from config, if set: keep only files under those subfolders.
  # NOTE: TIF_DIR/ASC_DIR are flat (no per-source subfolders) since layers land
  # there after resampling, so this only has something to match against when
  # clipping straight from a nested raw folder. If it matches nothing, don't
  # silently zero out every raster -- just skip the filter.
  active_dirs <- tryCatch(ConfigManager$current$data_sources$active_layer_dirs,
                          error = function(e) NULL)
  if (!is.null(active_dirs) && length(active_dirs) > 0) {
    keep_patterns <- paste0("/", active_dirs, "/", collapse = "|")
    # Normalize path separators for regex matching
    norm_files <- gsub("\\\\", "/", raster_files)
    matched <- raster_files[grepl(keep_patterns, norm_files)]
    if (length(matched) > 0) {
      raster_files <- matched
      cat("Filtered to", length(raster_files), "files under active layer folders\n")
    } else {
      cat("  active_layer_dirs set but input folder is flat (", input_folder,
          ") - ignoring subfolder filter, using all files found\n")
    }
  }
  
  if (length(raster_files) == 0) stop("No raster files found")
  
  cat("Found", length(raster_files), "raster files\n")
  clip_extent <- terra::ext(left, right, south, north)
  
  successful <- 0
  for (i in seq_along(raster_files)) {
    input_file <- raster_files[i]
    file_name <- tools::file_path_sans_ext(basename(input_file))
    output_file <- file.path(output_folder, paste0(file_name, ".asc"))
    cat("  Processing", i, "of", length(raster_files), ":", basename(input_file), "\n")
    
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
      cat("    Error:", e$message, "\n")
    })
  }
  
  cat("Successfully clipped", successful, "of", length(raster_files), "rasters\n")
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
  
  cat("Reading points from:", points_path, "\n")
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
  cat("Clipped", nrow(result_df), "points (saved to:", output_path, ")\n")
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
  cat("Processing", nrow(input_data), "models for selection\n")
  
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
    cat("  Model", i, "of", nrow(input_data), "\n")
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
  cat("Model selection saved to:", output_path, "\n")
  return(results)
}