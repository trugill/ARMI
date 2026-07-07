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
    # Resolution is set once, at resampling time (Step 0) -- every layer in
    # ASC_DIR shares the same grid, so any one of them can define the cells
    # for deduping points. Preferring the master ASC_DIR (rather than the
    # per-species ClippedLayers) means this can run as soon as Step 0 has
    # produced layers, without waiting on Step 4's per-species clip.
    asc_dir <- if (exists("ASC_DIR", envir = .GlobalEnv)) get("ASC_DIR", envir = .GlobalEnv) else NULL
    asc_files <- if (!is.null(asc_dir) && dir.exists(asc_dir)) {
      sort(list.files(asc_dir, pattern = "\\.asc$", full.names = TRUE))
    } else {
      character(0)
    }
    
    if (length(asc_files) == 0) {
      # Fall back to the per-species clipped layers (e.g. if Step 0 hasn't
      # produced a master grid but Step 4 has already clipped some layers).
      asc_files <- sort(list.files(PathManager$getLayersPath(), pattern = "\\.asc$",
                                   full.names = TRUE))
    }
    if (length(asc_files) == 0) {
      stop("No .asc files found for mask (checked ASC/ and ClippedLayers/). ",
           "Run Step 0 to produce resampled layers first.")
    }
    mask_file_path <- asc_files[1]
    cat("Using raster for cell resolution:", mask_file_path, "\n")
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

# ---- TrimDupes with Range Polygon Clipping ----
# Clips occurrences to a range polygon AND resamples to one point per cell.
# The output .trimmed.csv contains ONLY points inside the polygon, deduped by cell.
RFunctions$trimDupesWithRange <- function(points_file_path = NULL,
                                          range_shp_path = NULL,
                                          mask_file_path = NULL,
                                          output_file_path = NULL) {
  if (is.null(points_file_path)) points_file_path <- PathManager$getSpeciesPath()
  
  if (is.null(range_shp_path) || !file.exists(range_shp_path)) {
    stop("Range shapefile not found: ", range_shp_path)
  }
  
  if (is.null(mask_file_path)) {
    asc_dir <- if (exists("ASC_DIR", envir = .GlobalEnv)) get("ASC_DIR", envir = .GlobalEnv) else NULL
    asc_files <- if (!is.null(asc_dir) && dir.exists(asc_dir)) {
      sort(list.files(asc_dir, pattern = "\\.asc$", full.names = TRUE))
    } else {
      character(0)
    }
    if (length(asc_files) == 0) {
      asc_files <- sort(list.files(PathManager$getLayersPath(),
                                   pattern = "\\.asc$", full.names = TRUE))
    }
    if (length(asc_files) == 0) {
      stop("No .asc files found for cell resolution reference.")
    }
    mask_file_path <- asc_files[1]
    cat("Using raster for cell resolution:", mask_file_path, "\n")
  }
  
  if (is.null(output_file_path)) {
    base_name <- tools::file_path_sans_ext(basename(points_file_path))
    output_file_path <- file.path(dirname(points_file_path),
                                  paste0(base_name, ".trimmed.csv"))
  }
  
  cat("Clipping occurrences to range polygon:\n")
  cat("  Points:    ", points_file_path, "\n")
  cat("  Range SHP: ", range_shp_path, "\n")
  
  # ---- Read points ----
  original_points <- read.csv(points_file_path, stringsAsFactors = FALSE)
  species_name_from_file <- tools::file_path_sans_ext(basename(points_file_path))
  cat("  Original point count:", nrow(original_points), "\n")
  
  # Identify lat/lon columns
  cn <- names(original_points)
  lat_col <- if ("latitude" %in% cn) "latitude"
  else if ("decimalLatitude" %in% cn) "decimalLatitude"
  else cn[2]
  lon_col <- if ("longitude" %in% cn) "longitude"
  else if ("decimalLongitude" %in% cn) "decimalLongitude"
  else cn[3]
  
  # Drop rows with missing coordinates
  original_points <- original_points[
    !is.na(original_points[[lat_col]]) & !is.na(original_points[[lon_col]]), ]
  
  if (nrow(original_points) == 0) {
    stop("No valid coordinates in points file.")
  }
  
  # ---- Clip to polygon ----
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(TRUE))
  
  points_sf <- sf::st_as_sf(original_points,
                            coords = c(lon_col, lat_col),
                            crs = 4326)
  
  polygons_sf <- sf::st_read(range_shp_path, quiet = TRUE)
  if (any(!sf::st_is_valid(polygons_sf))) {
    polygons_sf <- sf::st_make_valid(polygons_sf)
    polygons_sf <- polygons_sf[sf::st_is_valid(polygons_sf), ]
  }
  if (sf::st_crs(points_sf) != sf::st_crs(polygons_sf)) {
    points_sf <- sf::st_transform(points_sf, sf::st_crs(polygons_sf))
  }
  
  clipped_sf <- sf::st_filter(points_sf, polygons_sf)
  if (sf::st_crs(clipped_sf)$input != "EPSG:4326") {
    clipped_sf <- sf::st_transform(clipped_sf, 4326)
  }
  
  cat("  After polygon clip:", nrow(clipped_sf), "points\n")
  
  if (nrow(clipped_sf) == 0) {
    warning("No occurrences inside range polygon - writing empty file")
    empty_df <- data.frame(species = character(0),
                           latitude = numeric(0),
                           longitude = numeric(0))
    write.csv(empty_df, output_file_path, row.names = FALSE)
    return(output_file_path)
  }
  
  coords <- sf::st_coordinates(clipped_sf)
  clipped_df <- data.frame(
    species   = species_name_from_file,
    latitude  = coords[, 2],
    longitude = coords[, 1],
    stringsAsFactors = FALSE
  )
  
  # ---- Resample to one point per raster cell ----
  clipped_vect <- terra::vect(clipped_df,
                              geom = c("longitude", "latitude"),
                              crs = "EPSG:4326")
  mask_raster <- terra::rast(mask_file_path)
  
  if (has_enmtools) {
    trimmed_points <- ENMTools::trimdupes.by.raster(points = clipped_vect,
                                                    mask = mask_raster)
  } else {
    cells <- terra::cells(mask_raster, clipped_vect)[, "cell"]
    keep_idx <- !duplicated(cells) & !is.na(cells)
    trimmed_points <- clipped_vect[keep_idx, ]
  }
  
  final_coords <- terra::crds(trimmed_points)
  final_df <- data.frame(
    species   = species_name_from_file,
    latitude  = final_coords[, 2],
    longitude = final_coords[, 1],
    stringsAsFactors = FALSE
  )
  
  cat("  After cell dedup:  ", nrow(final_df), "points (final)\n")
  
  write.csv(final_df, output_file_path, row.names = FALSE)
  cat("  Saved to:", output_file_path, "\n")
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
  
  # Keep only files that actually exist, remembering the intended name
  exists_mask <- file.exists(list_file)
  if (!all(exists_mask)) {
    missing <- top_vars[!exists_mask]
    cat("  Missing .asc for:", paste(missing, collapse = ", "), "\n")
  }
  list_file <- list_file[exists_mask]
  kept_names <- top_vars[exists_mask]
  
  if (length(list_file) == 0) stop("No matching .asc files found")
  
  ac <- raster::stack(list_file)
  # Force layer names to match the .asc filename stems exactly.
  # raster::stack otherwise strips shared prefixes / munges names,
  # which is what caused the downstream combo/_variables.txt mismatch.
  names(ac) <- kept_names
  
  ac.brick <- raster::brick(ac)
  names(ac.brick) <- kept_names
  ac.brick[ac.brick < -100] <- NA
  
  rc2.corr <- raster::layerStats(ac.brick, "pearson", na.rm = TRUE)
  
  # layerStats returns a list; the matrix is in $`pearson correlation coefficient`.
  # Force its row/col names too, in case raster still mangled them.
  corr_mat <- rc2.corr[[1]]
  rownames(corr_mat) <- kept_names
  colnames(corr_mat) <- kept_names
  rc2.corr[[1]] <- corr_mat
  
  save_file <- PathManager$getCorrelationPath()
  write.csv(rc2.corr, save_file)
  cat("Correlation matrix saved to:", save_file, "\n")
  return(save_file)
}

RFunctions$combos <- function(threshold = 0.8, csv_path = NULL, model_path = NULL) {
  if (is.null(csv_path)) csv_path <- PathManager$getCorrelationPath()
  if (is.null(model_path)) model_path <- PathManager$getModelsPath()
  
  cat("Generating variable combinations from:", csv_path, "\n")
  
  # check.names=FALSE preserves dots and special chars in variable names
  # (e.g. "wc2.1_bio_15_60m" won't become "wc2.1_bio_15_60m" -> "wc2_1_bio_15_60m")
  corr <- read.csv(csv_path, header = TRUE, check.names = FALSE)
  corr <- corr[, -c(1, ncol(corr))]
  colnames(corr) <- gsub("^pearson correlation coefficient\\.", "",
                         colnames(corr))
  colnames(corr) <- gsub("^pearson\\.correlation\\.coefficient\\.", "",
                         colnames(corr))
  
  cat("  Correlation matrix variables:\n")
  for (nm in colnames(corr)) cat("    ", nm, "\n")
  
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
        combo_vars <- colnames(corr)[combo]
        combo_name <- paste(combo_vars, collapse = "_")
        
        total_created <- total_created + 1
        combo_dir_name <- sprintf("combo_%05d", total_created)
        combo_dir <- file.path(model_path, combo_dir_name)
        if (!dir.exists(combo_dir)) dir.create(combo_dir, recursive = TRUE)
        
        writeLines(combo_vars, file.path(combo_dir, "_variables.txt"))
        
        combo_csv <- file.path(model_path, paste0("combinations_", size, ".csv"))
        write(paste(combo_dir_name, combo_name, sep = ","),
              file = combo_csv, append = TRUE)
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

# =============================================================================
# ENMTools-style model selection (auto-detects CSV format, computes AIC/AICc/BIC)
# =============================================================================

# ---- Read a points CSV and detect lat/lon columns -----------------------------
# ---- Read a points CSV and detect lat/lon columns -----------------------------
RFunctions$read_points_csv <- function(infile, verbose = TRUE) {
  
  # Read as raw bytes so we can detect and strip a UTF-8 BOM cleanly.
  # write.csv() on Windows sometimes writes one; if present, the first
  # column name becomes literally "\ufefflatitude" and every regex/name
  # match against "latitude" silently fails.
  fsize <- file.info(infile)$size
  con <- file(infile, "rb")
  raw_all <- readBin(con, "raw", n = fsize)
  close(con)
  
  if (length(raw_all) >= 3 &&
      raw_all[1] == as.raw(0xEF) &&
      raw_all[2] == as.raw(0xBB) &&
      raw_all[3] == as.raw(0xBF)) {
    raw_all <- raw_all[-(1:3)]
    if (verbose) message("  Stripped UTF-8 BOM from CSV")
  }
  
  # Normalize line endings and split
  text <- rawToChar(raw_all)
  text <- gsub("\r\n", "\n", text, fixed = TRUE)
  text <- gsub("\r",   "\n", text, fixed = TRUE)
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  
  if (length(lines) < 2) {
    stop(sprintf("CSV has fewer than 2 lines: %s", infile))
  }
  
  header <- trimws(lines[1])
  header_fields <- strsplit(header, ",", fixed = TRUE)[[1]]
  
  # Aggressive cleaning: whitespace, surrounding quotes, non-alphanumeric junk
  clean_name <- function(x) {
    x <- trimws(x)
    x <- gsub('^["\']+|["\']+$', "", x)
    x <- gsub("[^A-Za-z0-9_]", "", x)
    tolower(x)
  }
  header_fields_clean <- vapply(header_fields, clean_name, character(1))
  
  if (verbose) {
    message(sprintf("  CSV: %d columns [%s]",
                    length(header_fields),
                    paste(header_fields, collapse = ", ")))
    message(sprintf("  Cleaned:         [%s]",
                    paste(header_fields_clean, collapse = ", ")))
  }
  
  # ---- Name-based detection ----
  lat_idx <- NA_integer_
  lon_idx <- NA_integer_
  
  exact_lat <- which(header_fields_clean %in%
                       c("latitude", "lat", "decimallatitude", "y"))
  exact_lon <- which(header_fields_clean %in%
                       c("longitude", "lon", "long", "decimallongitude", "x"))
  if (length(exact_lat) > 0) lat_idx <- exact_lat[1]
  if (length(exact_lon) > 0) lon_idx <- exact_lon[1]
  
  if (is.na(lat_idx)) {
    p <- which(grepl("^lat", header_fields_clean))
    if (length(p) > 0) lat_idx <- p[1]
  }
  if (is.na(lon_idx)) {
    p <- which(grepl("^lon", header_fields_clean))
    if (length(p) > 0) lon_idx <- p[1]
  }
  
  # ---- Parse data rows ----
  raw_rows <- vector("list", length(lines) - 1)
  keep <- logical(length(lines) - 1)
  for (i in seq_along(raw_rows)) {
    ln <- trimws(lines[i + 1], which = "right")
    if (nchar(ln) == 0) { keep[i] <- FALSE; next }
    fields <- strsplit(ln, ",", fixed = TRUE)[[1]]
    fields <- trimws(fields)
    fields <- gsub('^["\']+|["\']+$', "", fields)
    raw_rows[[i]] <- fields
    keep[i] <- TRUE
  }
  raw_rows <- raw_rows[keep]
  
  if (length(raw_rows) == 0) {
    stop("CSV has no data rows: ", infile)
  }
  
  # ---- Value-plausibility fallback (used if name detection failed) ----
  #
  # Latitude is always in [-90, 90]; longitude is in [-180, 180]. A column
  # containing any value with |v|>90 CANNOT be latitude. This catches the
  # (species, lat, lon) vs (species, lon, lat) ambiguity that the old code
  # was hard-coding to the wrong choice.
  ncol_data <- length(raw_rows[[1]])
  n_check <- min(50, length(raw_rows))
  plausible_lat <- integer(ncol_data)
  plausible_lon <- integer(ncol_data)
  for (i in seq_len(ncol_data)) {
    for (j in seq_len(n_check)) {
      if (i > length(raw_rows[[j]])) next
      v <- suppressWarnings(as.numeric(raw_rows[[j]][i]))
      if (is.na(v)) next
      if (abs(v) <= 90)  plausible_lat[i] <- plausible_lat[i] + 1
      if (abs(v) <= 180) plausible_lon[i] <- plausible_lon[i] + 1
    }
  }
  
  if (is.na(lat_idx) || is.na(lon_idx)) {
    if (verbose) {
      message("  WARNING: lat/lon not detected by name -- probing values")
      for (i in seq_len(ncol_data)) {
        message(sprintf("    col %d: %d/%d in [-90,90], %d/%d in [-180,180]",
                        i, plausible_lat[i], n_check,
                        plausible_lon[i], n_check))
      }
    }
    
    # lon: highest plausible_lon, prefer cols that have OOB-for-lat values
    lon_cands <- order(-plausible_lon, plausible_lat)
    lat_cands <- order(-plausible_lat, -plausible_lon)
    
    if (is.na(lon_idx)) lon_idx <- lon_cands[1]
    if (is.na(lat_idx)) {
      for (c in lat_cands) if (c != lon_idx) { lat_idx <- c; break }
    }
    
    if (verbose) {
      message(sprintf("  Value-probe result: lat_idx=%d, lon_idx=%d",
                      lat_idx, lon_idx))
    }
  }
  
  # ---- Sanity check: if lat_idx has values outside [-90,90], SWAP ----
  if (plausible_lat[lat_idx] < n_check * 0.9 &&
      plausible_lat[lon_idx] >= n_check * 0.9) {
    if (verbose) {
      message(sprintf(
        "  !! Detected column swap: lat_idx=%d has %d/%d valid lats, ",
        lat_idx, plausible_lat[lat_idx], n_check),
        sprintf("but col %d has %d/%d. Swapping.",
                lon_idx, plausible_lat[lon_idx], n_check))
    }
    tmp <- lat_idx; lat_idx <- lon_idx; lon_idx <- tmp
  }
  
  has_species <- (length(header_fields) >= 3) &&
    !(1 %in% c(lat_idx, lon_idx))
  
  if (verbose) {
    lat_name <- if (lat_idx <= length(header_fields_clean))
      header_fields_clean[lat_idx] else "?"
    lon_name <- if (lon_idx <= length(header_fields_clean))
      header_fields_clean[lon_idx] else "?"
    message(sprintf("  Detected: lat_idx=%d (%s), lon_idx=%d (%s), has_species=%s",
                    lat_idx, lat_name, lon_idx, lon_name, has_species))
    message(sprintf("  Loaded %d data rows", length(raw_rows)))
  }
  
  list(
    header      = header,
    has_species = has_species,
    lat_idx     = lat_idx,
    lon_idx     = lon_idx,
    rows        = raw_rows
  )
}

# ---- Compute AIC/AICc/BIC for a single model ---------------------------------
RFunctions$compute_AIC <- function(datafile, csv_info, lambdasfile,
                                   verbose = TRUE) {
  loglikelihood <- 0
  
  # Count non-zero lambdas
  lambda_lines <- readLines(lambdasfile, warn = FALSE)
  nparams <- 0
  for (ln in lambda_lines) {
    parts <- strsplit(ln, ",", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    weight <- sub("\\s+", "", parts[2])
    if (!identical(weight, "0.0")) nparams <- nparams + 1
  }
  nparams <- nparams - 4  # Subtract MaxEnt metadata lines
  
  # Parse ASCII raster
  asc_lines <- readLines(datafile, warn = FALSE)
  fileparams <- list()
  data_lines <- character(0)
  
  for (ln in asc_lines) {
    if (grepl("^\\s*[0-9-]", ln)) {
      data_lines <- c(data_lines, trimws(ln, which = "right"))
    } else {
      parts <- strsplit(trimws(ln), "\\s+")[[1]]
      if (length(parts) >= 2) {
        fileparams[[tolower(parts[1])]] <- parts[2]
      }
    }
  }
  
  xll      <- as.numeric(fileparams[["xllcorner"]])
  yll      <- as.numeric(fileparams[["yllcorner"]])
  cellsize <- as.numeric(fileparams[["cellsize"]])
  ncols    <- as.integer(fileparams[["ncols"]])
  nrows    <- as.integer(fileparams[["nrows"]])
  
  if (verbose) {
    message(sprintf("  Raster geometry: xll=%g, yll=%g, cellsize=%g, ncols=%s, nrows=%s",
                    xll, yll, cellsize,
                    ifelse(is.na(ncols), "?", ncols),
                    ifelse(is.na(nrows), "?", nrows)))
    message(sprintf("  Data rows found: %d", length(data_lines)))
  }
  
  # Row 0 = bottom of map
  env_data <- rev(data_lines)
  
  # Pre-split for speed
  env_cells <- lapply(env_data, function(r) {
    v <- strsplit(r, "\\s+")[[1]]
    suppressWarnings(as.numeric(v[v != ""]))
  })
  
  probsum <- sum(sapply(env_cells, function(v) {
    ok <- !is.na(v) & v != -9999
    sum(v[ok])
  }))
  if (verbose) message(sprintf("  probsum = %g", probsum))
  
  npoints              <- 0
  points_processed     <- 0
  points_out_of_bounds <- 0
  points_no_data       <- 0
  
  for (fields in csv_info$rows) {
    points_processed <- points_processed + 1
    
    thisx <- suppressWarnings(as.numeric(fields[csv_info$lon_idx]))
    thisy <- suppressWarnings(as.numeric(fields[csv_info$lat_idx]))
    if (is.na(thisx) || is.na(thisy)) next
    
    row <- as.integer(floor((thisy - yll) / cellsize))
    col <- as.integer(floor((thisx - xll) / cellsize))
    
    if (row < 0 || row >= length(env_cells)) {
      points_out_of_bounds <- points_out_of_bounds + 1
      next
    }
    cells <- env_cells[[row + 1]]
    if (col < 0 || col >= length(cells)) {
      points_out_of_bounds <- points_out_of_bounds + 1
      next
    }
    
    layer_value <- cells[col + 1]
    
    if (!is.na(layer_value) && layer_value > 0) {
      loglikelihood <- loglikelihood + log(layer_value / probsum)
      npoints <- npoints + 1
    } else {
      points_no_data <- points_no_data + 1
    }
  }
  
  if (verbose) {
    message(sprintf("  === Summary: processed=%d, matched=%d, oob=%d, no-data/zero=%d ===",
                    points_processed, npoints,
                    points_out_of_bounds, points_no_data))
  }
  
  if (nparams >= npoints - 1) {
    AICscore <- NA_real_; AICcscore <- NA_real_; BICscore <- NA_real_
  } else {
    AICscore  <- 2 * nparams - 2 * loglikelihood
    AICcscore <- AICscore + (2 * nparams * (nparams + 1) / (npoints - nparams - 1))
    BICscore  <- nparams * log(npoints) - 2 * loglikelihood
  }
  
  list(
    loglikelihood = loglikelihood,
    nparams       = nparams,
    npoints       = npoints,
    AIC           = AICscore,
    AICc          = AICcscore,
    BIC           = BICscore
  )
}

# ---- Extract stats for one model and write to output connection --------------
RFunctions$modsel_extract_data <- function(ascfile, csv_info, lambdasfile,
                                           outcon, csvfile_label,
                                           verbose = TRUE) {
  if (verbose) {
    message(sprintf("\nExtracting data from %s using %s...", ascfile, csvfile_label))
  }
  res <- RFunctions$compute_AIC(ascfile, csv_info, lambdasfile, verbose = verbose)
  outline <- paste(
    csvfile_label, ascfile,
    res$loglikelihood, res$nparams, res$npoints,
    res$AIC, res$AICc, res$BIC,
    sep = ","
  )
  writeLines(outline, con = outcon)
  return(res)
}

# ---- Main entry point: process a modelSelection.csv control file -------------
# Control file format (no header, comma-separated):
#   <points_csv>,<ascii_file>,<lambdas_file>
# Writes: <control_file_stem>_model_selection.csv (with header)
RFunctions$modsel_execute <- function(modselfile, output_path = NULL,
                                      verbose = TRUE) {
  if (!file.exists(modselfile)) {
    stop(sprintf("Control file not found: %s", modselfile))
  }
  
  if (is.null(output_path)) {
    output_path <- sub("\\.csv$", "_model_selection.csv", modselfile,
                       ignore.case = TRUE)
    if (identical(output_path, modselfile)) {
      output_path <- paste0(modselfile, "_model_selection.csv")
    }
  }
  
  outcon <- file(output_path, open = "w")
  on.exit(close(outcon), add = TRUE)
  
  writeLines(
    "Points,ASCII file,Log Likelihood,Parameters,Sample Size,AIC score,AICc score,BIC score",
    con = outcon
  )
  
  control_lines <- readLines(modselfile, warn = FALSE)
  csv_cache <- list()
  processed <- 0
  skipped   <- 0
  
  for (raw in control_lines) {
    line <- trimws(raw, which = "right")
    line <- gsub("\"", "", line, fixed = TRUE)
    if (nchar(line) == 0) next
    
    fields <- strsplit(line, ",", fixed = TRUE)[[1]]
    if (length(fields) < 3) {
      message(sprintf("Skipping malformed line: %s", raw))
      skipped <- skipped + 1
      next
    }
    
    points_csv   <- fields[1]
    ascii_file   <- fields[2]
    lambdas_file <- fields[3]
    
    ready_to_go <- TRUE
    if (!file.exists(points_csv)) {
      message(sprintf("  MISSING points CSV:   %s", points_csv))
      ready_to_go <- FALSE
    }
    if (!file.exists(ascii_file)) {
      message(sprintf("  MISSING ASCII output: %s", ascii_file))
      ready_to_go <- FALSE
    }
    if (!file.exists(lambdas_file)) {
      message(sprintf("  MISSING lambdas:      %s", lambdas_file))
      ready_to_go <- FALSE
    }
    
    if (ready_to_go) {
      if (is.null(csv_cache[[points_csv]])) {
        csv_cache[[points_csv]] <- RFunctions$read_points_csv(points_csv,
                                                              verbose = verbose)
      }
      RFunctions$modsel_extract_data(
        ascii_file, csv_cache[[points_csv]], lambdas_file,
        outcon, csvfile_label = points_csv, verbose = verbose
      )
      processed <- processed + 1
    } else {
      skipped <- skipped + 1
    }
  }
  
  message(sprintf("\nModel selection done: %d processed, %d skipped", processed, skipped))
  message(sprintf("Output: %s", output_path))
  invisible(output_path)
}

# ---- Layer name resolution --------------------------------------------------
# Maps MaxEnt's shortened variable names (e.g. "BIO05") back to the real
# filenames on disk (e.g. "wc2.1_bio_05_60m.asc"). Handles WorldClim naming
# and falls back to raw filenames for custom layers.
RFunctions$buildLayerNameMap <- function(layer_dir) {
  if (!dir.exists(layer_dir)) stop("Layer directory does not exist: ", layer_dir)
  
  asc_files <- list.files(layer_dir, pattern = "\\.asc$", full.names = TRUE)
  if (length(asc_files) == 0) stop("No .asc files in: ", layer_dir)
  
  bases <- tools::file_path_sans_ext(basename(asc_files))
  name_map <- list()
  
  for (i in seq_along(asc_files)) {
    b    <- bases[i]
    path <- asc_files[i]
    
    # Extract the bioclim number if this is a WorldClim file
    match <- regmatches(
      b,
      regexec("^wc2\\.[0-9]+_bio_0*([0-9]+)_.*$", b, ignore.case = TRUE)
    )[[1]]
    
    if (length(match) == 2) {
      num_str <- match[2]                       # "5"  (no padding)
      num_int <- as.integer(num_str)            # 5
      padded  <- sprintf("BIO%02d", num_int)    # "BIO05"
      unpad   <- sprintf("BIO%d",  num_int)     # "BIO5"
      
      name_map[[padded]] <- path
      name_map[[unpad]]  <- path
      
      # Also case-insensitive variants
      name_map[[tolower(padded)]] <- path
      name_map[[tolower(unpad)]]  <- path
    }
    
    # Register the raw filename base too (for non-WorldClim custom layers)
    name_map[[b]]           <- path
    name_map[[toupper(b)]]  <- path
    name_map[[tolower(b)]]  <- path
  }
  
  return(name_map)
}


# Look up a single short/raw name -> full .asc path (NULL if not found)
RFunctions$resolveLayerFile <- function(name, name_map, layer_dir) {
  # Direct match
  hit <- name_map[[name]]
  if (!is.null(hit) && file.exists(hit)) return(hit)
  
  # Uppercase
  hit <- name_map[[toupper(name)]]
  if (!is.null(hit) && file.exists(hit)) return(hit)
  
  # Lowercase
  hit <- name_map[[tolower(name)]]
  if (!is.null(hit) && file.exists(hit)) return(hit)
  
  # Strip leading zeros from any BIO## and retry (BIO05 -> BIO5)
  if (grepl("^BIO0*[0-9]+$", name, ignore.case = TRUE)) {
    n <- as.integer(sub("^BIO0*", "", name, ignore.case = TRUE))
    for (candidate in c(sprintf("BIO%d", n), sprintf("BIO%02d", n))) {
      hit <- name_map[[candidate]]
      if (!is.null(hit) && file.exists(hit)) return(hit)
    }
  }
  
  # Last resort: try <name>.asc directly in the layer dir
  raw <- file.path(layer_dir, paste0(name, ".asc"))
  if (file.exists(raw)) return(raw)
  
  return(NULL)
}

# Recover the exact variable names that make up a combo. combos() writes a
# "_variables.txt" manifest into each combo directory containing exactly
# those names, one per line -- this is the safe way to read them back.
#
# Do NOT reconstruct this list by splitting the combo folder name on "_":
# the folder name is built by JOINING variable names with "_", but real
# variable names can themselves contain underscores (a SPECTRE layer like
# "spectre_paituli__hy_spectre_1_9_human_biomes_60m" has a dozen), so that
# join is not reversible by splitting. Doing so shreds multi-word variable
# names into meaningless single-word fragments that fail to resolve, and
# the whole combo gets skipped.
#
# Falls back to the legacy split-on-"_" behavior only for combo folders
# that predate this manifest (i.e. no "_variables.txt" present) -- re-run
# Step 9 (combos) to regenerate manifests for old combo folders.
RFunctions$getComboVariables <- function(combo_dir, folder_name) {
  manifest <- file.path(combo_dir, "_variables.txt")
  if (file.exists(manifest)) {
    vars <- readLines(manifest, warn = FALSE)
    vars <- vars[nzchar(trimws(vars))]
    if (length(vars) > 0) return(vars)
  }
  strsplit(folder_name, "_")[[1]]
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
  cat("Running ENMTools-style model selection...\n")
  cat("  Control file:", csv_path, "\n")
  cat("  Output:      ", output_path, "\n")
  
  RFunctions$modsel_execute(csv_path,
                            output_path = output_path,
                            verbose = TRUE)
  
  # Read back and reshape into the format Step 11 expects
  results <- read.csv(output_path, stringsAsFactors = FALSE)
  
  # Rename columns to match what findLowestScoreVariables looks for
  colnames(results) <- c("points_path", "ascii_file_path",
                         "loglikelihood", "parameter_count", "sample_size",
                         "aic_score", "aicc_score", "bic_score")
  
  # Coerce numeric columns (they may contain "x" from ENMTools convention)
  for (col in c("loglikelihood", "parameter_count", "sample_size",
                "aic_score", "aicc_score", "bic_score")) {
    results[[col]] <- suppressWarnings(as.numeric(results[[col]]))
  }
  
  # Overwrite the output file with the reshaped/normalized version
  write.csv(results, output_path, row.names = FALSE)
  
  cat("Model selection saved to:", output_path, "\n")
  return(results)
}