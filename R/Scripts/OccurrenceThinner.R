# =============================================================================
# OccurrenceThinner.R - Spatial thinning + range clipping (from 02_thin_occurrences.R)
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
})

OccurrenceThinner <- new.env()

# ── Find a reference raster to define the master grid ───────────────────────
OccurrenceThinner$findReferenceRaster <- function() {
  if (dir.exists(TIF_DIR)) {
    cands <- list.files(TIF_DIR, pattern = "\\.tif$", full.names = TRUE)
    if (length(cands) > 0) return(cands[1])
  }
  if (dir.exists(ASC_DIR)) {
    cands <- list.files(ASC_DIR, pattern = "\\.asc$", full.names = TRUE)
    if (length(cands) > 0) return(cands[1])
  }
  return(NULL)
}

OccurrenceThinner$loadReferenceTemplate <- function() {
  ref_path <- OccurrenceThinner$findReferenceRaster()
  if (!is.null(ref_path)) {
    cat("  Reference raster:", ref_path, "\n")
    r <- terra::rast(ref_path)
    if (terra::nlyr(r) > 1) r <- r[[1]]
    return(r)
  }
  cat("  No reference raster - synthesising from TARGET_RES_ARCMIN\n")
  res_deg <- TARGET_RES_ARCMIN / 60
  terra::rast(xmin = -180, xmax = 180, ymin = -90, ymax = 90,
              resolution = res_deg, crs = "EPSG:4326")
}

# ── Thin: one record per grid cell ──────────────────────────────────────────
OccurrenceThinner$thinToGrid <- function(df, template) {
  # Support both "decimalLongitude/Latitude" and "longitude/latitude" column names
  lon_col <- if ("decimalLongitude" %in% names(df)) "decimalLongitude" else "longitude"
  lat_col <- if ("decimalLatitude" %in% names(df)) "decimalLatitude" else "latitude"
  
  cells <- terra::cellFromXY(template, cbind(df[[lon_col]], df[[lat_col]]))
  out_of_extent <- is.na(cells)
  
  if (any(out_of_extent)) {
    cat("    Dropping", sum(out_of_extent), "records outside grid extent\n")
    df <- df[!out_of_extent, , drop = FALSE]
    cells <- cells[!out_of_extent]
  }
  
  keep_idx <- !duplicated(cells)
  thinned <- df[keep_idx, , drop = FALSE]
  
  list(df = thinned, n_in = nrow(df), n_out = nrow(thinned),
       n_dropped = nrow(df) - nrow(thinned),
       cell_res_arcmin = terra::res(template)[1] * 60)
}

# ── Clip to range polygon ───────────────────────────────────────────────────
OccurrenceThinner$clipToRange <- function(df, shp_path) {
  if (!file.exists(shp_path)) return(NULL)
  
  poly <- tryCatch(sf::st_read(shp_path, quiet = TRUE), error = function(e) NULL)
  if (is.null(poly) || nrow(poly) == 0) return(NULL)
  
  if (is.na(sf::st_crs(poly))) sf::st_crs(poly) <- 4326
  poly <- sf::st_transform(poly, 4326)
  poly <- tryCatch(sf::st_make_valid(poly), error = function(e) poly)
  
  lon_col <- if ("decimalLongitude" %in% names(df)) "decimalLongitude" else "longitude"
  lat_col <- if ("decimalLatitude" %in% names(df)) "decimalLatitude" else "latitude"
  
  pts <- sf::st_as_sf(df, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE)
  hits <- lengths(sf::st_intersects(pts, poly)) > 0
  kept <- df[hits, , drop = FALSE]
  
  list(df = kept, n_in = nrow(df), n_out = nrow(kept),
       n_dropped = nrow(df) - nrow(kept))
}

# ── Main thinning workflow ──────────────────────────────────────────────────
OccurrenceThinner$processSpecies <- function(species_tag) {
  cat("\n=== Thinning:", species_tag, "===\n")
  
  in_csv <- PathManager$speciesCsv
  out_csv <- PathManager$trimmedCsv
  clipped_csv <- PathManager$clippedCsv
  range_shp <- PathManager$getRangeShp()
  
  if (!file.exists(in_csv)) {
    cat("  Input CSV not found:", in_csv, "\n")
    return(NULL)
  }
  
  df <- read.csv(in_csv, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0) {
    write.csv(df, out_csv, row.names = FALSE)
    return(list(n_in = 0, n_thinned = 0))
  }
  
  cat("  Input records:", nrow(df), "\n")
  
  template <- OccurrenceThinner$loadReferenceTemplate()
  
  # Thin
  th <- OccurrenceThinner$thinToGrid(df, template)
  cat("  Thinned:", th$n_in, "->", th$n_out, "(dropped", th$n_dropped, ")\n")
  write.csv(th$df, out_csv, row.names = FALSE)
  cat("  Saved:", out_csv, "\n")
  
  # Clip to range
  n_clipped <- NA
  if (file.exists(range_shp) && EXTRACT_IUCN_RANGE) {
    cl <- OccurrenceThinner$clipToRange(th$df, range_shp)
    if (!is.null(cl)) {
      cat("  Range-clipped:", cl$n_in, "->", cl$n_out, "\n")
      write.csv(cl$df, clipped_csv, row.names = FALSE)
      n_clipped <- cl$n_out
    }
  }
  
  list(n_in = th$n_in, n_thinned = th$n_out, n_clipped = n_clipped)
}