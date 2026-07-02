# =============================================================================
# SHPExtractor.R  --  Extract per-species range polygons from IUCN-style
#                    multi-part shapefiles.
#
# Expects a source directory layout like:
#   <MAIN_PATH>/SHP/IUCN/
#     AMPHIBIANS/{AMPHIBIANS_PART1.shp, AMPHIBIANS_PART2.shp, ...}
#     REPTILES/  {REPTILES_PART1.shp,   REPTILES_PART2.shp,   ...}
#     MAMMALS/   {MAMMALS_PART1.shp,    MAMMALS_PART2.shp,    ...}
#
# Ported from 01_download_gbif.R
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
})

SHPExtractor <- new.env()

# Map GBIF class group -> IUCN subfolder name
SHPExtractor$GROUP_TO_SHP_DIR <- list(
  Amphibian = "AMPHIBIANS",
  Reptile   = "REPTILES",
  Mammal    = "MAMMALS"
  # Bird / Fish / Plant: no source folder configured -> skipped
)

SHPExtractor$SPECIES_FIELD_CANDIDATES <- c(
  "binomial", "BINOMIAL",
  "sci_name", "SCI_NAME",
  "scientific", "SCIENTIFIC",
  "species",  "SPECIES",
  "name",     "NAME"
)

# ----------------------------------------------------------------------------
# Detect which attribute column holds the scientific name
# ----------------------------------------------------------------------------
SHPExtractor$detectSpeciesField <- function(shp_path) {
  dbf_path <- sub("\\.shp$", ".dbf", shp_path, ignore.case = TRUE)
  
  if (file.exists(dbf_path) && requireNamespace("foreign", quietly = TRUE)) {
    hdr <- tryCatch(foreign::read.dbf(dbf_path, as.is = TRUE),
                    error = function(e) NULL)
    if (!is.null(hdr)) {
      hits <- intersect(SHPExtractor$SPECIES_FIELD_CANDIDATES, names(hdr))
      if (length(hits) > 0) return(hits[1])
    }
  }
  
  layer <- tools::file_path_sans_ext(basename(shp_path))
  sql_variants <- c(
    sprintf('SELECT * FROM "%s" LIMIT 1', layer),
    sprintf("SELECT * FROM %s LIMIT 1", layer)
  )
  for (sql in sql_variants) {
    hdr <- tryCatch(
      sf::st_read(shp_path, query = sql, quiet = TRUE),
      error = function(e) NULL, warning = function(w) NULL
    )
    if (!is.null(hdr)) {
      hits <- intersect(SHPExtractor$SPECIES_FIELD_CANDIDATES, names(hdr))
      if (length(hits) > 0) return(hits[1])
    }
  }
  NULL
}

# ----------------------------------------------------------------------------
# Read matching features from ONE PART shapefile.
#
# Strategy: DBF presence check first (skip file entirely if species absent),
# then wkt_filter global bbox read to bypass ExecuteSQL/.atx quirks on Windows.
# ----------------------------------------------------------------------------
SHPExtractor$readMatchingFeatures <- function(shp_path, sci_name) {
  field <- SHPExtractor$detectSpeciesField(shp_path)
  if (is.null(field)) {
    cat("    ! No species field found in", basename(shp_path), "\n")
    return(NULL)
  }
  
  dbf_path <- sub("\\.shp$", ".dbf", shp_path, ignore.case = TRUE)
  
  # ---- Step 1: DBF presence check (fast reject) ----
  if (file.exists(dbf_path) && requireNamespace("foreign", quietly = TRUE)) {
    dbf <- tryCatch(foreign::read.dbf(dbf_path, as.is = TRUE),
                    error = function(e) NULL)
    if (!is.null(dbf) && field %in% names(dbf)) {
      present <- any(toupper(as.character(dbf[[field]])) == toupper(sci_name))
      if (!present) return(NULL)
    }
  }
  
  # ---- Step 2: wkt_filter global read + R-side filter ----
  global_wkt <- "POLYGON ((-180 -90, 180 -90, 180 90, -180 90, -180 -90))"
  full <- tryCatch(
    sf::st_read(shp_path, wkt_filter = global_wkt, quiet = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(full)) {
    cat("    ! wkt_filter read failed on", basename(shp_path),
        "-- falling back to full read.\n")
    full <- tryCatch(sf::st_read(shp_path, quiet = TRUE),
                     error = function(e) NULL)
    if (is.null(full)) return(NULL)
  }
  
  keep <- toupper(as.character(full[[field]])) == toupper(sci_name)
  keep[is.na(keep)] <- FALSE
  feats <- full[keep, , drop = FALSE]
  
  if (nrow(feats) == 0) return(NULL)
  
  feats$.source_file   <- basename(shp_path)
  feats$.species_field <- field
  feats
}

# ----------------------------------------------------------------------------
# Extract range polygons for a species and write to output SHP dir.
#
# Args:
#   species_tag - "Genus_species"
#   sci_name    - "Genus species" (scientific name to match in shapefile)
#   group       - "Reptile" / "Amphibian" / "Mammal" / etc.
#   shp_source_dir - root path containing GROUP subfolders
#   out_shp_dir - where to write "<species_tag>_range.shp"
#
# Returns list with $status ("copied"|"not_found"|"source_missing"|"skipped")
#                   $out_shp (path if written, else NA)
# ----------------------------------------------------------------------------
SHPExtractor$extractSpeciesRange <- function(species_tag,
                                             sci_name,
                                             group,
                                             shp_source_dir,
                                             out_shp_dir) {
  
  if (is.null(group) || is.na(group) ||
      is.null(SHPExtractor$GROUP_TO_SHP_DIR[[group]])) {
    cat("  -> SHP: no source folder configured for group '",
        group %||% "NA", "'. Skipping.\n", sep = "")
    return(list(status = "skipped", out_shp = NA_character_))
  }
  
  group_dir <- file.path(shp_source_dir, SHPExtractor$GROUP_TO_SHP_DIR[[group]])
  if (!dir.exists(group_dir)) {
    cat("  ! SHP source folder missing:", group_dir, "\n")
    return(list(status = "source_missing", out_shp = NA_character_))
  }
  
  parts <- list.files(group_dir, pattern = "\\.shp$",
                      ignore.case = TRUE, full.names = TRUE)
  if (length(parts) == 0) {
    cat("  ! No .shp files inside", group_dir, "\n")
    return(list(status = "source_missing", out_shp = NA_character_))
  }
  
  cat(sprintf("  -> Searching %d range file(s) for '%s' ...\n",
              length(parts), sci_name))
  
  matches <- list()
  for (p in parts) {
    feats <- SHPExtractor$readMatchingFeatures(p, sci_name)
    if (!is.null(feats) && nrow(feats) > 0) {
      cat(sprintf("     match: %d feature(s) in %s\n",
                  nrow(feats), basename(p)))
      matches[[length(matches) + 1]] <- feats
    }
  }
  
  if (length(matches) == 0) {
    cat("  ! No range polygons found for '", sci_name, "' in ",
        SHPExtractor$GROUP_TO_SHP_DIR[[group]], "/.\n", sep = "")
    return(list(status = "not_found", out_shp = NA_character_))
  }
  
  combined <- tryCatch(
    do.call(rbind, matches),
    error = function(e) {
      cat("    ! rbind failed (", conditionMessage(e),
          ") -- reducing to common columns.\n", sep = "")
      common <- Reduce(intersect, lapply(matches, names))
      do.call(rbind, lapply(matches, function(x) x[, common, drop = FALSE]))
    }
  )
  
  if (!dir.exists(out_shp_dir)) {
    dir.create(out_shp_dir, recursive = TRUE, showWarnings = FALSE)
  }
  out_shp <- file.path(out_shp_dir, paste0(species_tag, "_range.shp"))
  sf::st_write(combined, out_shp, delete_dsn = TRUE, quiet = TRUE)
  cat("  OK Wrote", nrow(combined), "polygon(s) ->", out_shp, "\n")
  
  list(status = "copied", out_shp = out_shp, feature_count = nrow(combined))
}

# ----------------------------------------------------------------------------
# Get bounding-box extent from a range shapefile.
# Returns c(xmin, xmax, ymin, ymax) with a small buffer, or NULL on failure.
# ----------------------------------------------------------------------------
SHPExtractor$getExtentFromShapefile <- function(shp_path, buffer_deg = 1.0) {
  if (!file.exists(shp_path)) return(NULL)
  
  shp <- tryCatch(sf::st_read(shp_path, quiet = TRUE),
                  error = function(e) NULL)
  if (is.null(shp) || nrow(shp) == 0) return(NULL)
  
  bbox <- sf::st_bbox(shp)
  c(
    as.numeric(bbox["xmin"]) - buffer_deg,
    as.numeric(bbox["xmax"]) + buffer_deg,
    as.numeric(bbox["ymin"]) - buffer_deg,
    as.numeric(bbox["ymax"]) + buffer_deg
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a