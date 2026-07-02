# =============================================================================
# DataDownloader.R - GBIF occurrence + environmental layer downloads
# Wraps the standalone 01_download_gbif.R into the system architecture
# =============================================================================

suppressPackageStartupMessages({
  library(rgbif)
  library(dplyr)
  library(terra)
  library(geodata)
  library(sf)
})

DataDownloader <- new.env()

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

# ── Group classification (from 01_download_gbif.R) ──────────────────────────
DataDownloader$classifyGroup <- function(gbif_class = NA, gbif_order = NA, 
                                         gbif_kingdom = NA, gbif_phylum = NA) {
  cls <- if (is.null(gbif_class) || is.na(gbif_class)) "" else tolower(gbif_class)
  ord <- if (is.null(gbif_order) || is.na(gbif_order)) "" else tolower(gbif_order)
  kng <- if (is.null(gbif_kingdom) || is.na(gbif_kingdom)) "" else tolower(gbif_kingdom)
  phy <- if (is.null(gbif_phylum) || is.na(gbif_phylum)) "" else tolower(gbif_phylum)
  
  if (cls == "reptilia") return("Reptile")
  if (cls == "amphibia") return("Amphibian")
  if (cls == "mammalia") return("Mammal")
  if (cls == "aves") return("Bird")
  if (cls %in% c("actinopterygii", "chondrichthyes", "sarcopterygii",
                 "myxini", "petromyzonti", "agnatha")) return("Fish")
  
  reptile_orders <- c("squamata", "testudines", "crocodylia", "rhynchocephalia")
  amphibian_orders <- c("anura", "caudata", "urodela", "gymnophiona")
  mammal_orders <- c("rodentia", "chiroptera", "carnivora", "primates",
                     "artiodactyla", "cetartiodactyla", "perissodactyla",
                     "lagomorpha", "eulipotyphla", "cetacea")
  
  if (cls %in% reptile_orders   || ord %in% reptile_orders)   return("Reptile")
  if (cls %in% amphibian_orders || ord %in% amphibian_orders) return("Amphibian")
  if (cls %in% mammal_orders    || ord %in% mammal_orders)    return("Mammal")
  
  if (kng == "plantae") return("Plant")
  NA_character_
}

# ── GBIF taxon resolution ───────────────────────────────────────────────────
DataDownloader$resolveTaxonKey <- function(species_underscore, min_confidence = 90) {
  sci_name <- gsub("_", " ", species_underscore)
  m <- tryCatch(
    rgbif::name_backbone(name = sci_name, rank = "species", strict = FALSE),
    error = function(e) NULL
  )
  if (is.null(m) || is.null(m$usageKey)) return(NULL)
  
  status <- m$matchType %||% ""
  conf <- as.integer(m$confidence %||% 0)
  rank <- tolower(m$rank %||% "")
  if (rank != "species" || status %in% c("NONE", "HIGHERRANK") || 
      conf < min_confidence) return(NULL)
  
  list(
    key = m$usageKey, sci_name = sci_name,
    canonical_name = m$canonicalName %||% NA_character_,
    matched_name = m$scientificName %||% NA_character_,
    confidence = conf, class = m$class %||% NA_character_,
    order = m$order %||% NA_character_, kingdom = m$kingdom %||% NA_character_,
    phylum = m$phylum %||% NA_character_, family = m$family %||% NA_character_,
    genus = m$genus %||% NA_character_,
    group = DataDownloader$classifyGroup(m$class, m$order, m$kingdom, m$phylum)
  )
}

# ── Download GBIF occurrences ───────────────────────────────────────────────
DataDownloader$downloadGBIF <- function(species_tag, max_records = NULL, 
                                        page_size = NULL) {
  if (is.null(max_records)) max_records <- GBIF_MAX_RECORDS
  if (is.null(page_size)) page_size <- GBIF_PAGE_SIZE
  
  cat("Resolving:", species_tag, "\n")
  tk <- DataDownloader$resolveTaxonKey(species_tag, GBIF_MIN_CONFIDENCE)
  if (is.null(tk)) {
    cat("  Species not resolved against GBIF backbone\n")
    return(NULL)
  }
  
  cat("  Matched:", tk$matched_name, "| group:", tk$group, "\n")
  cat("  Downloading occurrences...\n")
  
  all_rows <- list(); start <- 0; fetched <- 0
  repeat {
    res <- rgbif::occ_search(taxonKey = tk$key, hasCoordinate = TRUE,
                             hasGeospatialIssue = FALSE,
                             limit = page_size, start = start)
    if (is.null(res$data) || nrow(res$data) == 0) break
    all_rows[[length(all_rows) + 1]] <- res$data
    fetched <- fetched + nrow(res$data)
    cat("    fetched", fetched, "records\n")
    if (nrow(res$data) < page_size || fetched >= max_records) break
    start <- start + page_size
  }
  
  df <- if (length(all_rows) == 0) tibble::tibble() else bind_rows(all_rows)
  
  list(data = df, taxon = tk, n_downloaded = fetched)
}

# ── Tidy and save occurrences ───────────────────────────────────────────────
DataDownloader$tidyAndSave <- function(df, species_tag, tk, csv_path) {
  if (nrow(df) == 0) {
    write.csv(data.frame(species_tag = character(),
                         decimalLongitude = numeric(),
                         decimalLatitude = numeric()),
              csv_path, row.names = FALSE)
    return(0)
  }
  
  keep_cols <- intersect(
    c("species", "scientificName", "decimalLongitude", "decimalLatitude",
      "countryCode", "year", "basisOfRecord", "kingdom", "phylum", "class",
      "order", "family", "genus", "gbifID"),
    names(df))
  
  out <- df[, keep_cols, drop = FALSE] %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
    filter(decimalLongitude >= -180, decimalLongitude <= 180,
           decimalLatitude >= -90, decimalLatitude <= 90) %>%
    distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)
  
  out$species_tag <- species_tag
  out$class_group <- tk$group %||% NA_character_
  
  write.csv(out, csv_path, row.names = FALSE)
  cat("  Saved", nrow(out), "records to:", csv_path, "\n")
  return(nrow(out))
}

# ── WorldClim download ──────────────────────────────────────────────────────
DataDownloader$downloadWorldClim <- function(target_res_arcmin = 5) {
  if (!DOWNLOAD_WORLDCLIM) {
    cat("WorldClim download disabled in config\n")
    return(NULL)
  }
  
  cat("\n=== Downloading WorldClim ===\n")
  res_tag <- sprintf("%dm", round(target_res_arcmin * 10))
  
  # Download raw to EnvironmentalLayers/WorldClim/
  bioclim_raw <- tryCatch(
    geodata::worldclim_global(var = "bio", res = 5, path = WORLDCLIM_RAW_DIR),
    error = function(e) {
      cat("  Failed:", e$message, "\n"); NULL
    }
  )
  if (is.null(bioclim_raw)) return(NULL)
  
  # Resample and write to TIF/ at master resolution
  cat("  Resampling to", target_res_arcmin, "arc-min and writing to TIF/\n")
  template <- DataDownloader$makeTargetTemplate(target_res_arcmin)
  
  out_paths <- character(0)
  for (i in seq_len(terra::nlyr(bioclim_raw))) {
    out_path <- file.path(TIF_DIR, sprintf("wc2.1_bio_%02d_%s.tif", i, res_tag))
    if (file.exists(out_path)) {
      out_paths <- c(out_paths, out_path); next
    }
    lyr <- bioclim_raw[[i]]
    names(lyr) <- sprintf("BIO%02d", i)
    lyr_res <- terra::resample(lyr, template, method = "bilinear")
    terra::writeRaster(lyr_res, out_path, overwrite = TRUE)
    cat("    BIO", sprintf("%02d", i), "->", basename(out_path), "\n")
    out_paths <- c(out_paths, out_path)
  }
  
  return(out_paths)
}

# ── Landcover & Footprint ───────────────────────────────────────────────────
DataDownloader$downloadLandcover <- function(target_res_arcmin = 5) {
  if (!DOWNLOAD_LANDCOVER) return(NULL)
  
  cat("\n=== Downloading Landcover ===\n")
  res_tag <- sprintf("%dm", round(target_res_arcmin * 10))
  template <- DataDownloader$makeTargetTemplate(target_res_arcmin)
  
  vars <- c("trees", "built")
  out_paths <- character(0)
  
  for (var in vars) {
    cat("  Downloading", var, "...\n")
    src <- tryCatch(geodata::landcover(var = var, path = LANDCOVER_RAW_DIR),
                    error = function(e) NULL)
    if (is.null(src)) next
    
    out_path <- file.path(TIF_DIR, sprintf("%s_%s.tif", var, res_tag))
    r <- terra::resample(src, template, method = "bilinear")
    names(r) <- var
    terra::writeRaster(r, out_path, overwrite = TRUE)
    out_paths <- c(out_paths, out_path)
    cat("    ->", basename(out_path), "\n")
  }
  
  return(out_paths)
}

DataDownloader$downloadFootprint <- function(target_res_arcmin = 5, year = 2009) {
  if (!DOWNLOAD_FOOTPRINT) return(NULL)
  
  cat("\n=== Downloading Human Footprint ===\n")
  res_tag <- sprintf("%dm", round(target_res_arcmin * 10))
  template <- DataDownloader$makeTargetTemplate(target_res_arcmin)
  
  src <- tryCatch(geodata::footprint(year = year, path = FOOTPRINT_RAW_DIR),
                  error = function(e) NULL)
  if (is.null(src)) return(NULL)
  
  out_path <- file.path(TIF_DIR, sprintf("footprint_%d_%s.tif", year, res_tag))
  r <- terra::resample(src, template, method = "bilinear")
  names(r) <- sprintf("footprint_%d", year)
  terra::writeRaster(r, out_path, overwrite = TRUE)
  cat("  ->", basename(out_path), "\n")
  return(out_path)
}

# ── SPECTRE download via WCS ────────────────────────────────────────────────
DataDownloader$downloadSPECTRE <- function(target_res_arcmin = 5) {
  if (!DOWNLOAD_SPECTRE) return(NULL)
  
  cat("\n=== Downloading SPECTRE (WCS) ===\n")
  
  if (!requireNamespace("gecko", quietly = TRUE)) {
    cat("  gecko package not installed - skipping SPECTRE\n")
    return(NULL)
  }
  
  res_tag <- sprintf("%dm", round(target_res_arcmin * 10))
  cache_path <- file.path(TIF_DIR, paste0("spectre_", res_tag, ".tif"))
  if (file.exists(cache_path)) {
    cat("  Cache hit:", cache_path, "\n")
    return(cache_path)
  }
  
  wcs_base <- "https://paituli.csc.fi/geoserver/wcs"
  res_deg <- target_res_arcmin / 60
  width <- round(360 / res_deg)
  height <- round(150 / res_deg)
  
  meta_fn <- tryCatch(getAnywhere("spectre.metadata")$objs[[1]],
                      error = function(e) NULL)
  if (is.null(meta_fn)) return(NULL)
  
  meta <- meta_fn()
  per_layer_dir <- file.path(SPECTRE_RAW_DIR, "wcs_layers")
  if (!dir.exists(per_layer_dir)) dir.create(per_layer_dir, recursive = TRUE)
  
  options(timeout = 600)
  layer_paths <- character(0)
  
  for (i in 1:min(21, nrow(meta))) {
    layer_name <- as.character(meta[i, 2])
    coverage_id <- as.character(meta[i, 3])
    out_file <- file.path(per_layer_dir, paste0(coverage_id, "_", res_tag, ".tif"))
    
    if (file.exists(out_file) && file.size(out_file) > 1024) {
      layer_paths <- c(layer_paths, out_file); next
    }
    
    url <- sprintf(paste0("%s?version=2.0.1&request=GetCoverage&coverageId=%s",
                          "&subset=Long(-180,180)&subset=Lat(-60,90)",
                          "&format=image/tiff&SCALESIZE=i(%d),j(%d)"),
                   wcs_base, coverage_id, width, height)
    
    ok <- tryCatch({
      utils::download.file(url, out_file, mode = "wb", quiet = TRUE); TRUE
    }, error = function(e) FALSE)
    
    if (ok && file.size(out_file) > 1024) {
      layer_paths <- c(layer_paths, out_file)
      cat("  Got:", layer_name, "\n")
    }
  }
  
  if (length(layer_paths) == 0) return(NULL)
  
  spectre <- terra::rast(layer_paths)
  template <- DataDownloader$makeTargetTemplate(target_res_arcmin)
  if (terra::crs(spectre) == "" || is.na(terra::crs(spectre))) {
    terra::crs(spectre) <- "EPSG:4326"
  }
  spectre <- terra::resample(spectre, template, method = "bilinear")
  terra::writeRaster(spectre, cache_path, overwrite = TRUE)
  cat("  ->", cache_path, "\n")
  return(cache_path)
}

# ── Process custom user layers ──────────────────────────────────────────────
DataDownloader$processCustomLayers <- function(target_res_arcmin = 5) {
  if (!USE_CUSTOM_LAYERS) return(NULL)
  
  cat("\n=== Processing Custom Layers ===\n")
  res_tag <- sprintf("%dm", round(target_res_arcmin * 10))
  template <- DataDownloader$makeTargetTemplate(target_res_arcmin)
  
  files <- list.files(CUSTOM_RAW_DIR, pattern = "\\.(tif|tiff|asc)$",
                      ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
  
  if (length(files) == 0) {
    cat("  No custom rasters found in:", CUSTOM_RAW_DIR, "\n")
    return(NULL)
  }
  
  out_paths <- character(0)
  for (f in files) {
    bn <- tools::file_path_sans_ext(basename(f))
    out_path <- file.path(TIF_DIR, paste0(bn, "_", res_tag, ".tif"))
    if (file.exists(out_path)) {
      out_paths <- c(out_paths, out_path); next
    }
    
    src <- tryCatch(terra::rast(f), error = function(e) NULL)
    if (is.null(src)) next
    
    method <- if (grepl("_categorical", tolower(bn))) "near" else "bilinear"
    r <- terra::resample(src, template, method = method)
    terra::writeRaster(r, out_path, overwrite = TRUE)
    cat("  ->", basename(out_path), "(", method, ")\n")
    out_paths <- c(out_paths, out_path)
  }
  
  return(out_paths)
}

# ── Master template builder ─────────────────────────────────────────────────
DataDownloader$makeTargetTemplate <- function(target_res_arcmin) {
  if (!is.null(USER_TEMPLATE_RASTER) && file.exists(USER_TEMPLATE_RASTER)) {
    return(terra::rast(USER_TEMPLATE_RASTER))
  }
  res_deg <- target_res_arcmin / 60
  terra::rast(xmin = -180, xmax = 180, ymin = -90, ymax = 90,
              resolution = res_deg, crs = "EPSG:4326")
}

# ── TIF -> ASC conversion ───────────────────────────────────────────────────
DataDownloader$convertTifToAsc <- function() {
  cat("\n=== Converting TIF -> ASC ===\n")
  tif_files <- list.files(TIF_DIR, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(tif_files) == 0) {
    cat("  No TIF files found in:", TIF_DIR, "\n")
    return(NULL)
  }
  
  for (tif in tif_files) {
    bn <- tools::file_path_sans_ext(basename(tif))
    asc_path <- file.path(ASC_DIR, paste0(bn, ".asc"))
    if (file.exists(asc_path)) next
    
    r <- terra::rast(tif)
    terra::NAflag(r) <- -9999
    terra::writeRaster(r, asc_path, filetype = "AAIGrid",
                       overwrite = TRUE, NAflag = -9999)
    cat("  ->", basename(asc_path), "\n")
  }
  
  cat("  ASC files in:", ASC_DIR, "\n")
}

# ── IUCN range extraction ───────────────────────────────────────────────────
DataDownloader$extractSpeciesRange <- function(species_tag, tk) {
  if (!EXTRACT_IUCN_RANGE) return("disabled")
  
  group_to_dir <- list(Amphibian = "AMPHIBIANS", Reptile = "REPTILES",
                       Mammal = "MAMMALS")
  
  group <- tk$group %||% NA_character_
  if (is.na(group) || is.null(group_to_dir[[group]])) {
    cat("  No IUCN source configured for group:", group, "\n")
    return("no_source")
  }
  
  group_dir <- file.path(SHP_SOURCE_DIR, group_to_dir[[group]])
  if (!dir.exists(group_dir)) return("source_missing")
  
  parts <- list.files(group_dir, pattern = "\\.shp$", 
                      ignore.case = TRUE, full.names = TRUE)
  if (length(parts) == 0) return("source_missing")
  
  sci_name <- tk$canonical_name %||% tk$sci_name
  cat("  Searching", length(parts), "shapefiles for:", sci_name, "\n")
  
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(TRUE))
  
  matches <- list()
  for (p in parts) {
    full <- tryCatch(sf::st_read(p, quiet = TRUE), error = function(e) NULL)
    if (is.null(full)) next
    
    field <- intersect(c("binomial", "BINOMIAL", "sci_name", "SCI_NAME"),
                       names(full))[1]
    if (is.na(field)) next
    
    keep <- toupper(as.character(full[[field]])) == toupper(sci_name)
    keep[is.na(keep)] <- FALSE
    if (any(keep)) matches[[length(matches) + 1]] <- full[keep, , drop = FALSE]
  }
  
  if (length(matches) == 0) return("not_found")
  
  combined <- do.call(rbind, matches)
  out_shp <- PathManager$getRangeShp()
  sf::st_write(combined, out_shp, delete_dsn = TRUE, quiet = TRUE)
  cat("  Wrote", nrow(combined), "polygon(s) to:", out_shp, "\n")
  
  return("copied")
}