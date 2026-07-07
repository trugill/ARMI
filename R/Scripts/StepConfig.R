# =============================================================================
# StepConfig.R - Per-step configuration constructors
# =============================================================================
StepConfig <- new.env()

# --- Step 0: Data download ---------------------------------------------------
StepConfig$Step0Config <- function(enabled = TRUE,
                                   download_worldclim = TRUE,
                                   download_spectre = FALSE,
                                   download_landcover = FALSE,
                                   download_footprint = FALSE,
                                   use_custom_layers = FALSE,
                                   download_gbif = TRUE,
                                   extract_iucn_range = FALSE,
                                   target_res_arcmin = 5) {
  list(
    type = "Step0Config",
    enabled = enabled,
    download_worldclim = download_worldclim,
    download_spectre = download_spectre,
    download_landcover = download_landcover,
    download_footprint = download_footprint,
    use_custom_layers = use_custom_layers,
    download_gbif = download_gbif,
    extract_iucn_range = extract_iucn_range,
    target_res_arcmin = target_res_arcmin
  )
}

# --- Step 1: GBIF / occurrence loading ---------------------------------------
StepConfig$Step1Config <- function(enabled = TRUE,
                                   occurrence_file = "",
                                   use_raw_gbif = TRUE,
                                   use_iucn_shp = FALSE,
                                   shp_path = "") {
  list(
    type = "Step1Config",
    enabled = enabled,
    occurrence_file = occurrence_file,
    use_raw_gbif = use_raw_gbif,
    use_iucn_shp = use_iucn_shp,
    shp_path = shp_path
  )
}

# --- Step 3: Trim duplicates -------------------------------------------------
StepConfig$Step3Config <- function(enabled = TRUE) {
  list(type = "Step3Config", enabled = enabled)
}

# --- Step 4: Clip rasters ----------------------------------------------------
StepConfig$Step4Config <- function(enabled = TRUE,
                                   clip_to_extent = TRUE,
                                   max_lat = NA, max_lon = NA,
                                   min_lat = NA, min_lon = NA,
                                   use_mask = FALSE,
                                   mask_max_lat = NA, mask_max_lon = NA,
                                   mask_min_lat = NA, mask_min_lon = NA,
                                   reuse_existing_layers = FALSE,
                                   reuse_existing_mask = FALSE) {
  list(
    type = "Step4Config",
    enabled = enabled,
    clip_to_extent = clip_to_extent,
    max_lat = max_lat, max_lon = max_lon,
    min_lat = min_lat, min_lon = min_lon,
    use_mask = use_mask,
    mask_max_lat = mask_max_lat, mask_max_lon = mask_max_lon,
    mask_min_lat = mask_min_lat, mask_min_lon = mask_min_lon,
    reuse_existing_layers = reuse_existing_layers,
    reuse_existing_mask = reuse_existing_mask
  )
}

# --- Step 6: Global Maxent model ---------------------------------------------
StepConfig$Step6Config <- function(enabled = TRUE,
                                   extra_args = character(0),
                                   skip_if_exists = FALSE) {
  list(
    type = "Step6Config",
    enabled = enabled,
    extra_args = extra_args,
    skip_if_exists = skip_if_exists
  )
}

# --- Step 7: Top variables ---------------------------------------------------
StepConfig$Step7Config <- function(enabled = TRUE,
                                   threshold = 0.05,
                                   top_count = 10) {
  list(
    type = "Step7Config",
    enabled = enabled,
    threshold = threshold,
    top_count = top_count
  )
}

# --- Step 8: Correlation -----------------------------------------------------
StepConfig$Step8Config <- function(enabled = TRUE) {
  list(type = "Step8Config", enabled = enabled)
}

# --- Step 9: Variable combinations -------------------------------------------
StepConfig$Step9Config <- function(enabled = TRUE, correlation_threshold = 0.8) {
  list(
    type = "Step9Config",
    enabled = enabled,
    correlation_threshold = correlation_threshold
  )
}

# --- Step 10: Run all permutations -------------------------------------------
StepConfig$Step10Config <- function(enabled = TRUE, required_tifs = character(0)) {
  list(
    type = "Step10Config",
    enabled = enabled,
    required_tifs = required_tifs
  )
}

# --- Step 11: Model selection ------------------------------------------------
StepConfig$Step11Config <- function(enabled = TRUE, selection_criterion = "AICc") {
  list(
    type = "Step11Config",
    enabled = enabled,
    selection_criterion = selection_criterion
  )
}

# --- Step 12: Beta optimization ----------------------------------------------
StepConfig$Step12Config <- function(enabled = TRUE,
                                    max_beta = 5.0,
                                    beta_increment = 0.5,
                                    selection_criterion = "AICc") {
  list(
    type = "Step12Config",
    enabled = enabled,
    max_beta = max_beta,
    beta_increment = beta_increment,
    selection_criterion = selection_criterion
  )
}

# --- Step 13: Final model ----------------------------------------------------
StepConfig$Step13Config <- function(enabled = TRUE,
                                    replicates = 10,
                                    projection_layers = character(0)) {
  list(
    type = "Step13Config",
    enabled = enabled,
    replicates = replicates,
    projection_layers = projection_layers
  )
}

# --- Validator (tolerant) ----------------------------------------------------
StepConfig$validate <- function(config) {
  # If not a proper step config, skip silently
  if (is.null(config) || !is.list(config)) return(TRUE)
  
  type <- config$type
  if (is.null(type) || length(type) != 1 || is.na(type) || !nzchar(type)) {
    return(TRUE)  # unknown/missing type: nothing to validate
  }
  
  switch(type,
         Step7Config = {
           if (config$threshold < 0 || config$threshold > 1)
             stop("Step 7: threshold must be 0-1")
           if (config$top_count < 1)
             stop("Step 7: top_count must be >= 1")
         },
         Step9Config = {
           if (config$correlation_threshold < 0 || config$correlation_threshold > 1)
             stop("Step 9: correlation threshold must be 0-1")
         },
         Step11Config = {
           if (!(config$selection_criterion %in% c("AIC", "AICc", "BIC")))
             stop("Step 11: selection criterion must be AIC, AICc, or BIC")
         },
         Step12Config = {
           if (config$max_beta < 1)
             stop("Step 12: max_beta must be >= 1")
           if (config$beta_increment <= 0)
             stop("Step 12: beta_increment must be > 0")
         },
         Step13Config = {
           if (config$replicates < 1)
             stop("Step 13: replicates must be >= 1")
         },
         # Default branch: no validation for unlisted types
         TRUE
  )
  
  return(TRUE)
}