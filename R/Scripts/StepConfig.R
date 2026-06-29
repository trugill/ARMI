# =============================================================================
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

