# =============================================================================
# MaxentRSystemGUI.R - Shiny GUI for ARMI (Updated UI)
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)

# Null-coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a
})

# Null-coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

# Optional: leaflet for the map preview. Falls back to plot if unavailable.
HAS_LEAFLET <- requireNamespace("leaflet", quietly = TRUE)
if (HAS_LEAFLET) library(leaflet)

# ============================================================================
# HELPERS
# ============================================================================

# Discover folders within EnvironmentalLayers/ that already contain raster files
discover_downloaded_layers <- function() {
  if (!exists("ENV_LAYERS_DIR") || !dir.exists(ENV_LAYERS_DIR)) {
    return(character(0))
  }
  
  subdirs <- list.dirs(ENV_LAYERS_DIR, full.names = FALSE, recursive = FALSE)
  with_content <- character(0)
  
  for (d in subdirs) {
    full <- file.path(ENV_LAYERS_DIR, d)
    files <- list.files(full, pattern = "\\.(tif|tiff|asc|nc|img|grd)$",
                        ignore.case = TRUE, recursive = TRUE)
    if (length(files) > 0) {
      with_content <- c(with_content, d)
    }
  }
  return(with_content)
}

# Common bioclimatic variable list for autocomplete
COMMON_BIOCLIM_VARS <- c(
  "bio01_AnnualMeanTemp", "bio02_MeanDiurnalRange", "bio03_Isothermality",
  "bio04_TempSeasonality", "bio05_MaxTempWarmestMonth", "bio06_MinTempColdestMonth",
  "bio07_TempAnnualRange", "bio08_MeanTempWettestQuarter", "bio09_MeanTempDriestQuarter",
  "bio10_MeanTempWarmestQuarter", "bio11_MeanTempColdestQuarter",
  "bio12_AnnualPrecipitation", "bio13_PrecipWettestMonth", "bio14_PrecipDriestMonth",
  "bio15_PrecipSeasonality", "bio16_PrecipWettestQuarter", "bio17_PrecipDriestQuarter",
  "bio18_PrecipWarmestQuarter", "bio19_PrecipColdestQuarter",
  "elevation", "slope", "aspect", "trees", "built", "footprint"
)

DOWNLOAD_SOURCES <- c(
  "WorldClim Bioclim (19 vars)" = "worldclim",
  "SPECTRE Threat Layers"       = "spectre",
  "Land Cover (trees, built)"   = "landcover",
  "Human Footprint (2009)"      = "footprint",
  "Custom user rasters"         = "custom"
)

# ============================================================================
# CONFIG SAVE/LOAD HELPERS
# ============================================================================

save_gui_to_config_shiny <- function(input) {
  cfg <- ConfigManager$current
  
  # Species
  cfg$species$species_auto   <- if (is.null(input$species_auto)) character(0) else input$species_auto
  cfg$species$species_manual <- if (is.null(input$species_manual)) character(0) else input$species_manual
  cfg$species$species_list   <- unique(c(cfg$species$species_auto, cfg$species$species_manual))
  cfg$species$download_mode  <- input$species_dl_mode
  cfg$species$use_raw_gbif   <- TRUE   # Mandatory now
  cfg$species$trim_occurrences <- TRUE # Mandatory now
  cfg$species$use_iucn_range <- input$use_iucn_range
  
  # Data sources - multi-select tags
  selected_sources <- if (is.null(input$dl_sources)) character(0) else input$dl_sources
  cfg$data_sources$download_worldclim <- "worldclim" %in% selected_sources
  cfg$data_sources$download_spectre   <- "spectre"   %in% selected_sources
  cfg$data_sources$download_landcover <- "landcover" %in% selected_sources
  cfg$data_sources$download_footprint <- "footprint" %in% selected_sources
  cfg$data_sources$use_custom_layers  <- "custom"    %in% selected_sources
  cfg$data_sources$download_gbif      <- TRUE  # Always on; species list controls behavior
  cfg$data_sources$extract_iucn_range <- input$use_iucn_range
  cfg$data_sources$override_existing  <- input$override_downloads
  cfg$data_sources$active_layer_dirs  <- if (is.null(input$active_layers)) character(0) else input$active_layers
  
  # GBIF
  cfg$gbif$max_records    <- as.integer(input$gbif_max_records)
  cfg$gbif$page_size      <- as.integer(input$gbif_page_size)
  cfg$gbif$min_confidence <- as.integer(input$gbif_min_confidence)
  cfg$gbif$skip_existing  <- isTRUE(input$skip_existing_gbif)
  
  # Resolution
  cfg$resolution$mode             <- input$res_mode
  cfg$resolution$target_arcmin    <- as.numeric(input$target_res)
  cfg$resolution$auto_align       <- input$auto_align
  cfg$resolution$user_template_raster <-
    if (!is.null(input$template_file) && nchar(input$template_file$datapath) > 0)
      input$template_file$datapath else NULL
  
  # Extent
  cfg$extent$clipping_strategy <- input$clipping_strategy
  if (input$clipping_strategy == "manual_bbox") {
    cfg$extent$max_lat <- as.numeric(input$max_lat)
    cfg$extent$max_lon <- as.numeric(input$max_lon)
    cfg$extent$min_lat <- as.numeric(input$min_lat)
    cfg$extent$min_lon <- as.numeric(input$min_lon)
  }
  cfg$extent$clip_to_extent <- (input$clipping_strategy == "auto_extent")
  
  # Mask
  cfg$mask$use_mask <- input$use_mask
  if (input$use_mask) {
    cfg$mask$max_lat <- as.numeric(input$mask_max_lat)
    cfg$mask$max_lon <- as.numeric(input$mask_max_lon)
    cfg$mask$min_lat <- as.numeric(input$mask_min_lat)
    cfg$mask$min_lon <- as.numeric(input$mask_min_lon)
  }
  
  # Model parameters
  cfg$model_parameters$perm_importance_threshold <- as.numeric(input$var_imp)
  cfg$model_parameters$top_variable_count        <- as.integer(input$top_count)
  cfg$model_parameters$correlation_threshold     <- as.numeric(input$corr_val)
  cfg$model_parameters$selection_criterion       <- input$criterion
  cfg$model_parameters$max_combo_size            <- as.integer(input$max_combo_size)
  
  # Beta
  cfg$beta_optimization$min  <- as.numeric(input$min_beta)
  cfg$beta_optimization$max  <- as.numeric(input$max_beta)
  cfg$beta_optimization$step <- as.numeric(input$beta_step)
  
  # Final model
  cfg$final_model$replicates     <- as.integer(input$replicates)
  cfg$final_model$replicate_type <- "crossvalidate"  # Hardcoded
  
  # Advanced - required bioclim vars (selectize returns a character vector)
  cfg$advanced$required_bioclim <-
    if (is.null(input$required_bioclim)) character(0) else input$required_bioclim
  cfg$advanced$java_memory_mb       <- as.integer(input$java_memory)
  cfg$advanced$skip_existing_models <- input$skip_existing
  
  # Cleanup
  cfg$cleanup$keep_clipped_layers <- input$keep_clipped_layers
  cfg$cleanup$keep_clipped_mask   <- input$keep_clipped_mask
  cfg$cleanup$use_shp <- TRUE  # Always on, controlled by extent strategy now
  
  ConfigManager$current <- cfg
  ConfigManager$save()
  ConfigManager$applyToGlobals()
}

load_gui_from_config_shiny <- function(session) {
  cfg <- ConfigManager$current
  
  # Species - split into auto/manual
  sp_auto   <- if (!is.null(cfg$species$species_auto))   cfg$species$species_auto   else character(0)
  sp_manual <- if (!is.null(cfg$species$species_manual)) cfg$species$species_manual else character(0)
  
  # Legacy: if no split, treat all as auto
  if (length(sp_auto) == 0 && length(sp_manual) == 0 && 
      length(cfg$species$species_list) > 0) {
    sp_auto <- cfg$species$species_list
  }
  
  updateSelectizeInput(session, "species_auto",   selected = sp_auto)
  updateSelectizeInput(session, "species_manual", selected = sp_manual)
  
  dl_mode <- if (!is.null(cfg$species$download_mode)) cfg$species$download_mode else "automatic"
  updateSelectInput(session, "species_dl_mode", selected = dl_mode)
  
  updateCheckboxInput(session, "use_iucn_range",
                      value = if (!is.null(cfg$species$use_iucn_range)) 
                        cfg$species$use_iucn_range else cfg$data_sources$extract_iucn_range)
  
  # Data sources - rebuild multi-select
  selected_sources <- character(0)
  if (isTRUE(cfg$data_sources$download_worldclim)) selected_sources <- c(selected_sources, "worldclim")
  if (isTRUE(cfg$data_sources$download_spectre))   selected_sources <- c(selected_sources, "spectre")
  if (isTRUE(cfg$data_sources$download_landcover)) selected_sources <- c(selected_sources, "landcover")
  if (isTRUE(cfg$data_sources$download_footprint)) selected_sources <- c(selected_sources, "footprint")
  if (isTRUE(cfg$data_sources$use_custom_layers))  selected_sources <- c(selected_sources, "custom")
  
  updateSelectizeInput(session, "dl_sources", selected = selected_sources)
  
  updateCheckboxInput(session, "override_downloads",
                      value = isTRUE(cfg$data_sources$override_existing))
  
  active_layers <- if (!is.null(cfg$data_sources$active_layer_dirs))
    cfg$data_sources$active_layer_dirs else character(0)
  updateSelectizeInput(session, "active_layers",
                       choices = discover_downloaded_layers(),
                       selected = active_layers)
  
  # GBIF
  updateNumericInput(session, "gbif_max_records",    value = cfg$gbif$max_records)
  updateNumericInput(session, "gbif_page_size",      value = cfg$gbif$page_size)
  updateNumericInput(session, "gbif_min_confidence", value = cfg$gbif$min_confidence)
  updateCheckboxInput(session, "skip_existing_gbif",
                      value = if (!is.null(cfg$gbif$skip_existing))
                        cfg$gbif$skip_existing else TRUE)
  
  # Resolution
  updateSelectInput(session, "res_mode",  selected = cfg$resolution$mode)
  updateNumericInput(session, "target_res", value = cfg$resolution$target_arcmin)
  updateCheckboxInput(session, "auto_align",
                      value = if (!is.null(cfg$resolution$auto_align)) 
                        cfg$resolution$auto_align else TRUE)
  
  # Extent
  strategy <- if (!is.null(cfg$extent$clipping_strategy)) cfg$extent$clipping_strategy
  else if (isTRUE(cfg$extent$clip_to_extent)) "auto_extent" else "manual_bbox"
  updateSelectInput(session, "clipping_strategy", selected = strategy)
  
  updateNumericInput(session, "max_lat", value = cfg$extent$max_lat)
  updateNumericInput(session, "max_lon", value = cfg$extent$max_lon)
  updateNumericInput(session, "min_lat", value = cfg$extent$min_lat)
  updateNumericInput(session, "min_lon", value = cfg$extent$min_lon)
  
  # Mask
  updateCheckboxInput(session, "use_mask", value = cfg$mask$use_mask)
  updateNumericInput(session, "mask_max_lat", value = cfg$mask$max_lat)
  updateNumericInput(session, "mask_max_lon", value = cfg$mask$max_lon)
  updateNumericInput(session, "mask_min_lat", value = cfg$mask$min_lat)
  updateNumericInput(session, "mask_min_lon", value = cfg$mask$min_lon)
  
  # Model parameters
  updateSliderInput(session, "var_imp",     value = cfg$model_parameters$perm_importance_threshold)
  updateNumericInput(session, "top_count",  value = cfg$model_parameters$top_variable_count)
  updateSliderInput(session, "corr_val",    value = cfg$model_parameters$correlation_threshold)
  updateRadioButtons(session, "criterion",  selected = cfg$model_parameters$selection_criterion)
  updateNumericInput(session, "max_combo_size", value = cfg$model_parameters$max_combo_size)
  
  # Beta
  updateNumericInput(session, "min_beta",  value = cfg$beta_optimization$min)
  updateNumericInput(session, "max_beta",  value = cfg$beta_optimization$max)
  updateNumericInput(session, "beta_step", value = cfg$beta_optimization$step)
  
  # Final
  updateNumericInput(session, "replicates", value = cfg$final_model$replicates)
  
  # Advanced (moved required vars)
  req_bio <- if (!is.null(cfg$advanced$required_bioclim)) cfg$advanced$required_bioclim
  else if (!is.null(cfg$advanced$required_tifs)) cfg$advanced$required_tifs
  else character(0)
  updateSelectizeInput(session, "required_bioclim", selected = req_bio)
  
  updateNumericInput(session, "java_memory", value = cfg$advanced$java_memory_mb)
  updateCheckboxInput(session, "skip_existing", value = cfg$advanced$skip_existing_models)
  
  # Cleanup
  updateCheckboxInput(session, "keep_clipped_layers", value = cfg$cleanup$keep_clipped_layers)
  updateCheckboxInput(session, "keep_clipped_mask",   value = cfg$cleanup$keep_clipped_mask)
}

# ============================================================================
# STYLES
# ============================================================================

mini_style_css <- "
html, body {
  height: 100%;
  margin: 0;
  padding: 0;
  background-color: #f5f5f5;
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
  overflow: hidden;
}
.gadget-container {
  display: flex;
  flex-direction: column;
  height: 100vh;
  max-width: 1100px;
  margin: 0 auto;
  background: #fff;
  box-shadow: 0 0 10px rgba(0,0,0,0.1);
}
.gadget-title-bar {
  background: linear-gradient(to bottom, #fafafa, #e8e8e8);
  border-bottom: 1px solid #c8c8c8;
  padding: 10px 15px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  flex-shrink: 0;
}
.gadget-title-bar h4 {
  margin: 0;
  font-size: 14px;
  font-weight: 600;
  color: #333;
}
.gadget-title-bar .btn { font-size: 12px; padding: 4px 12px; }
.nav-tabs-custom {
  display: flex;
  background: #f0f0f0;
  border-bottom: 1px solid #c8c8c8;
  flex-shrink: 0;
  overflow-x: auto;
}
.nav-tabs-custom .nav-link {
  flex: 1;
  min-width: 90px;
  text-align: center;
  padding: 10px 5px;
  font-size: 11px;
  color: #555;
  border: none;
  border-right: 1px solid #d8d8d8;
  background: transparent;
  cursor: pointer;
  text-decoration: none;
}
.nav-tabs-custom .nav-link:hover { background: #e8e8e8; }
.nav-tabs-custom .nav-link.active {
  background: #fff;
  color: #428bca;
  border-bottom: 2px solid #428bca;
  font-weight: 600;
}
.nav-tabs-custom .nav-link i { display: block; font-size: 16px; margin-bottom: 3px; }
.gadget-content {
  flex: 1;
  overflow-y: auto;
  padding: 15px 20px;
}
.gadget-content h4 {
  font-size: 13px;
  font-weight: 600;
  color: #444;
  text-transform: uppercase;
  letter-spacing: 0.5px;
  margin-top: 10px;
  margin-bottom: 10px;
  border-bottom: 1px solid #eee;
  padding-bottom: 5px;
}
.gadget-content hr { margin: 15px 0; border-color: #eee; }
.form-group { margin-bottom: 10px; }
.gadget-content label { font-size: 12px; font-weight: 500; }
.status-bar {
  flex-shrink: 0;
  padding: 8px 15px;
  background: #f8f8f8;
  border-top: 1px solid #ddd;
  font-size: 11px;
  color: #555;
}
.status-bar progress {
  width: 100%;
  height: 14px;
  margin-top: 4px;
}
.tab-content > .tab-pane { display: none; }
.tab-content > .tab-pane.active { display: block; }
.help-text {
  font-size: 11px;
  color: #888;
  margin-top: 2px;
  margin-bottom: 8px;
}
.disclaimer {
  background: #fff8e1;
  border-left: 3px solid #ffc107;
  padding: 6px 10px;
  font-size: 11px;
  color: #5d4037;
  margin: 5px 0;
}
.map-container {
  border: 1px solid #ddd;
  border-radius: 4px;
  margin: 10px 0;
}
"

# ============================================================================
# UI
# ============================================================================

ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(HTML(mini_style_css)),
    tags$link(rel = "stylesheet",
              href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css")
  ),
  
  div(class = "gadget-container",
      
      # ===== Title bar =====
      div(class = "gadget-title-bar",
          actionButton("cancel", "Cancel", class = "btn-default btn-sm"),
          h4("ARMI - Automated R Maxent Integration"),
          actionButton("done", "Process", class = "btn-primary btn-sm")
      ),
      
      # ===== Tab strip =====
      div(class = "nav-tabs-custom",
          actionLink("tab_data",     HTML('<i class="fa fa-database"></i>Data'),       class = "nav-link active"),
          actionLink("tab_species",  HTML('<i class="fa fa-paw"></i>Species'),         class = "nav-link"),
          actionLink("tab_extent",   HTML('<i class="fa fa-globe"></i>Extent'),        class = "nav-link"),
          actionLink("tab_model",    HTML('<i class="fa fa-sliders-h"></i>Model'),     class = "nav-link"),
          actionLink("tab_steps",    HTML('<i class="fa fa-list-check"></i>Steps'),    class = "nav-link"),
          actionLink("tab_advanced", HTML('<i class="fa fa-cogs"></i>Advanced'),       class = "nav-link"),
          actionLink("tab_config",   HTML('<i class="fa fa-save"></i>Config'),         class = "nav-link")
      ),
      
      # ===== Content =====
      div(class = "gadget-content",
          div(class = "tab-content",
              
              # ─────── DATA TAB ───────
              div(class = "tab-pane active", id = "pane_data",
                  
                  h4("Environmental Layer Sources"),
                  selectizeInput("dl_sources", 
                                 "Layers to download:",
                                 choices = DOWNLOAD_SOURCES,
                                 selected = "worldclim",
                                 multiple = TRUE,
                                 options = list(plugins = list("remove_button"),
                                                placeholder = "Click to add sources...")),
                  div(class = "help-text",
                      "Type to filter. Click X to remove a tag."),
                  
                  fluidRow(
                    column(6,
                           actionButton("override_btn",
                                        HTML('<i class="fa fa-redo"></i> Override Existing Downloads'),
                                        class = "btn-warning btn-sm",
                                        style = "width: 100%; margin-top: 5px;"),
                           checkboxInput("override_downloads",
                                         "Re-download even if cached", FALSE)
                    ),
                    column(6,
                           actionButton("refresh_layers",
                                        HTML('<i class="fa fa-sync"></i> Refresh Layer List'),
                                        class = "btn-default btn-sm",
                                        style = "width: 100%; margin-top: 5px;")
                    )
                  ),
                  
                  br(),
                  selectizeInput("active_layers",
                                 "Already-downloaded layers to use:",
                                 choices = NULL,
                                 multiple = TRUE,
                                 options = list(plugins = list("remove_button"),
                                                placeholder = "Click Refresh to load...")),
                  div(class = "help-text",
                      "Lists folders in EnvironmentalLayers/ that contain raster files."),
                  
                  hr(),
                  
                  h4("Spatial Resolution"),
                  fluidRow(
                    column(6,
                           selectInput("res_mode", "Resolution mode:",
                                       choices = c(
                                         "Auto (lowest resolution of all rasters)" = "auto",
                                         "Manual (specify arc-minutes)"            = "manual"
                                       ),
                                       selected = "manual")
                    ),
                    column(6,
                           conditionalPanel(
                             condition = "input.res_mode == 'manual'",
                             numericInput("target_res", "Target (arc-min):",
                                          5, min = 0.5, max = 60, step = 0.5)
                           )
                    )
                  ),
                  
                  checkboxInput("auto_align", 
                                HTML("<b>Automatically align rasters</b> &mdash; resamples all rasters to a common grid"),
                                value = TRUE),
                  div(class = "help-text",
                      "When enabled, every raster is reprojected and resampled onto a single ",
                      "master template so MaxEnt sees a uniform stack."),
                  
                  fileInput("template_file", "Master template raster (optional):",
                            accept = c(".tif", ".tiff", ".asc")),
                  
                  hr(),
                  
                  h4("Advanced"),
                  selectizeInput("required_bioclim",
                                 "Required bioclimatic variables in final model:",
                                 choices = COMMON_BIOCLIM_VARS,
                                 multiple = TRUE,
                                 options = list(plugins = list("remove_button"),
                                                placeholder = "Optional - leave empty to allow any combination",
                                                create = TRUE)),
                  div(class = "help-text",
                      "Only permutations including ALL of these variables will be considered.",
                      " Type to add custom variable names.")
              ),
              
              # ─────── SPECIES TAB ───────
              div(class = "tab-pane", id = "pane_species",
                  
                  h4("Species List"),
                  
                  selectInput("species_dl_mode", "Data download option:",
                              choices = c(
                                "Automatic download from GBIF only" = "automatic",
                                "Manual CSV loading only"           = "manual",
                                "Both (use both sets)"              = "both"
                              ),
                              selected = "automatic"),
                  
                  conditionalPanel(
                    condition = "input.species_dl_mode == 'automatic' || input.species_dl_mode == 'both'",
                    selectizeInput("species_auto",
                                   "Species to auto-download from GBIF:",
                                   choices = NULL,
                                   multiple = TRUE,
                                   options = list(
                                     plugins = list("remove_button"),
                                     placeholder = "Type \"Genus species\" and press Enter...",
                                     create = TRUE,
                                     createOnBlur = TRUE
                                   )),
                    div(class = "help-text",
                        "Format: Genus species (will be converted to Genus_species internally).")
                  ),
                  
                  conditionalPanel(
                    condition = "input.species_dl_mode == 'manual' || input.species_dl_mode == 'both'",
                    selectizeInput("species_manual",
                                   "Species with manually-loaded CSVs:",
                                   choices = NULL,
                                   multiple = TRUE,
                                   options = list(
                                     plugins = list("remove_button"),
                                     placeholder = "Type species name as it appears in CSV...",
                                     create = TRUE,
                                     createOnBlur = TRUE
                                   )),
                    div(class = "help-text",
                        "These species must have a CSV at: ",
                        tags$code("output/<Genus_species>/CSV/<Genus_species>.csv"))
                  ),
                  
                  hr(),
                  
                  h4("Occurrence File Handling"),
                  div(class = "disclaimer",
                      icon("info-circle"), " ",
                      "All input CSVs are treated as raw GBIF tab-delimited format and ",
                      "will be rewritten to species, latitude, longitude format. ",
                      "Duplicate trimming is always performed."),
                  
                  hr(),
                  
                  h4("Range Shapefile"),
                  checkboxInput("use_iucn_range",
                                "Clip occurrences outside of range polygon",
                                value = FALSE),
                  div(class = "disclaimer",
                      icon("exclamation-triangle"), " ",
                      strong("Range polygon must be present in species SHP folder."),
                      tags$br(),
                      "Expected location: ",
                      tags$code("output/<Genus_species>/SHP/<Genus_species>_range.shp"),
                      tags$br(),
                      "This can be auto-extracted from IUCN data if you place shapefiles in ",
                      tags$code("SHP/IUCN/"), " and enable ",
                      strong("Extract IUCN range"), " on the Data tab.",
                      tags$br(),
                      "When enabled, the ", tags$code(".trimmed.csv"),
                      " will contain ONLY occurrences within the range polygon, ",
                      "with one point per raster cell.",
                      tags$br(),
                      "IUCN data: ",
                      tags$a(href = "https://www.iucnredlist.org/resources/spatial-data-download",
                             "IUCN Red List Spatial Data", target = "_blank")),
                  
                  hr(),
                  
                  h4("GBIF Settings"),
                  checkboxInput("skip_existing_gbif",
                                "Do not re-download GBIF data if CSV already exists",
                                value = TRUE),
                  div(class = "help-text",
                      "When enabled, species with an existing ",
                      tags$code("output/<Genus_species>/CSV/<Genus_species>.csv"),
                      " will be skipped during the GBIF download step."),
                  fluidRow(
                    column(4, numericInput("gbif_max_records", "Max records:",
                                           100000, min = 100, step = 1000)),
                    column(4, numericInput("gbif_page_size", "Page size:",
                                           300, min = 50, max = 300, step = 50)),
                    column(4, numericInput("gbif_min_confidence", "Min confidence:",
                                           90, min = 50, max = 100))
                  )
              ),
              
              # ─────── EXTENT TAB ───────
              div(class = "tab-pane", id = "pane_extent",
                  
                  h4("Clipping Strategy"),
                  selectInput("clipping_strategy", "Clip rasters to:",
                              choices = c(
                                "Extent of occurrence points"      = "auto_extent",
                                "Manual bounding box"              = "manual_bbox",
                                "IUCN range polygons"              = "iucn_range"
                              ),
                              selected = "auto_extent"),
                  
                  conditionalPanel(
                    condition = "input.clipping_strategy == 'manual_bbox'",
                    div(class = "help-text",
                        icon("map-marker-alt"), " ",
                        "Adjust the coordinates below to define your study area."),
                    fluidRow(
                      column(6, numericInput("max_lat", "Max Lat (North):", 90,
                                             min = -90, max = 90, step = 0.5)),
                      column(6, numericInput("max_lon", "Max Lon (East):", 180,
                                             min = -180, max = 180, step = 0.5))
                    ),
                    fluidRow(
                      column(6, numericInput("min_lat", "Min Lat (South):", -90,
                                             min = -90, max = 90, step = 0.5)),
                      column(6, numericInput("min_lon", "Min Lon (West):", -180,
                                             min = -180, max = 180, step = 0.5))
                    ),
                    div(class = "map-container",
                        if (HAS_LEAFLET) leafletOutput("extent_map", height = "300px")
                        else plotOutput("extent_map_static", height = "300px"))
                  ),
                  
                  conditionalPanel(
                    condition = "input.clipping_strategy == 'iucn_range'",
                    div(class = "disclaimer",
                        icon("info-circle"), " ",
                        "Each species' MaxEnt model will use its IUCN range polygon as the ",
                        "clipping boundary. Requires IUCN data in ", 
                        tags$code("SHP/IUCN/"), ".")
                  ),
                  
                  hr(),
                  
                  h4("Projection Mask"),
                  checkboxInput("use_mask", 
                                "Use separate mask for projection extent", FALSE),
                  conditionalPanel(
                    condition = "input.use_mask == true",
                    div(class = "help-text",
                        "Final model predictions will be projected onto this extent ",
                        "(useful for future climate projections or larger study regions)."),
                    fluidRow(
                      column(6, numericInput("mask_max_lat", "Max Lat:", 90,
                                             min = -90, max = 90, step = 0.5)),
                      column(6, numericInput("mask_max_lon", "Max Lon:", 180,
                                             min = -180, max = 180, step = 0.5))
                    ),
                    fluidRow(
                      column(6, numericInput("mask_min_lat", "Min Lat:", -90,
                                             min = -90, max = 90, step = 0.5)),
                      column(6, numericInput("mask_min_lon", "Min Lon:", -180,
                                             min = -180, max = 180, step = 0.5))
                    ),
                    div(class = "map-container",
                        if (HAS_LEAFLET) leafletOutput("mask_map", height = "300px")
                        else plotOutput("mask_map_static", height = "300px"))
                  )
              ),
              
              # ─────── MODEL TAB ───────
              div(class = "tab-pane", id = "pane_model",
                  h4("Variable Selection"),
                  sliderInput("var_imp", "Permutation Importance Threshold:",
                              min = 0, max = 1, value = 0.05, step = 0.01),
                  div(class = "help-text",
                      "Variables with importance below this fraction are dropped."),
                  fluidRow(
                    column(6, numericInput("top_count", "Number of top variables:",
                                           10, min = 1, max = 50)),
                    column(6, numericInput("max_combo_size", "Max combination size:",
                                           7, min = 2, max = 15))
                  ),
                  hr(),
                  h4("Correlation Filtering"),
                  sliderInput("corr_val", "Correlation Threshold (|r|):",
                              min = 0, max = 1, value = 0.80, step = 0.01),
                  div(class = "help-text",
                      "Combinations with any pair exceeding this |r| are excluded."),
                  hr(),
                  h4("Model Selection"),
                  radioButtons("criterion", "Information Criterion:",
                               choices = c("AIC" = "AIC", "AICc" = "AICc", "BIC" = "BIC"),
                               selected = "AICc", inline = TRUE),
                  hr(),
                  h4("Regularization (Beta Multiplier)"),
                  fluidRow(
                    column(4, numericInput("min_beta", "Min Beta:",
                                           1.0, min = 0.1, step = 0.1)),
                    column(4, numericInput("max_beta", "Max Beta:",
                                           5.0, min = 1.0, step = 0.5)),
                    column(4, numericInput("beta_step", "Increment:",
                                           0.5, min = 0.1, step = 0.1))
                  ),
                  hr(),
                  h4("Final Model"),
                  numericInput("replicates", "Cross-validation replicates:",
                               10, min = 1, max = 100),
                  div(class = "help-text",
                      "Replicate type: cross-validation (fixed).")
              ),
              
              # ─────── STEPS TAB ───────
              div(class = "tab-pane", id = "pane_steps",
                  div(style = "color:#428bca; font-size: 12px; margin-bottom: 10px;",
                      icon("info-circle"),
                      " Toggle individual workflow steps. Steps depend on previous outputs."),
                  fluidRow(
                    column(4, actionButton("select_all_steps",   "Select All",
                                           class = "btn-sm", style = "width: 100%;")),
                    column(4, actionButton("invert_steps",       "Invert Selection",
                                           class = "btn-sm btn-info", style = "width: 100%;")),
                    column(4, actionButton("deselect_all_steps", "Deselect All",
                                           class = "btn-sm", style = "width: 100%;"))
                  ),
                  br(),
                  checkboxGroupInput(
                    "steps", NULL,
                    choiceNames = c(
                      "Step 0: Download environmental data",
                      "Step 1-2: GBIF download + clip to shapefile",
                      "Step 3: Remove duplicates (spatial thinning)",
                      "Step 4-5: Clip rasters to study extent",
                      "Step 6: Run global model",
                      "Step 7: Identify top variables",
                      "Step 8: Calculate correlation matrix",
                      "Step 9: Generate variable permutations",
                      "Step 10: Run all permutations",
                      "Step 11: Identify best model (AIC/AICc/BIC)",
                      "Step 12: Optimize regularization (beta)",
                      "Step 13: Run final optimized model"
                    ),
                    choiceValues = c("step0","step1","step3","step4","step6",
                                     "step7","step8","step9","step10",
                                     "step11","step12","step13"),
                    selected = c("step0","step1","step3","step4","step6",
                                 "step7","step8","step9","step10",
                                 "step11","step12","step13")
                  )
              ),
              
              # ─────── ADVANCED TAB ───────
              div(class = "tab-pane", id = "pane_advanced",
                  h4("MaxEnt Runtime"),
                  fluidRow(
                    column(6, numericInput("java_memory", "Java heap size (MB):",
                                           2048, min = 512, step = 512)),
                    column(6, checkboxInput("skip_existing", 
                                            "Skip models already completed", TRUE))
                  ),
                  hr(),
                  h4("Cleanup"),
                  checkboxInput("keep_clipped_layers", 
                                "Keep ClippedLayers/ folder after completion", TRUE),
                  checkboxInput("keep_clipped_mask", 
                                "Keep ClippedMask/ folder after completion", TRUE),
                  div(class = "help-text",
                      "Disable to save disk space on multi-species runs."),
                  hr(),
                  h4("Diagnostic Info"),
                  verbatimTextOutput("diag_info")
              ),
              
              # ─────── CONFIG TAB ───────
              div(class = "tab-pane", id = "pane_config",
                  h4("Configuration Management"),
                  p("Settings are saved as JSON in docs/config/user_config.json.",
                    style = "font-size: 12px; color: #666;"),
                  actionButton("save_cfg", "Save Current Config",
                               icon = icon("save"), width = "100%",
                               class = "btn-primary"),
                  br(), br(),
                  fileInput("cfg_file", "Load config from JSON file:", 
                            accept = ".json"),
                  actionButton("reset_cfg", "Reset to Factory Defaults",
                               icon = icon("undo"), width = "100%",
                               class = "btn-warning"),
                  hr(),
                  h4("Current Config Path"),
                  verbatimTextOutput("current_config_path")
              )
          )
      ),
      
      # ===== Status / Progress bar =====
      div(class = "status-bar",
          div(textOutput("progress_label", inline = TRUE)),
          tags$progress(id = "progbar", value = "0", max = "14")
      )
  ),
  
  # JavaScript for tab switching and progress
  tags$script(HTML("
    $(document).on('click', '.nav-link', function(e) {
      e.preventDefault();
      var id = $(this).attr('id');
      var paneMap = {
        'tab_data':     'pane_data',
        'tab_species':  'pane_species',
        'tab_extent':   'pane_extent',
        'tab_model':    'pane_model',
        'tab_steps':    'pane_steps',
        'tab_advanced': 'pane_advanced',
        'tab_config':   'pane_config'
      };
      $('.nav-link').removeClass('active');
      $(this).addClass('active');
      $('.tab-pane').removeClass('active');
      $('#' + paneMap[id]).addClass('active');
    });
    Shiny.addCustomMessageHandler('setProg', function(m){
      var p = document.getElementById('progbar');
      if (p) p.value = m.val;
    });
  "))
)

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {
  
  # ===== Initialize selectize choices on startup =====
  observe({
    updateSelectizeInput(session, "active_layers",
                         choices = discover_downloaded_layers(),
                         server = FALSE)
  })
  
  observeEvent(input$refresh_layers, {
    updateSelectizeInput(session, "active_layers",
                         choices = discover_downloaded_layers(),
                         server = FALSE)
    showNotification(
      paste("Found", length(discover_downloaded_layers()), "layer source folder(s)"),
      type = "message", duration = 3)
  })
  
  # ===== Initialize from config =====
  observe({
    if (!is.null(ConfigManager$current)) {
      load_gui_from_config_shiny(session)
    }
  })
  
  # ===== Override existing downloads button =====
  observeEvent(input$override_btn, {
    showModal(modalDialog(
      tagList(
        p("This will mark all existing downloads to be re-downloaded on next run."),
        p(strong("Warning:"), " this can take significant time and bandwidth."),
        p("Affected sources: ", paste(input$dl_sources, collapse = ", "))
      ),
      title = "Override Existing Downloads",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_override", "Enable Override", class = "btn-warning")
      )
    ))
  })
  
  observeEvent(input$confirm_override, {
    updateCheckboxInput(session, "override_downloads", value = TRUE)
    removeModal()
    showNotification("Override mode enabled - next run will re-download",
                     type = "warning", duration = 6)
  })
  
  # ===== Progress =====
  progress_state <- reactiveValues(step = 0, name = "Ready")
  
  output$progress_label <- renderText({
    if (progress_state$step == 0) "Ready"
    else sprintf("Step %d: %s", progress_state$step, progress_state$name)
  })
  
  update_progress <- function(step, name) {
    progress_state$step <- step
    progress_state$name <- name
    session$sendCustomMessage("setProg", list(val = step))
  }
  
  # ===== Maps (extent + mask) =====
  if (HAS_LEAFLET) {
    output$extent_map <- renderLeaflet({
      leaflet() %>%
        addProviderTiles("OpenStreetMap.Mapnik") %>%
        setView(lng = 0, lat = 20, zoom = 2)
    })
    
    observe({
      req(input$clipping_strategy == "manual_bbox")
      req(input$min_lat, input$max_lat, input$min_lon, input$max_lon)
      
      leafletProxy("extent_map") %>%
        clearShapes() %>%
        addRectangles(
          lng1 = input$min_lon, lat1 = input$min_lat,
          lng2 = input$max_lon, lat2 = input$max_lat,
          fillColor = "red", fillOpacity = 0.2,
          color = "red", weight = 2
        )
    })
    
    output$mask_map <- renderLeaflet({
      leaflet() %>%
        addProviderTiles("OpenStreetMap.Mapnik") %>%
        setView(lng = 0, lat = 20, zoom = 2)
    })
    
    observe({
      req(input$use_mask)
      req(input$mask_min_lat, input$mask_max_lat, 
          input$mask_min_lon, input$mask_max_lon)
      
      leafletProxy("mask_map") %>%
        clearShapes() %>%
        addRectangles(
          lng1 = input$mask_min_lon, lat1 = input$mask_min_lat,
          lng2 = input$mask_max_lon, lat2 = input$mask_max_lat,
          fillColor = "blue", fillOpacity = 0.2,
          color = "blue", weight = 2
        )
    })
  } else {
    # Fallback to base R plot if leaflet unavailable
    output$extent_map_static <- renderPlot({
      req(input$clipping_strategy == "manual_bbox")
      plot(c(-180, 180), c(-90, 90), type = "n", 
           xlab = "Longitude", ylab = "Latitude",
           main = "Manual Bounding Box (install 'leaflet' for interactive map)")
      rect(-180, -90, 180, 90, col = "lightgray", border = "gray")
      abline(h = 0, v = 0, col = "white")
      rect(input$min_lon, input$min_lat, input$max_lon, input$max_lat,
           col = rgb(1, 0, 0, 0.3), border = "red", lwd = 2)
    })
    
    output$mask_map_static <- renderPlot({
      req(input$use_mask)
      plot(c(-180, 180), c(-90, 90), type = "n", 
           xlab = "Longitude", ylab = "Latitude",
           main = "Projection Mask")
      rect(-180, -90, 180, 90, col = "lightgray", border = "gray")
      abline(h = 0, v = 0, col = "white")
      rect(input$mask_min_lon, input$mask_min_lat, 
           input$mask_max_lon, input$mask_max_lat,
           col = rgb(0, 0, 1, 0.3), border = "blue", lwd = 2)
    })
  }
  
  # ===== Step selection buttons =====
  all_steps <- c("step0","step1","step3","step4","step6",
                 "step7","step8","step9","step10",
                 "step11","step12","step13")
  
  observeEvent(input$select_all_steps, {
    updateCheckboxGroupInput(session, "steps", selected = all_steps)
  })
  observeEvent(input$deselect_all_steps, {
    updateCheckboxGroupInput(session, "steps", selected = character(0))
  })
  observeEvent(input$invert_steps, {
    current <- input$steps
    inverted <- setdiff(all_steps, current)
    updateCheckboxGroupInput(session, "steps", selected = inverted)
  })
  
  # ===== Config save/load/reset =====
  observeEvent(input$save_cfg, {
    tryCatch({
      save_gui_to_config_shiny(input)
      showNotification(paste0("Saved to: ", ConfigManager$current_path),
                       type = "message", duration = 4)
    }, error = function(e) {
      showNotification(paste("Save failed:", e$message),
                       type = "error", duration = 8)
    })
  })
  
  observeEvent(input$cfg_file, {
    req(input$cfg_file)
    tryCatch({
      ConfigManager$load(input$cfg_file$datapath)
      load_gui_from_config_shiny(session)
      showNotification("Configuration loaded", type = "message")
    }, error = function(e) {
      showNotification(paste("Load failed:", e$message),
                       type = "error", duration = 8)
    })
  })
  
  observeEvent(input$reset_cfg, {
    showModal(modalDialog(
      "Reset all settings to factory defaults?",
      title = "Confirm Reset",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_reset", "Reset", class = "btn-warning")
      )
    ))
  })
  
  observeEvent(input$confirm_reset, {
    ConfigManager$resetToDefaults()
    load_gui_from_config_shiny(session)
    removeModal()
    showNotification("Reset to factory defaults", type = "warning")
  })
  
  # ===== Diagnostic info =====
  output$diag_info <- renderText({
    paste0(
      "R version:    ", R.version.string, "\n",
      "MAIN_PATH:    ", get0("MAIN_PATH", ifnotfound = "(not set)"), "\n",
      "MaxEnt jar:   ", get0("MAXENT_JAR", ifnotfound = "(not set)"), "\n",
      "MaxEnt found: ", file.exists(get0("MAXENT_JAR", ifnotfound = "")), "\n",
      "Leaflet:      ", if (HAS_LEAFLET) "available" else "NOT INSTALLED (run install.packages('leaflet'))"
    )
  })
  
  output$current_config_path <- renderText({
    if (!is.null(ConfigManager$current_path)) ConfigManager$current_path
    else "(not loaded)"
  })
  
  # ===== Process button =====
  observeEvent(input$done, {
    # Gather species from both lists
    sp_auto   <- if (is.null(input$species_auto))   character(0) else input$species_auto
    sp_manual <- if (is.null(input$species_manual)) character(0) else input$species_manual
    
    all_species <- switch(input$species_dl_mode,
                          "automatic" = sp_auto,
                          "manual"    = sp_manual,
                          "both"      = unique(c(sp_auto, sp_manual)),
                          character(0)
    )
    
    if (length(all_species) == 0) {
      showNotification(
        "Add at least one species to the relevant list.",
        type = "error", duration = 6)
      return()
    }
    
    showModal(modalDialog(
      tagList(
        p("Begin processing with these settings?"),
        p(strong("Species count: "), length(all_species)),
        p(strong("  Auto-download: "), length(sp_auto)),
        p(strong("  Manual CSV: "),    length(sp_manual)),
        p(strong("Steps enabled: "),   length(input$steps)),
        p(strong("Selection criterion: "), input$criterion),
        p(strong("Clipping strategy: "),   input$clipping_strategy)
      ),
      title = "Confirm Processing",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_process", "Start", class = "btn-primary")
      ),
      easyClose = TRUE
    ))
  })
  
  observeEvent(input$confirm_process, {
    removeModal()
    
    tryCatch({
      cat("=== [1/4] Saving config ===\n"); flush.console()
      save_gui_to_config_shiny(input)
      
      selected_steps <- input$steps %||% character(0)
      has_step <- function(k) k %in% selected_steps
      
      normalize_species <- function(s) gsub("\\s+", "_", trimws(s))
      sp_auto <- if (length(input$species_auto) > 0)
        vapply(input$species_auto, normalize_species, character(1))
      else character(0)
      sp_manual <- if (length(input$species_manual) > 0)
        vapply(input$species_manual, normalize_species, character(1))
      else character(0)
      
      species_to_process <- switch(input$species_dl_mode %||% "automatic",
                                   "automatic" = sp_auto,
                                   "manual"    = sp_manual,
                                   "both"      = unique(c(sp_auto, sp_manual)),
                                   character(0)
      )
      
      if (length(species_to_process) == 0) {
        showNotification("No species to process.", type = "error", duration = 6)
        return()
      }
      
      selected_sources <- input$dl_sources %||% character(0)
      cat("=== [2/4] Building workflow config ===\n"); flush.console()
      config <- WorkflowConfiguration()

      
      config$step0 <- StepConfig$Step0Config(
        enabled            = has_step("step0"),
        download_worldclim = "worldclim" %in% selected_sources,
        download_spectre   = "spectre"   %in% selected_sources,
        download_landcover = "landcover" %in% selected_sources,
        download_footprint = "footprint" %in% selected_sources,
        use_custom_layers  = "custom"    %in% selected_sources,
        download_gbif      = (input$species_dl_mode %||% "automatic") %in% c("automatic", "both"),
        extract_iucn_range = isTRUE(input$use_iucn_range),
        target_res_arcmin  = as.numeric(input$target_res %||% 5)
      )
      
      strategy <- input$clipping_strategy %||% "auto_extent"
      config$step4 <- StepConfig$Step4Config(
        enabled        = has_step("step4"),
        clip_to_extent = (strategy == "auto_extent"),
        max_lat = if (strategy == "manual_bbox") input$max_lat else NA,
        max_lon = if (strategy == "manual_bbox") input$max_lon else NA,
        min_lat = if (strategy == "manual_bbox") input$min_lat else NA,
        min_lon = if (strategy == "manual_bbox") input$min_lon else NA,
        use_mask = isTRUE(input$use_mask),
        mask_max_lat = if (isTRUE(input$use_mask)) input$mask_max_lat else NA,
        mask_max_lon = if (isTRUE(input$use_mask)) input$mask_max_lon else NA,
        mask_min_lat = if (isTRUE(input$use_mask)) input$mask_min_lat else NA,
        mask_min_lon = if (isTRUE(input$use_mask)) input$mask_min_lon else NA
      )
      
      config$step3 <- StepConfig$Step3Config(enabled = has_step("step3"))
      config$step6 <- StepConfig$Step6Config(enabled = has_step("step6"))
      config$step7 <- StepConfig$Step7Config(
        enabled = has_step("step7"),
        threshold = as.numeric(input$var_imp %||% 0.05),
        top_count = as.integer(input$top_count %||% 10)
      )
      config$step8 <- StepConfig$Step8Config(enabled = has_step("step8"))
      config$step9 <- StepConfig$Step9Config(
        enabled = has_step("step9"),
        correlation_threshold = as.numeric(input$corr_val %||% 0.8)
      )
      
      req_bio <- input$required_bioclim %||% character(0)
      
      config$step10 <- StepConfig$Step10Config(
        enabled = has_step("step10"),
        required_tifs = req_bio
      )
      config$step11 <- StepConfig$Step11Config(
        enabled = has_step("step11"),
        selection_criterion = input$criterion %||% "AICc"
      )
      config$step12 <- StepConfig$Step12Config(
        enabled = has_step("step12"),
        max_beta = as.numeric(input$max_beta %||% 5),
        beta_increment = as.numeric(input$beta_step %||% 0.5),
        selection_criterion = input$criterion %||% "AICc"
      )
      config$step13 <- StepConfig$Step13Config(
        enabled = has_step("step13"),
        replicates = as.integer(input$replicates %||% 10),
        projection_layers = character(0)
      )
      
      cat("=== [3/4] Species loop:", length(species_to_process), "species ===\n"); flush.console()
      
      withProgress(message = "Processing species...", value = 0, {
        for (i in seq_along(species_to_process)) {
          sp <- species_to_process[i]
          cat("\n--- [", i, "/", length(species_to_process), "] ", sp, " ---\n", sep = "")
          flush.console()
          incProgress(1 / length(species_to_process), detail = sp)
          
          config$step1 <- StepConfig$Step1Config(
            enabled = has_step("step1"),
            occurrence_file = sp,
            use_raw_gbif = TRUE,
            use_iucn_shp = isTRUE(input$use_iucn_range),
            shp_path = ""
          )
          
          result <- tryCatch({
            ProjectManager$runMethodology(config, update_progress)
            "success"
          }, error = function(e) {
            cat("ERROR inside runMethodology:", e$message, "\n")
            cat("--- Traceback ---\n")
            print(sys.calls())
            flush.console()
            paste("Error for", sp, ":", e$message)
          })
          
          if (!identical(result, "success")) {
            showNotification(as.character(result), type = "error", duration = 10)
          }
        }
      })
      
      cat("=== [4/4] Workflow completed ===\n"); flush.console()
      showNotification("Workflow completed!", type = "message", duration = 6)
      update_progress(14, "Complete")
      
    }, error = function(e) {
      cat("=== FATAL ERROR ===\n")
      cat("Message:", e$message, "\n")
      cat("Call:", deparse(e$call), "\n")
      flush.console()
      showNotification(paste("Processing error:", e$message),
                       type = "error", duration = 15)
    })
  })
  
  observeEvent(input$cancel, { stopApp() })
}

# ============================================================================
# LAUNCH
# ============================================================================
create_maxent_gui <- function() {
  app <- shinyApp(ui, server)
  
  # runGadget opens the app in a dedicated RStudio window/pane
  # instead of the default browser.
  viewer <- shiny::dialogViewer(
    dialogName = "ARMI - Automated R Maxent Integration",
    width  = 1200,
    height = 850
  )
  
  shiny::runGadget(app, viewer = viewer)
}
