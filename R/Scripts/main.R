# =============================================================================
# main.R - ARMI Application entry point (Shiny edition)
# =============================================================================

# Detect script directory
script_dir <- tryCatch({
  if (!is.null(sys.frames()[[1]]$ofile)) {
    dirname(normalizePath(sys.frames()[[1]]$ofile))
  } else if (requireNamespace("rstudioapi", quietly = TRUE) && 
             rstudioapi::isAvailable()) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    getwd()
  }
}, error = function(e) getwd())

cat("Script directory:", script_dir, "\n")

# Install required packages if missing
required_pkgs <- c("shiny", "shinyjs", "jsonlite", "terra", "raster", 
                   "sp", "sf", "dplyr", "tools")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs)
}

# Source modules in dependency order
source(file.path(script_dir, "ConfigManager.R"))
source(file.path(script_dir, "PathManager.R"))
source(file.path(script_dir, "CSVPermutationAnalyzer.R"))
source(file.path(script_dir, "CSVSpeciesExtractor.R"))
source(file.path(script_dir, "RFunctions.R"))
source(file.path(script_dir, "DataDownloader.R"))
source(file.path(script_dir, "OccurrenceThinner.R"))
source(file.path(script_dir, "MaxentCaller.R"))
source(file.path(script_dir, "StepConfig.R"))
source(file.path(script_dir, "ProjectManager.R"))
source(file.path(script_dir, "MaxentRSystemGUI.R"))
source(file.path(script_dir, "SHPExtractor.R"))


# Load configuration
ConfigManager$load()
ConfigManager$applyToGlobals()

# Verify MaxEnt
if (file.exists(MAXENT_JAR)) {
  MaxentCaller$setMaxentJar(MAXENT_JAR)
  cat("MaxEnt jar OK:", MAXENT_JAR, "\n")
} else {
  cat("WARNING: MaxEnt jar not found at:", MAXENT_JAR, "\n")
  cat("Place maxent.jar in:", BIN_DIR, "\n")
}

cat("\n========================================\n")
cat("  ARMI - Automated R Maxent Integration\n")
cat("========================================\n\n")
cat("Launching Shiny GUI in browser...\n")
cat("(Press ESC in console to stop the app)\n\n")

create_maxent_gui()