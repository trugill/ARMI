# =============================================================================
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
  cat("WARNING: MaxEnt jar not found at:", MAXENT_JAR_PATH, "\n")
  cat("Set the correct path using: MaxentCaller$setMaxentJar(\"path/to/maxent.jar\")\n")
}

cat("\n========================================\n")
cat("R Automated Maxent-R System\n")
cat("========================================\n\n")
cat("Launching GUI...\n")

create_maxent_gui()

