# =============================================================================
# PathManager.R - Centralized path management
# =============================================================================

PathManager <- new.env()

PathManager$initialize <- function(species_file_path) {
  species_file_path <- normalizePath(species_file_path, mustWork = FALSE)
  PathManager$speciesPath <- species_file_path
  PathManager$speciesName <- tools::file_path_sans_ext(basename(species_file_path))
  
  # Walk up the tree to find the species root folder.
  # Structure is: <MAIN_PATH>/output/<Genus_species>/CSV[/Species]/<file>.csv
  # We want projectPath = <MAIN_PATH>/output/<Genus_species>
  #
  # Strategy: start from the file's directory and walk up until we find a
  # directory whose parent is "output" (or hit a sensible fallback).
  candidate <- dirname(species_file_path)
  project_path <- NULL
  for (i in 1:5) {  # walk up at most 5 levels
    parent <- dirname(candidate)
    if (basename(parent) == "output") {
      project_path <- candidate
      break
    }
    if (candidate == parent) break  # reached filesystem root
    candidate <- parent
  }
  
  # Fallback: assume CSV file is directly under <projectPath>/CSV/
  if (is.null(project_path)) {
    if (basename(dirname(species_file_path)) == "CSV") {
      project_path <- dirname(dirname(species_file_path))
    } else if (basename(dirname(dirname(species_file_path))) == "CSV") {
      project_path <- dirname(dirname(dirname(species_file_path)))
    } else {
      project_path <- dirname(species_file_path)  # last resort
    }
  }
  
  PathManager$projectPath <- project_path
  
  # All paths now hang cleanly off projectPath = <MAIN_PATH>/output/<Genus_species>/
  PathManager$csvPath        <- file.path(project_path, "CSV")
  PathManager$layersPath     <- file.path(project_path, "ClippedLayers")
  PathManager$clippedMaskPath <- file.path(project_path, "ClippedMask")
  PathManager$tifPath        <- file.path(project_path, "TIF")
  PathManager$modelsPath     <- file.path(project_path, "Models")
  PathManager$shapefilePath  <- file.path(project_path, "SHP")
  PathManager$maskPath       <- file.path(project_path, "Mask")
  PathManager$speciesDir     <- file.path(project_path, "Species")
  
  PathManager$correlationPath    <- file.path(PathManager$csvPath, "correlation.csv")
  PathManager$modelSelectionPath <- file.path(PathManager$csvPath, "modelSelection.csv")
  
  # Create essential directories
  dirs_to_create <- c(PathManager$csvPath,
                      PathManager$layersPath,
                      PathManager$clippedMaskPath,
                      PathManager$modelsPath,
                      PathManager$shapefilePath,
                      PathManager$speciesDir)
  for (d in dirs_to_create) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  cat("PathManager initialized:\n")
  cat("  projectPath:   ", PathManager$projectPath, "\n")
  cat("  speciesName:   ", PathManager$speciesName, "\n")
  cat("  speciesPath:   ", PathManager$speciesPath, "\n")
  cat("  csvPath:       ", PathManager$csvPath, "\n")
  cat("  modelsPath:    ", PathManager$modelsPath, "\n")
  cat("  layersPath:    ", PathManager$layersPath, "\n")
  cat("  clippedMaskPath:", PathManager$clippedMaskPath, "\n")
  cat("  shapefilePath: ", PathManager$shapefilePath, "\n")
  
  PathManager$savePaths()
}

PathManager$savePaths <- function() {
  paths_file <- file.path(PathManager$projectPath, "paths.txt")
  paths_content <- c(
    paste0("speciesPath=", PathManager$speciesPath),
    paste0("speciesName=", PathManager$speciesName),
    paste0("projectPath=", PathManager$projectPath),
    paste0("csvPath=", PathManager$csvPath),
    paste0("layersPath=", PathManager$layersPath),
    paste0("clippedMaskPath=", PathManager$clippedMaskPath),
    paste0("tifPath=", PathManager$tifPath),
    paste0("modelsPath=", PathManager$modelsPath),
    paste0("correlationPath=", PathManager$correlationPath),
    paste0("modelSelectionPath=", PathManager$modelSelectionPath),
    paste0("shapefilePath=", PathManager$shapefilePath),
    paste0("maskPath=", PathManager$maskPath)
  )
  writeLines(paths_content, paths_file)
}

PathManager$trimDupePath <- function() {
  base_name <- tools::file_path_sans_ext(basename(PathManager$speciesPath))
  new_name <- paste0(base_name, ".trimmed.csv")
  PathManager$speciesPath <- file.path(dirname(PathManager$speciesPath), new_name)
  PathManager$speciesName <- tools::file_path_sans_ext(new_name)
  PathManager$savePaths()
  return(PathManager$speciesPath)
}

PathManager$getSpeciesPath <- function() PathManager$speciesPath
PathManager$getSpeciesName <- function() PathManager$speciesName
PathManager$getLayersPath <- function() PathManager$layersPath
PathManager$getClippedMaskPath <- function() PathManager$clippedMaskPath
PathManager$getTifPath <- function() PathManager$tifPath
PathManager$getModelsPath <- function() PathManager$modelsPath
PathManager$getCSVPath <- function() PathManager$csvPath
PathManager$getCorrelationPath <- function() PathManager$correlationPath
PathManager$getModelSelectionPath <- function() PathManager$modelSelectionPath
PathManager$getShapefilePath <- function() PathManager$shapefilePath
PathManager$getMaskPath <- function() PathManager$maskPath
PathManager$getProjectPath <- function() PathManager$projectPath