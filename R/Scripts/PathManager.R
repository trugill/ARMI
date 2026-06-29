# =============================================================================
# PathManager.R - Centralized path management
# =============================================================================

PathManager <- new.env()

PathManager$initialize <- function(species_file_path) {
  species_file_path <- normalizePath(species_file_path, mustWork = FALSE)
  
  PathManager$speciesPath <- species_file_path
  PathManager$speciesName <- tools::file_path_sans_ext(basename(species_file_path))
  PathManager$projectPath <- dirname(dirname(species_file_path))
  
  PathManager$csvPath <- file.path(PathManager$projectPath, "CSV")
  PathManager$layersPath <- file.path(PathManager$projectPath, "ClippedLayers")
  PathManager$tifPath <- file.path(PathManager$projectPath, "TIF")
  PathManager$modelsPath <- file.path(PathManager$projectPath, "Models", 
                                      PathManager$speciesName)
  PathManager$correlationPath <- file.path(PathManager$csvPath, "correlation.csv")
  PathManager$modelSelectionPath <- file.path(PathManager$csvPath, "modelSelection.csv")
  PathManager$shapefilePath <- file.path(PathManager$projectPath, "Shapefiles")
  PathManager$maskPath <- file.path(PathManager$projectPath, "Mask")
  PathManager$speciesDir <- file.path(PathManager$projectPath, "Species")
  
  dirs_to_create <- c(PathManager$csvPath, PathManager$layersPath, 
                      PathManager$modelsPath, PathManager$speciesDir)
  for (d in dirs_to_create) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
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
PathManager$getTifPath <- function() PathManager$tifPath
PathManager$getModelsPath <- function() PathManager$modelsPath
PathManager$getCSVPath <- function() PathManager$csvPath
PathManager$getCorrelationPath <- function() PathManager$correlationPath
PathManager$getModelSelectionPath <- function() PathManager$modelSelectionPath
PathManager$getShapefilePath <- function() PathManager$shapefilePath
PathManager$getMaskPath <- function() PathManager$maskPath
PathManager$getProjectPath <- function() PathManager$projectPath

