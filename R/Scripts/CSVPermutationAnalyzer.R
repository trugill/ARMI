# =============================================================================
# CSVPermutationAnalyzer.R - CSV analysis utilities
# =============================================================================

CSVPermutationAnalyzer <- new.env()

CSVPermutationAnalyzer$getTopPermutationImportance <- function(csv_file_path, 
                                                                threshold, 
                                                                top_count) {
  if (!file.exists(csv_file_path)) stop("CSV file not found: ", csv_file_path)
  
  data <- read.csv(csv_file_path, stringsAsFactors = FALSE)
  perm_cols <- grep("permutation.importance", names(data), ignore.case = TRUE, value = TRUE)
  
  if (length(perm_cols) == 0) {
    warning("No permutation importance columns found")
    return(character(0))
  }
  
  values <- as.numeric(data[1, perm_cols])
  names(values) <- gsub("\\.permutation\\.importance", "", perm_cols, ignore.case = TRUE)
  
  threshold_pct <- threshold * 100
  filtered <- values[!is.na(values) & values >= threshold_pct]
  sorted <- sort(filtered, decreasing = TRUE)
  top_n <- head(sorted, top_count)
  return(names(top_n))
}

CSVPermutationAnalyzer$findLowestScoreVariables <- function(csv_file_path, selection_type) {
  if (!file.exists(csv_file_path)) stop("CSV file not found: ", csv_file_path)
  
  data <- read.csv(csv_file_path, stringsAsFactors = FALSE)
  
  score_col <- switch(toupper(selection_type),
                      "AIC" = "aic_score",
                      "AICC" = "aicc_score",
                      "BIC" = "bic_score",
                      "aicc_score")
  
  if (!(score_col %in% names(data))) stop("Score column not found: ", score_col)
  
  scores <- as.numeric(data[[score_col]])
  valid_idx <- which(!is.na(scores))
  
  if (length(valid_idx) == 0) return(NULL)
  
  min_idx <- valid_idx[which.min(scores[valid_idx])]
  path <- data$ascii_file_path[min_idx]
  folder_name <- basename(dirname(path))
  
  return(list(variables = folder_name, score = scores[min_idx], path = path))
}

CSVPermutationAnalyzer$addToFirstEmptyRow <- function(file_path, data) {
  if (!file.exists(file_path)) {
    df <- as.data.frame(t(data), stringsAsFactors = FALSE)
    write.csv(df, file_path, row.names = FALSE)
  } else {
    line <- paste(data, collapse = ",")
    write(line, file = file_path, append = TRUE)
  }
}

CSVPermutationAnalyzer$clearDirectory <- function(directory_path) {
  if (!dir.exists(directory_path)) return(TRUE)
  files <- list.files(directory_path, full.names = TRUE, recursive = FALSE)
  for (f in files) {
    if (dir.exists(f)) unlink(f, recursive = TRUE, force = TRUE)
    else file.remove(f)
  }
  return(TRUE)
}

CSVPermutationAnalyzer$refreshCSV <- function(path) {
  if (file.exists(path)) file.remove(path)
  file.create(path)
}

CSVPermutationAnalyzer$deleteDirectory <- function(path) {
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
}

CSVPermutationAnalyzer$createDirectory <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

CSVPermutationAnalyzer$getFolderNames <- function(directory) {
  if (!dir.exists(directory)) return(character(0))
  list.dirs(directory, full.names = FALSE, recursive = FALSE)
}

CSVPermutationAnalyzer$getFileNames <- function(directory) {
  if (!dir.exists(directory)) return(character(0))
  files <- list.files(directory, full.names = FALSE, recursive = FALSE)
  files[!file.info(file.path(directory, files))$isdir]
}

