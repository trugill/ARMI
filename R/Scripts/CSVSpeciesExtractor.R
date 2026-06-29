# =============================================================================
# CSVSpeciesExtractor.R - Extract species from GBIF tab-delimited files
# =============================================================================

CSVSpeciesExtractor <- new.env()

CSVSpeciesExtractor$processCSVFile <- function(input_file_path, is_temp = FALSE) {
  if (!file.exists(input_file_path)) stop("Input file does not exist: ", input_file_path)
  
  cat("Processing GBIF file:", input_file_path, "\n")
  
  data <- tryCatch({
    read.delim(input_file_path, sep = "\t", stringsAsFactors = FALSE,
               quote = "", fill = TRUE, header = TRUE)
  }, error = function(e) stop("Error reading file: ", e$message))
  
  col_names <- tolower(names(data))
  species_col <- which(grepl("species", col_names))[1]
  lat_col <- which(grepl("decimallatitude|latitude", col_names))[1]
  lon_col <- which(grepl("decimallongitude|longitude", col_names))[1]
  
  if (is.na(species_col) || is.na(lat_col) || is.na(lon_col)) {
    stop("Required columns not found. Need: species, latitude, longitude")
  }
  
  output_dir <- if (is_temp) {
    file.path(dirname(input_file_path), "temp_species")
  } else {
    file.path(dirname(input_file_path), "Species")
  }
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  species_names <- unique(data[[species_col]])
  species_names <- species_names[!is.na(species_names) & nchar(species_names) > 0]
  
  cat("Found", length(species_names), "unique species\n")
  
  last_output <- NULL
  for (sp in species_names) {
    clean_name <- gsub("[^a-zA-Z0-9]", "_", sp)
    sp_rows <- data[data[[species_col]] == sp & 
                    !is.na(data[[lat_col]]) & 
                    !is.na(data[[lon_col]]), ]
    if (nrow(sp_rows) == 0) next
    
    out_df <- data.frame(
      species = clean_name,
      latitude = as.numeric(sp_rows[[lat_col]]),
      longitude = as.numeric(sp_rows[[lon_col]]),
      stringsAsFactors = FALSE
    )
    out_df <- out_df[complete.cases(out_df), ]
    if (nrow(out_df) == 0) next
    
    output_file <- file.path(output_dir, paste0(clean_name, ".csv"))
    write.csv(out_df, output_file, row.names = FALSE)
    last_output <- output_file
    cat("  Saved:", clean_name, "(", nrow(out_df), "records)\n")
  }
  
  return(last_output)
}

