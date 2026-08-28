# Load libraries
library(ballgown)
library(rtracklayer)

# Simple test to check if file exists
test_path <- "\\\\wsl.localhost\\Ubuntu\\home\\brendan\\majiq_quant\\deltapsi_for_R.tsv"
cat("Testing path:", test_path, "\n")
cat("File exists:", file.exists(test_path), "\n")

if(file.exists(test_path)) {
  # Try to read first few lines
  cat("\nReading first 3 lines:\n")
  lines <- readLines(test_path, n = 3)
  for(i in 1:length(lines)) {
    cat("Line", i, ":", lines[i], "\n")
  }
  
  # Try to load the data
  cat("\nLoading full file...\n")
  majiq_data <- read.delim(test_path, stringsAsFactors = FALSE)
  cat("Success! Loaded", nrow(majiq_data), "rows and", ncol(majiq_data), "columns\n")
  cat("Column names:", paste(colnames(majiq_data), collapse = ", "), "\n")
} else {
  cat("\nFile not found. Listing files in directory:\n")
  dir_path <- "\\\\wsl.localhost\\Ubuntu\\home\\brendan\\majiq_quant"
  if(file.exists(dir_path)) {
    files <- list.files(dir_path)
    cat("Files in", dir_path, ":\n")
    for(f in files) {
      cat("  -", f, "\n")
    }
  } else {
    cat("Directory not found either.\n")
  }
}