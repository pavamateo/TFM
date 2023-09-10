#### READ THE COMPARABLE

# Define the path
path <- "/Users/mateopava/Desktop/"

# Get all file names ending with .txt or .tsv
files <- list.files(path, pattern = "segmentjoin.summary.txt$|segmentjoin.events.tsv$")

# Read files
for (x in files) {
  full_path <- file.path(path, x)
  
  # Extract the number (XXX) from the filename
  number <- gsub("segmentjoin.(\\d+).summary.txt|segmentjoin.(\\d+).events.tsv", "\\1\\2", x)
  
  # Check if the file is a summary or events file
  if (grepl("segmentjoin.summary.txt$", x)) {
    # Read summary file
    var_name <- paste0("joinseg", number, "summary")
    assign(var_name, read.table(full_path, "\t", header = TRUE), envir = .GlobalEnv)
  } else if (grepl("segmentjoin.events.tsv$", x)) {
    # Read events file
    var_name <- paste0("joinseg", number)
    assign(var_name, read.csv(full_path, "\t", header = TRUE), envir = .GlobalEnv)
  }
}
