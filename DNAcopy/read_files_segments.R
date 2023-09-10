# Define the path
path <- "/Users/mateopava/Desktop/segment_files"

# Get a list of all ".txt" and ".tsv" files
files_summary <- list.files(path, pattern = "\\.summary\\.txt$", full.names = TRUE)
files_events <- list.files(path, pattern = "\\.events\\.tsv$", full.names = TRUE)

# Initialize empty lists to store the data
data_summary <- list()
data_events <- list()

# Read the ".txt" files
for (i in 1:length(files_summary)) {
  data_summary[[i]] <- read.table(files_summary[i], sep = "\t", header = TRUE)
  names(data_summary)[i] <- basename(files_summary[i])
}

# Read the ".tsv" files
for (i in 1:length(files_events)) {
  data_events[[i]] <- read.csv(files_events[i], sep = "\t", header = TRUE)
  names(data_events)[i] <- basename(files_events[i])
}

# Remove the 7th and 8th elements from the list which are PRO 16
data_events <- data_events[-c(7, 8)]
