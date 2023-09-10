# read correspondecny
correspond<- read.csv("/Users/mateopava/Desktop/sample_correspondencies_CNA_FINAL.tsv.completed.tsv.vf.tsv", sep = "\t")

#filter correspond to only wes and wgs available
correspond_filtered <- correspond[!is.na(correspond$WES_fileID) & !is.na(correspond$sWGS_fileID), ]

# counts for panel, WES, and sWGS

wes_count <- sum(!is.na(correspond$WES_fileID))
wgs_count <- sum(!is.na(correspond$sWGS_fileID))
p300_count <- sum(!is.na(correspond$P300_fileID))

print(paste("Number of WES_fileID: ", wes_count))
print(paste("Number of sWGS_fileID: ", wgs_count))
print(paste("Number of P300_fileID: ", p300_count))

# From 219 patients, 171 have sWGS data, 82 WES and 126 Panel.
# From those 219, 70 patients are comparable between sWGS and WES methodologies.

# Specify the path to directory
path_to_directory <- "/Users/mateopava/Desktop/Mateo_master_genecalls"

# Get a list of all csv files in the directory
wgs_files <- list.files(path_to_directory, pattern = "\\.tsv$")

# Use lapply to read all csv files into a list of data frames
list_of_dataframes <- lapply(paste0(path_to_directory, "/", wgs_files), function(x) read.csv(x, header=TRUE, sep="\t", skip = 4))

# this has all data
list_of_dataframes <- setNames(list_of_dataframes, wgs_files)
df_names <- names(list_of_dataframes)
length(df_names) # total number of entries
length(unique(df_names)) # number of unique entries

# now subset to only the 70 comparable
library(stringr)

# Extract 'PROxxx' from 'WES_fileID' column and store it in a character vector
PRO_ids <- str_extract(correspond_filtered$WES_fileID, "PRO\\d+")
length(PRO_ids) # total number of entries
length(unique(PRO_ids)) # number of unique entries

# filter to exract the 70
wgs_files_no_ext <- sub("\\.gene_calls\\.tsv$", "", wgs_files)

# list with 70 comparable
selected_dataframes <- list_of_dataframes[wgs_files_no_ext %in% correspond_filtered$sWGS_fileID]



