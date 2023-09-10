setwd("/Users/mateopava/Desktop/CNA_Mateo_MSC/allpanel")
all_files <- list.files()

# get all the files to read them in a df
panel_files <- all_files[grepl(".cna.tsv$", all_files)]

list_of_panel <- lapply(panel_files, function(file) {
  tryCatch({
    df <- read.csv(file, sep = "\t")
    
    # Extract the unique identifier from the filename without the "_cna.tsv" part
    identifier <- gsub("_cna\\.tsv$", "", basename(file))
    
    # Add a new column to store the identifier
    df$Origin <- identifier
    return(df)
  }, error = function(e) {
    message(paste("Error reading file:", file))
    return(NULL)
  })
})

# Assuming genes_of_interest is already defined
genes_of_interest
# Assuming genes_of_interest is already defined

# Function to expand based on genes and then filter
filter_genes_expanded <- function(df) {
  # Expand the genes
  expanded_df <- df[rep(seq_len(nrow(df)), times = sapply(strsplit(as.character(df$gene), ","), length)), ]
  expanded_df$gene <- unlist(strsplit(as.character(df$gene), ","))
  
  # Filter based on genes of interest
  filtered_df <- expanded_df[expanded_df$gene %in% genes_of_interest, ]
  
  return(filtered_df)
}

# Apply the function to each dataframe in the list
list_of_panel_filtered <- lapply(list_of_panel, filter_genes_expanded)

# Display the first few rows of the first filtered dataframe
head(list_of_panel_filtered[[1]])


# Merge log_df with compare_tech to add P300_fileID column
updated_log_df <- merge(log_df, compare_tech[, c("PRO_ID", "P300_fileID")], 
                        by.x = "PRO_simple_id", by.y = "PRO_ID", all.x = TRUE)

####

# Remove the 'length' column from each dataframe in the list if it exists
list_of_panel_filtered_cleaned <- lapply(list_of_panel_filtered, function(df) {
  if ("length" %in% colnames(df)) {
    df$length <- NULL
  }
  return(df)
})

list_of_panel_standardized <- lapply(list_of_panel_filtered_cleaned, function(df) {
  if ("length..in.MB." %in% colnames(df)) {
    df$length..in.MB. <- NULL
  }
  return(df)
})

panel_df <- do.call(rbind, list_of_panel_standardized)


mostmix <- merge(updated_log_df, panel_df, 
                   by.x = c("P300_fileID", "Gene"), 
                   by.y = c("Origin", "gene"), 
                   all.x = TRUE)
# Subset dataframe to retain only the desired columns
subset_mostmix <- mostmix[, c("PRO_id", "Gene", "CNVkit_cn", "ASCAT_cn", "ichorcna_log2", "qdnaseq_log2", "log2", "X..copias..normalizado.por.TA.")]

# Rename columns
colnames(subset_mostmix) <- c("PRO_id", "Gene", "CNVkit_cn", "ASCAT_cn", "ichorcna_log2", "qdnaseq_log2", "panel_log2", "panel_cn")

# Reorder columns so PRO_id and Gene are the first two columns
subset_mostmix <- subset_mostmix[, c("PRO_id", "Gene", "CNVkit_cn", "ASCAT_cn", "ichorcna_log2", "qdnaseq_log2", "panel_log2", "panel_cn")]





