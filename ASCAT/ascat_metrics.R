# GENE annotate for ASCAT
setwd("/Users/mateopava/Desktop/ascat/ascat_analyze")
all_files <- list.files()

# get all the files to read them in a df
ascat_metrics <- all_files[grepl(".metrics.txt$", all_files)]

#list it
list_of_metrics <- lapply(ascat_metrics, function(file) {
  df <- read.table(file, header = TRUE, sep = "\t") # assuming they are tab-delimited
  
  # Extract the unique identifier from the filename
  identifier <- gsub("^(.*)\\_vs.*\\.metrics\\.txt$", "\\1", file)
  
  # Add a new column to store the identifier
  df$Origin <- identifier
  
  return(df)
})

# Combine dataframes vertically using do.call with rbind
metrics <- do.call(rbind, list_of_metrics)

# Set Origin as row names and drop the Origin column
rownames(metrics) <- metrics$Origin
metrics$Origin <- NULL

# Save the table
write.table(metrics, file = "metrics.txt", row.names = F, col.names = T, quote = F, sep = "\t")


############

# Define a function to calculate the expected neutral copy number
calculate_expected_neutral_copy_number <- function(purity, ploidy) {
  return(ploidy * purity + 2 * (1 - purity))
}

# Iterate over each dataframe in ASCAT_annotated and calculate the expected neutral copy number
for (i in 1:length(ASCAT_annotated)) {
  origin <- unique(ASCAT_annotated[[i]]$Origin)
  
  # Check if the origin exists in the metrics dataframe
  if (origin %in% metrics$Origin) {
    current_purity <- metrics[metrics$Origin == origin, "purity"]
    current_ploidy <- metrics[metrics$Origin == origin, "ploidy"]
    
    expected_neutral_copy_number <- calculate_expected_neutral_copy_number(as.numeric(current_purity), as.numeric(current_ploidy))
    ASCAT_annotated[[i]]$ExpectedNeutralCopyNumber <- rep(expected_neutral_copy_number, nrow(ASCAT_annotated[[i]]))
  } else {
    # If the origin does not exist in the metrics dataframe, assign NA
    ASCAT_annotated[[i]]$ExpectedNeutralCopyNumber <- rep(NA, nrow(ASCAT_annotated[[i]]))
  }
}

# Return the head of the first dataframe in the list to check the changes
head(ASCAT_annotated[[1]])

## COPY NUMBER AFTER PLOIDY AND PURITY
# Iterate over each dataframe in ASCAT_annotated and calculate the rounded neutral copy number and inferred call
for (i in 1:length(ASCAT_annotated)) {
  
  # Round the ExpectedNeutral value to the nearest whole number
  ASCAT_annotated[[i]]$RoundedNeutral <- round(ASCAT_annotated[[i]]$ExpectedNeutralCopyNumber)
  
  # Infer the call based on the RoundedNeutral value and observed copy number (nSum)
  ASCAT_annotated[[i]]$InferredCall <- ifelse(
    ASCAT_annotated[[i]]$nSum == ASCAT_annotated[[i]]$RoundedNeutral, "neutral",
    ifelse(ASCAT_annotated[[i]]$nSum < ASCAT_annotated[[i]]$RoundedNeutral, "loss", "gain")
  )
}

# Return the head of the first dataframe in the list to check the changes
head(ASCAT_annotated[[1]])


#Introducing a Buffer: Given that the expected neutral value is 2.6 (close to 3), you might consider introducing a buffer. For instance:
#Segments with copy numbers of 2 or 3 can be considered "neutral" or "ambiguous."
#Segments with copy numbers < 2 would be considered "loss."
#Segments with copy numbers > 3 would be considered "gain."
#This approach acknowledges the uncertainty introduced by the non-integer expected neutral copy number.



