#grab the old environment
new_env <- new.env()

load("junio_28_2023.RData", new_env)


# genes<-new_env$filtered_genes153_all

# use Update the chromosome column of genesanalyze
genesanalyze$chromosome_name <- paste0("chr", genesanalyze$chromosome_name) # remember that before you were using without chr


# Create a function for the procedure
apply_procedure <- function(df, genesanalyze) {
  
  df$genesanalyze <- NA  # initialize an empty vector
  
  matching_genesanalyze <- vector("list", nrow(df))  # create a list to store matching genesanalyze for each row
  
  for (j in 1:nrow(df)) {
    # Find matching genesanalyze
    matching_genesanalyze[[j]] <- genesanalyze$hgnc_symbol[
      between(genesanalyze[, 3], df[j, 2], df[j, 3]) &
        genesanalyze[, 5] == df[j, 1] # df column 1 should be chromosome, gene column 5 should be chromosome
    ]
  }
  
  # If there are matching genesanalyze, paste them together separated by comma
  df$genesanalyze <- sapply(matching_genesanalyze, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
  
  df <- na.omit(df)  # remove rows with NA in 'genesanalyze' column
  
  # Split genesanalyze into separate rows
  df <- separate_rows(df, genesanalyze, sep = ",")
  
  return(df)
}

# Apply the procedure to each data frame in the list
data_events <- lapply(data_events, apply_procedure, genesanalyze = genesanalyze)  

# Add column to each dataframe
for (i in seq_along(data_events)) {
  data_events[[i]]$ID <- names(data_events)[i]
}

# Loop through each dataframe
for (i in seq_along(data_events)) {
  # Move the ID column to the first position
  data_events[[i]] <- data_events[[i]][, c("ID", setdiff(names(data_events[[i]]), "ID"))]
}

# below is for AR
ARresultsn_df$ID <- paste0("PRO", str_extract(ARresultsn_df$ID, "\\d+"))

# Loop through each dataframe in the list
for (i in seq_along(data_events)) {
  # Update the 'ID' column
  data_events[[i]]$ID <- paste0("PRO", stringr::str_extract(data_events[[i]]$ID, "\\d+"))
}


# now i have all the WES files with the 422 gene annotations

