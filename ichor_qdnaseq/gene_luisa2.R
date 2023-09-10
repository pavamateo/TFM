
genesqands <- list.files(path = "/Users/mateopava/Desktop/Luisa_DATOS/oncogenecalls", full.names = TRUE)

ichorCNA_files <- grep("ichorCNA.tsv$", genesqands, value = TRUE)
qdnaseq_files <- grep("qdnaseq.tsv$", genesqands, value = TRUE)

# Adjusted function to modify column names
adjust_columns <- function(df_list, file_names) {
  # Adjust the column names for the first dataframe
  prefix <- tools::file_path_sans_ext(basename(file_names[1]))
  colnames(df_list[[1]]) <- paste(prefix, colnames(df_list[[1]]), sep="_")
  
  # Adjust the column names and drop the genes column for the remaining dataframes
  for(i in 2:length(df_list)) {
    prefix <- tools::file_path_sans_ext(basename(file_names[i]))
    colnames(df_list[[i]]) <- paste(prefix, colnames(df_list[[i]]), sep="_")
    df_list[[i]] <- df_list[[i]][,-1] # Drop the genes column
  }
  
  return(df_list)
}

# Adjust the columns for qdnaseq_list
adjust_qdnaseq_columns <- function(df_list, file_names) {
  # Rename the log2 column of the first dataframe
  prefix <- tools::file_path_sans_ext(basename(file_names[1]))
  colnames(df_list[[1]])[2] <- paste(prefix, "log2R_qdnaseq", sep="_")
  
  # For the subsequent dataframes, drop the genes column and rename the log2 column
  for(i in 2:length(df_list)) {
    prefix <- tools::file_path_sans_ext(basename(file_names[i]))
    colnames(df_list[[i]]) <- paste(prefix, "log2R_qdnaseq", sep="_")
  }
  
  return(df_list)
}

ichorCNA_list <- lapply(ichorCNA_files, read.table, header = TRUE, sep = "\t")
qdnaseq_list <- lapply(qdnaseq_files, read.table, header = TRUE, sep = "\t")

# Adjust the columns for both ichorCNA_list and qdnaseq_list
ichorCNA_list <- adjust_columns(ichorCNA_list, ichorCNA_files)
qdnaseq_list <- adjust_qdnaseq_columns(qdnaseq_list, qdnaseq_files)

# Column bind the adjusted dataframes
ichorCNA_df <- do.call(cbind, ichorCNA_list)
colnames(ichorCNA_df)[1] <- "gene"
qdnaseq_df <- do.call(cbind, qdnaseq_list)
colnames(qdnaseq_df)[1] <- "gene"

genes<- qdnaseq_df[1]
genes<- as.vector(genes)

## use dnacopy varscan WES
filtered_data_events <- lapply(data_events, function(df) {
  df[df$genes %in% genes$gene, ]
})
#
filtered_ordered_data_events <- lapply(data_events, function(df) {
  filtered_df <- df[df$genes %in% genes$gene, ]
  ordered_df <- filtered_df[order(filtered_df$genes), ]
  return(ordered_df)
})

#
transformed_data_events <- lapply(filtered_ordered_data_events, function(df) {
  df[, c("ID", "genes", "seg_mean", "event_type")]
})

### intento de df
# Create template dataframe
template_df <- data.frame(gene = genes$gene)

# Function to process each dataframe in the list
process_fun <- function(df, name) {
  # Merge with template
  merged_df <- merge(template_df, df, by.x = "gene", by.y = "genes", all.x = TRUE)
  
  # Fill missing values
  merged_df$ID[is.na(merged_df$ID)] <- name
  merged_df$seg_mean[is.na(merged_df$seg_mean)] <- "NA"
  merged_df$event_type[is.na(merged_df$event_type)] <- "NA"
  
  # Rename columns
  colnames(merged_df)[2:4] <- paste(name, colnames(merged_df)[2:4], sep = "_")
  
  return(merged_df)
}

# Apply the function to each dataframe in the list ACA ES DNACOPY
processed_dfs <- lapply(seq_along(transformed_data_events), function(i) {
  process_fun(transformed_data_events[[i]], names(transformed_data_events)[i])
})

# Combine all processed dataframes into a single dataframe
combined_df <- Reduce(function(...) merge(..., by = "gene", all = TRUE), processed_dfs)

# Order by the original gene list
combined_df <- combined_df[order(match(combined_df$gene, genes$gene)), ]

# Identify columns that have 'ID' in their names
id_columns <- grep("ID", colnames(combined_df))

# Remove these columns
combined_df <- combined_df[,-id_columns]

# now have to change the name from joinxxx to PRO_xxx_

# Diagnostic step 1: Print extracted segment numbers
extracted_segment_numbers <- sapply(colnames(combined_df), function(seg_col) {
  if (grepl("segment", seg_col)) {
    return(gsub("segment([0-9]{3})join\\.events\\.tsv", "\\1", seg_col))
  } else {
    return(NA)
  }
})
print(extracted_segment_numbers)

# Diagnostic step 2: Print matched PROxxx IDs
matched_pro_ids <- sapply(extracted_segment_numbers, function(segment_number) {
  if (!is.na(segment_number)) {
    return(paste0("PRO", sprintf("%03d", as.numeric(segment_number))))
  } else {
    return(NA)
  }
})
print(matched_pro_ids)

# Diagnostic step 3: Print final matching column names
final_matching_colnames <- sapply(matched_pro_ids, function(matching_pro_id) {
  if (!is.na(matching_pro_id)) {
    return(colnames(ichorCNA_df)[grepl(paste0(matching_pro_id, ".*_ichorCNA_log2R_ichorCNA"), colnames(ichorCNA_df))])
  } else {
    return(NA)
  }
})
print(final_matching_colnames)

###################

# Extract and rename columns based on the segment number
new_colnames <- sapply(colnames(combined_df), function(col) {
  if (grepl("segment", col)) {
    return(gsub("segment([0-9]+).*", "PRO\\1", col))
  } else {
    return(col)
  }
})

# Assign the new column names to combined_df
colnames(combined_df) <- new_colnames

# Check the updated column names
head(colnames(combined_df))

##

# Exclude the 'gene' column for renaming
pro_cols <- colnames(combined_df)[grepl("^PRO[0-9]+", colnames(combined_df)) & colnames(combined_df) != "gene"]

# Create a sequence for renaming
rename_sequence <- rep(c("_log2", "_call"), length(pro_cols)/2)

# Rename the columns
colnames(combined_df)[colnames(combined_df) %in% pro_cols] <- paste0(pro_cols, rename_sequence)

# Check the updated column names
head(colnames(combined_df))

###

# Extract PROxxx patterns from combined_df
pro_patterns_combined <- unique(gsub("(PRO[0-9]+)_.*", "\\1", colnames(combined_df)[grepl("^PRO[0-9]+_", colnames(combined_df))]))

log2R_cols_ichor <- colnames(ichorCNA_df)[grepl("log2R", colnames(ichorCNA_df))]
print(log2R_cols_ichor)

merged_data <- merge(ichorCNA_df, combined_df, by = "gene")

#### try corrleation

# 1. Identify columns
ichor_columns <- grep("*.ichorCNA_log2R_ichorCNA", names(merged_data), value=TRUE)
log2_columns <- grep("PRO[0-9]+_log2", names(merged_data), value=TRUE)

# Ensure that the columns are in the same order for correlation
ichor_columns <- sort(ichor_columns)
log2_columns <- sort(log2_columns)

# Convert the character columns in log2_columns to numeric
merged_data[log2_columns] <- lapply(merged_data[log2_columns], as.numeric)

# Now, you can re-run the correlation
correlations <- sapply(1:length(ichor_columns), function(i) {
  cor(merged_data[[ichor_columns[i]]], merged_data[[log2_columns[i]]], method = "pearson", use = "complete.obs")
})

correlations

correlations_df <- data.frame(ichor_column = ichor_columns, log2_column = log2_columns, correlation = correlations)

correlations_df

####try to merge them by id
add_suffix_to_duplicated <- function(ids) {
  # Find duplicated IDs
  dups <- duplicated(ids)
  
  # Add ".1" suffix to the second occurrence
  ids[dups] <- paste0(ids[dups], ".1")
  
  return(ids)
}


ichor_ids <- gsub("(.+?)_.+\\.ichorCNA_log2R_ichorCNA", "\\1", ichor_columns)
ichor_ids <- add_suffix_to_duplicated(ichor_ids)
log2_ids <- gsub("(.+?)(\\.1)?_log2", "\\1\\2", log2_columns)

common_ids <- intersect(ichor_ids, log2_ids)

matched_ichor_columns <- ichor_columns[ichor_ids %in% common_ids]
matched_ichor_columns <- setdiff(matched_ichor_columns, "PRO054_19-9250_DNAFF_NA_PC_HN00154307.genecalls.ichorCNA_log2R_ichorCNA")
matched_log2_columns <- log2_columns[log2_ids %in% common_ids]

##
# Create a new dataframe with the "gene" column and the desired columns
subset_ichorCNA_df <- ichorCNA_df[,matched_ichor_columns]

# View the first few rows of the subset dataframe
head(subset_ichorCNA_df)


# Initialize a vector to store the correlations
correlation_values <- numeric(length(matched_ichor_columns))

# Calculate correlations for each pair
for (i in 1:length(matched_ichor_columns)) {
  ichor_col <- merged_data[[matched_ichor_columns[i]]]
  log2_col <- merged_data[[matched_log2_columns[i]]]
  
  correlation_values[i] <- cor(ichor_col, log2_col, method = "pearson", use = "complete.obs")
}
# Create a named list with the column pairs as names and correlations as values
correlations <- setNames(as.list(correlation_values), paste(matched_ichor_columns, matched_log2_columns, sep = " vs. "))

# Display the correlations
correlations

#########3

get_pro_ids <- function(column_names) {
  # Use regex to extract the unique identifier
  return(gsub("(PRO[0-9]+).*", "\\1", column_names))
}

log2_ids <- get_pro_ids(matched_log2_columns)
ichor_ids <- get_pro_ids(matched_ichor_columns)


######3

# Check which genes are in the dataframe
genes_present <- genes$gene[genes$gene %in% new_env$filtered_genes153_all$hgnc_symbol]

# Print the genes that are present
print(genes_present)

gene_analyze <- new_env$filtered_genes153_all[new_env$filtered_genes153_all$hgnc_symbol %in% genes_present, ]

# View the subsetted dataframe
View(gene_analyze)

write.table(gene_analyze, file = "genes.csv", row.names = F, col.names = T, quote = F, sep = "\t")

# subsetting to compare log2R
swgs_ichorv <- read_xlsx("/Users/mateopava/Desktop/Luisa_DATOS/logR_qdnaseqVsichorCNA.xlsx")
# Extract unique PRO identifiers from matched_ichor_columns
pro_full <- gsub("(PRO\\d+).*", "\\1", matched_ichor_columns)



#############3visaje
# Extract column names from subset_ichorCNA_df without the "Gene" column
ichor_columns <- colnames(subset_ichorCNA_df)

# Extract the identifier portion before "genecalls.ichorCNA_log2R_"
identifiers <- gsub("(.*)genecalls.ichorCNA_log2R_.*", "\\1", ichor_columns)

identifiers <- gsub("(_[^-]+)-", "\\1.", identifiers)

# Replace the identifier with an empty string to get corresponding qdnaseq columns
qdnaseq_columns <- paste0(identifiers, "qdnaseq")

# Extract the required columns from swgs_ichorv
qdnaseq_dfs <- swgs_ichorv %>% select(c("Gene", qdnaseq_columns))

# Merge the two data frames based on the "Gene" column
merged_df <- merge(subset_ichorCNA_df, qdnaseq_df, by = "Gene")



############ ultimo
# Extract column names from subsetichocna without the "Gene" column
ichor_columns <- colnames(subset_ichorCNA_df)

# Replace ".ichorCNA" with ".qdnaseq" to get corresponding qdnaseq columns
qdnaseq_columns <- gsub(".ichorCNA", ".qdnaseq", ichor_columns)

# Extract the required columns from swgs_ichorv
qdnaseq_dfs <- swgs_ichorv %>% select(c("Gene", qdnaseq_columns))

# Merge the two data frames based on the "Gene" column
merged_df <- merge(subsetichocna, qdnaseq_df, by = "Gene")


###### CALLS FOR DNACOPY WITH THE SELECTED GENES

allcalls <- read_xlsx("/Users/mateopava/Desktop/Luisa_DATOS/logR_qdnaseqVsichorCNA_withcall.xlsx")

# Identify columns with "calls" in their name
call_cols <- grepl("call", colnames(combined_df))

# Extract those columns
copy_calls <- combined_df[, call_cols]

# Handle duplicate column names
duplicated_cols <- duplicated(colnames(copy_calls))
if(any(duplicated_cols)){
  # If there are duplicated columns, rename them with a suffix
  new_colnames <- colnames(copy_calls)
  new_colnames[duplicated_cols] <- paste0(new_colnames[duplicated_cols], "_dup")
  colnames(copy_calls) <- new_colnames
}

# Iterate through all columns and replace 'amplification' with 'gain' and 'deletion' with 'loss'
copy_calls[] <- apply(copy_calls, 2, function(x) {
  ifelse(x == 'amplification', 'gain', ifelse(x == 'deletion', 'loss', x))
})

allcopycalls<- cbind(allcalls, copy_calls) #gives me all copycalls

# Extract all unique PRO patterns
pro_patterns_copy <- unique(gsub("(PRO\\d+).*", "\\1", colnames(allcopycalls)))

# For each PRO pattern, get the matching columns and combine
ordered_columns <- unlist(sapply(pro_patterns_copy, function(pattern) {
  grep(pattern, colnames(allcopycalls), value = TRUE)
}))

# Rearrange columns based on the ordered columns
allcopycalls2 <- allcopycalls[, ordered_columns]

# try to fix to only ones we have
# Extract PRO patterns from matched_ichor_columns
pro_patterns_to_keep <- unique(gsub("(PRO\\d+).*", "\\1", matched_ichor_columns))

# Identify columns in allcopycalls2 that match these patterns
matching_cols <- grep(paste(pro_patterns_to_keep, collapse="|"), colnames(allcopycalls2), value = TRUE)

# Subset the data frame to keep only the matching columns
filtered_allcopycalls2 <- allcopycalls2[, matching_cols]

filtered_allcopycalls2[filtered_allcopycalls2 == "NA"] <- NA
filtered_allcopycalls2 <- cbind(allcopycalls[, 1], filtered_allcopycalls2)
colnames(filtered_allcopycalls2)[1] <- "Gene"

###### DEAL WITH DUPLICATION
# Create a dataframe to store duplicated columns
duplicated_df <- data.frame()

# Extract the pattern from the column names
pattern <- "(PRO\\d+)_(.*?)_DNA"
replacement <- "\\1_\\2_call"

# Loop over the column names to check for duplicates
new_colnames <- character(length = ncol(filtered_allcopycalls2))
for (i in seq_along(filtered_allcopycalls2)) {
  colname <- names(filtered_allcopycalls2)[i]
  new_colname <- gsub(pattern, replacement, colname)
  
  # Check if the new column name is a duplicate
  if (new_colname %in% new_colnames) {
    duplicated_df[[colname]] <- filtered_allcopycalls2[[colname]]
  } else {
    new_colnames[i] <- new_colname
  }
}

# Update column names
names(filtered_allcopycalls2) <- new_colnames


# make counts

call_counts <- colSums(!is.na(filtered_allcopycalls2[,-1]))

unique_calls_per_gene <- apply(filtered_allcopycalls2, 1, function(row) {
  unique(na.omit(row))
}) # this is unneccesary



