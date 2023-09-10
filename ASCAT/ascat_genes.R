# GENE annotate for ASCAT
setwd("/Users/mateopava/Desktop/ascat/ascat_analyze")
all_files <- list.files()

# get all the files to read them in a df
ascat_files <- all_files[grepl(".cnvs.txt$", all_files)]

#list it
list_of_ascat_dfs <- lapply(ascat_files, function(file) {
  df <- read.table(file, header = TRUE, sep = "\t") # assuming they are tab-delimited
  
  # Extract the unique identifier from the filename
  identifier <- gsub("^(.*)\\_vs.*\\.cnvs\\.txt$", "\\1", file)
  
  # Add a new column to store the identifier
  df$Origin <- identifier
  
  return(df)
})

# ascat_df <- do.call(rbind, list_of_ascat_dfs) not neccesary

# read
genesanalyze<- read.csv("/Users/mateopava/Desktop/genes.csv", sep = "\t")

# ANNOTATION

# Read the genes file
genesanalyze <- read.csv("/Users/mateopava/Desktop/genes.csv", sep = "\t")
# Remove the "chr" prefix from the chromosome column (assuming the column name is 'chromosome')
genesanalyze$chromosome_name <- gsub("chr", "", genesanalyze$chromosome_name)


apply_procedure <- function(df, genes) {
  
  df$genes <- NA  # initialize an empty vector
  
  matching_genes <- vector("list", nrow(df))  # create a list to store matching genes for each row
  
  for (j in 1:nrow(df)) {
    # Find matching genes
    matching_genes[[j]] <- genesanalyze$hgnc_symbol[
      between(genesanalyze[, 3], df[j, 2], df[j, 3]) &
        genesanalyze[, 5] == df[j, 1] # assuming column 1 in df is chromosome and column 5 in genes is chromosome
    ]
  }
  
  # If there are matching genes, paste them together separated by comma
  df$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
  
  df <- na.omit(df)  # remove rows with NA in 'genes' column
  
  # Split genes into separate rows
  df <- separate_rows(df, genes, sep = ",")
  
  return(df)
}

# Apply the procedure to each data frame in the list_of_ascat_dfs using genesanalyze
list_of_annotated_dfs <- lapply(list_of_ascat_dfs, apply_procedure, genes = genesanalyze)

# Apply the call
ASCAT_annotated <- lapply(list_of_annotated_dfs, function(df) {
  # Summing up nMajor and nMinor
  df$nSum <- df$nMajor + df$nMinor
  
  # Creating call column based on nSum and chromosome
  df$call <- ifelse(df$chr == "X", 
                    ifelse(df$nSum == 1, "neutral", ifelse(df$nSum < 1, "loss", "gain")), 
                    ifelse(df$nSum == 2, "neutral", ifelse(df$nSum < 2, "loss", "gain"))
  )
  
  # Creating InferredCall column based on nSum and chromosome
  df$AbsCall <- ifelse(df$chr == "X", 
                       ifelse(df$nSum == 1, "neutral",
                              ifelse(df$nSum > 1, "gain", "loss")),
                       ifelse(df$nSum == 0, "homdel",
                              ifelse(df$nSum == 1, "hetdel",
                                     ifelse(df$nSum == 2, "neutral",
                                            ifelse(df$nSum == 3, "gain",
                                                   ifelse(df$nSum == 4, "amp",
                                                          ifelse(df$nSum == 5, "amp",
                                                                 ifelse(df$nSum >= 6, "hlamp", NA))))))))
  return(df)
})

# create the calls df
reshape_fun_ascat <- function(df) {
  df %>%
    dplyr::select(genes, Origin, call) %>%
    dplyr::group_by(genes, Origin) %>%
    dplyr::summarise(call = paste(unique(call), collapse = "/")) %>%
    tidyr::spread(key = Origin, value = call) %>%
    dplyr::right_join(data.frame(genes = genes_of_interest), by = "genes") %>%
    dplyr::arrange(match(genes, genes_of_interest))
}

# Apply the reshaping function to each df in the ASCAT_annotated list and merge them together
reshaped_ascat_list <- lapply(ASCAT_annotated, reshape_fun_ascat)
reshaped_ascat_df <- Reduce(function(x, y) dplyr::left_join(x, y, by = "genes"), reshaped_ascat_list)

reshaped_ascat_df[-1] <- lapply(reshaped_ascat_df[-1], function(col) {
  ifelse(col == "neutral", NA, col)
})

# seeing if the logR are needed
#logR001<- read.csv("/Users/mateopava/Desktop/PRO001_18-27110_DNAFFPE_NA_C_HN00131668_vs_PRO001_NA_DNASAL_NA_C_HN00131668.tumour_tumourLogR.txt", sep = "\t")


# have to consider that the call is being assumed as the sum of the cn of the alleles, nMinor and nMajor without taking into account ploidy or purity
cols_to_rename_ascat <- colnames(reshaped_ascat_df)[-1] # Exclude the genes column
new_colnames_ascat <- paste("ascat", cols_to_rename_ascat, sep = "_")
names(reshaped_ascat_df)[names(reshaped_ascat_df) %in% cols_to_rename_ascat] <- new_colnames_ascat


#Expected neutral copy number=Ploidy×Tumor Purity+2×(1−Tumor Purity)

#AFTER PURITY PLOIDY

reshape_fun_ascat_inferred <- function(df) {
  df %>%
    dplyr::select(genes, Origin, InferredCall) %>%
    dplyr::group_by(genes, Origin) %>%
    dplyr::summarise(call = paste(unique(InferredCall), collapse = "/")) %>%
    tidyr::spread(key = Origin, value = call) %>%
    dplyr::right_join(data.frame(genes = genes_of_interest), by = "genes") %>%
    dplyr::arrange(match(genes, genes_of_interest))
}

# Apply the reshaping function to each df in the ASCAT_annotated list based on InferredCall and merge them together
reshaped_ascat_inferred_list <- lapply(ASCAT_annotated, reshape_fun_ascat_inferred)
reshaped_ascat_inferred_df <- Reduce(function(x, y) dplyr::left_join(x, y, by = "genes"), reshaped_ascat_inferred_list)

reshaped_ascat_inferred_df[-1] <- lapply(reshaped_ascat_inferred_df[-1], function(col) {
  ifelse(col == "neutral", NA, col)
})

cols_to_rename_ascat_ <- colnames(reshaped_ascat_inferred_df)[-1] # Exclude the genes column
new_colnames_ascat_ <- paste("ascat", cols_to_rename_ascat_, sep = "_")
names(reshaped_ascat_inferred_df)[names(reshaped_ascat_inferred_df) %in% cols_to_rename_ascat_] <- new_colnames_ascat_


WES_techs_inferred <- dplyr::left_join(reshaped_df, reshaped_ascat_inferred_df, by = c("gene" = "genes"))

# Extract unique sample names (excluding prefixes and 'gene' column)
sample_names <- unique(gsub("^cnvkit_|^ascat_", "", colnames(WES_techs_inferred)[-1]))

# Create the ordered column sequence
ordered_cols <- c("gene", unlist(lapply(sample_names, function(sample) {
  c(paste("cnvkit", sample, sep = "_"), paste("ascat", sample, sep = "_"))
})))

# Subset the dataframe by the ordered columns
WES_techs_inferred <- WES_techs_inferred[, ordered_cols] # well done
WES_techs_inferred <- as.data.frame(WES_techs_inferred) #change from ttbl to df

