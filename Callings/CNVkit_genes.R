#CNVKIT WES
setwd("/Users/mateopava/Desktop/example/cnvkit")
cnvkit_files <- list.files()

# get all the files to read them in a df
cnvkit <- cnvkit_files[grepl("\\.call\\.cns$", cnvkit_files)]

#list it
list_of_cnvkit_dfs <- lapply(cnvkit, function(file) {
  df <- read.table(file, header = TRUE, sep = "\t") # assuming they are tab-delimited
  
  # Extract the unique identifier from the filename
  identifier <- gsub("^(.*)(\\.call\\.cns)$", "\\1", file)
  
  # Add a new column to store the identifier
  df$Origin <- identifier
  return(df)
})

# Filter genes
swgs_ichorcna <- read_xlsx("/Users/mateopava/Desktop/Luisa_DATOS/logR_qdnaseqVsichorCNA_withcall.xlsx")
genes<- swgs_ichorcna$Gene
# Define the genes of interest
genes_of_interest <- genes # the genes list is obtained from another script

# Function to expand, filter the genes, and then get unique rows
filter_genes_expanded <- function(df) {
  # Expand the genes
  expanded_df <- df[rep(seq_len(nrow(df)), times = sapply(strsplit(df$gene, ","), length)), ]
  expanded_df$gene <- unlist(strsplit(df$gene, ","))
  
  # Filter based on genes of interest
  filtered_df <- expanded_df[expanded_df$gene %in% genes_of_interest, ]
  
  # Get unique rows
  unique_df <- unique(filtered_df)
  
  return(unique_df)
}

# Apply the function to each dataframe in the list
cnvkit_annotated <- lapply(list_of_cnvkit_dfs, filter_genes_expanded)

# assuming the calls
cnvkit_annotated <- lapply(cnvkit_annotated, function(df) {
  
  # Creating call column based on cn and chromosome
  df$call <- ifelse(df$chromosome == "chrX",
                    ifelse(df$cn == 1, 'neutral', ifelse(df$cn < 1, 'loss', 'gain')),
                    ifelse(df$cn == 2, 'neutral', ifelse(df$cn < 2, 'loss', 'gain'))
  )
  
  return(df)
})

# create the call df for all
# Define the reshaping function
reshape_fun <- function(df) {
  df %>%
    dplyr::select(gene, Origin, call) %>%
    dplyr::group_by(gene, Origin) %>%
    dplyr::summarise(call = paste(unique(call), collapse = "/")) %>%
    tidyr::spread(key = Origin, value = call) %>%
    dplyr::right_join(data.frame(gene = genes_of_interest), by = "gene") %>%
    dplyr::arrange(match(gene, genes_of_interest))
}


# Apply the reshaping function to each df in the list and merge them together
reshaped_list <- lapply(cnvkit_annotated, reshape_fun)
reshaped_df <- Reduce(function(x, y) dplyr::left_join(x, y, by = "gene"), reshaped_list)

reshaped_df[-1] <- lapply(reshaped_df[-1], function(col) {
  ifelse(col == "neutral", NA, col)
})



####### TECH COMPARISON ########

# Rename columns in reshaped_df
cols_to_rename <- colnames(reshaped_df)[-1] # Exclude the gene column
new_colnames <- paste("cnvkit", cols_to_rename, sep = "_")
names(reshaped_df)[names(reshaped_df) %in% cols_to_rename] <- new_colnames

# Rename columns in reshaped_ascat_df
cols_to_rename_ascat <- colnames(reshaped_ascat_df)[-1] # Exclude the genes column
new_colnames_ascat <- paste("ascat", cols_to_rename_ascat, sep = "_")
names(reshaped_ascat_df)[names(reshaped_ascat_df) %in% cols_to_rename_ascat] <- new_colnames_ascat

# Join the dataframes
WES_techs <- dplyr::left_join(reshaped_df, reshaped_ascat_df, by = c("gene" = "genes"))

# Extract unique sample names (excluding prefixes and 'gene' column)
sample_names <- unique(gsub("^cnvkit_|^ascat_", "", colnames(WES_techs)[-1]))

# Create the ordered column sequence
ordered_cols <- c("gene", unlist(lapply(sample_names, function(sample) {
  c(paste("cnvkit", sample, sep = "_"), paste("ascat", sample, sep = "_"))
})))

# Subset the dataframe by the ordered columns
WES_techs <- WES_techs[, ordered_cols] # well done
WES_techs <- as.data.frame(WES_techs) #change from ttbl to df

total_ascat_samples <- ncol(WES_techs[grep("^ascat_", colnames(WES_techs))])
total_cnvkit_samples <- ncol(WES_techs[grep("^cnvkit_", colnames(WES_techs))])

ascat_gains <- sum(WES_techs[grep("^ascat_", colnames(WES_techs))] == "gain", na.rm = TRUE)
ascat_losses <- sum(WES_techs[grep("^ascat_", colnames(WES_techs))] == "loss", na.rm = TRUE)

cnvkit_gains <- sum(WES_techs[grep("^cnvkit_", colnames(WES_techs))] == "gain", na.rm = TRUE)
cnvkit_losses <- sum(WES_techs[grep("^cnvkit_", colnames(WES_techs))] == "loss", na.rm = TRUE)

genes_to_check <- c("BRCA1","BRCA2", "AR", "RAD51", "ATM")

gene_counts <- lapply(genes_to_check, function(gene) {
  gene_row <- WES_techs[WES_techs$gene == gene, ]
  
  ascat_gain_count <- sum(gene_row[grep("^ascat_", colnames(gene_row))] == "gain", na.rm = TRUE)
  ascat_loss_count <- sum(gene_row[grep("^ascat_", colnames(gene_row))] == "loss", na.rm = TRUE)
  total_ascat_for_gene <- sum(!is.na(gene_row[grep("^ascat_", colnames(gene_row))]))
  
  cnvkit_gain_count <- sum(gene_row[grep("^cnvkit_", colnames(gene_row))] == "gain", na.rm = TRUE)
  cnvkit_loss_count <- sum(gene_row[grep("^cnvkit_", colnames(gene_row))] == "loss", na.rm = TRUE)
  total_cnvkit_for_gene <- sum(!is.na(gene_row[grep("^cnvkit_", colnames(gene_row))]))
  
  return(list(
    gene = gene,
    ascat_gain = ascat_gain_count,
    ascat_loss = ascat_loss_count,
    total_ascat = total_ascat_for_gene,
    cnvkit_gain = cnvkit_gain_count,
    cnvkit_loss = cnvkit_loss_count,
    total_cnvkit = total_cnvkit_for_gene
  ))
})

gene_counts_df <- do.call(rbind.data.frame, gene_counts) # total means the genes that were called as gain/loss

########################### above works
library(dplyr)
library(tidyr)

# Extracting column names
wes_cols <- colnames(WES_techs)

# Get columns related to each technology
cnvkit_cols <- wes_cols[grep("^cnvkit_", wes_cols)]
ascat_cols <- wes_cols[grep("^ascat_", wes_cols)]

# Create a function to compute the metrics for one sample
compare_calls <- function(cnvkit_col, ascat_col) {
  df <- dplyr::ungroup(WES_techs) %>%
    dplyr::select(gene, cnvkit = cnvkit_col, ascat = ascat_col)
  
  both_gain = sum(df$cnvkit == 'gain' & df$ascat == 'gain', na.rm = TRUE)
  both_loss = sum(df$cnvkit == 'loss' & df$ascat == 'loss', na.rm = TRUE)
  cnvkit_gain_ascat_loss = sum(df$cnvkit == 'gain' & df$ascat == 'loss', na.rm = TRUE)
  cnvkit_loss_ascat_gain = sum(df$cnvkit == 'loss' & df$ascat == 'gain', na.rm = TRUE)
  cnvkit_call_ascat_na = sum((df$cnvkit %in% c('gain', 'loss')) & is.na(df$ascat))
  ascat_call_cnvkit_na = sum((df$ascat %in% c('gain', 'loss')) & is.na(df$cnvkit))
  
  return(c(both_gain, both_loss, cnvkit_gain_ascat_loss, cnvkit_loss_ascat_gain, cnvkit_call_ascat_na, ascat_call_cnvkit_na))
}


# Apply the function to all samples
results <- mapply(compare_calls, cnvkit_cols, ascat_cols)

# Convert to a dataframe
results_df <- as.data.frame(t(results))
colnames(results_df) <- c("both_gain", "both_loss", "cnvkit_gain_ascat_loss", "cnvkit_loss_ascat_gain", "cnvkit_call_ascat_na", "ascat_call_cnvkit_na")
results_df$sample <- cnvkit_cols

library(dplyr)
library(tidyr)
library(pheatmap)

# Convert 'gain', 'loss', and 'NA' to numeric values for heatmap
heatmap_data <- WES_techs %>%
  mutate(across(starts_with("cnvkit"), ~case_when(
    . == "gain" ~ 1,
    . == "loss" ~ -1,
    TRUE ~ 0
  ))) %>%
  mutate(across(starts_with("ascat"), ~case_when(
    . == "gain" ~ 1,
    . == "loss" ~ -1,
    TRUE ~ 0
  )))

heatmap_matrix <- as.matrix(heatmap_data)
rownames(heatmap_matrix) <- heatmap_data$gene
heatmap_matrix <- heatmap_matrix[,-1]

# Define colors for heatmap
colors <- colorRampPalette(c("red", "white", "green"))(3)

# Plot heatmap
# Plot heatmap with rotated x-axis labels and adjusted font size
pheatmap(heatmap_matrix, 
         color = colors, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "none", 
         fontsize_row = 8, 
         fontsize_col = 5, # Adjust this value as needed
         angle_col = 90)   # Rotate x-axis labels to be vertical


########## POSSIBLE
ascat_cols <- grep("ascat", names(WES_techs_inferred), value = TRUE)
ascat_data <- unlist(WES_techs_inferred[ascat_cols])
ascat_data[is.na(ascat_data)] <- "neutral"
ascat_counts <- table(ascat_data)

cnvkit_cols <- grep("cnvkit", names(WES_techs_inferred), value = TRUE)
cnvkit_data <- unlist(WES_techs_inferred[cnvkit_cols])
cnvkit_data[is.na(cnvkit_data)] <- "neutral"
cnvkit_counts <- table(cnvkit_data)


################ WORKS
genes_to_check <- c("BRCA1", "BRCA2", "AR", "RAD51", "ATM")

# Filter the dataframe for the specified genes
filtered_df <- WES_techs_inferred[WES_techs_inferred$gene %in% genes_to_check, ]

# Initialize lists to store the counts for each gene
ascat_gene_counts <- list()
cnvkit_gene_counts <- list()

for (gene in genes_to_check) {
  # Extract data for the current gene
  gene_data <- filtered_df[filtered_df$gene == gene, ]
  
  # Count ASCAT Calls
  ascat_cols <- grep("ascat", names(gene_data), value = TRUE)
  ascat_data <- unlist(gene_data[ascat_cols])
  ascat_data[is.na(ascat_data)] <- "neutral"
  ascat_gene_counts[[gene]] <- table(ascat_data)
  
  # Count CNVkit Calls
  cnvkit_cols <- grep("cnvkit", names(gene_data), value = TRUE)
  cnvkit_data <- unlist(gene_data[cnvkit_cols])
  cnvkit_data[is.na(cnvkit_data)] <- "neutral"
  cnvkit_gene_counts[[gene]] <- table(cnvkit_data)
}

list(ASCAT = ascat_gene_counts, CNVkit = cnvkit_gene_counts)

# Remove neutral counts from the lists
ascat_gene_counts <- lapply(ascat_gene_counts, function(x) x[names(x) != "neutral"])
cnvkit_gene_counts <- lapply(cnvkit_gene_counts, function(x) x[names(x) != "neutral"])

# Combine ASCAT and CNVkit counts into one dataframe for ggplot
combined_list <- list(ASCAT = ascat_gene_counts, CNVkit = cnvkit_gene_counts)
combined_df <- melt(combined_list, id.vars = c("gene", "call"))

######3
# Reshaping the dataframe
ascat_melted <- data.frame(
  gene = combined_df$L2,
  technology = "ASCAT",
  call = combined_df$ascat_data,
  count = combined_df$value
)

cnvkit_melted <- data.frame(
  gene = combined_df$L2,
  technology = "CNVkit",
  call = combined_df$cnvkit_data,
  count = combined_df$value
)

final_df <- rbind(ascat_melted, cnvkit_melted)
final_df <- final_df[!is.na(final_df$call), ]
final_df <- subset(final_df, !(call %in% c("loss/gain", "neutral/loss")))

# Plotting with explicit colors and side-by-side bars
colors <- c("gain" = "green", "loss" = "red")

ggplot(final_df, aes(x = gene, y = count, fill = call, group = interaction(technology, call))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors) + 
  labs(title = "Comparison of Calls between ASCAT and CNVkit excluding Neutral", 
       x = "Gene", 
       y = "Count", 
       fill = "Call Type") +
  facet_grid(. ~ technology)

# Convert the list counts to matrices for correlation computation
# Filter out unwanted call labels from the list counts
ascat_gene_counts <- lapply(ascat_gene_counts, function(x) x[names(x) %in% c("gain", "loss")])
cnvkit_gene_counts <- lapply(cnvkit_gene_counts, function(x) x[names(x) %in% c("gain", "loss")])

# Convert the filtered list counts to matrices for correlation computation
ascat_matrix <- do.call(rbind, ascat_gene_counts)
cnvkit_matrix <- do.call(rbind, cnvkit_gene_counts)

# Compute the correlation
cor(as.vector(ascat_matrix), as.vector(cnvkit_matrix), method = "spearman")

genes <- rownames(ascat_matrix)
correlations <- sapply(genes, function(gene) {
  row_ascat <- ascat_matrix[gene, , drop = TRUE]
  row_cnvkit <- cnvkit_matrix[gene, , drop = TRUE]
  
  # Ensure vectors have the same length before computing correlation
  if(length(row_ascat) == length(row_cnvkit)) {
    return(cor(row_ascat, row_cnvkit, method = "spearman"))
  } else {
    return(NA)
  }
})

names(correlations) <- genes
correlations


######

# Function to process each dataframe
process_df <- function(df) {
  # Merge the dataframe with the compare_tech dataframe based on matching columns
  merged_df <- merge(df, compare_tech, by.x = "Origin", by.y = "WES_fileID", all.x = TRUE)
  
  # Apply the formula to calculate copynumber
  merged_df$infer_cn <- (2^(merged_df$log2 + 1) - 2 + 2 * merged_df$HandE) / merged_df$HandE
  
  # Reorder the columns to place infer_cn at column 7
  column_order <- c(names(df)[1:6], "infer_cn", names(df)[7:length(names(df))])
  return(merged_df[, column_order])
}

# Apply the function to all dataframes in the cnvkit_annotated list
final_list <- lapply(cnvkit_annotated, process_df)

# Check the first dataframe in the final_list
head(final_list[[1]])

## corrplot package
cor_matrix <- cor(compare_tech[, c("HandE", "TumorFractionIchorCNA", "TumorFractionASCAT")], use = "pairwise.complete.obs")
# Compute p-values for correlations
cor_pvalues <- cor.mtest(compare_tech[, c("HandE", "TumorFractionIchorCNA", "TumorFractionASCAT")])$p

corrplot(cor_matrix, method = "circle")

corrplot(cor_matrix, method = "circle", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 45, 
         addCoef.col = "black",
         p.mat = cor_pvalues, 
         sig.level = 0.05, 
         insig = "blank", 
         diag = FALSE)

### ttest summary
# Extract relevant values from t-test objects
t_values <- c(t_test_ASCAT_Ichor$statistic, 
              t_test_ASCAT_HandE$statistic, 
              t_test_Ichor_HandE$statistic)

df_values <- c(t_test_ASCAT_Ichor$parameter, 
               t_test_ASCAT_HandE$parameter, 
               t_test_Ichor_HandE$parameter)

p_values <- c(t_test_ASCAT_Ichor$p.value, 
              t_test_ASCAT_HandE$p.value, 
              t_test_Ichor_HandE$p.value)

lower_ci <- c(t_test_ASCAT_Ichor$conf.int[1], 
              t_test_ASCAT_HandE$conf.int[1], 
              t_test_Ichor_HandE$conf.int[1])

upper_ci <- c(t_test_ASCAT_Ichor$conf.int[2], 
              t_test_ASCAT_HandE$conf.int[2], 
              t_test_Ichor_HandE$conf.int[2])

mean_x <- c(t_test_ASCAT_Ichor$estimate[1], 
            t_test_ASCAT_HandE$estimate[1], 
            t_test_Ichor_HandE$estimate[1])

mean_y <- c(t_test_ASCAT_Ichor$estimate[2], 
            t_test_ASCAT_HandE$estimate[2], 
            t_test_Ichor_HandE$estimate[2])

# Organize the values into a data frame
ttest_summary <- data.frame(
  Comparison = c("ASCAT vs Ichor", "ASCAT vs HandE", "Ichor vs HandE"),
  t_value = t_values,
  df = df_values,
  p_value = p_values,
  lower_ci = lower_ci,
  upper_ci = upper_ci,
  mean_x = mean_x,
  mean_y = mean_y
)

ttest_summary




