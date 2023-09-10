#READ ALL PANEL FILES
path_to_panel <- "/Users/mateopava/Desktop/CNA_Mateo_MSC/allpanel"

# Get a list of all tsv files in the directory
panel_files <- list.files(path_to_panel, pattern = "\\.tsv$")

# Use lapply to read all csv files into a list of data frames
list_of_panel <- lapply(paste0(path_to_panel, "/", panel_files), function(x) read.csv(x, header=TRUE, sep="\t"))
list_of_panel <- setNames(list_of_panel, panel_files)
names(list_of_panel) <- gsub("_cna.tsv", "", names(list_of_panel))
panel_names <- names(list_of_panel)

# COUNT HOW MANY PANELS HAVE NO SIGNIFICANT INFORMATION
sum(sapply(list_of_panel, nrow) == 1) #95
# ONLY 35 HAVE SIGNIFICANT CALLS

# make the genes as rows, not separated by commas
list_of_panel_separated <- lapply(list_of_panel, function(df) {
  df %>% 
    separate_rows(gene, sep = ",")
})

# now get rid of repeating genes
list_of_panel_unique <- lapply(list_of_panel_separated, function(df) {
  df %>% 
    distinct(gene, .keep_all = TRUE)
})

#now only significant panels
list_of_panel_unique <- list_of_panel_unique[sapply(list_of_panel_unique, nrow) > 1]

#rename the 6th column having the cnv info for easier access
list_of_panel_unique <- lapply(list_of_panel_unique, function(df) {
  names(df)[6] <- "CNV"
  return(df)
})

# rearranging the list_of_dataframes #chech script shallows files to see any problems
# Use gsub() to replace ".gene_calls.tsv" with ""
wgs_names <- gsub(".gene_calls.tsv", "", df_names)

# Assign the new names back to the list of dataframes
names(list_of_dataframes) <- wgs_names


#function to see comparison #THIS WORKS BUT I CANT SEEM TO USE IT
list_of_common_genes <- list()
list_of_common_cnv <- list()

for (i in 1:nrow(correspond_panel_s)) {
  panel_name <- correspond_panel_s$P300_fileID[i]
  wgs_name <- correspond_panel_s$sWGS_fileID[i]
  
  # Make sure both names exist in their respective lists
  if (panel_name %in% names(list_of_panel_unique) && wgs_name %in% names(list_of_dataframes)) {
    df1 <- list_of_panel_unique[[panel_name]]
    df2 <- list_of_dataframes[[wgs_name]]
    
    # common 'gene' values
    common_genes <- intersect(df1$gene, df2$gene)
    list_of_common_genes[[paste(panel_name, wgs_name, sep = "_")]] <- common_genes
    
    # common 'cnv' values
    # replace 'cnv' with the actual column name in df2 that represents 'cnv'
    common_cnv <- intersect(df1$CNV, df2$cn_ichorCNA)
    list_of_common_cnv[[paste(panel_name, wgs_name, sep = "_")]] <- common_cnv
  }
}

# NOW CREATE A TABLE LIKE THE ONE FOR WES
# Initialize a list to store merged dataframes
list_of_merged_data <- list()

for (i in 1:nrow(correspond_panel_s)) {
  panel_name <- correspond_panel_s$P300_fileID[i]
  wgs_name <- correspond_panel_s$sWGS_fileID[i]
  
  # Make sure both names exist in their respective lists
  if (panel_name %in% names(list_of_panel_unique) && wgs_name %in% names(list_of_dataframes)) {
    df1 <- list_of_panel_unique[[panel_name]]
    df2 <- list_of_dataframes[[wgs_name]]
    
    # Merge dataframes on 'gene'
    merged_data <- merge(df1, df2, by = "gene")
    list_of_merged_data[[paste(panel_name, wgs_name, sep = "_")]] <- merged_data
  }
}

### see if this works/ IT WORKS

library(reshape2)

# Generate plots for each merged dataframe
list_of_plots <- lapply(list_of_merged_data, function(df) {
  
  # Reshape to long format
  df_long <- melt(df, id.vars = "gene", measure.vars = c("CNV", "cn_ichorCNA"))
  
  # Generate plot
  p <- ggplot(df_long, aes(x = gene, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Gene", y = "Value", title = "CNV and cn_ichorCNA by Gene", fill = "Source") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(fill = guide_legend(title = "Source"))
  
  return(p)
})

# Print a specific plot as an example
print(list_of_plots[[1]])

# correlation see if it works/ IT WORKS
# Generate correlation coefficients for each merged dataframe
list_of_correlations <- lapply(list_of_merged_data, function(df) {
  
  # Calculate correlation
  cor_coeff <- cor(df$CNV, df$cn_ichorCNA, method = "pearson", use = "pairwise.complete.obs")
  
  return(cor_coeff)
})

### NEXT TRY
# Create an empty data frame to store the results
results <- data.frame()

# Loop through each dataframe in the list
for (gene_name in names(list_of_merged_data)) {
  df <- list_of_merged_data[[gene_name]]
  
  # Ensure the dataframe has both sources (panel and WGS)
  if(length(unique(df$source)) == 2){
    # Apply t-test
    t_test_result <- t.test(df$CNV ~ df$source)
    
    # Store the results
    results <- rbind(results, data.frame(
      gene = gene_name,
      p_value = t_test_result$p.value
    ))
    
    # Uncomment the following lines if you want to use Mann-Whitney U test
    # wilcox_test_result <- wilcox.test(df$CNV ~ df$source)
    # results <- rbind(results, data.frame(
    #   gene = gene_name,
    #   wilcox_p_value = wilcox_test_result$p.value
    # ))
  }
}

# Print the results
print(results)

# make table with only shallow and panel
correspond_panel_s <- correspond[!is.na(correspond$P300_fileID) & !is.na(correspond$sWGS_fileID), ]

