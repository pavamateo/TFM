# First read wes_wgsfinal.csv

wes_wgsfinal <- read.csv("/Users/mateopava/wes_wgsfinal.csv", sep = "\t")
write.table(heatmapwes, file = "heatmapall.csv", row.names = F, col.names = T, quote = F, sep = "\t")


# Merge the two dataframes based on Gene symbol
heatmapwes <- wes_wgsfinal %>%
  left_join(genesanalyze, by = c("Gene" = "hgnc_symbol"))

create_heatmap <- function(df) {
  # Mapping between call types and numeric values
  call_mapping <- c("DEEPDEL" = 1, "LOSS" = 2, "NEUTRAL" = 3, "GAIN" = 4, "AMP" = 5)
  
  # Convert the call column
  df$call_numeric <- as.numeric(call_mapping[df$call])
  
  # Create the heatmap
  ggplot(df, aes(x = technology, y = Gene, fill = call_numeric)) + 
    geom_tile() +
    scale_fill_gradient2(
      name = "Calls",
      low = "blue",
      mid = "lightgreen",
      high = "red",
      midpoint = 3,
      breaks = c(1, 3, 5),
      labels = c("DEEPDEL", "NEUTRAL", "AMP")
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} # should show genes/ WORKS


# Convert the data from wide to long format
heatmap_data <- heatmapwes %>%
  select(PRO_id, Gene, chromosome_name, cnvkit_call, ichorcna_call, qdnaseq_call, panel_call, DNAcopy_call, ASCAT_call) %>%
  gather(technology, call, -PRO_id, -Gene, -chromosome_name)

# Filter the data to include only the technologies you are interested in
filtered_data <- heatmap_data %>%
  filter(technology %in% c("cnvkit_call","ASCAT_call", "DNAcopy_call", "ichorcna_call", "qdnaseq_call", "panel_call"))

# Keep only the first occurrence of each combination of PRO_id, Gene, and technology
distinct_data <- filtered_data %>%
  distinct(PRO_id, Gene, technology, .keep_all = TRUE)

# Convert to long format
long_datas <- distinct_data %>%
  select(Gene, PRO_id, technology, call)

long_datas$sample_technology <- paste0(long_datas$PRO_id, "-", long_datas$technology) # this is for ascat

oncoprint_matrix <- long_datas %>%
  select(Gene, sample_technology, call) %>%
  spread(key = sample_technology, value = call) #ASCAT

rownames(oncoprint_matrix) <- oncoprint_matrix$Gene
oncoprint_matrix$Gene <- NULL  # remove the Gene column since it's now the row names

oncoprint_matrix <- as.matrix(oncoprint_matrix) 
oncoprint_matrix[is.na(oncoprint_matrix)] <- ""

# List of technologies
techs <- c("cnvkit_call", "DNAcopy_call","ASCAT_call", "ichorcna_call", "qdnaseq_call","panel_call")

# Create a list of data frames for each technology
dfs_by_tech <- lapply(techs, function(tech) {
  oncoprint_matrix[, grepl(tech, colnames(oncoprint_matrix))]
})

# Convert the list to a named list for easy access
names(dfs_by_tech) <- techs

# You can access each data frame using the technology name
cnv_onco <- dfs_by_tech$cnvkit_call
dna_onco <- dfs_by_tech$DNAcopy_call
qdna_onco <- dfs_by_tech$qdnaseq_call
ichor_onco <- dfs_by_tech$ichorcna_call
panel_onco <- dfs_by_tech$panel_call
ascat_onco <- dfs_by_tech$ASCAT_call

# Melt the data for each technology using R's reshape2 package
cnv_onco_melted <- melt(cnv_onco, variable.name = "Sample", value.name = "Alteration")
dna_onco_melted <- melt(dna_onco, variable.name = "Sample", value.name = "Alteration")
ichor_onco_melted <- melt(ichor_onco, variable.name = "Sample", value.name = "Alteration")
panel_onco_melted <- melt(panel_onco, variable.name = "Sample", value.name = "Alteration")
qdna_onco_melted <- melt(qdna_onco, variable.name = "Sample", value.name = "Alteration")
ascat_onco_melted <- melt(ascat_onco, variable.name = "Sample", value.name = "Alteration")

# Combine all melted data
combined_melted <- rbind(cnv_onco_melted, dna_onco_melted,ascat_onco_melted, ichor_onco_melted, qdna_onco_melted, panel_onco_melted) #works fine/ original

### combined_melted y ponerlo cada 48 samples
combined_melted$Tech <- gsub("^.+-(.+)$", "\\1", combined_melted$Var2)
# Ensure "panel" appears last in the Tech levels
tech_levels <- setdiff(unique(combined_melted$Tech), "panel")
tech_levels <- c(tech_levels, "panel")
combined_melted$Tech <- factor(combined_melted$Tech, levels = tech_levels)


# Extract the base sample names
combined_melted$BaseSample <- gsub("(.+?)(-.+)$", "\\1", combined_melted$Var2)

combined_melted <- combined_melted[order(combined_melted$BaseSample, combined_melted$Var1, combined_melted$Tech), ]

combined_melted$BaseSample <- factor(combined_melted$BaseSample, levels = unique(combined_melted$BaseSample))
combined_melted$Tech <- factor(combined_melted$Tech, levels = unique(combined_melted$Tech))

combined_melted$Var2 <- factor(combined_melted$Var2, levels = unique(combined_melted$Var2))

# Plot works but shows all the techs
ggplot(combined_melted, aes(x = Var2 , y = Var1, fill = Alteration)) +
  geom_tile() +
  scale_fill_manual(values = c("NEUTRAL" = "lightgreen", 
                               "GAIN" = "pink", 
                               "LOSS" = "lightblue1", 
                               "AMP" = "red", 
                               "DEEPDEL" = "blue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##### 4 plots

# Split the data into four groups
unique_sample <- unique(combined_melted$BaseSample)
quarter_length <- length(unique_sample) %/% 4

group1_samples <- unique_sample[1:quarter_length]
group2_samples <- unique_sample[(quarter_length + 1):(2*quarter_length)]
group3_samples <- unique_sample[(2*quarter_length + 1):(3*quarter_length)]
group4_samples <- unique_sample[(3*quarter_length + 1):length(unique_sample)]

group1_data <- combined_melted[combined_melted$BaseSample %in% group1_samples, ]
group2_data <- combined_melted[combined_melted$BaseSample %in% group2_samples, ]
group3_data <- combined_melted[combined_melted$BaseSample %in% group3_samples, ]
group4_data <- combined_melted[combined_melted$BaseSample %in% group4_samples, ]

# Order each group by BaseSample and Tech
group1_data <- group1_data[order(group1_data$BaseSample, group1_data$Tech), ]
group2_data <- group2_data[order(group2_data$BaseSample, group2_data$Tech), ]
group3_data <- group3_data[order(group3_data$BaseSample, group3_data$Tech), ]
group4_data <- group4_data[order(group4_data$BaseSample, group4_data$Tech), ]

combined_melted$Alteration <- factor(combined_melted$Alteration, 
                                     levels = c("AMP", "GAIN", "NEUTRAL", "LOSS", "DEEPDEL"))
# Create a plotting function
create_plot <- function(data) {
  ggplot(data, aes(x = Var2, y = Var1, fill = Alteration)) +
    geom_tile() +
    scale_fill_manual(values = c("NEUTRAL" = "lightgreen", 
                                 "GAIN" = "pink", 
                                 "LOSS" = "lightblue1", 
                                 "AMP" = "red", 
                                 "DEEPDEL" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = function(x) {
      sapply(x, function(label) gsub("_call", "", label))
    })
}

# Generate the four plots
plot1 <- create_plot(group1_data)
plot2 <- create_plot(group2_data)
plot3 <- create_plot(group3_data)
plot4 <- create_plot(group4_data)

# Modify plot1 with new axis labels
plot1 <- plot1 + labs(x = "Sample", y = "Gene")

# Modify plot2 with new axis labels
plot2 <- plot2 + labs(x = "Sample", y = "Gene")

# Modify plot3 with new axis labels
plot3 <- plot3 + labs(x = "Sample", y = "Gene")

# Modify plot4 with new axis labels
plot4 <- plot4 + labs(x = "Sample", y = "Gene")

# Display the plots
plot1
plot2
plot3
plot4
