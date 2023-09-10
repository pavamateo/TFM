westable # file read it if not in environment

install.packages("ComplexUpset")
install.packages("ComplexHeatmap")
BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(ggplot2)

library(tidyr)
# hande duplicates
duplicates_for_same_gene <- filtered_data %>%
  group_by(Gene, PRO_id) %>%
  filter(n() > 1) %>%
  ungroup()



# Filter the data to include only the technologies you are interested in
filtered_data <- heatmap_data %>%
  filter(technology %in% c("cnvkit_call","ASCAT_call", "DNAcopy_call", "ichorcna_call", "qdnaseq_call", "panel_call"))

# Keep only the first occurrence of each combination of PRO_id, Gene, and technology
distinct_data <- filtered_data %>%
  distinct(PRO_id, Gene, technology, .keep_all = TRUE)

# Transform the data into the desired matrix format
oncoprint_data <- distinct_data %>%
  select(Gene, PRO_id, call) %>%
  spread(key = PRO_id, value = call)


# Extract the matrix from the dataframe
mat <- as.matrix(oncoprint_data[,-1])
rownames(mat) <- oncoprint_data$Gene


col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "#008000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  # big red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  # small green
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["MUT"], col = NA))
  }
)

# Now, call oncoPrint with the matrix
oncoPrint(mat, alter_fun = alter_fun, col = col)


# Find duplicated rows
duplicated_rows <- distinct_data[duplicated(distinct_data[, c("Gene", "PRO_id")]), ]

# Print the duplicated rows
head(duplicated_rows)

#### try

# Convert to long format
long_datas <- distinct_data %>%
  select(Gene, PRO_id, technology, call)

# Determine the most frequent call for each gene-sample pair
aggregated_data <- long_datas %>%
  group_by(Gene, PRO_id) %>%
  count(call) %>%
  arrange(-n) %>%
  slice(1) %>%
  ungroup() %>%
  select(-n)

# Convert to wide format for oncoprint
data_wide <- aggregated_data %>%
  spread(key = PRO_id, value = call)

data_wide <- as.data.frame(data_wide)

# Replace NA or empty values with "NO_MUT"
data_wide[is.na(data_wide) | data_wide == ""] <- "NO_MUT"


# Define the colors for different alterations
alter_col = c(
  "GAIN" = "lightcoral",
  "AMP" = "red",
  "LOSS" = "lightblue",
  "DEEPDEL" = "blue",
  "NEUTRAL" = "lightgreen"
)
# Update the colors to include "NO_MUT"
alter_col["NO_MUT"] = "#FFFFFF"  # White for no mutation


# Generate the oncoprint
# Define the colors and alter functions
alter_fun = list(
  AMP = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["AMP"], col = NA)),
  GAIN = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["GAIN"], col = NA)),
  NEUTRAL = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["NEUTRAL"], col = NA)),
  LOSS = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["LOSS"], col = NA)),
  DEEPDEL = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["DEEPDEL"], col = NA))
)
alter_fun$NO_MUT = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["NO_MUT"], col = NA))


# Create the oncoprint
oncoPrint(data_wide, alter_fun = alter_fun, col = alter_col)

###
rownames(data_wide) <- data_wide$Gene
data_wide$Gene <- NULL  # remove the Gene column since it's now the row names


# Define the colors for different alterations
alter_col = c(
  "GAIN" = "lightcoral",
  "AMP" = "red",
  "LOSS" = "lightblue",
  "DEEPDEL" = "blue",
  "NEUTRAL" = "lightgreen",
  "NO_MUT" = "white"
)

# Define the alterations functions for the oncoprint
alter_fun = list(
  GAIN = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["GAIN"], col = NA)),
  AMP = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["AMP"], col = NA)),
  LOSS = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["LOSS"], col = NA)),
  DEEPDEL = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["DEEPDEL"], col = NA)),
  NEUTRAL = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["NEUTRAL"], col = NA)),
  NO_MUT = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = alter_col["NO_MUT"], col = NA))
)

# Generate the oncoprint
oncoPrint(data_wide, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)

### try to fit everything in 1 single plot

# Create a new column that combines sample and technology
long_datasa$sample_technology <- paste0(long_datasa$PRO_id, "-", long_datasa$technology)
long_datas$sample_technology <- paste0(long_datas$PRO_id, "-", long_datas$technology) # this is for ascat

# Convert the data into a wide matrix format for oncoprint
oncoprint_matrix <- long_datasa %>%
  select(Gene, sample_technology, call) %>%
  spread(key = sample_technology, value = call) # no ascat

oncoprint_matrix <- long_datas %>%
  select(Gene, sample_technology, call) %>%
  spread(key = sample_technology, value = call) #ASCAT

# Convert to a traditional matrix
rownames(oncoprint_matrix) <- oncoprint_matrix$Gene
oncoprint_matrix$Gene <- NULL  # remove the Gene column since it's now the row names

oncoprint_matrix <- as.matrix(oncoprint_matrix) 
oncoprint_matrix[is.na(oncoprint_matrix)] <- ""

# Now, plot the oncoprint
oncoPrint(oncoprint_matrix, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)


########## DEAL WITH COLUMNS

# Create a sequence
secu <- 1:ncol(oncoprint_matrix)

# Insert NA after every 5th number
sep_seq <- as.vector(rbind(matrix(secu, ncol=5), ""))
sep_seq <- sep_seq[!is.na(sep_seq)]

# Plot the oncoprint with separators
oncoPrint(oncoprint_matrix, alter_fun = alter_fun, col = alter_col, column_order = sep_seq)

oncoPrint(...) + oncoPrint(...) + rowAnnotation(bar1 = anno_barplot(freq_table1), bar2 = anno_barplot(freq_table2))

######### summing oncoprints

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

oncoPrint(cnv_onco, alter_fun = alter_fun, col = alter_col,pct_side = FALSE)
oncoPrint(dna_onco, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)
oncoPrint(qdna_onco, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)
oncoPrint(ichor_onco, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)
oncoPrint(panel_onco, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)
oncoPrint(ascat_onco, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)


# save the oncoprints
# Save cnv_onco to PNG
png("cnv_onco.png", width=10, height=8, units="in", res=300)
oncoPrint(cnv_onco, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)
dev.off()

# Save dna_onco to PNG
png("dna_onco.png", width=10, height=8, units="in", res=300)
oncoPrint(dna_onco, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)
dev.off()

# Save qdna_onco to PNG
png("qdna_onco.png", width=10, height=8, units="in", res=300)
oncoPrint(qdna_onco, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)
dev.off()

# Save ichor_onco to PNG
png("ichor_onco.png", width=10, height=8, units="in", res=300)
oncoPrint(ichor_onco, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)
dev.off()

# Save panel_onco to PNG
png("panel_onco.png", width=10, height=8, units="in", res=300)
oncoPrint(panel_onco, alter_fun = alter_fun, col = alter_col, pct_side = FALSE)
dev.off()



(oncoPrint(cnv_onco, alter_fun = alter_fun, col = alter_col) + oncoPrint(dna_onco, alter_fun = alter_fun, col = alter_col) + oncoPrint(qdna_onco,alter_fun = alter_fun, col = alter_col) 
+ oncoPrint(ichor_onco, alter_fun = alter_fun, col = alter_col) + oncoPrint(panel_onco, alter_fun = alter_fun, col = alter_col))

# ggplot type heatmap for 1 tech
# Melt the data
cnv_onco_melted <- melt(cnv_onco, id.vars = NULL)

# Plot
ggplot(cnv_onco_melted, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("NEUTRAL" = "lightgreen", 
                               "GAIN" = "pink", 
                               "LOSS" = "lightblue", 
                               "AMP" = "red", 
                               "DEEPDEL" = "blue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# try for all techs

# Melt the data for each technology using R's reshape2 package
cnv_onco_melted <- melt(cnv_onco, variable.name = "Sample", value.name = "Alteration")
dna_onco_melted <- melt(dna_onco, variable.name = "Sample", value.name = "Alteration")
ichor_onco_melted <- melt(ichor_onco, variable.name = "Sample", value.name = "Alteration")
panel_onco_melted <- melt(panel_onco, variable.name = "Sample", value.name = "Alteration")
qdna_onco_melted <- melt(qdna_onco, variable.name = "Sample", value.name = "Alteration")
ascat_onco_melted <- melt(ascat_onco, variable.name = "Sample", value.name = "Alteration")

# Combine all melted data
combined_melted <- rbind(cnv_onco_melted, dna_onco_melted,ascat_onco_melted, ichor_onco_melted, qdna_onco_melted, panel_onco_melted) #works fine/ original

# Plot
ggplot(combined_melted, aes(x = Var2, y = Var1, fill = Alteration)) +
  geom_tile() +
  scale_fill_manual(values = c("NEUTRAL" = "lightgreen", 
                               "GAIN" = "pink", 
                               "LOSS" = "lightblue1", 
                               "AMP" = "red", 
                               "DEEPDEL" = "blue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# organize

# Assign a factor level for the sample order
combined_melted$SampleFactor <- as.numeric(factor(combined_melted$Var2, levels = unique(combined_melted$Var2)))

# Create a new column that combines the sample name and its factor level
combined_melted$Order <- with(combined_melted, interaction(SampleFactor, Var2, drop=TRUE))

# Sort the data frame by this new column
combined_melted <- combined_melted[order(combined_melted$SampleFactor),]

# Extract the sample name (assuming the format is always like 'PROxxxxx-technique')
combined_melted$plot <- gsub("-.*", "", combined_melted$Var2)

# Plot works perfect
ggplot(combined_melted, aes(x = Var2, y = Var1, fill = Alteration)) +
  geom_tile() +
  scale_fill_manual(values = c("NEUTRAL" = "lightgreen", 
                               "GAIN" = "pink", 
                               "LOSS" = "lightblue", 
                               "AMP" = "red", 
                               "DEEPDEL" = "blue")) +
  scale_x_discrete(labels = function(x) {
    s <- unique(substr(x, 1, regexpr("-", x)-1))
    rep_len(s, length(x))
  }) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), # Adjusting the font size here
        aspect.ratio = 0.5) # Adjusting the height-to-width ratio

# dont know
ggplot(combined_melted, aes(x = Var2, y = Var1, fill = Alteration)) +
  geom_tile() +
  scale_fill_manual(values = c("NEUTRAL" = "lightgreen", 
                               "GAIN" = "pink", 
                               "LOSS" = "lightblue", 
                               "AMP" = "red", 
                               "DEEPDEL" = "blue")) +
  scale_x_discrete(labels = function(x) {
    s <- unique(substr(x, 1, regexpr("-", x)-1))
    rep_len(s, length(x))
  }) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), # Adjusting the font size here
        aspect.ratio = 0.5) # Adjusting the height-to-width ratio


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

# ggplot(combined_melted, aes(x = BaseSample, fill = Alteration

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

######### splitting into 4

# Split the data into four quarters
unique_sample <- unique(combined_melted$BaseSample)
quarter_length <- length(unique_sample) %/% 4

first_quarter_samples <- unique_sample[1:quarter_length]
second_quarter_samples <- unique_sample[(quarter_length + 1):(2*quarter_length)]
third_quarter_samples <- unique_sample[(2*quarter_length + 1):(3*quarter_length)]
fourth_quarter_samples <- unique_sample[(3*quarter_length + 1):length(unique_sample)]

first_quarter_data <- combined_melted[combined_melted$BaseSample %in% first_quarter_samples, ]
second_quarter_data <- combined_melted[combined_melted$BaseSample %in% second_quarter_samples, ]
third_quarter_data <- combined_melted[combined_melted$BaseSample %in% third_quarter_samples, ]
fourth_quarter_data <- combined_melted[combined_melted$BaseSample %in% fourth_quarter_samples, ]

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
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Generate the plots
plot1 <- create_plot(first_quarter_data)
plot2 <- create_plot(second_quarter_data)
plot3 <- create_plot(third_quarter_data)
plot4 <- create_plot(fourth_quarter_data)

# Display the plots
plot1
plot2
plot3
plot4


## FINAL YES
# Combine all melted data
combined_melted <- rbind(cnv_onco_melted, dna_onco_melted,ascat_onco_melted, ichor_onco_melted, qdna_onco_melted, panel_onco_melted)
# extract the technology
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

write.table(combined_melted, file="combined_melted.csv", row.names=F, col.names=T, quote=F, sep="\t")


# Save plot1
ggsave(filename = "samples1.png", plot = plot1, width = 10, height = 8, dpi = 300)

# Save plot2
ggsave(filename = "samples2.png", plot = plot2, width = 10, height = 8, dpi = 300)

# Save plot3
ggsave(filename = "samples3.png", plot = plot3, width = 10, height = 8, dpi = 300)

# Save plot4
ggsave(filename = "samples4.png", plot = plot4, width = 10, height = 8, dpi = 300)

### all samples

library(ggplot2)
library(dplyr)
library(purrr)

# Function to create the plot for a given sample's data
plot_for_sample <- function(sample_data) {
  ggplot(sample_data, aes(x = Tech, y = Var1, fill = Alteration)) +
    geom_tile() +
    scale_fill_manual(values = c("NEUTRAL" = "lightgreen", 
                                 "GAIN" = "pink", 
                                 "LOSS" = "lightblue1", 
                                 "AMP" = "red", 
                                 "DEEPDEL" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Split the data by sample and create a list of plots
plots_list <- combined_melted %>%
  split(.$BaseSample) %>%
  map(plot_for_sample)

# The plots_list now contains a ggplot for each of the samples

