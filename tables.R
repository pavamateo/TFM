install.packages(c("officer", "flextable"))
library(officer)
library(flextable)

# Create the data frame for table thresholds
cnv_thresholds <- data.frame(
  Condition = c(rep("Low TF ≤ 0.3", 5), rep("High/Normal TF > 0.3", 5), rep("chr X", 4)),
  Call = c("AMP", "GAIN", "NEUTRAL", "LOSS", "DEEPDEL", 
           "AMP", "GAIN", "NEUTRAL", "LOSS", "DEEPDEL", 
           "AMP", "GAIN", "NEUTRAL", "DEEPDEL"),
  Threshold_from = c("0.8", "0.25", "-0.25", "-0.7", "< -0.7", 
                     "0.7", "0.2", "-0.2", "-0.58", "< -0.58", 
                     "1.6", "0.6", "-1.3", "-1.3"),
  Threshold_to = c("> 0,8", "0.8", "0.25", "-0.25", "-0.7", 
                   "< -0.7", "0.7", "0.2", "-0.2", "-0.58", 
                   "> 1.6", "1.6", "0.6", "< -1.3")
)

# save it 
table_grob <- tableGrob(cnv_thresholds)
png(filename = "thresh.png", width = 800, height = 600)
grid.draw(table_grob)
dev.off()

tablegenes <- genesanalyze[,c(1,2,3,4,5)]
# Convert data frame to flextable
ft <- flextable::flextable(head(tablegenes))

# Create a new Word document
doc <- officer::read_docx()

# Add the flextable to the document and save it
doc <- flextable::save_as_docx(ft, path = "my_table.docx")

wes<- as_data_frame(cnvkit_vs_DNAcopy)

ws <- flextable(wes)

tumortable<-tumorfraction[,c(1,2,3,4)]

tf <- flextable(head(tumortable))

doc <- officer::read_docx()

# Add the flextable to the document and save it
doc <- flextable::save_as_docx(mt, path = "my_table.docx")

mt<- flextable(modetable)


modetable <- mode_data %>%
  filter(Gene %in% genes_pro)

print(modetable)

westable <- wes_wgsfinal
westable$HandE <- NULL

write.table(westable, file = "finaltechs.csv", row.names = F, col.names = T, quote = F, sep = "\t")

westable <- read.csv("/Users/mateopava/Desktop/CNA_Mateo_MSC/allpanel/finaltechs.csv", sep = "\t")
tumorfraction <- read.csv("/Users/mateopava/Desktop/CNA_Mateo_MSC/allpanel/TumorFractionFinal.csv", sep = "\t")

methods_df <- data.frame(
  Method = c("CNVKit", "DNAcopy", "ASCAT", "ichorCNA", "QDNAseq"),
  Technology = c("WES/Panel", "WES", "WES", "sWGS", "sWGS"),
  Ploidy = c("NO", "NO", "YES", "YES", "NO"),
  TF = c("NO", "NO", "YES", "YES", "NO"),
  Implementation = c("Python", "Varscan/R", "nf-core/sarek", "*", "R"),
  Segmentation_method = c("*","*","*","*","*"),
  Version = c("*", "*", "*", "*", "*")

)

# CONTINGENCY TABLES

# Merge the table with itself to get pairs of technologies for each Gene and BaseSample
contingency <- merge(combined_melted, combined_melted, by = c("Var1", "BaseSample"))

# Filter the merged data to get the desired pairs of technologies
cnvkit_vs_dnacopy <- contingency[contingency$Tech.x == "cnvkit_call" & contingency$Tech.y == "DNAcopy_call", ]
cnvkit_vs_ascat <- contingency[contingency$Tech.x == "cnvkit_call" & contingency$Tech.y == "ASCAT_call", ]
ascat_vs_dnacopy <- contingency[contingency$Tech.x == "ASCAT_call" & contingency$Tech.y == "DNAcopy_call", ]
ichorcna_vs_qdnaseq <- contingency[contingency$Tech.x == "ichorcna_call" & contingency$Tech.y == "qdnaseq_call", ]
# adding panel
cnvkit_vs_panel <- contingency[contingency$Tech.x == "cnvkit_call" & contingency$Tech.y == "panel_call", ]
ascat_vs_panel <- contingency[contingency$Tech.x == "ASCAT_call" & contingency$Tech.y == "panel_call", ]
ichorcna_vs_panel <- contingency[contingency$Tech.x == "ichorcna_call" & contingency$Tech.y == "panel_call", ]
qdnaseq_vs_panel <- contingency[contingency$Tech.x == "qdnaseq_call" & contingency$Tech.y == "panel_call", ]
dnacopy_vs_panel <- contingency[contingency$Tech.x == "DNAcopy_call" & contingency$Tech.y == "panel_call", ]
# adding cross methods
cnvkit_vs_ichorcna <- contingency[contingency$Tech.x == "cnvkit_call" & contingency$Tech.y == "ichorcna_call", ]
cnvkit_vs_qdnaseq <- contingency[contingency$Tech.x == "cnvkit_call" & contingency$Tech.y == "qdnaseq_call", ]
ascat_vs_qdnaseq <- contingency[contingency$Tech.x == "ASCAT_call" & contingency$Tech.y == "qdnaseq_call", ]
ascat_vs_ichorcna <- contingency[contingency$Tech.x == "ASCAT_call" & contingency$Tech.y == "ichorcna_call", ]
dnacopy_vs_ichorcna <- contingency[contingency$Tech.x == "DNAcopy_call" & contingency$Tech.y == "ichorcna_call", ]
dnacopy_vs_qdnaseq <- contingency[contingency$Tech.x == "DNAcopy_call" & contingency$Tech.y == "qdnaseq_call", ]

# Create contingency tables
table_cnvkit_vs_dnacopy <- table(cnvkit_vs_dnacopy$Alteration.x, cnvkit_vs_dnacopy$Alteration.y)
table_cnvkit_vs_ascat <- table(cnvkit_vs_ascat$Alteration.x, cnvkit_vs_ascat$Alteration.y)
table_ascat_vs_dnacopy <- table(ascat_vs_dnacopy$Alteration.x, ascat_vs_dnacopy$Alteration.y)
table_ichorcna_vs_qdnaseq <- table(ichorcna_vs_qdnaseq$Alteration.x, ichorcna_vs_qdnaseq$Alteration.y)

# add panel
table_cnvkit_vs_panel <- table(cnvkit_vs_panel$Alteration.x, cnvkit_vs_panel$Alteration.y)
table_ascat_vs_panel <- table(ascat_vs_panel$Alteration.x, ascat_vs_panel$Alteration.y)
table_ichorcna_vs_panel <- table(ichorcna_vs_panel$Alteration.x, ichorcna_vs_panel$Alteration.y)
table_qdnaseq_vs_panel <- table(qdnaseq_vs_panel$Alteration.x, qdnaseq_vs_panel$Alteration.y)
table_dnacopy_vs_panel <- table(dnacopy_vs_panel$Alteration.x, dnacopy_vs_panel$Alteration.y)

# adding cross methods
table_cnvkit_vs_ichorcna <- table(cnvkit_vs_ichorcna$Alteration.x, cnvkit_vs_ichorcna$Alteration.y)
table_cnvkit_vs_qdnaseq <- table(cnvkit_vs_qdnaseq$Alteration.x, cnvkit_vs_qdnaseq$Alteration.y)
table_ascat_vs_qdnaseq <- table(ascat_vs_qdnaseq$Alteration.x, ascat_vs_qdnaseq$Alteration.y)
table_ascat_vs_ichorcna <- table(ascat_vs_ichorcna$Alteration.x, ascat_vs_ichorcna$Alteration.y)
table_dnacopy_vs_ichorcna <- table(dnacopy_vs_ichorcna$Alteration.x, dnacopy_vs_ichorcna$Alteration.y)
table_dnacopy_vs_qdnaseq <- table(dnacopy_vs_qdnaseq$Alteration.x, dnacopy_vs_qdnaseq$Alteration.y)


tables_to_save <- list(
  cnvkit_vs_dnacopy = table_cnvkit_vs_dnacopy,
  cnvkit_vs_ascat = table_cnvkit_vs_ascat,
  ascat_vs_dnacopy = table_ascat_vs_dnacopy,
  ichorcna_vs_qdnaseq = table_ichorcna_vs_qdnaseq,
  cnvkit_vs_panel = table_cnvkit_vs_panel,
  ascat_vs_panel = table_ascat_vs_panel,
  qdnaseq_vs_panel = table_qdnaseq_vs_panel,
  dnacopy_vs_panel = table_dnacopy_vs_panel,
  ichorcna_vs_panel = table_ichorcna_vs_panel,
  cnvkit_vs_ichorcna = table_cnvkit_vs_ichorcna,
  cnvkit_vs_qdnaseq = table_cnvkit_vs_qdnaseq,
  ascat_vs_qdnaseq = table_ascat_vs_qdnaseq,
  ascat_vs_ichorcna = table_ascat_vs_ichorcna,
  dnacopy_vs_ichorcna = table_dnacopy_vs_ichorcna,
  dnacopy_vs_qdnaseq = table_dnacopy_vs_qdnaseq
)


# Iterate over each table
for (table_name in names(tables_to_save)) {
  # Define a filename based on table name
  filename <- paste0(table_name, ".png")
  # Save the table as an image
  save_table_as_image(tables_to_save[[table_name]], filename)
}

# Display the tables
table_cnvkit_vs_dnacopy
table_cnvkit_vs_ascat
table_ascat_vs_dnacopy
table_ichorcna_vs_qdnaseq
table_cnvkit_vs_panel
table_ascat_vs_panel
table_qdnase_vs_panel
table_dnacopy_vs_panel
table_ichorcna_vs_panel
table_cnvkit_vs_ichorcna
table_cnvkit_vs_qdnaseq
table_ascat_vs_qdnaseq
table_ascat_vs_ichorcna
table_dnacopy_vs_ichorcna
table_dnacopy_vs_qdnaseq

# create 1 single data
summary_data <- data.frame(
  Alteration = c("AMP", "GAIN", "NEUTRAL", "LOSS", "DEEPDEL"),
  cnvkit_vs_dnacopy = rowSums(table_cnvkit_vs_dnacopy),
  cnvkit_vs_ascat = rowSums(table_cnvkit_vs_ascat),
  ascat_vs_dnacopy = rowSums(table_ascat_vs_dnacopy),
  ichorcna_vs_qdnaseq = rowSums(table_ichorcna_vs_qdnaseq),
  cnvkit_vs_panel = rowSums(table_cnvkit_vs_panel),
  ascat_vs_panel = rowSums(table_ascat_vs_panel),
  qdnase_vs_panel = rowSums(table_qdnase_vs_panel),
  dnacopy_vs_panel = rowSums(table_dnacopy_vs_panel),
  ichorcna_vs_panel = rowSums(table_ichorcna_vs_panel)
)
  
  
  
  
# below not necessary/ dont reme,ber what it is
sample_data <- combined_melted[combined_melted$BaseSample == "PRO001_18", ]
genes_with_panel_calls <- unique(sample_data[sample_data$Tech == "panel", "Var1"])
sample_data <- sample_data[sample_data$Var1 %in% genes_with_panel_calls, ]
compare_to_panel <- function(tech_name, data) {
  panel_calls <- data[data$Tech == "panel", "Alteration"]
  tech_calls <- data[data$Tech == "cnvkit", "Alteration"]
  table(panel_calls, tech_calls)
}
table_cnvkit_vs_panel <- compare_to_panel("cnvkit", sample_data)
table_dnacopy_vs_panel <- compare_to_panel("DNAcopy", sample_data)
table_ascat_vs_panel <- compare_to_panel("ASCAT", sample_data)


## corr matrixx THIS IS USEFUL
log2_columns <- wes_wgsfinal[, grepl("_log2", colnames(wes_wgsfinal))]
log2_columns$ASCAT_log2_diploid<- NULL

# Compute correlation matrix
cor_matrix <- cor(log2_columns, use = "pairwise.complete.obs")

print(cor_matrix)

# Compute p-values for correlations
cor_pvalues <- cor.mtest(log2_columns)$p

# Plot correlation matrix
corrplot(cor_matrix, method = "circle", p.mat = cor_pvalues, sig.level = 0.05, insig = "blank")


cor_matrix <- cor(tumorfraction[, c("HandE", "TumorFractionIchorCNA", "TumorFractionASCAT")], use = "pairwise.complete.obs")
# Compute p-values for correlations
cor_pvalues <- cor.mtest(tumorfraction[, c("HandE", "TumorFractionIchorCNA", "TumorFractionASCAT")])$p
corrplot(cor_matrix, method = "circle")

# Concordance index

# Given table as matrix
table_ichorcna_vs_qdnaseq

# Calculate concordance index
co <- nrow(table_ichorcna_vs_qdnaseq)
concordance_index <- sum(diag(table_ichorcna_vs_qdnaseq)) - sum(table_ichorcna_vs_qdnaseq) + sum(diag(table_ichorcna_vs_qdnaseq))
concordance_index

### all concordance for techs

# Function to compute the concordance index for a contingency table
concordance_index_func <- function(contingency_table) {
  # Calculate concordance using the sum of the diagonal values
  concordance <- sum(diag(contingency_table))
  
  # Calculate the total number of observations in the table
  total_observations <- sum(contingency_table)
  
  # Calculate and return the concordance index
  return(concordance / total_observations)
}

# List of tables
concor_list <- list(
  cnvkit_vs_dnacopy = table_cnvkit_vs_dnacopy,
  cnvkit_vs_ascat = table_cnvkit_vs_ascat,
  ascat_vs_dnacopy = table_ascat_vs_dnacopy,
  ichorcna_vs_qdnaseq = table_ichorcna_vs_qdnaseq,
  cnvkit_vs_panel = table_cnvkit_vs_panel,
  ascat_vs_panel = table_ascat_vs_panel,
  qdnaseq_vs_panel = table_qdnaseq_vs_panel,
  dnacopy_vs_panel = table_dnacopy_vs_panel,
  ichorcna_vs_panel = table_ichorcna_vs_panel,
  cnvkit_vs_ichorcna = table_cnvkit_vs_ichorcna,
  cnvkit_vs_qdnaseq = table_cnvkit_vs_qdnaseq,
  ascat_vs_qdnaseq = table_ascat_vs_qdnaseq,
  ascat_vs_ichorcna = table_ascat_vs_ichorcna,
  dnacopy_vs_ichorcna = table_dnacopy_vs_ichorcna,
  dnacopy_vs_qdnaseq = table_dnacopy_vs_qdnaseq
  
)

# Calculate the concordance index for each table
concordance_indices <- lapply(concor_list, concordance_index_func)

# Print the results
print(concordance_indices)

summary_table <- data.frame(
  Comparison = names(concordance_indices),
  Concordance_Index = unlist(concordance_indices)
)

summary_table$Concordance_Index <- round(summary_table$Concordance_Index, 3)

# Renaming the methods
names_to_replace <- c("cnvkit", "ascat", "dnacopy", "ichorcna", "qdnaseq")
correct_names <- c("CNVkit", "ASCAT", "DNAcopy", "ichorCNA", "QDNAseq")

for (i in 1:length(names_to_replace)) {
  summary_table$Comparison <- gsub(names_to_replace[i], correct_names[i], summary_table$Comparison)
}
# save it 
table_grob <- tableGrob(summary_table)
png(filename = "concordance.png", width = 800, height = 600)
grid.draw(table_grob)
dev.off()


# download gene table

genefordownload <-genesanalyze
genefordownload$entrezgene_id <- NULL
genefordownload
# Rename the columns
colnames(genefordownload) <- c("Gene", "Start", "End", "Chromosome")
table_grob <- tableGrob(head(genefordownload))
png(filename = "table_image.png", width = 800, height = 600)
grid.draw(table_grob)
dev.off()



# genes contingency

# Filter data for BRCA1
brca1_data <- wes_wgsfinal[wes_wgsfinal$Gene == "BRCA1", ]

# Ensure unique PRO_id
brca1_data <- brca1_data[!duplicated(brca1_data$PRO_id), ]

# Create the specified contingency tables
table_cnvkit_vs_ascat <- table(brca1_data$cnvkit_call, brca1_data$ASCAT_call)
table_cnvkit_vs_DNAcopy <- table(brca1_data$cnvkit_call, brca1_data$DNAcopy_call)
table_DNAcopy_vs_ascat <- table(brca1_data$DNAcopy_call, brca1_data$ASCAT_call)
table_ichorcna_vs_qdnaseq <- table(brca1_data$ichorcna_call, brca1_data$qdnaseq_call)
table_panel_vs_cnvkit <- table(brca1_data$panel_call, brca1_data$cnvkit_call)
table_panel_vs_ichorcna <- table(brca1_data$panel_call, brca1_data$ichorcna_call)
table_panel_vs_qdnaseq <- table(brca1_data$panel_call, brca1_data$qdnaseq_call)
table_panel_vs_DNAcopy <- table(brca1_data$panel_call, brca1_data$DNAcopy_call)
table_panel_vs_ascat <- table(brca1_data$panel_call, brca1_data$ASCAT_call)

# You can print out these tables to inspect them
print(table_cnvkit_vs_ascat)
print(table_cnvkit_vs_DNAcopy)
print(table_DNAcopy_vs_ascat)
print(table_ichorcna_vs_qdnaseq)
print(table_panel_vs_cnvkit)
print(table_panel_vs_ichorcna)
print(table_panel_vs_qdnaseq)
print(table_panel_vs_DNAcopy)
print(table_panel_vs_ascat)

genes_of_interest <- c("AR", "RAD51", "ATM")
tables_list <- list()

for (gene in genes_of_interest) {
  gene_data <- wes_wgsfinal[wes_wgsfinal$Gene == gene, ]
  
  # Ensure unique PRO_id
  gene_data <- gene_data[!duplicated(gene_data$PRO_id), ]
  
  # Create the specified contingency tables
  tables_list[[gene]] <- list(
    cnvkit_vs_ascat = table(gene_data$cnvkit_call, gene_data$ASCAT_call),
    cnvkit_vs_DNAcopy = table(gene_data$cnvkit_call, gene_data$DNAcopy_call),
    DNAcopy_vs_ascat = table(gene_data$DNAcopy_call, gene_data$ASCAT_call),
    ichorcna_vs_qdnaseq = table(gene_data$ichorcna_call, gene_data$qdnaseq_call),
    panel_vs_cnvkit = table(gene_data$panel_call, gene_data$cnvkit_call),
    panel_vs_ichorcna = table(gene_data$panel_call, gene_data$ichorcna_call),
    panel_vs_qdnaseq = table(gene_data$panel_call, gene_data$qdnaseq_call),
    panel_vs_DNAcopy = table(gene_data$panel_call, gene_data$DNAcopy_call),
    panel_vs_ascat = table(gene_data$panel_call, gene_data$ASCAT_call)
  )
}
# You can access them using:
# tables_list$BRCA2$cnvkit_vs_ascat and so on

save_table_as_image <- function(tbl, filename) {
  if (nrow(tbl) == 0 || ncol(tbl) == 0) {
    message(paste("Skipping empty table for:", filename))
    return()
  }
  png(filename, width=800, height=800)
  grid.table(tbl)
  dev.off()
}

# Iterate over each gene and its associated tables
for (gene in names(tables_list)) {
  for (comparison in names(tables_list[[gene]])) {
    # Define a filename based on gene and comparison
    filename <- paste0(gene, "_", comparison, ".png")
    # Save the table as an image
    save_table_as_image(tables_list[[gene]][[comparison]], filename)
  }
}

# This will save each contingency table as a PNG image with filenames indicating the gene and comparison. Adjust the width and height in the save_table_as_image function if needed.








