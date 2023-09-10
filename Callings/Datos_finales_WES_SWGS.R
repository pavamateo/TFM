formato_tech <- read_xlsx("/Users/mateopava/Desktop/Luisa_DATOS/FORMATEOsample_correspondencies_withTFandHyE.xlsx")

swgs_ichorcna <- read_xlsx("/Users/mateopava/Desktop/Luisa_DATOS/logR_qdnaseqVsichorCNA_withcall.xlsx")

swgs_ichorv <- read_xlsx("/Users/mateopava/Desktop/Luisa_DATOS/logR_qdnaseqVsichorCNA.xlsx")

# this to make it logical TRUE/FALSE for upset
#formato_tech$P300_assayDone <- as.logical(formato_tech$P300_assayDone)
#formato_tech$WES_assayDone <- as.logical(formato_tech$WES_assayDone)
#formato_tech$sWGS_assayDone <- as.logical(formato_tech$sWGS_assayDone)
formato_tech <-as.data.frame(formato_tech)
swgs_ichorcna <-as.data.frame(swgs_ichorcna)
swgs_ichorv <-as.data.frame(swgs_ichorv)

subset_data <- formato_tech[, c("P300_assayDone", "WES_assayDone", "sWGS_assayDone")]
upset(subset_data, nsets = 3)
# Convert data to long format
long_data <- subset_data %>% 
  pivot_longer(everything(), names_to = "assay", values_to = "done")

# Count occurrences of each assay being done or not
count_data <- long_data %>% 
  group_by(assay, done) %>% 
  tally()

# Convert back to a wide format
wide_data <- count_data %>% 
  pivot_wider(names_from = "assay", values_from = "n")

wide_data<-as.data.frame(wide_data)

# Use UpSetR to visualize
upset(wide_data)


#### FROM HERE ITS NOT GOOD

# Store gene names as row names
# Set row names to gene names
rownames(swgs_ichorcna) <- swgs_ichorcna$Gene
swgs_ichorcna$Gene <- NULL # Remove the gene column after setting it as row names

# ASCAT
# I WANT TO READ ALL FILES

# PRO001_18-27110_DNAFFPE_NA_C_HN00131668_vs_PRO001_NA_DNASAL_NA_C_HN00131668.metrics.txt
# PRO001_18-27110_DNAFFPE_NA_C_HN00131668_vs_PRO001_NA_DNASAL_NA_C_HN00131668.purityploidy.txt

# set wd to the directory of the files
setwd("/Users/mateopava/Desktop/ascat/ascat_analyze")
all_files <- list.files()

# get all the files to read them in a df
metrics_files <- all_files[grepl(".metrics.txt$", all_files)]
purityploidy_files <- all_files[grepl(".purityploidy.txt$", all_files)]


#list it
list_of_metrics_dfs <- lapply(metrics_files, function(file) {
  df <- read.table(file, header = TRUE, sep = "\t") # assuming they are tab-delimited
  
  # Extract the unique identifier from the filename
  identifier <- gsub("(^PRO[0-9]+).*", "\\1", file)
  
  # Add a new column to store the identifier
  df$Origin <- identifier
  return(df)
})

metrics_df <- do.call(rbind, list_of_metrics_dfs)

#write.csv(metrics_df, file = "metrics.csv", row.names = FALSE)

# purity only
list_of_purityploidy_dfs <- lapply(purityploidy_files, function(file) {
  df <- read.table(file, header = TRUE, sep = "\t") # assuming they are tab-delimited
  
  # Extract the unique identifier from the filename
  identifier <- gsub(".*(PRO[0-9]+_[0-9]+_[A-Z]+_[A-Z]+_[A-Z]+_HN[0-9]+).*", "\\1", file)
  
  # Add a new column to store the identifier
  df$Origin <- identifier
  return(df)
})

purityploidy_df <- do.call(rbind, list_of_purityploidy_dfs)
purityploidy_df$Origin <- gsub("_vs.*", "", purityploidy_df$Origin)


# Convert character "NA" to actual NA
formato_tech$TumorFractionIchorCNA[formato_tech$TumorFractionIchorCNA == "NA"] <- NA

formato_tech[formato_tech == "NA"] <- NA #convert all

# Convert the column to numeric
formato_tech$TumorFractionIchorCNA <- as.numeric(formato_tech$TumorFractionIchorCNA)

# Now divide by 100 to convert to percentages
formato_tech$TumorFractionIchorCNA <- formato_tech$TumorFractionIchorCNA / 100

# comparison from ichorcna and ascat, swgs/wes

formato_tech_filtered <- formato_tech[complete.cases(formato_tech$TumorFractionIchorCNA), ]
formato_tech_filtered <- formato_tech_filtered[formato_tech_filtered$WES_fileID != "<NA>", ]
formato_tech_filtered <- formato_tech_filtered[!is.na(formato_tech_filtered$TumorFractionIchorCNA), ]


compare_tech <- merge(formato_tech_filtered, purityploidy_df, by.x="WES_fileID", by.y="Origin", all.x=TRUE)
# this below works
compare_tech <- compare_tech[!is.na(compare_tech$AberrantCellFraction), ]

compare_tech$HandE <- as.numeric(compare_tech$HandE)

compare_tech$HandE <- compare_tech$HandE / 100

rownames(compare_tech)<- NULL

# correcting NA in HandE
mean_HandE <- mean(compare_tech$HandE, na.rm = TRUE)
compare_tech$HandE[is.na(compare_tech$HandE)] <- mean_HandE





#comparing
summary(compare_tech[, c("TumorFractionIchorCNA", "AberrantCellFraction", "HandE")])


# Summary statistics for TumorFractionIchorCNA
summary(compare_tech$TumorFractionIchorCNA)

cor(compare_tech$TumorFractionIchorCNA, compare_tech$AberrantCellFraction, method="pearson") # 0.4855477 


library(ggplot2)

# Create a data frame for plotting
dfplot <- data.frame(
  Index = 1:nrow(compare_tech),
  TumorFraction = compare_tech$TumorFractionIchorCNA,
  AberrantCell = compare_tech$AberrantCellFraction
)

# Scatter plot with regression lines
ggplot(dfplot, aes(x=Index)) +
  geom_point(aes(y=TumorFraction, color="Tumor Fraction IchorCNA")) +
  geom_point(aes(y=AberrantCell, color="Aberrant Cell Fraction")) +
  geom_smooth(aes(y=TumorFraction, color="Tumor Fraction IchorCNA"), method='lm', se=FALSE) +
  geom_smooth(aes(y=AberrantCell, color="Aberrant Cell Fraction"), method='lm', se=FALSE) +
  scale_color_manual(values=c("Tumor Fraction IchorCNA"="blue", "Aberrant Cell Fraction"="red")) +
  ggtitle("Scatter Plot of Tumor Fraction and Aberrant Cell Fraction with Regression Lines") +
  xlab("Sample Index") +
  ylab("Value") +
  theme(legend.title=element_blank())

#Histogram
# Determine x-axis limits
x_limits <- range(c(compare_tech$TumorFractionIchorCNA, compare_tech$AberrantCellFraction))

# Plot histograms with consistent scales
par(mfrow=c(2,1))

hist(compare_tech$TumorFractionIchorCNA, 
     main="Histogram of Tumor Fraction sWGS IchorCNA", 
     xlab="Tumor Fraction IchorCNA", 
     col="lightblue", 
     border="black", 
     xlim=x_limits, 
     freq=FALSE)

hist(compare_tech$AberrantCellFraction, 
     main="Histogram of Tumor Fraction WES ASCAT", 
     xlab="Tumor Fraction ASCAT", 
     col="lightblue", 
     border="black", 
     xlim=x_limits, 
     freq=FALSE)

#boxplot
boxplot(compare_tech$TumorFractionIchorCNA, compare_tech$AberrantCellFraction, 
        names=c("Tumor Fraction IchorCNA", "Aberrant Cell Fraction"),
        main="Boxplot of Tumor Fraction IchorCNA and Aberrant Cell Fraction", 
        col=c("pink", "lightblue"))

# Paired t-test
t.test(compare_tech$TumorFractionIchorCNA, compare_tech$AberrantCellFraction, paired=TRUE)
# Wilcoc
wilcox.test(compare_tech$TumorFractionIchorCNA, compare_tech$AberrantCellFraction, paired=TRUE)

shapiro.test(compare_tech$TumorFractionIchorCNA)
shapiro.test(compare_tech$AberrantCellFraction)

############### histology, swgs-ichorna, wes-ascat
# Pairwise scatter plots
pairs(compare_tech[, c("TumorFractionIchorCNA", "AberrantCellFraction", "HandE")], 
      pch = 19, col = "blue",
      main = "Pairwise Scatter Plots")

# Histograms
# Determine the overall range across the three variables
overall_range <- range(c(compare_tech$TumorFractionIchorCNA, compare_tech$AberrantCellFraction, compare_tech$HandE))

# Plot histograms with the same x-axis scale
par(mfrow=c(3,1))
hist(compare_tech$TumorFractionIchorCNA, main="Histogram of Tumor Fraction IchorCNA sWGS", xlab="Tumor Fraction IchorCNA", col="lightblue", border="black", xlim=overall_range)
hist(compare_tech$AberrantCellFraction, main="Histogram of Tumor Fraction ASCAT WES", xlab="Tumor Fraction ASCAT", col="lightgreen", border="black", xlim=overall_range)
hist(compare_tech$HandE, main="Histogram of Pathology HandE", xlab="Tumor Fraction HandE", col="lightpink", border="black", xlim=overall_range)

# Wilcoxon signed-rank test between TumorFractionIchorCNA and AberrantCellFraction
wilcox1 <- wilcox.test(compare_tech$TumorFractionIchorCNA, compare_tech$AberrantCellFraction, paired=TRUE)

# Wilcoxon signed-rank test between TumorFractionIchorCNA and HandE
wilcox2 <- wilcox.test(compare_tech$TumorFractionIchorCNA, compare_tech$HandE, paired=TRUE)

# Wilcoxon signed-rank test between AberrantCellFraction and HandE
wilcox3 <- wilcox.test(compare_tech$AberrantCellFraction, compare_tech$HandE, paired=TRUE)

wilcox1
wilcox2
wilcox3


# Scatter plot of TumorFractionIchorCNA vs HandE
par(mfrow=c(2,1))
plot(compare_tech$HandE, compare_tech$TumorFractionIchorCNA, main="TumorFractionIchorCNA vs HandE", xlab="HandE", ylab="TumorFractionIchorCNA", col="blue", pch=19)
abline(a=0, b=1, col="red", lty=2) # line of perfect agreement

# Compute and display R-squared for TumorFractionIchorCNA vs HandE
model1 <- lm(TumorFractionIchorCNA ~ HandE, data=compare_tech)
rsq1 <- summary(model1)$r.squared
text(min(compare_tech$HandE) + 0.05, max(compare_tech$TumorFractionIchorCNA) - 0.05, paste("R^2 =", round(rsq1, 3)), pos=4, col="blue")

# Scatter plot of AberrantCellFraction vs HandE
plot(compare_tech$HandE, compare_tech$AberrantCellFraction, main="AberrantCellFraction vs HandE", xlab="HandE", ylab="AberrantCellFraction", col="green", pch=19)
abline(a=0, b=1, col="red", lty=2) # line of perfect agreement

# Compute and display R-squared for AberrantCellFraction vs HandE
model2 <- lm(AberrantCellFraction ~ HandE, data=compare_tech)
rsq2 <- summary(model2)$r.squared
text(0.1, 0.9, paste("R^2 =", round(rsq2, 3)), pos=4, col="green")

### compare test

# Compute differences
compare_tech$Diff_TumorichorCNA_HandE <- compare_tech$TumorFractionIchorCNA - compare_tech$HandE
compare_tech$Diff_TumorASCAT_HandE <- compare_tech$TumorFractionASCAT - compare_tech$HandE

# Summary of differences
summary(compare_tech[, c("Diff_TumorichorCNA_HandE", "Diff_TumorASCAT_HandE")])

cor_TumorichorCNA_HandE <- cor(compare_tech$HandE, compare_tech$TumorFractionIchorCNA, method="pearson")
cor_TumorASCAT_HandE <- cor(compare_tech$HandE, compare_tech$TumorFractionASCAT, method="pearson")

cor_TumorichorCNA_HandE
cor_TumorASCAT_HandE

# Bland-Altman plot for TumorFractionIchorCNA vs HandE
plot((compare_tech$HandE + compare_tech$TumorFractionIchorCNA)/2, compare_tech$Diff_Tumor_HandE, 
     xlab="Average of HandE and TumorFractionIchorCNA", ylab="Difference (Tumor - HandE)",
     main="Bland-Altman Plot: TumorFractionIchorCNA vs HandE", pch=19, col="blue")
abline(h=mean(compare_tech$Diff_Tumor_HandE), col="red", lty=2) # mean difference line

# Bland-Altman plot for AberrantCellFraction vs HandE
plot((compare_tech$HandE + compare_tech$AberrantCellFraction)/2, compare_tech$Diff_Aberrant_HandE, 
     xlab="Average of HandE and AberrantCellFraction", ylab="Difference (Aberrant - HandE)",
     main="Bland-Altman Plot: AberrantCellFraction vs HandE", pch=19, col="green")
abline(h=mean(compare_tech$Diff_Aberrant_HandE), col="red", lty=2) # mean difference line

View(compare_tech)
# ORGANIZING RESUME DATA remember the above code has old name
# Assuming df is your dataframe and you want to rename the column "old_name" to "new_name"
colnames(compare_tech)[colnames(compare_tech) == "AberrantCellFraction"] <- "TumorFractionASCAT"
colnames(compare_tech)[colnames(compare_tech) == "Ploidy"] <- "ASCAT_Ploidy"
colnames(compare_tech)[colnames(compare_tech) == "Diff_Tumor_HandE"] <- "Diff_TumorichorCNA_HandE"
colnames(compare_tech)[colnames(compare_tech) == "Diff_Aberrant_HandE"] <- "Diff_TumorASCAT_HandE"

# Reordering the columns
compare_tech <- compare_tech[, c(2:11, 1, 12:ncol(compare_tech))]
compare_tech$Perc_Deviation_TumorichorCNA <- (compare_tech$Diff_TumorichorCNA_HandE / compare_tech$HandE) * 100
compare_tech$Perc_Deviation_TumorASCAT <- (compare_tech$Diff_TumorASCAT_HandE / compare_tech$HandE) * 100

# Summary of percentage deviations
summary(compare_tech[, c("Perc_Deviation_TumorichorCNA", "Perc_Deviation_TumorASCAT")])

# Subset the dataframe to keep only the relevant columns
subset_df <- compare_tech[, c("TumorFractionASCAT", "TumorFractionIchorCNA", "HandE")]

# Descriptive Statistics
desc_stats <- data.frame(
  Mean = colMeans(subset_df, na.rm = TRUE),
  Median = apply(subset_df, 2, median, na.rm = TRUE),
  Std_Dev = apply(subset_df, 2, sd, na.rm = TRUE)
)
print(desc_stats)

library(ggplot2)

# Scatter plot comparing TumorFractionASCAT vs. TumorFractionIchorCNA
ggplot(subset_df, aes(x = TumorFractionASCAT, y = TumorFractionIchorCNA)) +
  geom_point() +
  labs(title = "Comparison between TumorFractionASCAT and TumorFractionIchorCNA",
       x = "Tumor Fraction ASCAT", y = "Tumor Fraction IchorCNA") +
  geom_smooth(method = "lm", se = FALSE, color = "red")

# Scatter plot comparing TumorFractionASCAT vs. HandE
ggplot(subset_df, aes(x = TumorFractionASCAT, y = HandE)) +
  geom_point() +
  labs(title = "Comparison between TumorFractionASCAT and HandE",
       x = "Tumor Fraction ASCAT", y = "HandE") +
  geom_smooth(method = "lm", se = FALSE, color = "red")

# Scatter plot comparing TumorFractionIchorCNA vs. HandE
ggplot(subset_df, aes(x = TumorFractionIchorCNA, y = HandE)) +
  geom_point() +
  labs(title = "Comparison between TumorFractionIchorCNA and HandE",
       x = "Tumor Fraction IchorCNA", y = "HandE") +
  geom_smooth(method = "lm", se = FALSE, color = "red")

# Compute correlation coefficients
cor_ASCAT_Ichor <- cor(subset_df$TumorFractionASCAT, subset_df$TumorFractionIchorCNA, method = "spearman", use = "complete.obs")
cor_ASCAT_HandE <- cor(subset_df$TumorFractionASCAT, subset_df$HandE, method = "spearman", use = "complete.obs")
cor_Ichor_HandE <- cor(subset_df$TumorFractionIchorCNA, subset_df$HandE, method = "spearman", use = "complete.obs")

correlations <- data.frame(
  Pair = c("ASCAT vs Ichor", "ASCAT vs HandE", "Ichor vs HandE"),
  Correlation = c(cor_ASCAT_Ichor, cor_ASCAT_HandE, cor_Ichor_HandE)
)
print(correlations)


####
library(ggplot2)

# Melt the data for easier plotting
library(reshape2)
# Melt the data for easier plotting
melted_data <- melt(subset_df, id.vars = NULL)

# Boxplots with outliers as dots
ggplot(melted_data, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA) +  # this removes the default outliers
  geom_jitter(width = 0.2, shape = 16, color = "red", size = 2) +  # this adds custom jittered points for all data, highlighting potential outliers
  labs(title = "Comparison of Tumor Fraction Measurements", 
       x = "Method", y = "Tumor Fraction Value") +
  theme_minimal()


# T-test comparisons
t_test_ASCAT_Ichor <- t.test(subset_df$TumorFractionASCAT, subset_df$TumorFractionIchorCNA)
t_test_ASCAT_HandE <- t.test(subset_df$TumorFractionASCAT, subset_df$HandE)
t_test_Ichor_HandE <- t.test(subset_df$TumorFractionIchorCNA, subset_df$HandE)

# Print results
print(t_test_ASCAT_Ichor)
print(t_test_ASCAT_HandE)
print(t_test_Ichor_HandE)


