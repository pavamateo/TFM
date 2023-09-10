# possible plot using PRO153 after merge segmentation
joinseg153geneC <- joinseg153gene
joinseg153geneC <- joinseg153geneC %>% 
separate_rows(genes, sep = ",")
colnames(joinseg153geneC)[14] <- "genes"
# C is with the panel
panel153C<- left_join(`L21-3701_cna`, joinseg153geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel153A <- panel153C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]


#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data <- panel153A[!is.na(panel153A$log2) & !is.na(panel153A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long <- gene_data %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long$variable <- factor(gene_data_long$variable, 
                                  levels = c("log2", "seg_mean"), 
                                  labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO153 Comparison Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data <- gene_data %>% distinct(genes, .keep_all = TRUE)
correlation153 <- cor.test(unique_gene_data$log2, unique_gene_data$seg_mean)
# see correlation with a plot
plot(unique_gene_data$log2, unique_gene_data$seg_mean)
# cor.test(unique_gene_data$log2, unique_gene_data$seg_mean, method = "spearman") gives 68% because there are repeated values for genes

##### NEED THE PLOTS FOR THE FOLLOWINGS IDS
#[1] "PRO18"    "PRO28"    "PRO35"    "PRO43"    "PRO55"   
#[6] "PRO57"    "PRO59"    "PRO68"    "PRO81"    "PRO111"  
#[11] "PRO113"   "PRO149.2" "PRO151"   "PRO153"   "PRO157"  
#[16] "PRO179"  
# out of this
################################# 17 	L19-0565
joinseg17geneC <- joinseg17segmentjoin.events.tsv
joinseg17geneC <- joinseg17geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg17geneC)[14] <- "genes"
# C is with the panel
panel17C<- left_join(`L19-0565_cna`, joinseg17geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel17A <- panel17C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data17 <- panel17A[!is.na(panel17A$log2) & !is.na(panel17A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long17 <- gene_data17 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long17$variable <- factor(gene_data_long17$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long17, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO017 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data17 <- gene_data17 %>% distinct(genes, .keep_all = TRUE)
correlation17 <- cor.test(unique_gene_data17$log2, unique_gene_data17$seg_mean) #cant correlate with 1 value

################################## 18 mutant
joinseg18geneC <- joinseg18segmentjoin.events.tsv
joinseg18geneC <- joinseg18geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg18geneC)[14] <- "genes"
# C is with the panel
panel18C<- left_join(`L21-2508_cna`, joinseg153geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel18A <- panel18C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data18 <- panel18A[!is.na(panel18A$log2) & !is.na(panel18A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long18 <- gene_data18 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long18$variable <- factor(gene_data_long18$variable, 
                                  levels = c("log2", "seg_mean"), 
                                  labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long18, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO018 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data18 <- gene_data18 %>% distinct(genes, .keep_all = TRUE)
correlation18 <- cor.test(unique_gene_data18$log2, unique_gene_data18$seg_mean)
# see correlation with a plot
plot(unique_gene_data18$log2, unique_gene_data18$seg_mean)



################################## 18.2 mutant
joinseg18.2geneC <- joinseg18.2
joinseg18.2geneC <- joinseg18.2geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg18.2geneC)[14] <- "genes"
# C is with the panel
panel18.2C<- left_join(`L20-4110_cna`, joinseg18.2geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel18.2A <- panel18.2C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data18.2 <- panel18.2A[!is.na(panel18.2A$log2) & !is.na(panel18.2A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long18.2 <- gene_data18.2 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long18.2$variable <- factor(gene_data_long18.2$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long18.2, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO018.2 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data18.2 <- gene_data18.2 %>% distinct(genes, .keep_all = TRUE)
plot(unique_gene_data18.2$log2, unique_gene_data18.2$seg_mean)

############# 28 L20-3102 mutant
joinseg28geneC <- joinseg28segmentjoin.events.tsv
joinseg28geneC <- joinseg28geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg28geneC)[14] <- "genes"
# C is with the panel
panel28C<- left_join(`L20-3102_cna`, joinseg28geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel28A <- panel28C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data28 <- panel28A[!is.na(panel28A$log2) & !is.na(panel28A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long28 <- gene_data28 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long28$variable <- factor(gene_data_long28$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long28, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO018 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data28 <- gene_data28 %>% distinct(genes, .keep_all = TRUE)
correlation28 <- cor.test(unique_gene_data28$log2, unique_gene_data28$seg_mean) #cant correlate with 1 value
############# 35 L21-2751 mutant

joinseg35geneC <- joinseg35segmentjoin.events.tsv
joinseg35geneC <- joinseg35geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg35geneC)[14] <- "genes"
# C is with the panel
panel35C<- left_join(`L21-2751_cna`, joinseg35geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel35A <- panel35C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data35 <- panel35A[!is.na(panel35A$log2) & !is.na(panel35A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long35 <- gene_data35 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long35$variable <- factor(gene_data_long35$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long35, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO035 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data35 <- gene_data35 %>% distinct(genes, .keep_all = TRUE)
correlation35 <- cor.test(unique_gene_data35$log2, unique_gene_data35$seg_mean) #not enough
############# 43 L19-0648 WT
joinseg43geneC <- joinseg43segmentjoin.events.tsv
joinseg43geneC <- joinseg43geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg43geneC)[14] <- "genes"
# C is with the panel
panel43C<- left_join(`L19-0648_cna`, joinseg43geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel43A <- panel43C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data43 <- panel43A[!is.na(panel43A$log2) & !is.na(panel43A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long43 <- gene_data43 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long43$variable <- factor(gene_data_long43$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long43, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO043 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data43 <- gene_data43 %>% distinct(genes, .keep_all = TRUE)
correlation43 <- cor.test(unique_gene_data43$log2, unique_gene_data43$seg_mean) #not enough
############# 55 L19-1139 mutant
joinseg55geneC <- joinseg55segmentjoin.events.tsv
joinseg55geneC <- joinseg55geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg55geneC)[14] <- "genes"
# C is with the panel
panel55C<- left_join(`L19-1139_cna`, joinseg55geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel55A <- panel55C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data55 <- panel55A[!is.na(panel55A$log2) & !is.na(panel55A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long55 <- gene_data55 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long55$variable <- factor(gene_data_long55$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long55, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO055 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data55 <- gene_data55 %>% distinct(genes, .keep_all = TRUE)
correlation55 <- cor.test(unique_gene_data55$log2, unique_gene_data55$seg_mean) #cant
############# 57 L19-1367 mutant
joinseg57geneC <- joinseg57segmentjoin.events.tsv
joinseg57geneC <- joinseg57geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg57geneC)[14] <- "genes"
# C is with the panel
panel57C<- left_join(`L19-1367_cna`, joinseg57geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel57A <- panel57C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data57 <- panel57A[!is.na(panel57A$log2) & !is.na(panel57A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long57 <- gene_data57 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long57$variable <- factor(gene_data_long57$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long57, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO057 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data57 <- gene_data57 %>% distinct(genes, .keep_all = TRUE)
correlation57 <- cor.test(unique_gene_data57$log2, unique_gene_data57$seg_mean)
############# 59 L21-3688 mutant
joinseg59geneC <- joinseg59segmentjoin.events.tsv
joinseg59geneC <- joinseg59geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg59geneC)[14] <- "genes"
# C is with the panel
panel59C<- left_join(`L21-3688_cna`, joinseg59geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel59A <- panel59C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data59 <- panel59A[!is.na(panel59A$log2) & !is.na(panel59A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long59 <- gene_data59 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long59$variable <- factor(gene_data_long59$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long59, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data59 <- gene_data59 %>% distinct(genes, .keep_all = TRUE)
correlation59 <- cor.test(unique_gene_data59$log2, unique_gene_data59$seg_mean)
############# 68 L19-3826 Mutant
joinseg68geneC <- joinseg68segmentjoin.events.tsv
joinseg68geneC <- joinseg68geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg68geneC)[14] <- "genes"
# C is with the panel
panel68C<- left_join(`L19-3826_cna`, joinseg68geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel68A <- panel68C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data68 <- panel68A[!is.na(panel68A$log2) & !is.na(panel68A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long68 <- gene_data68 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long68$variable <- factor(gene_data_long68$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long68, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO068 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data68 <- gene_data68 %>% distinct(genes, .keep_all = TRUE)
correlation68 <- cor.test(unique_gene_data68$log2, unique_gene_data68$seg_mean) #cant
############# 81 L20-0350 Mutant
joinseg81geneC <- joinseg81segmentjoin.events.tsv
joinseg81geneC <- joinseg81geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg81geneC)[14] <- "genes"
# C is with the panel
panel81C<- left_join(`L20-0350_cna`, joinseg81geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel81A <- panel81C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data81 <- panel81A[!is.na(panel81A$log2) & !is.na(panel81A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long81 <- gene_data81 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long81$variable <- factor(gene_data_long81$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long81, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO081 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data81 <- gene_data81 %>% distinct(genes, .keep_all = TRUE)
correlation81 <- cor.test(unique_gene_data81$log2, unique_gene_data81$seg_mean)
############# 111 L20-4112 Mutant
joinseg111geneC <- joinseg111segmentjoin.events.tsv
joinseg111geneC <- joinseg111geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg111geneC)[14] <- "genes"
# C is with the panel
panel111C<- left_join(`L20-4112_cna`, joinseg111geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel111A <- panel111C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data111 <- panel111A[!is.na(panel111A$log2) & !is.na(panel111A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long111 <- gene_data111 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long111$variable <- factor(gene_data_long111$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long111, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO111 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data111 <- gene_data111 %>% distinct(genes, .keep_all = TRUE)
correlation111 <- cor.test(unique_gene_data111$log2, unique_gene_data111$seg_mean)
############# 113 L20-6010 mutant
joinseg113geneC <- joinseg113segmentjoin.events.tsv
joinseg113geneC <- joinseg113geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg113geneC)[14] <- "genes"
# C is with the panel
panel113C<- left_join(`L20-6010_cna`, joinseg113geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel113A <- panel113C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data113 <- panel113A[!is.na(panel113A$log2) & !is.na(panel113A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long113 <- gene_data113 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long113$variable <- factor(gene_data_long113$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long113, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO113 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data113 <- gene_data113 %>% distinct(genes, .keep_all = TRUE)
correlation113 <- cor.test(unique_gene_data113$log2, unique_gene_data113$seg_mean)
############ 131 L21-8381
joinseg131geneC <- joinseg131
joinseg131geneC <- joinseg131geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg131geneC)[14] <- "genes"
# C is with the panel
panel131C<- left_join(`L21-8381_cna`, joinseg131geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel131A <- panel131C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data131 <- panel131A[!is.na(panel131A$log2) & !is.na(panel131A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long131 <- gene_data131 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long131$variable <- factor(gene_data_long131$variable, 
                                      levels = c("log2", "seg_mean"), 
                                      labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long131, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO0131 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data131 <- gene_data131 %>% distinct(genes, .keep_all = TRUE)
correlation131 <- cor.test(unique_gene_data131$log2, unique_gene_data131$seg_mean)


############# 149 L21-1563 Mutant
joinseg149geneC <- joinseg149segmentjoin.events.tsv
joinseg149geneC <- joinseg149geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg149geneC)[14] <- "genes"
# C is with the panel
panel149C<- left_join(`L21-1563_cna`, joinseg149geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel149A <- panel149C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data149 <- panel149A[!is.na(panel149A$log2) & !is.na(panel149A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long149 <- gene_data149 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long149$variable <- factor(gene_data_long149$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long149, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO149 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data149 <- gene_data149 %>% distinct(genes, .keep_all = TRUE)
correlation149 <- cor.test(unique_gene_data149$log2, unique_gene_data149$seg_mean) #cant
############# 151 L21-1735 mutant
joinseg151geneC <- joinseg151segmentjoin.events.tsv
joinseg151geneC <- joinseg151geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg151geneC)[14] <- "genes"
# C is with the panel
panel151C<- left_join(`L21-1735_cna`, joinseg151geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel151A <- panel151C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data151 <- panel151A[!is.na(panel151A$log2) & !is.na(panel151A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long151 <- gene_data151 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long151$variable <- factor(gene_data_long151$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long151, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO151 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data151 <- gene_data151 %>% distinct(genes, .keep_all = TRUE)
correlation151 <- cor.test(unique_gene_data151$log2, unique_gene_data151$seg_mean)
############# 157 L21-2752 mutant
joinseg157geneC <- joinseg157segmentjoin.events.tsv
joinseg157geneC <- joinseg157geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg157geneC)[14] <- "genes"
# C is with the panel
panel157C<- left_join(`L21-2752_cna`, joinseg157geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel157A <- panel157C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data157 <- panel157A[!is.na(panel157A$log2) & !is.na(panel157A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long157 <- gene_data157 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long157$variable <- factor(gene_data_long157$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long157, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO157 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data157 <- gene_data157 %>% distinct(genes, .keep_all = TRUE)
correlation157 <- cor.test(unique_gene_data157$log2, unique_gene_data157$seg_mean) #nothing comparable from wes
############# 179 L21-5659 Mutant

joinseg179geneC <- joinseg179segmentjoin.events.tsv
joinseg179geneC <- joinseg179geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg179geneC)[14] <- "genes"
# C is with the panel
panel179C<- left_join(`L21-5659_cna`, joinseg179geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel179A <- panel179C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data179 <- panel179A[!is.na(panel179A$log2) & !is.na(panel179A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long179 <- gene_data179 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long179$variable <- factor(gene_data_long179$variable, 
                                    levels = c("log2", "seg_mean"), 
                                    labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long179, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO179 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data179 <- gene_data179 %>% distinct(genes, .keep_all = TRUE)
correlation179 <- cor.test(unique_gene_data179$log2, unique_gene_data179$seg_mean)

########### 179.2 L21-7137
joinseg179.2geneC <- joinseg179.2
joinseg179.2geneC <- joinseg179.2geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg179.2geneC)[14] <- "genes"
# C is with the panel
panel179.2C<- left_join(`L21-7137_cna`, joinseg179.2geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel179.2A <- panel179.2C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]

#create bar plot
# Remove rows with missing log2 or segment mean values
gene_data179.2 <- panel179.2A[!is.na(panel179.2A$log2) & !is.na(panel179.2A$seg_mean), ]

# reshape the data into a 'long' format
gene_data_long179.2 <- gene_data179.2 %>% 
  pivot_longer(cols = c("log2", "seg_mean"), 
               names_to = "variable", 
               values_to = "value")

# rename the levels of variable for a proper legend
gene_data_long179.2$variable <- factor(gene_data_long179.2$variable, 
                                     levels = c("log2", "seg_mean"), 
                                     labels = c("Log2 panel", "Log2 WES"))

# plot
ggplot(gene_data_long179.2, aes(x = genes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(x = "Gene", y = "Log2 value", title = "PRO179.2 Comparison of Log2 values for Each Gene") +
  scale_fill_manual(values = c("Log2 panel" = "steelblue", "Log2 WES" = "darkorange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "top")

# perform spearman correlation
unique_gene_data179.2 <- gene_data179.2 %>% distinct(genes, .keep_all = TRUE)
correlation179.2 <- cor.test(unique_gene_data179.2$log2, unique_gene_data179.2$seg_mean)


################### DO IT FOR PRO32 ######## Since PRO32 is no significant for the panel, we need to plot regardless
joinseg32$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg32))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg32)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg32[j, 2], joinseg32[j, 3]) &
      filtered_genes153_all[, 5] == joinseg32[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg32$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
# apply NA omit to get rid of the genomic regions without an annotated gene
joinseg32gene <- na.omit(joinseg32)
############## THE ORGANIZATION IS DONE TO COMPARE, IF THERE IS NOTHING TO COMPARE DONT DO IT
#now organize data
joinseg32geneC <- joinseg32gene
#this separates the genes that are joined in a same genomic region
joinseg32geneC <- joinseg32geneC %>% 
  separate_rows(genes, sep = ",")
colnames(joinseg32geneC)[14] <- "genes"
# C is with the panel
panel32C<- left_join(`L19-3813_cna`, joinseg32geneC, by = "genes")
# A is with the panel and only necessary columns according to me
panel32A <- panel32C[, c(1, 2, 3, 4, 5, 10, 11, 12,13, 17:22)]
############ WONT WORK WITH THE NON SIGNIFICATIVE