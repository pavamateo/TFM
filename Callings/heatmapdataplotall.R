library(dplyr)
library(tidyr)
library(ggplot2)

combined_melted<- read.csv("/Users/mateopava/combined_melted.csv",sep = "\t")

#plot with all technologies
p <- ggplot(combined_melted, aes(x = Var2 , y = Var1, fill = Alteration)) +
  geom_tile() +
  scale_fill_manual(values = c("NEUTRAL" = "lightgreen", 
                               "GAIN" = "pink", 
                               "LOSS" = "lightblue1", 
                               "AMP" = "red", 
                               "DEEPDEL" = "blue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# for putting in 4 different plots
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

####

p +  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####### tests
#y names
yNames <- NULL
for (i in 1:48){
  yNames=c(yNames,"","",samples[i],"","","")
}

p1 <- p + scale_x_discrete(labels=yNames) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1
