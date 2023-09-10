# MERGESEGMENT LIST
# 5. Merge adjacent segments of similar copy number and classify events by size (large-scale or focal)
# You can download the mergeSegements.pl utility script from the VarScan scripts folder on SourceForge. 

arm<- read.delim("/Users/mateopava/Desktop/arm.txt", "\t", header = FALSE)

arm <- read.table("/Users/mateopava/Desktop/arm.txt", sep = " ", header = FALSE)
colnames(arm) <- c("chromosome", "position_start", "position_end", "arm_name")
arm$chromosome <- paste("chr", arm$chromosome, sep = "")

write.table(arm, "arm_updated.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

# maybe use one example withouth the re center############ THIS ONLY TO SEE BEFORE
beforesegment153<- segments.p(PRO153seg)
beforesegment153<- beforesegment153 %>%
  filter(pval < 0.05)
colnames(beforesegment153)[3] <- "chr_start"
colnames(beforesegment153)[4] <- "chr_stop"
         
beforesegment153$sample <- beforesegment153$ID
beforesegment153 <- beforesegment153 %>%
select(ID, sample, everything())
#create the file to pass to the perl script         
write.table(beforesegment153, file = "segment153b", row.names = FALSE)
#read the file after the perl script
joinseg153summaryb<-read.table("/Users/mateopava/Desktop/153segmentjoinbefore.summary.txt", "\t", header = TRUE)
joinseg153b<- read.csv("/Users/mateopava/Desktop/153segmentjoinbefore.events.tsv", "\t", header = TRUE)


###read files for PRO018.2
segments18.2<-segments.p(try18.2seg)
segments18.2filter <- segments18.2 %>%
  filter(pval < 0.05) 
colnames(segments18.2filter)[3] <- "chr_start"
colnames(segments18.2filter)[4] <- "chr_stop"

segments18.2filter$sample <- segments18.2filter$ID
segments18.2filter <- segments18.2filter %>%
  select(ID, sample, everything())

write.table(segments18.2filter, file = "segment18tf", row.names = FALSE)
##read
joinseg18.2summary<-read.table("/Users/mateopava/Desktop/18tfsegmentjoin.summary.txt", "\t", header = TRUE)
joinseg18.2<- read.csv("/Users/mateopava/Desktop/18tfsegmentjoin.events.tsv", "\t", header = TRUE)

# read files for 131
segments131<-segments.p(try131seg)
segments131filter <- segments131 %>%
  filter(pval < 0.05) 
colnames(segments131filter)[3] <- "chr_start"
colnames(segments131filter)[4] <- "chr_stop"

segments131filter$sample <- segments131filter$ID
segments131filter <- segments131filter %>%
  select(ID, sample, everything())

write.table(segments131filter, file = "segment131f", row.names = FALSE)
##read
joinseg131summary<-read.table("/Users/mateopava/Desktop/131segmentjoin.summary.txt", "\t", header = TRUE)
joinseg131<- read.csv("/Users/mateopava/Desktop/131segmentjoin.events.tsv", "\t", header = TRUE)

####read files for 179.2
segments179.2<-segments.p(try179.2seg)
segments179.2filter <- segments179.2 %>%
  filter(pval < 0.05) 
colnames(segments179.2filter)[3] <- "chr_start"
colnames(segments179.2filter)[4] <- "chr_stop"

segments179.2filter$sample <- segments179.2filter$ID
segments179.2filter <- segments179.2filter %>%
  select(ID, sample, everything())

write.table(segments179.2filter, file = "segment179tf", row.names = FALSE)
##read
joinseg179.2summary<-read.table("/Users/mateopava/Desktop/179.2segmentjoin.summary.txt", "\t", header = TRUE)
joinseg179.2<- read.csv("/Users/mateopava/Desktop/179.2segmentjoin.events.tsv", "\t", header = TRUE)

##############LOOK THE CODE FROM HERE########################

#segment.p to obtain pvalue and differente statistics ALL THIS ARE AFTER VARSCAN CENTERED
# using this makes the mergeSegmentpl work
segments153<-segments.p(try153seg)
segments153filter <- segments153 %>%
  filter(pval < 0.05) 
colnames(segments153filter)[3] <- "chr_start"
colnames(segments153filter)[4] <- "chr_stop"

segments153filter$sample <- segments153filter$ID
segments153filter <- segments153filter %>%
  select(ID, sample, everything())


# write.table(segments153, file = "segment153", row.names = FALSE) dont use
write.table(segments153filter, file = "segment153f", row.names = FALSE)

##### 17
segments17<-segments.p(try17seg)
segments17filter <- segments17 %>%
  filter(pval < 0.05) 
colnames(segments17filter)[3] <- "chr_start"
colnames(segments17filter)[4] <- "chr_stop"

segments17filter$sample <- segments17filter$ID
segments17filter <- segments17filter %>%
  select(ID, sample, everything())


# write.table(segments153, file = "segment153", row.names = FALSE) dont use
write.table(segments17filter, file = "segment17f", row.names = FALSE)

segments32<-segments.p(try32seg)
segments32filter <- segments32 %>%
  filter(pval < 0.05) 
colnames(segments32filter)[3] <- "chr_start"
colnames(segments32filter)[4] <- "chr_stop"

segments32filter$sample <- segments32filter$ID
segments32filter <- segments32filter %>%
  select(ID, sample, everything())

write.table(segments32filter, file = "segment32f", row.names = FALSE)

joinseg32summary<-read.table("/Users/mateopava/Desktop/32segmentjoin.summary.txt", "\t", header = TRUE)
joinseg32<- read.csv("/Users/mateopava/Desktop/32segmentjoin.events.tsv", "\t", header = TRUE)

joinseg153summary<-read.table("/Users/mateopava/Desktop/153segmentjoin.summary.txt", "\t", header = TRUE)
joinseg153<- read.csv("/Users/mateopava/Desktop/153segmentjoin.events.tsv", "\t", header = TRUE)

joinseg32<- joinseg32 %>%
  filter(p_value < 0.05)
# for the before

###############THIS IS THE GENE ANNOTATION FOR THE MERGESEGMENTS AFTER THE MERGE SCRIPT use 153 as first example
joinseg153$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg153))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg153)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg153[j, 2], joinseg153[j, 3]) &
      filtered_genes153_all[, 5] == joinseg153[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg153$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg153gene <- na.omit(joinseg153)
########## 17
joinseg17segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg17segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg17segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg17segmentjoin.events.tsv[j, 2], joinseg17segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg17segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg17segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg17segmentjoin.events.tsv <- na.omit(joinseg17segmentjoin.events.tsv)


########## 18
joinseg18segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg18segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg18segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg18segmentjoin.events.tsv[j, 2], joinseg18segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg18segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg18segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg18segmentjoin.events.tsv <- na.omit(joinseg18segmentjoin.events.tsv)

############# 18.2
joinseg18.2$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg18.2))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg18.2)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg18.2[j, 2], joinseg18.2[j, 3]) &
      filtered_genes153_all[, 5] == joinseg18.2[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg18.2$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg18.2 <- na.omit(joinseg18.2)

############# 28
joinseg28segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg28segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg28segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg28segmentjoin.events.tsv[j, 2], joinseg28segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg28segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg28segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg28segmentjoin.events.tsv <- na.omit(joinseg28segmentjoin.events.tsv)
############# 35

joinseg35segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg35segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg35segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg35segmentjoin.events.tsv[j, 2], joinseg35segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg35segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg35segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg35segmentjoin.events.tsv <- na.omit(joinseg35segmentjoin.events.tsv)


############# 43

joinseg43segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg43segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg43segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg43segmentjoin.events.tsv[j, 2], joinseg43segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg43segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg43segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg43segmentjoin.events.tsv <- na.omit(joinseg43segmentjoin.events.tsv)


############# 55

joinseg55segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg55segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg55segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg55segmentjoin.events.tsv[j, 2], joinseg55segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg55segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg55segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg55segmentjoin.events.tsv <- na.omit(joinseg55segmentjoin.events.tsv)


############# 57

joinseg57segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg57segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg57segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg57segmentjoin.events.tsv[j, 2], joinseg57segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg57segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg57segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg57segmentjoin.events.tsv <- na.omit(joinseg57segmentjoin.events.tsv)


############# 59

joinseg59segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg59segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg59segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg59segmentjoin.events.tsv[j, 2], joinseg59segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg59segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg59segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg59segmentjoin.events.tsv <- na.omit(joinseg59segmentjoin.events.tsv)


############# 68

joinseg68segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg68segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg68segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg68segmentjoin.events.tsv[j, 2], joinseg68segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg68segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg68segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg68segmentjoin.events.tsv <- na.omit(joinseg68segmentjoin.events.tsv)


############# 81

joinseg81segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg81segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg81segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg81segmentjoin.events.tsv[j, 2], joinseg81segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg81segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg81segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg81segmentjoin.events.tsv <- na.omit(joinseg81segmentjoin.events.tsv)


############# 111

joinseg111segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg111segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg111segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg111segmentjoin.events.tsv[j, 2], joinseg111segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg111segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg111segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg111segmentjoin.events.tsv <- na.omit(joinseg111segmentjoin.events.tsv)


############# 113

joinseg113segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg113segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg113segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg113segmentjoin.events.tsv[j, 2], joinseg113segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg113segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg113segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg113segmentjoin.events.tsv <- na.omit(joinseg113segmentjoin.events.tsv)


############# 131 

joinseg131$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg131))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg131)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg131[j, 2], joinseg131[j, 3]) &
      filtered_genes153_all[, 5] == joinseg131[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg131$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg131 <- na.omit(joinseg131)

############# 149

joinseg149segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg149segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg149segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg149segmentjoin.events.tsv[j, 2], joinseg149segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg149segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg149segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg149segmentjoin.events.tsv <- na.omit(joinseg149segmentjoin.events.tsv)


############# 151

joinseg151segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg151segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg151segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg151segmentjoin.events.tsv[j, 2], joinseg151segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg151segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg151segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg151segmentjoin.events.tsv <- na.omit(joinseg151segmentjoin.events.tsv)


############# 157

joinseg157segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg157segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg157segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg157segmentjoin.events.tsv[j, 2], joinseg157segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg157segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg157segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg157segmentjoin.events.tsv <- na.omit(joinseg157segmentjoin.events.tsv)


############# 179

joinseg179segmentjoin.events.tsv$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg179segmentjoin.events.tsv))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg179segmentjoin.events.tsv)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg179segmentjoin.events.tsv[j, 2], joinseg179segmentjoin.events.tsv[j, 3]) &
      filtered_genes153_all[, 5] == joinseg179segmentjoin.events.tsv[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg179segmentjoin.events.tsv$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg179segmentjoin.events.tsv <- na.omit(joinseg179segmentjoin.events.tsv)

######### 179.2

joinseg179.2$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(joinseg179.2))  # create a list to store matching genes for each row

for (j in 1:nrow(joinseg179.2)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], joinseg179.2[j, 2], joinseg179.2[j, 3]) &
      filtered_genes153_all[, 5] == joinseg179.2[j, 1]
  ]
}

# If there are matching genes, paste them together separated by comma
joinseg179.2$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

joinseg179.2 <- na.omit(joinseg179.2)


# try172out_expanded <- try172out %>% #use this if neccesary to split
# separate_rows(genes, sep = ",")


