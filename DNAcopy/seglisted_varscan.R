#this will be the script for all the recentering
#first normalize the log2 ratio. adjusted or raw? use the raw as per varscan
#### When value is - use recenter down
#### When value is + use recenter up
scale(PRO016$adjusted_log_ratio)
try16 <- read.csv("/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/example16",sep = "\t")
### el varscan quita mucha info solo para la 16

try131 <- read.csv("/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/example131",sep = "\t")
colnames(try131)[1] <- "chr"
colnames(try131)[2] <- "start"
colnames(try131)[3] <- "end"
try131CNA.object <- CNA(genomdat = try131$adjusted_log_ratio, chrom = try131$chr, maploc = try131$start, data.type = 'logratio', sampleid = "PRO131")
try131CNA.smoothed <- smooth.CNA(try131CNA.object)
try131seg <- segment(try131CNA.smoothed, verbose = 1) #verbose 1 tells me which sample is being analyzed

try17 <- read.csv("/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/example17",sep = "\t")
colnames(try17)[1] <- "chr"
colnames(try17)[2] <- "start"
colnames(try17)[3] <- "end"
try17CNA.object <- CNA(genomdat = try17$adjusted_log_ratio, chrom = try17$chr, maploc = try17$start, data.type = 'logratio', sampleid = "PRO17")
try17CNA.smoothed <- smooth.CNA(try17CNA.object)
try17seg <- segment(try17CNA.smoothed, verbose = 1) #verbose 1 tells me which sample is being analyzed

plot(try17seg, plot.type = "w")
try17out <- try17seg$output
try17out$state <- ifelse(try17out$seg.mean > 0.25, "gain", 
                         ifelse(try17out$seg.mean < -0.25, "loss", "no_change"))

try17out$genes <- NA  # initialize an empty vector
for (j in 1:nrow(try17out)) {
  
  for (g in 1:nrow(filtered_genes153_all)) {
    if (between(filtered_genes153_all[g,3], try17out[j,3], try17out[j,4]) &
        filtered_genes153_all[g, 5] == try17out[j,2]) {
      try17out[j,8] <- filtered_genes153_all[g,2]  # assign gene symbol
      break  # end the inner loop
    }
  }
}
##
try18 <- read.csv("/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/example18", sep = "\t")
scale(PRO018$adjusted_log_ratio)
colnames(try18)[1] <- "chr"
colnames(try18)[2] <- "start"
colnames(try18)[3] <- "end"
try18CNA.object <- CNA(genomdat = try18$adjusted_log_ratio, chrom = try18$chr, maploc = try18$start, data.type = 'logratio', sampleid = "PRO18")
try18CNA.smoothed <- smooth.CNA(try18CNA.object)
try18seg <- segment(try18CNA.smoothed, verbose = 1) #verbose 1 tells me which sample is being analyzed

plot(try18seg, plot.type = "w")
try18out <- try18seg$output
try18out$state <- ifelse(try18out$seg.mean > 0.25, "gain",
                         ifelse(try18out$seg.mean < -0.25, "loss", "no_change"))

try18out$genes <- NA  # initialize an empty vector
for (j in 1:nrow(try18out)) {
  
  for (g in 1:nrow(filtered_genes153_all)) {
    if (between(filtered_genes153_all[g, 3], try18out[j, 3], try18out[j, 4]) &
        filtered_genes153_all[g, 5] == try18out[j, 2]) {
      try18out[j, 8] <- filtered_genes153_all[g, 2]  # assign gene symbol
      break  # end the inner loop
    }
  }
}

try18.2 <- read.csv("/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/example18.2", sep = "\t")
scale(PRO018.2$adjusted_log_ratio)
colnames(try18.2)[1] <- "chr"
colnames(try18.2)[2] <- "start"
colnames(try18.2)[3] <- "end"
try18.2CNA.object <- CNA(genomdat = try18.2$adjusted_log_ratio, chrom = try18.2$chr, maploc = try18.2$start, data.type = 'logratio', sampleid = "PRO18")
try18.2CNA.smoothed <- smooth.CNA(try18.2CNA.object)
try18.2seg <- segment(try18.2CNA.smoothed, verbose = 1) #verbose 1 tells me which sample is being analyzed

plot(try18.2seg, plot.type = "w")
try18.2out <- try18.2seg$output
try18.2out$state <- ifelse(try18.2out$seg.mean > 0.25, "gain",
                         ifelse(try18.2out$seg.mean < -0.25, "loss", "no_change"))

try18.2out$genes <- NA  # initialize an empty vector
for (j in 1:nrow(try18.2out)) {
  
  for (g in 1:nrow(filtered_genes153_all)) {
    if (between(filtered_genes153_all[g, 3], try18.2out[j, 3], try18.2out[j, 4]) &
        filtered_genes153_all[g, 5] == try18.2out[j, 2]) {
      try18.2out[j, 8] <- filtered_genes153_all[g, 2]  # assign gene symbol
      break  # end the inner loop
    }
  }
}

try179.2<- read.csv("/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/example179.2", sep = "\t")
scale(PRO179.2$adjusted_log_ratio)
colnames(try179.2)[1] <- "chr"
colnames(try179.2)[2] <- "start"
colnames(try179.2)[3] <- "end"
try179.2CNA.object <- CNA(genomdat = try179.2$adjusted_log_ratio, chrom = try179.2$chr, maploc = try179.2$start, data.type = 'logratio', sampleid = "PRO179.2")
try179.2CNA.smoothed <- smooth.CNA(try179.2CNA.object)
try179.2seg <- segment(try179.2CNA.smoothed, verbose = 1) #verbose 1 tells me which sample is being analyzed

plot(try179.2seg, plot.type = "w")

try32 <- read.csv("/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/example32",sep = "\t")
colnames(try32)[1] <- "chr"
colnames(try32)[2] <- "start"
colnames(try32)[3] <- "end"
try32CNA.object <- CNA(genomdat = try32$adjusted_log_ratio, chrom = try32$chr, maploc = try32$start, data.type = 'logratio', sampleid = "PRO32")
try32CNA.smoothed <- smooth.CNA(try32CNA.object)
try32seg <- segment(try32CNA.smoothed, verbose = 1) #verbose 1 tells me which sample is being analyzed

plot(try32seg, plot.type = "w")
try32out <- try32seg$output
try32out$state <- ifelse(try32out$seg.mean > 0.25, "gain", 
                         ifelse(try32out$seg.mean < -0.25, "loss", "no_change"))

try32out$genes <- NA  # initialize an empty vector
for (j in 1:nrow(try32out)) {
  
  for (g in 1:nrow(filtered_genes153_all)) {
    if (between(filtered_genes153_all[g,3], try32out[j,3], try32out[j,4]) &
        filtered_genes153_all[g, 5] == try32out[j,2]) {
      try32out[j,8] <- filtered_genes153_all[g,2]  # assign gene symbol
      break  # end the inner loop
    }
  }
}

try153 <- read.csv("/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/example153",sep = "\t")
colnames(try153)[1] <- "chr"
colnames(try153)[2] <- "start"
colnames(try153)[3] <- "end"
try153CNA.object <- CNA(genomdat = try153$adjusted_log_ratio, chrom = try153$chr, maploc = try153$start, data.type = 'logratio', sampleid = "PRO153")
try153CNA.smoothed <- smooth.CNA(try153CNA.object)
try153seg <- segment(try153CNA.smoothed, verbose = 1) #verbose 1 tells me which sample is being analyzed

plot(try153seg, plot.type = "w")
try153out <- try153seg$output
try153out$state <- ifelse(try153out$seg.mean > 0.25, "gain", 
                         ifelse(try153out$seg.mean < -0.25, "loss", "no_change"))

try153out$genes <- NA  # initialize an empty vector
for (j in 1:nrow(try153out)) {
  
  for (g in 1:nrow(filtered_genes153_all)) {
    if (between(filtered_genes153_all[g,3], try153out[j,3], try153out[j,4]) &
        filtered_genes153_all[g, 5] == try153out[j,2]) {
      try153out[j,8] <- filtered_genes153_all[g,2]  # assign gene symbol
      break  # end the inner loop
    }
  }
}
#########trying to loop for the new DNAcopy

file_path <- "/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/"

# Create empty lists to store the objects
seg_list <- list()
#using numers 17:206 because they are the files for patients
for (i in 17:206) {
  file_name <- paste0("example", i)
  try_name <- paste0("try", i)
  
  # Use tryCatch to handle the error if the file is missing
  tryCatch(
    {
      # Read the file
      file <- read.csv(paste0(file_path, file_name), sep = "\t")
      
      # Modify column names
      colnames(file)[1] <- "chr"
      colnames(file)[2] <- "start"
      colnames(file)[3] <- "end"
      
      # Create CNA object
      cna_object <- CNA(genomdat = file$adjusted_log_ratio, chrom = file$chr, maploc = file$start, data.type = 'logratio', sampleid = paste0("PRO", i))
      
      # Smooth the CNA object
      smoothed_cna <- smooth.CNA(cna_object)
      
      # Segment the smoothed CNA
      seg <- segment(smoothed_cna, verbose = 1)
      
      # Plot the segment
      plot(seg, plot.type = "w")
      
      # Save the objects to the respective lists
      assign(try_name, file)
      seg_list[[i]] <- seg
    },
    error = function(e) {
      message(paste("No file for", file_name, "- Skipping this file."))
    }
  )
}

##try 2
seg_list <- list()
file_path<- "~/Desktop/CNA_Mateo_MSC/allWES/"
# Using numbers 17:206 because they are the numbers of my files.
# I know my lowest file is 17
# and my largest file is 207. Adjust accordingly
for (i in 17:206) {
  file_name <- paste0("example", i)
  try_name <- paste0("try", i)
  
  # Use tryCatch to handle the error if the file is missing
  tryCatch(
    {
      # Read the file
      file <- read.csv(paste0(file_path, file_name), sep = "\t")
      
      # Modify column names
      colnames(file)[1] <- "chr"
      colnames(file)[2] <- "start"
      colnames(file)[3] <- "end"
      
      # Order the chr column based on chr_order
      file$chr <- factor(file$chr, levels = chr_order)
      
      # Create CNA object
      cna_object <- CNA(genomdat = file$adjusted_log_ratio, chrom = file$chr, maploc = file$start, data.type = 'logratio', sampleid = paste0("PRO", i))
      
      # Smooth the CNA object
      smoothed_cna <- smooth.CNA(cna_object)
      
      # Segment the smoothed CNA
      seg <- segment(smoothed_cna, verbose = 1)
      
      # Save the objects to the respective lists
      assign(try_name, file)
      seg_list[[i]] <- seg
    },
    error = function(e) {
      message(paste("No file for", file_name, "- Skipping this file."))
    }
  )
}


try18out<- try18seg$output
try17seg <- seg_list[[17]]
try18seg <-seg_list[[18]]
try28seg <-seg_list[[28]]
try32seg <-seg_list[[32]]
try34seg <- seg_list[[34]]
try35seg <- seg_list[[35]]
try40seg <- seg_list[[40]]
try43seg <- seg_list[[43]]
try44seg <- seg_list[[44]]
try55seg <- seg_list[[55]]
try57seg <- seg_list[[57]]
try58seg <- seg_list[[58]]
try59seg <- seg_list[[59]]
try60seg <- seg_list[[60]]
try68seg <- seg_list[[68]]
try81seg <- seg_list[[81]]
try85seg <- seg_list[[85]]
try93seg <- seg_list[[93]]
try111seg <- seg_list[[111]]
try112seg <- seg_list[[112]]
try113seg <- seg_list[[113]]
try114seg <- seg_list[[114]]
try120seg <- seg_list[[120]]
try122seg <- seg_list[[122]]
try136seg <- seg_list[[136]]
try149seg <- seg_list[[149]]
try151seg <- seg_list[[151]]
try153seg <- seg_list[[153]]
try157seg <- seg_list[[157]]
try172seg <- seg_list[[172]]
try176seg <- seg_list[[176]]
try179seg <- seg_list[[179]]
try185seg <- seg_list[[185]]
try206seg <- seg_list[[206]]

#redo 206 since it was not clear
solo206 <- read.csv("/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/example206",sep = "\t")
colnames(solo206)[1] <- "chr"
colnames(solo206)[2] <- "start"
colnames(solo206)[3] <- "end"
solo206$chr <- factor(solo206$chr, levels = chr_order)

try206_CNA.object <- CNA(genomdat = solo206$adjusted_log_ratio, chrom = solo206$chr, maploc = solo206$start, data.type = 'logratio', sampleid = "PRO206")
try206_CNA.smoothed <- smooth.CNA(try206_CNA.object)
try206_seg <- segment(try206_CNA.smoothed, verbose = 1) #verbose 1 tells me which sample is being analyzed

plot(try206_seg, plot.type = "w")


# GENE ASSIGNMENT and state assignment for the comparable
# this code gives multiple gene matching for a genomic region

#18
try18out <- try18seg$output
try18out$state <- ifelse(try18out$seg.mean > 0.25, "gain", 
                         ifelse(try18out$seg.mean < -0.25, "loss", "no_change"))

try18out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try18out))  # create a list to store matching genes for each row

for (j in 1:nrow(try18out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try18out[j, 3], try18out[j, 4]) &
      filtered_genes153_all[, 5] == try18out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try18out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

# now separate the genes column into multiple rows so i can left join after
try18out_expanded <- try18out %>%
  separate_rows(genes, sep = ",")
# TRY it again with another id eg 28

try28out <- try28seg$output
try28out$state <- ifelse(try28out$seg.mean > 0.25, "gain", 
                         ifelse(try28out$seg.mean < -0.25, "loss", "no_change"))

try28out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try28out))  # create a list to store matching genes for each row

for (j in 1:nrow(try28out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try28out[j, 3], try28out[j, 4]) &
      filtered_genes153_all[, 5] == try28out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try28out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

try28out_expanded <- try28out %>%
  separate_rows(genes, sep = ",")

#32
try32out <- try32seg$output
try32out$state <- ifelse(try32out$seg.mean > 0.25, "gain", 
                         ifelse(try32out$seg.mean < -0.25, "loss", "no_change"))

try32out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try32out))  # create a list to store matching genes for each row

for (j in 1:nrow(try32out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try32out[j, 3], try32out[j, 4]) &
      filtered_genes153_all[, 5] == try32out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try32out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

try32out_expanded <- try32out %>%
  separate_rows(genes, sep = ",")

# For try34out
try34out <- try34seg$output
try34out$state <- ifelse(try34out$seg.mean > 0.25, "gain", 
                         ifelse(try34out$seg.mean < -0.25, "loss", "no_change"))

try34out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try34out))  # create a list to store matching genes for each row

for (j in 1:nrow(try34out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try34out[j, 3], try34out[j, 4]) &
      filtered_genes153_all[, 5] == try34out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try34out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

try34out_expanded <- try34out %>%
  separate_rows(genes, sep = ",")
#35
try35out <- try35seg$output
try35out$state <- ifelse(try35out$seg.mean > 0.25, "gain", 
                         ifelse(try35out$seg.mean < -0.25, "loss", "no_change"))

try35out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try35out))  # create a list to store matching genes for each row

for (j in 1:nrow(try35out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try35out[j, 3], try35out[j, 4]) &
      filtered_genes153_all[, 5] == try35out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try35out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

try35out_expanded <- try35out %>%
  separate_rows(genes, sep = ",")

#40
try40out <- try40seg$output
try40out$state <- ifelse(try40out$seg.mean > 0.25, "gain", 
                         ifelse(try40out$seg.mean < -0.25, "loss", "no_change"))

try40out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try40out))  # create a list to store matching genes for each row

for (j in 1:nrow(try40out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try40out[j, 3], try40out[j, 4]) &
      filtered_genes153_all[, 5] == try40out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try40out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

try40out_expanded <- try40out %>%
  separate_rows(genes, sep = ",")
#43
try43out <- try43seg$output
try43out$state <- ifelse(try43out$seg.mean > 0.25, "gain", 
                         ifelse(try43out$seg.mean < -0.25, "loss", "no_change"))

try43out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try43out))  # create a list to store matching genes for each row

for (j in 1:nrow(try43out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try43out[j, 3], try43out[j, 4]) &
      filtered_genes153_all[, 5] == try43out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try43out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

try43out_expanded <- try43out %>%
  separate_rows(genes, sep = ",")
#44
try44out <- try44seg$output
try44out$state <- ifelse(try44out$seg.mean > 0.25, "gain", 
                         ifelse(try44out$seg.mean < -0.25, "loss", "no_change"))

try44out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try44out))  # create a list to store matching genes for each row

for (j in 1:nrow(try44out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try44out[j, 3], try44out[j, 4]) &
      filtered_genes153_all[, 5] == try44out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try44out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

try44out_expanded <- try44out %>%
  separate_rows(genes, sep = ",")
#55
try55out <- try55seg$output
try55out$state <- ifelse(try55out$seg.mean > 0.25, "gain", 
                         ifelse(try55out$seg.mean < -0.25, "loss", "no_change"))

try55out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try55out))  # create a list to store matching genes for each row

for (j in 1:nrow(try55out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try55out[j, 3], try55out[j, 4]) &
      filtered_genes153_all[, 5] == try55out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try55out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try55out_expanded <- try55out %>%
  separate_rows(genes, sep = ",")
#57
try57out <- try57seg$output
try57out$state <- ifelse(try57out$seg.mean > 0.25, "gain", 
                         ifelse(try57out$seg.mean < -0.25, "loss", "no_change"))

try57out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try57out))  # create a list to store matching genes for each row

for (j in 1:nrow(try57out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try57out[j, 3], try57out[j, 4]) &
      filtered_genes153_all[, 5] == try57out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try57out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try57out_expanded <- try57out %>%
  separate_rows(genes, sep = ",")
#58
try58out <- try58seg$output
try58out$state <- ifelse(try58out$seg.mean > 0.25, "gain", 
                         ifelse(try58out$seg.mean < -0.25, "loss", "no_change"))

try58out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try58out))  # create a list to store matching genes for each row

for (j in 1:nrow(try58out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try58out[j, 3], try58out[j, 4]) &
      filtered_genes153_all[, 5] == try58out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try58out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try58out_expanded <- try58out %>%
  separate_rows(genes, sep = ",")
#59
try59out <- try59seg$output
try59out$state <- ifelse(try59out$seg.mean > 0.25, "gain", 
                         ifelse(try59out$seg.mean < -0.25, "loss", "no_change"))

try59out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try59out))  # create a list to store matching genes for each row

for (j in 1:nrow(try59out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try59out[j, 3], try59out[j, 4]) &
      filtered_genes153_all[, 5] == try59out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try59out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

try59out_expanded <- try59out %>%
  separate_rows(genes, sep = ",")
#60
try60out <- try60seg$output
try60out$state <- ifelse(try60out$seg.mean > 0.25, "gain", 
                         ifelse(try60out$seg.mean < -0.25, "loss", "no_change"))

try60out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try60out))  # create a list to store matching genes for each row

for (j in 1:nrow(try60out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try60out[j, 3], try60out[j, 4]) &
      filtered_genes153_all[, 5] == try60out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try60out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try60out_expanded <- try60out %>%
  separate_rows(genes, sep = ",")

#68
try68out <- try68seg$output
try68out$state <- ifelse(try68out$seg.mean > 0.25, "gain", 
                         ifelse(try68out$seg.mean < -0.25, "loss", "no_change"))

try68out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try68out))  # create a list to store matching genes for each row

for (j in 1:nrow(try68out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try68out[j, 3], try68out[j, 4]) &
      filtered_genes153_all[, 5] == try68out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try68out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try68out_expanded <- try68out %>%
  separate_rows(genes, sep = ",")

#81
try81out <- try81seg$output
try81out$state <- ifelse(try81out$seg.mean > 0.25, "gain", 
                         ifelse(try81out$seg.mean < -0.25, "loss", "no_change"))

try81out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try81out))  # create a list to store matching genes for each row

for (j in 1:nrow(try81out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try81out[j, 3], try81out[j, 4]) &
      filtered_genes153_all[, 5] == try81out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try81out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try81out_expanded <- try81out %>%
  separate_rows(genes, sep = ",")

#85
try85out <- try85seg$output
try85out$state <- ifelse(try85out$seg.mean > 0.25, "gain", 
                         ifelse(try85out$seg.mean < -0.25, "loss", "no_change"))

try85out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try85out))  # create a list to store matching genes for each row

for (j in 1:nrow(try85out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try85out[j, 3], try85out[j, 4]) &
      filtered_genes153_all[, 5] == try85out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try85out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

try85out_expanded <- try85out %>%
  separate_rows(genes, sep = ",")
#93
try93out <- try93seg$output
try93out$state <- ifelse(try93out$seg.mean > 0.25, "gain", 
                         ifelse(try93out$seg.mean < -0.25, "loss", "no_change"))

try93out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try93out))  # create a list to store matching genes for each row

for (j in 1:nrow(try93out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try93out[j, 3], try93out[j, 4]) &
      filtered_genes153_all[, 5] == try93out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try93out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)

try93out_expanded <- try93out %>%
  separate_rows(genes, sep = ",")
#111
try111out <- try111seg$output
try111out$state <- ifelse(try111out$seg.mean > 0.25, "gain", 
                         ifelse(try111out$seg.mean < -0.25, "loss", "no_change"))

try111out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try111out))  # create a list to store matching genes for each row

for (j in 1:nrow(try111out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try111out[j, 3], try111out[j, 4]) &
      filtered_genes153_all[, 5] == try111out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try111out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try111out_expanded <- try111out %>%
  separate_rows(genes, sep = ",")
#112
try112out <- try112seg$output
try112out$state <- ifelse(try112out$seg.mean > 0.25, "gain", 
                         ifelse(try112out$seg.mean < -0.25, "loss", "no_change"))

try112out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try112out))  # create a list to store matching genes for each row

for (j in 1:nrow(try112out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try112out[j, 3], try112out[j, 4]) &
      filtered_genes153_all[, 5] == try112out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try112out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try112out_expanded <- try112out %>%
  separate_rows(genes, sep = ",")
#113
try113out <- try113seg$output
try113out$state <- ifelse(try113out$seg.mean > 0.25, "gain", 
                         ifelse(try113out$seg.mean < -0.25, "loss", "no_change"))

try113out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try113out))  # create a list to store matching genes for each row

for (j in 1:nrow(try113out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try113out[j, 3], try113out[j, 4]) &
      filtered_genes153_all[, 5] == try113out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try113out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try113out_expanded <- try113out %>%
  separate_rows(genes, sep = ",")
#114
try114out <- try114seg$output
try114out$state <- ifelse(try114out$seg.mean > 0.25, "gain", 
                         ifelse(try114out$seg.mean < -0.25, "loss", "no_change"))

try114out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try114out))  # create a list to store matching genes for each row

for (j in 1:nrow(try114out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try114out[j, 3], try114out[j, 4]) &
      filtered_genes153_all[, 5] == try114out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try114out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try114out_expanded <- try114out %>%
  separate_rows(genes, sep = ",")
#120
try120out <- try120seg$output
try120out$state <- ifelse(try120out$seg.mean > 0.25, "gain", 
                         ifelse(try120out$seg.mean < -0.25, "loss", "no_change"))

try120out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try120out))  # create a list to store matching genes for each row

for (j in 1:nrow(try120out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try120out[j, 3], try120out[j, 4]) &
      filtered_genes153_all[, 5] == try120out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try120out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try120out_expanded <- try120out %>%
  separate_rows(genes, sep = ",")
#122
try122out <- try122seg$output
try122out$state <- ifelse(try122out$seg.mean > 0.25, "gain", 
                         ifelse(try122out$seg.mean < -0.25, "loss", "no_change"))

try122out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try122out))  # create a list to store matching genes for each row

for (j in 1:nrow(try122out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try122out[j, 3], try122out[j, 4]) &
      filtered_genes153_all[, 5] == try122out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try122out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try122out_expanded <- try122out %>%
  separate_rows(genes, sep = ",")
#136
try136out <- try136seg$output
try136out$state <- ifelse(try136out$seg.mean > 0.25, "gain", 
                         ifelse(try136out$seg.mean < -0.25, "loss", "no_change"))

try136out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try136out))  # create a list to store matching genes for each row

for (j in 1:nrow(try136out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try136out[j, 3], try136out[j, 4]) &
      filtered_genes153_all[, 5] == try136out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try136out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try136out_expanded <- try136out %>%
  separate_rows(genes, sep = ",")
#151
try151out <- try151seg$output
try151out$state <- ifelse(try151out$seg.mean > 0.25, "gain", 
                         ifelse(try151out$seg.mean < -0.25, "loss", "no_change"))

try151out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try151out))  # create a list to store matching genes for each row

for (j in 1:nrow(try151out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try151out[j, 3], try151out[j, 4]) &
      filtered_genes153_all[, 5] == try151out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try151out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try151out_expanded <- try151out %>%
  separate_rows(genes, sep = ",")
#153
try153out <- try153seg$output
try153out$state <- ifelse(try153out$seg.mean > 0.25, "gain", 
                         ifelse(try153out$seg.mean < -0.25, "loss", "no_change"))

try153out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try153out))  # create a list to store matching genes for each row

for (j in 1:nrow(try153out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try153out[j, 3], try153out[j, 4]) &
      filtered_genes153_all[, 5] == try153out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try153out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try153out_expanded <- try153out %>%
  separate_rows(genes, sep = ",")
#157
try157out <- try157seg$output
try157out$state <- ifelse(try157out$seg.mean > 0.25, "gain", 
                         ifelse(try157out$seg.mean < -0.25, "loss", "no_change"))

try157out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try157out))  # create a list to store matching genes for each row

for (j in 1:nrow(try157out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try157out[j, 3], try157out[j, 4]) &
      filtered_genes153_all[, 5] == try157out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try157out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try157out_expanded <- try157out %>%
  separate_rows(genes, sep = ",")
#172
try172out <- try172seg$output
try172out$state <- ifelse(try172out$seg.mean > 0.25, "gain", 
                         ifelse(try172out$seg.mean < -0.25, "loss", "no_change"))

try172out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try172out))  # create a list to store matching genes for each row

for (j in 1:nrow(try172out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try172out[j, 3], try172out[j, 4]) &
      filtered_genes153_all[, 5] == try172out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try172out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try172out_expanded <- try172out %>%
  separate_rows(genes, sep = ",")
#176
try176out <- try176seg$output
try176out$state <- ifelse(try176out$seg.mean > 0.25, "gain", 
                         ifelse(try176out$seg.mean < -0.25, "loss", "no_change"))

try176out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try176out))  # create a list to store matching genes for each row

for (j in 1:nrow(try176out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try176out[j, 3], try176out[j, 4]) &
      filtered_genes153_all[, 5] == try176out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try176out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try176out_expanded <- try176out %>%
  separate_rows(genes, sep = ",")
#179
try179out <- try179seg$output
try179out$state <- ifelse(try179out$seg.mean > 0.25, "gain", 
                         ifelse(try179out$seg.mean < -0.25, "loss", "no_change"))

try179out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try179out))  # create a list to store matching genes for each row

for (j in 1:nrow(try179out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try179out[j, 3], try179out[j, 4]) &
      filtered_genes153_all[, 5] == try179out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try179out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try179out_expanded <- try179out %>%
  separate_rows(genes, sep = ",")
#185
try185out <- try185seg$output
try185out$state <- ifelse(try185out$seg.mean > 0.25, "gain", 
                         ifelse(try185out$seg.mean < -0.25, "loss", "no_change"))

try185out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try185out))  # create a list to store matching genes for each row

for (j in 1:nrow(try185out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try185out[j, 3], try185out[j, 4]) &
      filtered_genes153_all[, 5] == try185out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try185out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try185out_expanded <- try185out %>%
  separate_rows(genes, sep = ",")
#206
try206out <- try206seg$output
try206out$state <- ifelse(try206out$seg.mean > 0.25, "gain", 
                         ifelse(try206out$seg.mean < -0.25, "loss", "no_change"))

try206out$genes <- NA  # initialize an empty vector

matching_genes <- vector("list", nrow(try206out))  # create a list to store matching genes for each row

for (j in 1:nrow(try206out)) {
  # Find matching genes
  matching_genes[[j]] <- filtered_genes153_all$hgnc_symbol[
    between(filtered_genes153_all[, 3], try206out[j, 3], try206out[j, 4]) &
      filtered_genes153_all[, 5] == try206out[j, 2]
  ]
}

# If there are matching genes, paste them together separated by comma
try206out$genes <- sapply(matching_genes, function(x) if (length(x) > 0) paste(x, collapse = ",") else NA)
try206out_expanded <- try206out %>%
  separate_rows(genes, sep = ",")

# now want to intersect panel with the wes that i just got with genes, this script is modified from the WORKINGCODEscript
#WES 16 is faulty
#missing 18.2

panel18<- left_join(`L19-0565_cna`, try18out_expanded, by = "genes")
panel28<- left_join(`L20-3102_cna`, try28out_expanded, by = "genes")
panel32<- left_join(`L19-3813_cna`, try32out_expanded, by = "genes")
panel34<- left_join(`L19-2057_cna`, try34out_expanded, by = "genes")
panel35<- left_join(`L21-2751_cna`, try35out_expanded, by = "genes")
panel40<- left_join(`L20-3645_cna`, try40out_expanded, by = "genes")
panel43<- left_join(`L19-0648_cna`, try43out_expanded, by = "genes")
panel44<- left_join(`L19-4581_cna`, try44out_expanded, by = "genes")
panel55<- left_join(`L19-1139_cna`, try55out_expanded, by = "genes")
panel57<- left_join(`L19-1367_cna`, try57out_expanded, by = "genes")
panel58<- left_join(`L19-4580_cna`, try58out_expanded, by = "genes")
panel59<- left_join(`L21-3688_cna`, try59out_expanded, by = "genes")
panel60<- left_join(`L19-2689_cna`, try60out_expanded, by = "genes")
panel68<- left_join(`L19-3826_cna`, try68out_expanded, by = "genes")
panel81<- left_join(`L20-0350_cna`, try81out_expanded, by = "genes")
panel85<- left_join(`L20-6363_cna`, try85out_expanded, by = "genes")
panel93<- left_join(`L20-6275_cna`, try93out_expanded, by = "genes")
panel111<- left_join(`L20-4112_cna`, try111out_expanded, by = "genes")
panel112<- left_join(`L20-4109_cna`, try112out_expanded, by = "genes")
panel113<- left_join(`L20-6010_cna`, try113out_expanded, by = "genes")
panel114<- left_join(`L20-5347_cna`, try114out_expanded, by = "genes")
panel120<- left_join(`L20-5648_cna`, try120out_expanded, by = "genes")
panel122<- left_join(`L20-5639_cna`, try122out_expanded, by = "genes")
panel136<- left_join(`L21-0969_cna`, try136out_expanded, by = "genes")
panel151<- left_join(`L21-1735_cna`, try151out_expanded, by = "genes")
panel153<- left_join(`L21-3701_cna`, try153out_expanded, by = "genes")
panel157<- left_join(`L21-2752_cna`, try157out_expanded, by = "genes")
panel172<- left_join(`L21-4600_cna`, try172out_expanded, by = "genes")
panel176<- left_join(`L21-4814_cna`, try176out_expanded, by = "genes")
panel179<- left_join(`L21-5659_cna`, try179out_expanded, by = "genes")
panel185<- left_join(`L21-5861_cna`, try185out_expanded, by = "genes")
panel206<- left_join(`L21-5938_cna`, try206out_expanded, by = "genes")

#now combine everything to a single dataframe
centeredpanel <- rbind(panel18, panel28, panel32, panel34, panel35, panel40, panel43, panel44, panel55, panel57, panel58, panel59, panel60, panel68, panel81, panel85, panel93, panel111, panel112, panel113, panel114, panel120, panel122, panel136, panel149, panel151, panel153, panel157, panel166, panel172, panel176, panel179, panel185, panel206)
centeredpanel$chr <- factor(centeredpanel$chr, levels = x_order)
# rearrange columns, this dataframe has no significants for panel
centeredpanel <- centeredpanel[, c(10, 1:9, 11:ncol(centeredpanel))]
same_signalc <- ifelse(centeredpanel$log2 > 0 & centeredpanel$seg.mean > 0, "Match",
                       ifelse(centeredpanel$log2 < 0 & centeredpanel$seg.mean < 0, "Match", "Mixed"))

centeredpanel$same_sign <- same_signalc

# panels 32, 34, 40, 44, 60, 85, 93, 112, 120, 136, 172, 176, 185 and 206 have no significant panels
# However WES does have changes 

# make a dataframe without the no signficant
centeredpanel_wo<- centeredpanel[!grepl("No signi", centeredpanel$chr), ]
centeredpanel_wo$chr <- factor(centeredpanel_wo$chr, levels = x_order)
same_signal <- ifelse(centeredpanel_wo$log2 > 0 & centeredpanel_wo$seg.mean > 0, "Match",
                                             ifelse(centeredpanel_wo$log2 < 0 & centeredpanel_wo$seg.mean < 0, "Match", "Mixed"))
centeredpanel_wo$same_sign <- same_signal

# only columns of interest

newpanel <- centeredpanel_wo[, c(1,2,3,4,5,6,12,13,17)]

# columns 3,4,5,6 make reference to the values in panel, 7,8 are WES

newpanel_filter <- newpanel[complete.cases(newpanel), ]
x_order <- gsub("chr", "", chr_order)
newpanel_filter$chr <- factor(newpanel_filter$chr, levels = x_order)

percofcount <- sum(newpanel_filter$same_sign != "Mixed")

# Calculate the percentage of matches
total <- (percofcount / nrow(newpanel_filter)) * 100 #88% match

panelcount <- tableby(same_sign ~ ID, data = newpanel_filter)
summary(panelcount, text=TRUE)


# save the correspond table for only comparable WES and panel files
write.csv(full_table, file = "samples.csv", row.names = FALSE)

### BELOW ARE THE WES THAT HAVE COMPARABLE WITH WES
#17 docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO017_19-3303_DNAFFPE_NA_C_HN00151141.copynumber --recenter-up 3 --output-file /data/example17
#18 docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO018_17-41236_DNAFFPE_NA_C_HN00131668.copynumber --recenter-down 2.67 --output-file /data/example18
#18.2 docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO018_20-19722_DNAFFPE_NA_C_HN00151141.copynumber --recenter-up 2.96 --output-file /data/example18.2
#28 docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO028_04-13962_DNAFFPE_NA_C_HN00188109.copynumber --recenter-down .78 --output-file /data/example28
#32 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO032_M18-2884_DNAFFPE_NA_C_HN00131668.copynumber --recenter-down .46 --output-file /data/example32
#34 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO034_17-32279_DNAFFPE_NA_C_HN00168186.copynumber --recenter-up .43 --output-file /data/example34
#35 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO035_19-39081_DNAFFPE_NA_C_HN00151141.copynumber --recenter-up 1.54 --output-file /data/example35
#40 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO040_18-40640_DNAFFPE_NA_C_HN00151141.copynumber --recenter-up 1.20 --output-file /data/example40
#43 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO043_19-3611_DNAFFPE_NA_C_HN00151141.copynumber --recenter-up 2.75 --output-file /data/example43
#44 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO044_18-41762_DNAFFPE_NA_C_HN00168186.copynumber --recenter-up 1.42 --output-file /data/example44
#55 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO055_19-9615_DNAFFPE_NA_C_HN00151141.copynumber --recenter-up 1.81 --output-file /data/example55
#57 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO057_19-10550_DNAFFPE_NA_C_HN00151141.copynumber --recenter-up 2.45 --output-file /data/example57
#58 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO058_19-8627_DNAFFPE_NA_C_HN00168186.copynumber --recenter-up 1.13 --output-file /data/example58
#59 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO059_VH-21-B-017663_DNAFFPE_NA_C_HN00188109.copynumber --recenter-down .021 --output-file /data/example59
#60 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO060_M19-1182_DNAFFPE_NA_C_HN00168186.copynumber --recenter-up 0.43 --output-file /data/example60
#68 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO068_M19-1607_DNAFFPE_NA_C_HN00131668.copynumber --recenter-up 0.43 --output-file /data/example68
#81 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO081_19-33842A2_DNAFFPE_NA_C_HN00180257.copynumber --recenter-up 0.91 --output-file /data/example81
#85 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO085_18-14267A1_DNAFFPE_NA_C_HN00180257.copynumber --recenter-up 0.23 --output-file /data/example85
#93 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO093_14-1531B1_DNAFFPE_NA_C_HN00180257.copynumber --recenter-down .22 --output-file /data/example93
#111 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO111_VH20B016682_DNAFFPE_NA_C_HN00168186.copynumber --recenter-up 0.41 --output-file /data/example111
#112 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO112_M20-1730_DNAFFPE_NA_C_HN00188109.copynumber --recenter-up 3.29 --output-file /data/example112
#113 HECHO y todos los de arriba
#114 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO114_20-21567_DNAFFPE_NA_C_HN00168186.copynumber --recenter-up 2.69 --output-file /data/example114
#120 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO120_20-28233A1_DNAFFPE_NA_C_HN00151141.copynumber --recenter-up 2.50 --output-file /data/example120
#122 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO122_19-31982A1_DNAFFPE_NA_C_HN00180257.copynumber --recenter-up 0.21 --output-file /data/example122
#131 docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO131_21-54672_DNAFF_NA_C_HN00176021.copynumber --recenter-up 1.81 --output-file /data/example131
#136 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO136_21-2457_DNAFFPE_NA_C_HN00168186.copynumber --recenter-up 2.05 --output-file /data/example136
#151 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO151_20-31538A1_DNAFFPE_NA_C_HN00188109.copynumber --recenter-up 2.23 --output-file /data/example151
#153 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO153_21-7492_DNAFFPE_NA_C_HN00188109.copynumber --recenter-up 2.48 --output-file /data/example153
#157 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO157_21-12787_DNAFF_NA_C_HN00176021.copynumber --recenter-up 1.54 --output-file /data/example157
#172 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO172_M21-1908_DNAFFPE_NA_C_HN00188109.copynumber --recenter-down 1.06 --output-file /data/example172
#176 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO176_21-23725_DNAFFPE_NA_C_HN00188109.copynumber --recenter-down .43 --output-file /data/example176
#179 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO179_21-25869_DNAFFPE_NA_C_HN00176021.copynumber --recenter-up 1.32 --output-file /data/example179
#179.2 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO179_21-26303B1_DNAFFPE_NA_C_HN00188109.copynumber --recenter-up .055 --output-file /data/example179.2
#185 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO185_B21-21941B_DNAFFPE_NA_C_HN00176021.copynumber --recenter-up 0.75 --output-file /data/example185
#206 HECHO docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/PRO206_M21-2369_DNAFFPE_NA_C_HN00188109.copynumber --recenter-up 0.32 --output-file /data/example206

####
#library(gtools)
# get all dataframes
all_pros <- ls(pattern = "^PRO\\d{3}(\\.2)?$")
all_pros <- all_pros[mixedorder(all_pros)]

# Predefined levels
levels_region_call <- c("amp", "del", "neutral")

# Use sapply to apply a function to all data frames
result <- sapply(all_pros, function(all_pros) {
  df_ <- get(all_pros)
  
  # Count the frequency of each level
  counts <- table(factor(df_$region_call, levels = levels_region_call))
  
  # Return the counts
  return(counts)
})

# Convert the result to a data frame
df_result <- as.data.frame(t(result))

#make a column for the PROID
df_result <- df_result %>% rownames_to_column("PROID")

# now for the recentered
#this try to get table for the re centered that are comparable
all_try <- ls(pattern = "^try\\d{2,3}$")
all_try <- all_try[mixedorder(all_try)]

levels_region_call <- c("amp", "del", "neutral")

resulttry <- sapply(all_try, function(all_try) {
  dftry <- get(all_try)
  
  # Count the frequency of each level
  countstry <- table(factor(dftry$region_call, levels = levels_region_call))
  
  # Return the counts
  return(countstry)
})

# Convert the result to a data frame
dftry_result <- as.data.frame(t(resulttry))

dftry_result <- dftry_result %>% rownames_to_column("PROID")
dftry_result$PROID <- gsub("try", "PRO", dftry_result$PROID)


# simple graph showing amp, del and neutral
# Reshape data for ggplot
gene_summary_long <- gene_summary %>%
  pivot_longer(cols = c(Gain, Loss),
               names_to = "state",
               values_to = "percentage")

ggplot(gene_summary_long, aes(x = genes, y = percentage, fill = state)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Gain" = "blue", "Loss" = "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Percentage") +
  xlab("Genes")




#until here


#create the join columns
result_df <- left_join(df_result, tablefor, by = "PROID") #this has all statuses

#keep columns of interest
complete_table <- result_df[, c(1,2,3,4,6,7,8)]

# this is for including more information such as castration status

tablefor <- correspond_table

tablefor <- tablefor[tablefor$PROID %in% all_pros, ]

tablepanel <- left_join(df_result, correspond_table, by = "PROID")



#this code is to see how the PureCN works which is a 
#ABSOLUTE algorithm user to identidy ploidy and 
#copy number calling
retSegmented <- runAbsoluteCN(seg.file = seg.file,
                              interval.file = interval.file, vcf.file = vcf.file,
                              max.candidate.solutions = 1, genome = "hg19",
                              test.purity = seq(0.3,0.7,by = 0.05), verbose = FALSE,
                              plot.cnv = FALSE)



