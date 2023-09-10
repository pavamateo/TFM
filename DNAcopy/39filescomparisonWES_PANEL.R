# READING 39 WES FILES FOR COMPARISON
chr_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
weslist <- list()

# Specify the path where the files are
path_ <- "/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/"

for (f in wes_files) {
  if (endsWith(f, ".copynumber")) {
    # Include the full path of the file when reading it
    data_ <- read.csv(paste0(path_, f), sep = "\t")
    
    # Rename column names
    colnames(data_)[1] <- "chr"
    colnames(data_)[2] <- "start"
    colnames(data_)[3] <- "end"
    
    # Apply the level order to the "chr" column
    data_$chr <- factor(data_$chr, levels = chr_order)
    
    assign(f, data_)  # Assign the data to the environment with the f as the name
    weslist[[f]] <- data_  # Store the data in weslist with the f as the key
  }
}

################## NOW APPLY DNACOPY
for (filename in names(weslist)) {
  data_ <- get(filename)
  
  # Ensure the chr column is ordered according to chr_order
  data_$chr <- factor(data_$chr, levels = chr_order)
  
  # Extract the first 8 characters from filename
  sample_id <- substr(filename, start = 1, stop = 8)
  
  # Perform CNA analysis
  CNA.object <- CNA(genomdat = data_$log2_ratio, chrom = data_$chr, maploc = data_$start, data.type = 'logratio', sampleid = sample_id)
  CNA.smoothed <- smooth.CNA(CNA.object)
  seg <- segment(CNA.smoothed, verbose = 0) #verbose 1 tells me which sample is being analyzed
  
  # Store the segmented data frame
  assign(paste0(sample_id, "seg"), seg)
}

