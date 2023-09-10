library(DNAcopy)
# THIS GIVES ME THE DNACOPY OUTPUTS INTO THE ENV AS: tryxxxseg #can name them as you like
# IN ORDER FOR THIS SCRIPT TO WORK, THE INPUT .copynumber FILES ARE NECESSARY
# Specify directory path
dir_path <- "/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/"

# Get all file names in the directory
after_var <- list.files(dir_path, pattern = "^examplePRO\\d+(\\.\\d+)?$", full.names = TRUE)

chr_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")


# Loop over all files
for (ile in after_var) {
  # Extract the sample id from the ile name
  sampleid <- sub("examplePRO", "", basename(ile))
  
  # Read the data
  data <- read.csv(ile, sep = "\t")
  
  # Rename columns
  colnames(data)[1] <- "chr"
  colnames(data)[2] <- "start"
  colnames(data)[3] <- "end"
  
  # Apply the level order to the "chr" column
  data$chr <- factor(data$chr, levels = chr_order)
  
  # Create CNA object
  CNA.object <- CNA(genomdat = data$adjusted_log_ratio, chrom = data$chr, maploc = data$start, data.type = 'logratio', sampleid = paste0("PRO", sampleid))
  
  # Smooth the data
  CNA.smoothed <- smooth.CNA(CNA.object)
  
  # Perform segmentation
  seg <- segment(CNA.smoothed, verbose = 0) # If you want to see the command happening use 1 instead of 0
  
  # Assign the segmentation result to a dynamically created variable name
  assign(paste0("try", sampleid, "seg"), seg, envir = .GlobalEnv)
}
# THE OUTPUT IS DNAcopy type objects

# Now we need to store the seg output so we can perform segments.p
# List all objects in the environment
all_objects <- ls()

# Find the ones that match the pattern "tryxxxseg"
seg_objects <- grep("^try[0-9]+(\\.[0-9]+)?seg$", all_objects, value = TRUE)

# Loop over each matching object
for (obj_name in seg_objects) {
  # Get the object
  obj <- get(obj_name)
  
  # Check if 'output' is a valid field of the object
  if("output" %in% names(obj)) {
    # Extract the 'output' part
    output <- obj$output
    
    # Create a new variable name
    new_var_name <- gsub("seg$", "out", obj_name)
    
    # Assign it to a new variable
    assign(new_var_name, output, envir = .GlobalEnv)
  }
}
# THIS GIVES ME ALL THE OUTPUT DATAFRAMES FROM THE DNAcopy seg TYPE OBJECT