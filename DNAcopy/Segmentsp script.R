library(dplyr)
#THIS IS CHANGING A BIT, I DONT NEED TO DIRECTLY USE THE FILES FROM allWES DIRECTORY
# JUST NEED ALL THE TRYXXXSEG FILES AND PERFORM SEGMENTS.P FUNCTION
#AFTERWARDS, DO CREATE CVS FILES TO DO THE PERL SCRIPT AND HAVE ALL THE NECESSARY FILES TO READ
# Get all the files in the directory
all_objects <- ls()

# Find the ones that match the pattern "tryxxxout"
segmentpro <- seg_objects <- grep("^try[0-9]+(\\.[0-9]+)?seg$", all_objects, value = TRUE)


for (sample_id in segmentpro) {

  
  # Calculate the segments and filter them
  segments <- segments.p(segmentpro)
  segments_filtered <- segments %>%
    filter(pval < 0.05)
  
  # Rename the columns
  colnames(segments_filtered)[3] <- "chr_start"
  colnames(segments_filtered)[4] <- "chr_stop"
  
  # Add the sample column
  segments_filtered$sample <- segments_filtered$ID
  segments_filtered <- segments_filtered %>%
    select(ID, sample, everything())
  
  # Save the filtered segments
  write.table(segments_filtered, file = paste0("segment", sample_id, "f"), row.names = FALSE)
}


for (sample_id in segmentpro) {
  # Retrieve the segmented data
  seg_data <- get(sample_id)
  
  # Calculate the segments and filter them
  segments <- segments.p(seg_data)
  segments_filtered <- segments %>%
    filter(pval < 0.05)
  
  # Rename the columns
  colnames(segments_filtered)[3] <- "chr_start"
  colnames(segments_filtered)[4] <- "chr_stop"
  
  # Add the sample column
  segments_filtered$sample <- segments_filtered$ID
  segments_filtered <- segments_filtered %>%
    select(ID, sample, everything())
  
  # Save the filtered segments
  write.table(segments_filtered, file = paste0("segment", gsub("^try|out$", "", sample_id), "f"), row.names = FALSE)
}

#THIS GIVES US THE FILES READY FOR THE VARSCAN LAST STEP
