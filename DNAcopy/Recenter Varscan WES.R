# OBTAIN MEAN FROM THE RAW_RATIO TO DO RE CENTERING
# RUNNING THIS SCRIPT TAKES A LONG TIME, USING A DOCKER IMAGE VARSCAN TO PERFORM
# THIS GIVES ME THE FILES TO READ AND PERFORM DNAcopy AFTERWARDS

mean_list <- list()

for (pro in all_pros) {
  data_frame <- get(pro)
  mean_value <- mean(data_frame$raw_ratio)
  mean_list[[pro]] <- mean_value
}

# mean_list now has all the mean values from each patient regarding log2 raw_ratios

# now when to call the mean i would use mean_list$(the PRO i need)
#### When value is - use recenter down
#### When value is + use recenter up

## list of all WES with re center
# CREATE DF TO STORE THE NECESSARY INFO TO DO THE RE CENTER
# use all_pros
# use file_names
file_names <- list.files(pattern = "\\.copynumber\\.called$")
files_center <- file_names
# Use sub to replace ".called" with nothing
files_center <- sub(".called", "", files_center)
# use mean_list

mean_list <- lapply(mean_list, function(x) {
  # Round to two decimal places
  x <- round(x, 2)
  # If the value is positive, prepend "up ", else prepend "down "
  if (x > 0) {
    paste("up", abs(x))
  } else {
    paste("down", abs(x))
  }
})


# Create an empty data frame with specified column names
recentertable <- data.frame("PRO" = character(82),
                 "File name" = character(82),
                 "Recenter value" = character(82),
                 stringsAsFactors=FALSE)

recentertable$PRO <- all_pros
recentertable$File.name <- files_center
recentertable$Recenter.value <- mean_list

# now get the commands

for (i in 1:nrow(recentertable)) {
  pro_num <- gsub("PRO", "", recentertable$PRO[i])  # Extract the number from the PRO name
  
  command <- paste0("docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/",
                    recentertable$PRO[i],
                    "_",
                    recentertable$File.name[i],
                    ".copynumber --recenter-",
                    recentertable$Recenter.value[i],
                    " --output-file /data/example",
                    pro_num)
  
  print(command)
}

# WITH THE ABOVE I GET ALL THE COMMANDS TO PERFORM THE VARSCAN RECENTERING

# Open a connection to the file
file_conn <- file("docker_commands.sh")

for (i in 1:nrow(recentertable)) {
  command <- paste0("docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/",
                    recentertable$PRO[i],
                    "_",
                    recentertable$File.name[i],
                    ".copynumber --recenter-",
                    recentertable$Recenter.value[i],
                    " --output-file /data/example",
                    recentertable$PRO[i])
  
  writeLines(command, file_conn)
}

close(file_conn)

#############
docker_path <- "/Users/mateopava/Desktop/CNA_Mateo_MSC/allWES/docker_commands.sh"

# Clear out the file if it already exists
if (file.exists(docker_path)) {
  file.remove(docker_path)
}

for (i in 1:nrow(recentertable)) {
  command <- paste0("docker run --rm -it -v /Users/mateopava/Desktop/CNA_Mateo_MSC/allWES:/data alexcoppe/varscan copyCaller /data/",
      
                    recentertable$File.name[i],
                    " --recenter-",
                    recentertable$Recenter.value[i],
                    " --output-file /data/example",
                    recentertable$PRO[i])
  
  # Append to the file
  cat(command, "\n", file = docker_path, append = TRUE)
}

# NOW I HAVE A SCRIPT WITH ALL THE VARSCAN COMMANDS. FILENAME: docker_commands.sh
# chmod +x docker_commands.sh  # Make the script executable
# ./docker_commands.sh  # Run the script

# AFTER THIS IM OBTAINING ALL THE RE CENTERED FILES: new filenames: examplePRO001
# THRESHOLDS FOR THIS STEP ARE THE FOLLOWING, USING THE DEFAULTS, 
# USING MEAN VALUES FOR EACH PATIENT FOR RE CENTER:

# OPTIONS:
#--output-file	Output file to contain the calls
#--min-coverage	Minimum read depth at a position to make a call [8]
#--amp-threshold	Lower bound for log ratio to call amplification [0.25]
#--del-threshold	Upper bound for log ratio to call deletion (provide as positive number) [0.25]
#--min-region-size	Minimum size (in bases) for a region to be counted [10]
#--recenter-up	Recenter data around an adjusted baseline > 0 [0]
#--recenter-down	Recenter data around an adjusted baseline < 0 [0]


### BELOW ARE THE WES THAT HAVE COMPARABLE WITH WES

# 17,18,18.2,28,32,34,35,40,43,44,55,57,58,59,60,68,81,85,93,111,112,113,114,120,122
# 131,136,151,153,157,172,176,179,179.2,185,206