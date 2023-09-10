# normal calls

for (col in log2_cols) {
  new_col_name <- gsub("_log2", "_call_normal", col)
  wes_wgsfinal[[new_col_name]] <- mapply(make_call, wes_wgsfinal[[col]], wes_wgsfinal$HandE)
}

### this is when using low TF  >0.7 considered as gain,  < -0.58 loss
### this is using normal TF >0.8 gain, <-0.7 loss
# since i considered HandE i can say AMP or DEEPDEL
# Define a function to make the call based on the value and threshold
wes_wgsfinal <- merge(wes_wgsfinal, tumorfraction[, c("PRO_ID", "HandE")], 
                      by.x = "PRO_id", by.y = "PRO_ID", all.x = TRUE)

make_call <- function(value, tf) {
  thresholds <- list(
    high_amp = 0.8,    # high tf AMP x>0.8
    high_gain = 0.25, # high tf GAIN 0.8>x>.25
    high_neutral = -0.25, # high tf NEUTRAL -.25<x<.25
    high_loss = -0.7, # high tf LOSS -0.7 <x<-.25
    high_deepdel = -0.7, # high tf DEEPDEL x < -0.7
    low_amp = 0.7, # low tf AMP x> 0.7
    low_gain = 0.2, # low tf GAIN 0.7 >x> 0.2
    low_neutral = -0.2, # low tf NEUTRAL -0.2 <x< 0.2
    low_loss = -0.58, # low tf LOSS -0.58 <x< -0.20
    low_deepdel = -0.58 #low tf DEEPDEL x < -0.58
    
  )
  
  # Check for NA values first
  if (is.na(value) || is.na(tf)) {
    return(NA)
  }
  
  if (tf <= 0.3) {  # Low TF
    if (value > thresholds$low_gain) {
      return('AMP')
    } 
    if (value >= thresholds$gain_low) {
      return('GAIN')
    }
    else if (value <- thresholds$gain_low) {
      return('GAIN')
    }
    else if (value <- thresholds$low_loss) {
      return('DEEPDEL')
    } else if (value <- thresholds$neutral_low) {
      return('NEUTRAL')
    }
  } else {  # Normal/High TF
    if (value > thresholds$high_gain) {
      return('AMP')
    } else if (value < thresholds$high_loss) {
      return('DEEPDEL')
    } else {
      return('NEUTRAL')
    }
  }
}

for (col in log2_cols) {
  new_col_name <- gsub("_log2", "_call", col)
  wes_wgsfinal[[new_col_name]] <- mapply(make_call, wes_wgsfinal[[col]], wes_wgsfinal$HandE)
}

# 4. Remove the merged `HandE` column
wes_wgsfinal$HandE <- NULL

# the values are not NEUTRAL, should be gain or loss, i dont know what threshold to use

make_call <- function(value, gene, tf) {
  if (is.na(value)) {
    return(NA)
  }
  
  # Thresholds for diploid chromosomes
  thresholds <- list(
    low_gain = 0.7,
    low_loss = -0.58,
    high_gain = 0.8,
    high_loss = -0.7
  )
  
  # Adjusted thresholds for the X chromosome in males (haploid state)
  x_thresholds <- list(
    low_gain = 0.5,    # Adjust AMP should be >1.6
    low_loss = -1.3,   # Adjust GAIN should be 1.6>x>0.6
    high_gain = 1.6,   # Adjust NEUTRAL should be -1.3< x < 0.6
    high_loss = -1.3   # Adjust DEEPDEL should be -1.3 < x
  )
  
  # Determine which set of thresholds to use based on tumor fraction and gene
  if (gene %in% c("AR", "KDM6A")) {
    if (tf <= 0.3) {
      gain_thresh = x_thresholds$low_gain
      loss_thresh = x_thresholds$low_loss
    } else {
      gain_thresh = x_thresholds$high_gain
      loss_thresh = x_thresholds$high_loss
    }
  } else {
    if (tf <= 0.3) {
      gain_thresh = thresholds$low_gain
      loss_thresh = thresholds$low_loss
    } else {
      gain_thresh = thresholds$high_gain
      loss_thresh = thresholds$high_loss
    }
  }
  
  # Make the call based on determined thresholds
  if (value > gain_thresh) {
    return('AMP')
  } else if (value < loss_thresh) {
    return('DEEPDEL')
  } else {
    return('NEUTRAL')
  }
}

# Apply the function
for (col in log2_cols) {
  new_col_name <- gsub("_log2", "_call", col)
  wes_wgsfinal[[new_col_name]] <- mapply(make_call, wes_wgsfinal[[col]], wes_wgsfinal$Gene, wes_wgsfinal$HandE)
}

# still need to correct the values for AR and KDM6A genes which are the ones in chr X
# NEW CORRECTION
make_call <- function(value, tf, gene) {
  
  # Thresholds for Normal/High TF
  neutral_low <- -0.25
  neutral_high <- 0.25
  gain_thresh <- 0.25
  amp_thresh <- 0.8
  loss_thresh <- -0.25
  deepdel_thresh <- -0.7
  
  # Thresholds for Low TF
  low_tf_neutral_low <- -0.2
  low_tf_neutral_high <- 0.2
  low_tf_gain_thresh <- 0.2
  low_tf_amp_thresh <- 0.7
  low_tf_loss_thresh <- -0.2
  low_tf_deepdel_thresh <- -0.58
  
  # X chromosome thresholds
  x_gain_thresh = 0.6
  x_amp_thresh = 1.6
  x_loss_thresh = -1.3
  x_neutral_low = -1.3
  x_neutral_high = 0.6
  
  # Check for NA values first
  if (is.na(value) || is.na(tf)) {
    return(NA)
  }
  
  # Determine which set of thresholds to use based on gene
  if (gene %in% c("AR", "KDM6A")) {
    if (value > x_amp_thresh) {
      return('AMP')
    } else if (value > x_gain_thresh && value <= x_amp_thresh) {
      return('GAIN')
    } else if (value >= x_neutral_low && value <= x_neutral_high) {
      return('NEUTRAL')
    } else if (value <= x_loss_thresh) {
      return('DEEPDEL')
    }
  } else {
    if (tf <= 0.3) {  # Low TF
      if (value > low_tf_amp_thresh) {
        return('AMP')
      } else if (value > low_tf_gain_thresh && value <= low_tf_amp_thresh) {
        return('GAIN')
      } else if (value >= low_tf_neutral_low && value <= low_tf_neutral_high) {
        return('NEUTRAL')
      } else if (value > low_tf_deepdel_thresh && value <= low_tf_loss_thresh) {
        return('LOSS')
      } else if (value <= low_tf_deepdel_thresh) {
        return('DEEPDEL')
      }
    } else {  # Normal/High TF
      if (value >= neutral_low && value <= neutral_high) {
        return('NEUTRAL')
      } else if (value > gain_thresh && value < amp_thresh) {
        return('GAIN')
      } else if (value >= amp_thresh) {
        return('AMP')
      } else if (value < loss_thresh && value > deepdel_thresh) {
        return('LOSS')
      } else if (value <= deepdel_thresh) {
        return('DEEPDEL')
      }
    }
  }
}

# Apply the function
for (col in log2_cols) {
  new_col_name <- gsub("_log2", "_call", col)
  wes_wgsfinal[[new_col_name]] <- mapply(make_call, wes_wgsfinal[[col]], wes_wgsfinal$HandE, wes_wgsfinal$Gene)
}


# table for thresholds
# Create the data frame
cnv_thresholds <- data.frame(
  Condition = c(rep("Low TF", 5), rep("High/Normal TF ≤ 0.3", 5), rep("chrX genes AR, KDM6A", 4)),
  Call = c("AMP", "GAIN", "NEUTRAL", "LOSS", "DEEPDEL", 
           "AMP", "GAIN", "NEUTRAL", "LOSS", "DEEPDEL", 
           "AMP", "GAIN", "NEUTRAL", "DEEPDEL"),
  Threshold_from = c("0.8", "0.25", "-0.25", "-0.7", "-Inf", 
                     "0.7", "0.2", "-0.2", "-0.58", "-Inf", 
                     "1.6", "0.6", "-1.3", "-Inf"),
  Threshold_to = c("Inf", "0.8", "0.25", "-0.25", "-0.7", 
                   "Inf", "0.7", "0.2", "-0.2", "-0.58", 
                   "Inf", "1.6", "0.6", "-1.3")
)

library(gridExtra)
library(grid)
library(ggplot2)

# Convert table to grid
table_grob <- tableGrob(cnv_thresholds)

# Display the table on the plotting device
grid.newpage()
grid.draw(table_grob)

# Save the table as an image
ggsave("table_image.png", plot = table_grob, width = 10, height = 5)
write.table(cnv_thresholds, file = "thresholds.csv", row.names = F, col.names = T, quote = F, sep = "\t")

# Print the table
print(cnv_thresholds)


# Print the table
print(cnv_thresholds)

#order 

wes_wgsfinal<- completed

log2_cols <- c("cnvkit_log2", "ichorcna_log2", "qdnaseq_log2", "panel_log2", "DNAcopy_log2", "ASCAT_log2")

wes_wgsfinal$cnvkit_call<- NULL 
wes_wgsfinal$ichorcna_call<- NULL
wes_wgsfinal$qdnaseq_call<- NULL
wes_wgsfinal$panel_call<- NULL
wes_wgsfinal$DNAcopy_call<- NULL
