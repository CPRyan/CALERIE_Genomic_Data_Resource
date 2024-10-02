###########
# Load packages
###########
library(tidyverse)
library(parallel)
library(forcats)


###########
# Load Data
###########

load(here::here("Data", "CALERIE_EPIC_data.Rdata"))


pheno_data <-read_csv(pheno_data, here::here("Output/Data", "CALERIE_blood_basic_pheno_data_with_DNAm.csv"))  


remove_barcode <-pheno_data %>% filter(deidnum == "24270") %>% pull(Barcode)
remove_barcode
# 203548860043_R05C01

betas <-betas[,colnames(betas) !='203548860043_R05C01']


# Match rows to columns
pheno_data <-sample_file %>% 
  slice(match(colnames(betas), Basename)) 


identical(colnames(betas), pheno_data$Basename)


###########
# Function
###########

summarize_betas <- function(betas, visit_id, filename_suffix) {
  
  # Setup parallelization
  numCores <-detectCores()
  cl <- makeCluster(numCores)
  
  # Determine the barcodes based on the visit ID
  if (visit_id == "full") {
    selected_barcodes <- pheno_data %>% filter(Visit.ID %in% c("Baseline", "12 Month", "24 Month")) %>% pull(Basename)
  } else {
    selected_barcodes <- pheno_data %>% filter(Visit.ID == visit_id) %>% pull(Basename)
  }
  
  # Subset the betas matrix to include only the columns that match the selected barcodes
  selected_betas <- betas[, which(colnames(betas) %in% selected_barcodes)]
  
  # Remove rows with any NA values
  selected_betas_nona <- na.omit(selected_betas)
  
  # Initialize a summary data frame with the probe names as the row names
  summary_df <- data.frame(probe = rownames(selected_betas_nona), row.names = rownames(selected_betas_nona))
  
  # Calculate various summary statistics for each row (probe)
  summary_df$na_count <- parApply(cl, selected_betas_nona, 1, function(x) sum(is.na(x)))
  summary_df$mean <- parApply(cl, selected_betas_nona, 1, mean, na.rm = TRUE)
  summary_df$sd <- parApply(cl, selected_betas_nona, 1, sd, na.rm = TRUE)
  summary_df$median <- parApply(cl, selected_betas_nona, 1, median, na.rm = TRUE)
  summary_df$min <- parApply(cl, selected_betas_nona, 1, min, na.rm = TRUE)
  summary_df$max <- parApply(cl, selected_betas_nona, 1, max, na.rm = TRUE)
  

  # Write the summary data frame
  write_csv(summary_df, here::here("Output/Data/DNAm_Summary_Info", paste0("CALERIE_Blood_", filename_suffix, "_bval_summary_all_probes.csv")))
  
}


###########
# All 
###########

summarize_betas(betas = betas, visit_id = 'full', filename_suffix = 'full')
summarize_betas(betas = betas, visit_id = 'Baseline', filename_suffix = 'baseline')
summarize_betas(betas = betas, visit_id = '12 Month', filename_suffix = '12mo')
summarize_betas(betas = betas, visit_id = '24 Month', filename_suffix = "24mo")
