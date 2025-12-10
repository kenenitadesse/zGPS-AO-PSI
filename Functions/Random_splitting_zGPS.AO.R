
source("Random_Data_splitting.R")
#################################################
library(dplyr)
library(questionr)
library(pscl)
library(MASS)
library(ggplot2)
load("merge_list.rda")
load("drug_ae_veneto.rda")
load("dd.meddra.rda")

# ------------------------
# Run the analysis
# ------------------------
suppressMessages({
  suppressWarnings({

train_fracs <- c(0.5, 0.7, 0.8, 0.9)
seeds <- 123:132

# Function to generate a unique name
assign_name <- function(frac, seed) {
  sprintf("frac%.1f_seed%d", frac, seed)
}

# Create an empty list to store results
all_results_random_split <- list()

total_start <- Sys.time()
# Run zGPS.AO_tool for each combination and store in the list
for (frac in train_fracs) {
  for (s in seeds) {
    
    
    # Run the tool
    res_tmp <- zGPS.AO_random_split(dd.meddra, drug_ae_veneto,
                         merge_list,
                         min_freq = 15,
                         min_size = 15,
                         train_frac = frac,
                         seed = s)
    
    # Store in list with unique name
    obj_name <- assign_name(frac, s)
    all_results_random_split[[obj_name]] <- res_tmp
  }
}
# End total timer
total_end <- Sys.time()

  })
})


