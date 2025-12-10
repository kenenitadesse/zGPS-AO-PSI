# =============================================================================
# zGPS.AO Data Thinning Analysis
# =============================================================================
# Purpose: Perform data thinning, fit zero-inflated and negative binomial
#          models, and compute zGPS.AO parameters for pharmacovigilance analysis.
# =============================================================================


# Load functions and packages
source("Thinning_zGPS.AO_functions.R")
library(dplyr)
library(questionr)
library(datathin)
library(pscl)
library(MASS)
library(ggplot2)

# Load data
load("merge_list.rda")
load("drug_ae_veneto.rda")
load("dd.meddra.rda")

# ------------------------
# Thinning parameters
# ------------------------
suppressMessages({
  suppressWarnings({
thinning_ratio_list <- list(c(0.5, 0.5),
                     c(0.7, 0.3),
                     c(0.8, 0.2),
                     c(0.9, 0.1))
seeds <- 123:132

assign_name <- function(thinning_ratio, seed) {
  sprintf("res_epsilon%.1f_%.1f_seed%d",
          thinning_ratio[1]*100, thinning_ratio[2]*100, seed)
}

# ------------------------
# Run zGPS.AO on all thinning ratios and seeds
# ------------------------

all_results_thin <- list()
total_start <- Sys.time()
for (eps in thinning_ratio_list) {
  for (s in seeds) {
    
    obj_name <- assign_name(eps, s)

    t_start <- Sys.time()
    
    res_tmp <- tryCatch({
      
      zGPS.AO_thinning(dd.meddra, drug_ae_veneto,
                merge_list,
                min_freq = 15,
                min_size = 15,
                thinning_ratio = eps,
                seed = s)
      
    }, error = function(e) {
      cat("Error in", obj_name, ":", e$message, "\n")
      return(NULL)
    })
    
    all_results_thin[[obj_name]] <- res_tmp
    
    t_end <- Sys.time()
  }
}


# End total timer
total_end <- Sys.time()
  })
})