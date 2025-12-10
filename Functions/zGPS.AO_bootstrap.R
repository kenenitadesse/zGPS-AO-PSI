# ----------------------------
# Function: zGPS.AO_bootstrap
# ----------------------------
# Purpose:
#   Perform bootstrap to evaluate the accuracy and statistical significance 
#   of the zero–inflated Gamma–Poisson shrinker with AE ontology (zGPS.AO) analysis.
#
# Arguments:
#   zGPS_result : Output object from zGPS.AO_tool containing observed data and parameters
#   n_perm      : Number of bootstrap samples (default 20)
#   n_cores     : Number of CPU cores for parallel computation (default: half of detected cores)
#   seed        : Random seed for reproducibility (default 1234)
#
# Returns:
#   A zGPS.AO S3 object containing updated bootstrap results and calculated q-values
# ----------------------------

##########################

zGPS.AO_bootstrap = function(zGPS_result,
                          n_perm = 20,
                          n_cores = round(detectCores() / 2),seed = 1234) {
  
  list2env(zGPS_result, envir = environment())
  
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    set.seed(seed)
    new_resample_results =
      foreach(
        i = 1:n_perm,
        .packages = c('dplyr', 'pscl', 'MASS', 'questionr'),
        .export = c("drug_ae_veneto", "pairs2reports", "reports2pairs", "get_all_params", "make_table",'translate_coefs',
                    'update_lambda_hat'),
        .verbose = FALSE
      ) %dorng% {
        report_data_perm = pairs2reports(drug_ae_veneto)
        report_data_perm$AE_Set = sample(report_data_perm$AE_Set)
        pair_data_perm = reports2pairs(report_data_perm)
        
        get_all_params(grp_lst, pair_data_perm, merge_list = merge_list)
      }
    stopCluster(cl)
  } else {
    new_resample_results = list()
    for (i in 1:n_perm) {
      report_data_perm = pairs2reports(drug_ae_veneto)
      report_data_perm$AE_Set = sample(report_data_perm$AE_Set)
      pair_data_perm = reports2pairs(report_data_perm)
      new_resample_results[[i]] =
        get_all_params(grp_lst,
                       pair_data_perm,
                       merge_list = merge_list)
    }
  }
  
  zGPS_result$resample_results =
    c(zGPS_result$resample_results,
      new_resample_results)
  
  zGPS_result = get_pval(zGPS_result)
  return(zGPS_result)
}


