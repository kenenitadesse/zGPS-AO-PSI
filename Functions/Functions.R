##########################################
# ===== Function Definitions for zGPS.AO Analysis =====
# These functions support data preprocessing, model fitting, 
# and parameter extraction for the zero–inflated Gamma–Poisson shrinker 
# with AE ontology (zGPS.AO, for short) analysis.
##########################################

# ----------------------------
# Function: select_grp_size
# Purpose: Remove AE groups with fewer than `min_size` AEs
# Arguments:
#   grp_lst: List of data.frames, each representing an AE group
#   min_size: Minimum number of AEs required to keep the group (default 5)
# Returns: Filtered list of AE groups
# ----------------------------


select_grp_size = function(grp_lst, min_size = 5) {
  indi = rep(T, length(grp_lst))
  for (i in 1:length(grp_lst)) {
    if (nrow(grp_lst[[i]]) < min_size)
      indi[i] = F
  }
  grp_lst = grp_lst[indi]
  return(grp_lst)
}



# ----------------------------
# Function: delete_scarce_AEs
# Purpose: Remove rare AEs below a frequency threshold
# Arguments:
#   grp_lst: List of AE groups
#   dds: Data.frame with drug-AE pairs
#   min_freq: Minimum occurrence count for an AE to be retained (default 10)
# Returns: Updated list of AE groups
# ----------------------------




delete_scarse_AEs = function(grp_lst, dds, min_freq = 10) {
  tb = table(dds$DRUG_TYPE, dds$AE_NAME)
  AE_counts = apply(tb, 2, sum)
  AEs = names(which(AE_counts >= min_freq))
  for (i in 1:length(grp_lst)) {
    grp_lst[[i]] %>%
      filter(AE_NAME %in% AEs) ->
      grp_lst[[i]]
  }
  return(grp_lst)
}

# ----------------------------
# Function: pairs2reports
# Purpose: Collapse drug-AE pairs into unique reports per patient
# Arguments:
#   dds: Data.frame with drug-AE pairs
# Returns: Data.frame with one row per ID, listing all drugs and AEs
# ----------------------------

pairs2reports = function(dds) {
  dds %>%
    group_by(ID) %>%
    summarise(
      DRUG_Set = paste0(unique(DRUG_TYPE), collapse = ';'),
      AE_Set = paste0(unique(AE_NAME), collapse = ';')
    ) %>%
    ungroup() ->
    dds_reports
  return(dds_reports)
}

# ----------------------------
# Function: reports2pairs
# Purpose: Expand reports back into individual drug-AE pairs
# Arguments:
#   dds_reports: Data.frame with ID, DRUG_Set, AE_Set
# Returns: Data.frame of all drug-AE pairs per report
# ----------------------------

reports2pairs = function(dds_reports) {
  #browser()
  dds_pairs_lst = lapply(1:nrow(dds_reports), function(i) {
    DRUG_Set = dds_reports$DRUG_Set[i]
    AE_Set = dds_reports$AE_Set[i]
    ID = dds_reports$ID[i]
    data = data.frame(ID = ID,
                      expand.grid(
                        strsplit(DRUG_Set, split = ';')[[1]],
                        strsplit(AE_Set, split = ';')[[1]],
                        stringsAsFactors = F
                      ))
    return(data)
  })
  dds_pairs = bind_rows(dds_pairs_lst)
  names(dds_pairs)[2:3] = c("DRUG_TYPE", "AE_NAME")
  return(dds_pairs)
}


# ----------------------------
# Function: make_table
# Purpose: Compute weighted counts and expected counts for drug-AE pairs
# Arguments:
#   dds: Data.frame with drug-AE pairs
#   merge_list: List defining drug merging rules
# Returns: List containing observed table (tb) and expected table (E_tb)
# ----------------------------



make_table = function(dds,
                      merge_list) {
  dds %>%
    group_by(ID) %>%
    mutate(wt = 1 / length(unique(DRUG_TYPE))) %>%
    ungroup() ->
    dds
  
  for (i in 1:length(merge_list)) {
    indi = dds$DRUG_TYPE %in% merge_list[[i]][[1]]
    dds$DRUG_TYPE[indi] = merge_list[[i]][[2]]
  }
  
  indi_others = !dds$DRUG_TYPE %in% sapply(merge_list, function(x)
    x[[2]])
  dds$DRUG_TYPE[indi_others] = 'others'
  
  tb = wtd.table(dds$DRUG_TYPE, dds$AE_NAME, dds$wt)
  tb = round(tb)
  E_tb = apply(tb, 1, sum) %o% apply(tb, 2, sum) / sum(tb)
  tb = tb[rownames(tb) != 'others', ]
  E_tb = E_tb[rownames(E_tb) != 'others', ]
  
  tb = as.matrix(tb)
  E_tb = as.matrix(E_tb)
  
  tb_lst = list(tb = tb,
                E_tb = E_tb)
  
  return(tb_lst)
}

# ----------------------------
# Function: translate_coefs
# Purpose: Convert zero-inflated Gamma-Poisson model coefficients into interpretable parameters
# Arguments:
#   drug: Drug name
#   AE_grp: AE group name
#   model: Zero-inflated or negative binomial model
# Returns: Named vector of r, p, beta, and s
# ----------------------------

translate_coefs = function(drug, AE_grp, model) {
  if (is.list(model$coefficients)) {
    indi = names(model$coefficients$zero) %in% c('(Intercept)',
                                                 paste0('drug', drug),
                                                 paste0('AE_grp', AE_grp))
    
    a = sum(model$coefficients$zero[indi])
    p = exp(a) / (1 + exp(a))
    
    beta = exp(sum(model$coefficients$count[indi]))
  } else {
    indi = names(model$coefficients) %in% c('(Intercept)',
                                            paste0('drug', drug),
                                            paste0('AE_grp', AE_grp))
    p = 0
    beta = exp(sum(model$coefficients[indi]))
  }
  
  r = model$theta
  s = (1 - p) * beta
  
  return(c(
    r = r,
    p = p,
    beta = beta,
    s = s
  ))
}


# ----------------------------
# Function: update_lambda_hat
# Purpose: Compute empirical Bayes estimate for individual AE relative risk in zGPS.AO
# Arguments:
#   y: Observed count
#   E: Expected count
#   r, p, beta: Model parameters
# Returns: lambda_hat (Bayesian estimate)
# ----------------------------

update_lambda_hat = function(y, E, r, p, beta) {
  if (y == 0) {
    p_hat = p / (p + (1 - p) * (r / (r + E * beta)) ^ r)
    lambda_hat = (1 - p_hat) * r * beta / (r + E * beta)
  } else if (y > 0) {
    lambda_hat = beta * (r + y) / (r + E * beta)
  }
  
  return(lambda_hat)
}


# ----------------------------
# Function: get_all_params
# Purpose: Fit zero–inflated Gamma–Poisson shrinker (zGPS.AO) models for all AE groups 
#          and compute parameters for all drug-AE pairs
# Arguments:
#   grp_lst: List of AE groups
#   dds: Data.frame with drug-AE pairs
#   merge_list: List defining drug merging rules
# Returns: Data.frame with fitted s, r, p, beta, lambda_hat for all pairs
# ----------------------------

get_all_params = function(grp_lst,
                          dds,
                          merge_list) {
  tb_lst = make_table(dds, merge_list = merge_list)
  tb = tb_lst$tb
  E_tb = tb_lst$E_tb
  
  
  
  data_lst = lapply(grp_lst, function(x) {
    AEs = x$AE_NAME
    grp_tb = tb[, AEs]
    grp_E_tb = E_tb[, AEs]
    
    grp_data = data.frame(
      y = as.vector(t(grp_tb)),
      E = as.vector(t(grp_E_tb))+0.0001,
      drug = rep(rownames(grp_tb), each = ncol(grp_tb)),
      AE_grp = rep(x$GROUP_NAME[1], ncol(grp_tb) * nrow(grp_tb)),
      AE = rep(AEs, nrow(grp_tb))
    )
    
    if (any(grp_data$y == 0)) {
      model = zeroinfl(y ~ drug + offset(log(E)) | drug,
                       data = grp_data,
                       dist = "negbin")
      
    } else {
      model = glm.nb(y ~ drug + offset(log(E)),
                     data = grp_data)
    }
    grp_data %>%
      mutate(y_fit = fitted(model)) %>%
      group_by(drug) %>%
      mutate(s = mean(y_fit / E)) %>%
      rowwise() %>%
      mutate(
        r = translate_coefs(drug, AE_grp, model)['r'],
        p = translate_coefs(drug, AE_grp, model)['p'],
        beta = translate_coefs(drug, AE_grp, model)['beta'],
        lambda_hat = update_lambda_hat(y, E, r, p, beta)
      ) %>%
      ungroup() %>%
      dplyr::select(-y_fit) ->
      grp_data
    
    return(grp_data)
  })
  
  big_data = bind_rows(data_lst)
  
  return(big_data)
}

# ----------------------------
# Function: get_pval
# Purpose: Compute p-values for group-level s and individual AE lambda_hat 
#          using bootstrap results in zGPS.AO
# Arguments:
#   zGPS_result: zGPS.AO object containing big_data and resample_results
# Returns: zGPS_result updated with s_pval and lambda_pval
# ----------------------------


get_pval = function(zGPS_result) {
  list2env(zGPS_result, envir = environment())
  #browser()
  s = big_data$s
  max_s_seq = sapply(resample_results, function(big_data_perm) {
    max(big_data_perm$s)})
  max_s_seq = c(max(big_data$s), max_s_seq)
  
  lambda_mat = sapply(resample_results, function(big_data_perm) {
    big_data_perm$lambda_hat })
  
  lambda_mat = cbind(big_data$lambda_hat, lambda_mat)
  
  s_pval = sapply(big_data$s, function(s) {
    mean(s<=max_s_seq)})
  lambda_pval = apply(lambda_mat, 1, function(lambda)
    mean(lambda[1] <= lambda))
  zGPS_result$big_data$s_pval = s_pval
  zGPS_result$big_data$lambda_pval = lambda_pval
  return(zGPS_result)
}

# ----------------------------
# Global Variables for dplyr/ggplot
# Purpose: Avoid R CMD check notes for non-standard evaluation
# ----------------------------

utils::globalVariables(
  c(
    "AE",
    "AE_NAME",
    "AE_grp",
    "E",
    "ID",
    "VAX_TYPE",
    "lambda_hat",
    "lambda_pval",
    "p",
    "p_val",
    "r",
    "s",
    "s_pval",
    "vaccine",
    "y",
    "y_fit"
  )
)

