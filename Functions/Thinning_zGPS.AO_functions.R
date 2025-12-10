# ------------------------
# Filter groups by minimum size
# ------------------------

select_grp_size = function(grp_lst, min_size = 5) {
  indi = rep(T, length(grp_lst))
  for (i in 1:length(grp_lst)) {
    if (nrow(grp_lst[[i]]) < min_size)
      indi[i] = F
  }
  grp_lst = grp_lst[indi]
  return(grp_lst)
}



# ------------------------
# Remove AEs below a frequency threshold
# ------------------------

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



# ------------------------
# Create contingency tables and perform Poisson thinning
# ------------------------

make_table <- function(dds, 
                       merge_list, 
                       thinning_ratio = c(0.7, 0.3),  
                       seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  # Weight adjustment for multiple drugs per report
  
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
  
  # Create weighted contingency table
  tb <- wtd.table(dds$DRUG_TYPE, dds$AE_NAME, dds$wt) %>% round()
  tb <- tb[rownames(tb) != "others", , drop = FALSE]
  
  # Expected counts
  row_sums <- apply(tb, 1, sum)
  col_sums <- apply(tb, 2, sum)
  total_sum <- sum(tb)
  E_tb <- ((row_sums %o% col_sums) / total_sum)
  
  # Poisson thinning 
  
  thinned_data <- datathin(tb, family = "poisson", epsilon = thinning_ratio)
  train_counts <- as.matrix(thinned_data[,,1])
  v_counts <- as.matrix(thinned_data[,,2])
  
  # Compute expected counts for training and validation sets by scaling the full expected counts
  # according to the thinning ratios (e.g., 70% training, 30% validation).
  # A small constant (1e-5) is added to avoid zero expected counts.
  
  
  # Scale expected counts
  E_t <- (thinning_ratio[1] * E_tb)+1e-5 
  E_v <- (thinning_ratio[2] * E_tb)+1e-5

  
    tb_lst=list(
    tb = as.matrix(tb),
    E_tb = as.matrix(E_tb),
    train_counts = train_counts,
    v_counts = v_counts,
    E_t = as.matrix(E_t),
    E_v = as.matrix(E_v))
  return(tb_lst)
  
}


# ------------------------
# Translate model coefficients to parameters
# ------------------------

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

# ------------------------
# Compute empirical lambda_hat
# ------------------------

update_lambda_hat = function(y, E, r, p, beta) {
  if (y == 0) {
    p_hat = p / (p + (1 - p) * (r / (r + E * beta)) ^ r)
    lambda_hat = (1 - p_hat) * r * beta / (r + E * beta)
  } else if (y > 0) {
    lambda_hat = beta * (r + y) / (r + E * beta)
  }
  
  return(lambda_hat)
}


# ------------------------
# Fit models and compute parameters for all AE groups
# ------------------------

get_all_params = function(grp_lst,dds,merge_list, thinning_ratio = c(0.7, 0.3), seed = 123) {
  tb_lst <- make_table(dds, merge_list = merge_list, thinning_ratio = thinning_ratio, seed = seed)
  tb <- tb_lst$tb               
  E_tb <- tb_lst$E_tb           
  train_counts <- tb_lst$train_counts  
  v_counts <- tb_lst$v_counts          
  E_t <- tb_lst$E_t             
  E_v <- tb_lst$E_v             

  data_lst = lapply(grp_lst, function(x) {
    AEs = x$AE_NAME
    grp_tb <- tb[, AEs, drop = FALSE]
    grp_E_tb <- E_tb[, AEs, drop = FALSE]
    grp_train <- train_counts[, AEs, drop = FALSE]
    grp_E_t <- E_t[, AEs, drop = FALSE]
    grp_v <- v_counts[, AEs, drop = FALSE]
    grp_E_v <- E_v[, AEs, drop = FALSE]
    
    
    grp_data = data.frame(
      y = as.vector(t(grp_tb)),
      y_train = as.vector(t(grp_train)),
      y_val = as.vector(t(grp_v)),
      E = as.vector(t(grp_E_tb)),
      E_train = as.vector(t(grp_E_t)),
      E_val = as.vector(t(grp_E_v)),
      drug = rep(rownames(grp_tb), each = ncol(grp_tb)),
      AE_grp = rep(x$GROUP_NAME[1], ncol(grp_tb) * nrow(grp_tb)),
      AE = rep(AEs, nrow(grp_tb))
    )
    
    
    train_data <- grp_data %>%
      dplyr::select(y_train, E_train, drug, AE_grp, AE) %>%
      rename(y = y_train, E = E_train)
    
    
    if (any(train_data$y == 0)) {
      model <- zeroinfl(
        y ~ drug + offset(log(E)) | drug,
        data = train_data,
        dist = "negbin"
      )
    } else {
      model <- glm.nb(
        y ~ drug + offset(log(E)),
        data = train_data
      )
    }
    
    
    val_data <- grp_data %>%
      dplyr::select(y_val, E_val, drug, AE_grp, AE) %>%
      rename(y = y_val, E = E_val)
    
    val_pred <- predict(model, newdata = val_data, type = "response")
    
    
    grp_data$sq_error <- (val_data$y - val_pred)^2    # squared error 
    grp_data$predicted_values <- val_pred            #  predicted values from the model
    
    
    

    grp_data %>%
      mutate(y_fit = fitted(model)) %>%
      group_by(drug) %>%
      mutate(s = mean(y_fit / E_train)) %>%
      rowwise() %>%
      mutate(
        r = translate_coefs(drug, AE_grp, model)['r'],
        p = translate_coefs(drug, AE_grp, model)['p'],
        beta = translate_coefs(drug, AE_grp, model)['beta'],
        lambda_hat = update_lambda_hat(y_train, E_train, r, p, beta)
      ) %>%
      ungroup() %>%
      dplyr::select(-y_fit) -> grp_data
    
    
    return(grp_data)
  })
  
  big_data = bind_rows(data_lst)
  
  return(big_data)
}

# ------------------------
# Full zGPS.AO data thinning analysis tool
# ------------------------



zGPS.AO_thinning <- function(
    grp_data,
    pair_data,
    merge_list,
    min_freq = 15,      
    min_size = 20,      
    thinning_ratio = c(0.7, 0.3),  
    seed = 123      
) {
  set.seed(seed) 
  grp_lst = split(grp_data, grp_data$GROUP_NAME)
  grp_lst = delete_scarse_AEs(grp_lst, pair_data, min_freq = min_freq)
  grp_lst = select_grp_size(grp_lst, min_size = min_size)

  big_data <- suppressWarnings(
    get_all_params(
      grp_lst = grp_lst,
      dds = pair_data,
      merge_list = merge_list,
      thinning_ratio = thinning_ratio,
      seed = seed
    )
  )

  
  zGPS.AO_result = list(
    big_data = big_data,
    pair_data = pair_data,
    grp_lst = grp_lst,
    merge_list = merge_list
  )
  class(zGPS.AO_result) = 'zGPS.AO'
  
  return(zGPS.AO_result)
}



