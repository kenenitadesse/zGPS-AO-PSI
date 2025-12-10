
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
# Create contingency table and split into train/validation
# ------------------------

make_table <- function(dds, merge_list, train_frac = 0.7, seed = 123) {
  
  dds %>%
    group_by(ID) %>%
    mutate(wt = 1 / length(unique(DRUG_TYPE))) %>%
    ungroup() -> dds
  
  for (i in 1:length(merge_list)) {
    indi = dds$DRUG_TYPE %in% merge_list[[i]][[1]]
    dds$DRUG_TYPE[indi] = merge_list[[i]][[2]]
  }
  
  indi_others = !dds$DRUG_TYPE %in% sapply(merge_list, function(x) x[[2]])
  dds$DRUG_TYPE[indi_others] = 'others'
  
  tb = wtd.table(dds$DRUG_TYPE, dds$AE_NAME, dds$wt)
  tb = round(tb)
  
  E_tb = apply(tb, 1, sum) %o% apply(tb, 2, sum) / sum(tb)
  
  tb = tb[rownames(tb) != 'others', ]
  E_tb = E_tb[rownames(E_tb) != 'others', ]
  
  tb = as.matrix(tb)
  E_tb = as.matrix(E_tb)
  
  split_tb_include_zeros_na <- function(tb, E_tb, train_frac = 0.7, seed = 123) {
    set.seed(seed)
    all_indices <- which(!is.na(tb), arr.ind = TRUE)
    n_total <- nrow(all_indices)
    
    train_size <- floor(train_frac * n_total)
    sampled_indices <- sample(1:n_total, size = n_total, replace = FALSE)
    train_indices <- sampled_indices[1:train_size]
    test_indices <- sampled_indices[(train_size + 1):n_total]
    
    tb_train <- tb
    tb_test <- tb
    E_train <- E_tb
    E_test <- E_tb
    
    for (i in 1:n_total) {
      row <- all_indices[i, 1]
      col <- all_indices[i, 2]
      
      if (i %in% train_indices) {
        tb_test[row, col] <- NA
        E_test[row, col] <- NA
      } else {
        tb_train[row, col] <- NA
        E_train[row, col] <- NA
      }
    }
    
    return(list(
      train = tb_train,
      test = tb_test,
      E_train = E_train,
      E_test = E_test
    ))
  }
  
  split_result = split_tb_include_zeros_na(tb, E_tb, train_frac = train_frac, seed = seed)
  
  tb_lst = list(
    tb = tb,
    E_tb = E_tb,
    train = split_result$train,
    test = split_result$test,
    E_train = split_result$E_train,
    E_test = split_result$E_test
  )
  
  return(tb_lst)
}

# ------------------------
# Extract coefficients from model
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
# Fit models and compute parameters for all groups
# ------------------------


get_all_params = function(grp_lst, dds, merge_list, train_frac = 0.7, seed = 123) {
  
  tb_lst = make_table(dds, merge_list = merge_list, train_frac = train_frac, seed = seed)
  tb = tb_lst$train
  tb2 = tb_lst$test
  E_tb = tb_lst$E_train
  E_tb2 = tb_lst$E_test
  
  data_lst = lapply(grp_lst, function(x) {
    AEs = x$AE_NAME
    grp_tb = tb[, AEs]
    grp_tb2 = tb2[, AEs]
    grp_E_tb = E_tb[, AEs]
    grp_E_tb2 = E_tb2[, AEs]
    
    # Prepare data for model fitting
    grp_data = data.frame(
      y = as.vector(t(grp_tb)),
      E = as.vector(t(grp_E_tb))+0.00001,
      drug = rep(rownames(grp_tb), each = ncol(grp_tb)),
      AE_grp = rep(x$GROUP_NAME[1], ncol(grp_tb) * nrow(grp_tb)),
      AE = rep(AEs, nrow(grp_tb))
    )
    
    val_data = data.frame(
      y = as.vector(t(grp_tb2)),
      E = as.vector(t(grp_E_tb2))+0.00001,
      drug = rep(rownames(grp_tb2), each = ncol(grp_tb2)),
      AE_grp = rep(x$GROUP_NAME[1], ncol(grp_tb2) * nrow(grp_tb2)),
      AE = rep(AEs, nrow(grp_tb2))
    )
    
    # Remove NAs
    grp_data <- grp_data[!is.na(grp_data$y) & !is.na(grp_data$E), ]
    val_data <- val_data[!is.na(val_data$y) & !is.na(val_data$E), ]
    
    # Fit model
    if (any(grp_data$y == 0)) {
      model = zeroinfl(y ~ drug + offset(log(E)) | drug,
                       data = grp_data,
                       dist = "negbin")
    } else {
      model = glm.nb(y ~ drug + offset(log(E)),
                     data = grp_data)
    }
    
    # Validation predictions and metrics
    val_pred <- predict(model, newdata = val_data, type = "response")
    val_data$predicted_values <- val_pred
    val_data$sq_error <- (val_data$y - val_pred)^2
    
    # Post-model calculations for training data
    grp_data = grp_data %>%
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
      dplyr::select(-y_fit)
    
    return(list(train = grp_data, validation = val_data))
  })
  
  # Combine all data
  train_data = bind_rows(lapply(data_lst, function(x) x$train))
  validation_data = bind_rows(lapply(data_lst, function(x) x$validation))
  
  return(list(train = train_data, validation = validation_data))
}

# ------------------------
# Main zGPS.AO_tool function
# ------------------------


zGPS.AO_random_split = function(grp_data, pair_data, merge_list, min_freq = 15, min_size = 20,
                        train_frac = 0.8, seed = 42) {
  
  grp_lst = split(grp_data, grp_data$GROUP_NAME)
  grp_lst = delete_scarse_AEs(grp_lst, pair_data, min_freq = min_freq)
  grp_lst = select_grp_size(grp_lst, min_size = min_size)
  
  results = suppressWarnings(get_all_params(grp_lst, pair_data, merge_list,
                                            train_frac = train_frac, seed = seed))
  
  zGPS.AO_result = list(
    big_data = results,  
    pair_data = pair_data,
    grp_lst = grp_lst,
    merge_list = merge_list
  )
  
  class(zGPS.AO_result) = 'zGPS.AO'
  return(zGPS.AO_result)
}

