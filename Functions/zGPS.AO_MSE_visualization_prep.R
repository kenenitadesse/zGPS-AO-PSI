# source("Stratified_splitting_zGPS.AO.R")
# source("Random_splitting_zGPS.AO.R")
# source("Thinning_zGPS.AO.R")
# Helper to parse parameters from names
parse_params <- function(name, method = c("split", "thin")) {
  method <- match.arg(method)
  parts <- strsplit(name, "_")[[1]]
  if (method == "split") {
    frac <- as.numeric(sub("frac", "", parts[1]))
    seed <- as.integer(sub("seed", "", parts[2]))
    list(train_frac = frac, seed = seed)
  } else {
    eps_train <- as.numeric(sub("epsilon", "", parts[2])) / 100
    eps_val <- as.numeric(parts[3]) / 100
    seed <- as.integer(sub("seed", "", parts[4]))
    list(eps_train = eps_train, eps_val = eps_val, seed = seed)
  }
}

# Generic summarization function for MSE by group variables
summarize_mse <- function(results_list, group_vars, method = c("split", "thin"), validation = FALSE) {
  method <- match.arg(method)
  map_dfr(names(results_list), function(name) {
    params <- parse_params(name, method)
    
    # Extract big_data, or validation subset if needed
    df <- if (validation && method == "split") {
      results_list[[name]]$big_data$validation
    } else {
      results_list[[name]]$big_data
    }
    
    df %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(MSE = mean(sq_error, na.rm = TRUE), .groups = "drop") %>%
      bind_cols(as_tibble(params))
  })
}

# Summarize for all methods and groupings

# Stratified splitting: drug-level and drug+AE

mse_by_drug_stratified <- summarize_mse(all_results_stratified, group_vars = "drug", method = "split")
mse_summary_by_frac <- mse_by_drug_stratified %>%
  group_by(drug, train_frac) %>%
  summarise(mean_MSE_Sp_s = mean(MSE, na.rm = TRUE), .groups = "drop")


mse_by_drug_ae_stratified  <- summarize_mse(all_results_stratified, group_vars = c("drug", "AE_grp"), method = "split")
mse_summary_by_drug_ae <- mse_by_drug_ae_stratified %>%
  group_by(drug, AE_grp, train_frac) %>%
  summarise(mean_MSE_SP_s_bth = mean(MSE, na.rm = TRUE), .groups = "drop")


# Data thinning: drug-level and drug+AE
mse_by_drug_thin <- summarize_mse(all_results_thin, group_vars = "drug", method = "thin")
mse_summary_by_thin_frac <- mse_by_drug_thin %>%
  group_by(drug, eps_train, eps_val) %>%
  summarise(mean_MSE_th = mean(MSE, na.rm = TRUE), .groups = "drop")


mse_by_drug_ae_thin <- summarize_mse(all_results_thin, group_vars = c("drug", "AE_grp"), method = "thin")
mse_summary_by_drug_ae_thin <- mse_by_drug_ae_thin %>%
  group_by(drug, AE_grp, eps_train, eps_val) %>%
  summarise(mean_MSE_th_both = mean(MSE, na.rm = TRUE), .groups = "drop")


# Random splitting (validation subset): drug-level and drug+AE
mse_by_drug_random <- summarize_mse(all_results_random_split, group_vars = "drug", method = "split", validation = TRUE)
mse_summary_by_frac_random <- mse_by_drug_random %>%
  group_by(drug, train_frac) %>%
  summarise(mean_MSE_sp_r = mean(MSE, na.rm = TRUE), .groups = "drop")


mse_by_drug_ae_random <- summarize_mse(all_results_random_split, group_vars = c("drug", "AE_grp"), method = "split", validation = TRUE)
mse_summary_by_drug_ae_random <- mse_by_drug_ae_random %>%
  group_by(drug, AE_grp, train_frac) %>%
  summarise(mean_MSE_sp_r_bth = mean(MSE, na.rm = TRUE), .groups = "drop")




# Combine drug-level summaries (aligning by drug and train_frac / eps params as appropriate)
drug_level_mse_summary <- bind_cols(
  mse_summary_by_frac,
  mean_MSE_thin_frac = mse_summary_by_thin_frac[, 4],
  mean_MSE_frac_random = mse_summary_by_frac_random[, 3]
)

drug_level_mse_wide <- drug_level_mse_summary %>%
  pivot_wider(
    id_cols = drug,
    names_from = train_frac,
    values_from = c(mean_MSE_Sp_s, mean_MSE_th, mean_MSE_sp_r),
    names_sep = "_"
  )

# Combine drug-AE level summaries
drug_ae_level_mse_summary <- bind_cols(
  mse_summary_by_drug_ae,
  mean_MSE_thin_thin = mse_summary_by_drug_ae_thin[,5],
  mean_MSE_frac_random = mse_summary_by_drug_ae_random[,4]
)

drug_ae_level_mse_wide <- drug_ae_level_mse_summary %>%
  pivot_wider(
    id_cols = c(drug, AE_grp),
    names_from = train_frac,
    values_from = c(mean_MSE_SP_s_bth, mean_MSE_th_both, mean_MSE_sp_r_bth,),
    names_sep = "_"
  )



# ==== Generate plots  ====


plot_zGPS.AO_psi <- function(df_wide, exclude_pattern = NULL, title_expr = NULL) {
  if (!is.null(exclude_pattern)) {
    cols_to_keep <- !grepl(exclude_pattern, colnames(df_wide))
    df_wide <- df_wide[, cols_to_keep]
  }
  
  # Compute means across drugs (exclude first column 'drug')
  col_means <- colMeans(df_wide[, -1], na.rm = TRUE)
  
  tidy_df <- enframe(col_means, name = "metric_frac", value = "mean_MSE") %>%
    separate(metric_frac, into = c("metric", "train_frac"), sep = "_(?=\\d)", convert = TRUE) %>%
    mutate(metric = recode(metric,
                           mean_MSE_sp_r = "Data splitting (randomly)",
                           mean_MSE_Sp_s = "Data splitting (stratified)",
                           mean_MSE_th = "Data thinning"))
  
  # Define a palette with your preferred colors
  base_colors <- c("#F8766D", "#00BFC4")
  
  # Only assign colors to metrics that actually exist
  metrics_present <- unique(tidy_df$metric)
  colors_to_use <- setNames(base_colors[seq_along(metrics_present)], metrics_present)
  
  ggplot(tidy_df, aes(x = factor(train_frac), y = mean_MSE, fill = metric)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(
      x = expression(plain("Training proportion")~"(" * epsilon * ")"),
      y = "MSE",
      title = title_expr,
      fill = "Method"
    ) +
    scale_fill_manual(values = colors_to_use) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "bottom",
      panel.border = element_blank()
    )
}




###########Group-level MSE by epsilon################

# ===========================
# Drug-AE-level MSE Analysis & Plots
# ===========================



# Function to plot drug-AE-level MSE with optional exclusion of certain methods
plot_zGPS.AO_psi_group <- function(df_wide, exclude_pattern = NULL, title_expr = NULL) {
  # Filter out columns matching exclusion pattern if given
  if (!is.null(exclude_pattern)) {
    cols_to_keep <- !grepl(exclude_pattern, colnames(df_wide))
    df_wide <- df_wide[, cols_to_keep, drop = FALSE]
  }
  
  # Regex pattern to detect mean MSE columns with train fractions 0.5,0.7,0.8,0.9
  mse_cols_pattern <- "mean_MSE_.+_(0\\.5|0\\.7|0\\.8|0\\.9)$"
  
  # Select and pivot longer into method and train_frac columns
  df_long <- df_wide %>%
    select(matches(mse_cols_pattern)) %>%
    pivot_longer(
      cols = everything(),
      names_to = c("method", "train_frac"),
      names_pattern = "^mean_MSE_(.+)_((?:0\\.5|0\\.7|0\\.8|0\\.9))$",
      values_to = "mean_MSE"
    ) %>%
    mutate(
      train_frac = factor(train_frac, levels = c("0.5", "0.7", "0.8", "0.9")),
      method = recode(method,
                      "SP_s_bth" = "Data splitting (stratified)",
                      "th_both" = "Data thinning",
                      "sp_r_bth" = "Data splitting (randomly)",
                      .default = method
      )
    )
  
  # Plot boxplot with methods fill and train_frac x-axis
  ggplot(df_long, aes(x = train_frac, y = mean_MSE, fill = method)) +
    geom_boxplot(position = position_dodge(width = 0.7), width = 0.7, outlier.size = 0.8) +
    labs(
      title=title_expr,
      x = expression(plain("Training proportion")~"(" * epsilon * ")"),
      y = "MSE",
      fill = "Method"
    ) +
    theme_bw(base_size = 20) +
    theme(legend.position = "bottom",panel.border = element_blank()  )
}


############MSE across drugs by epsilon#########

plot_zGPS.AO_drug_eps <- function(drug_level_mse_wide, title_expr = "Drug-level MSE by Training Fraction") {
  
  # Prepare long-format data
  df_long <- drug_level_mse_wide %>%
    dplyr::select(drug, dplyr::matches("mean_MSE_.+_(0\\.5|0\\.7|0\\.8|0\\.9)$")) %>%
    tidyr::pivot_longer(
      cols = -drug,
      names_to = c("method", "train_frac"),
      names_pattern = "^mean_MSE_(.+)_((?:0\\.5|0\\.7|0\\.8|0\\.9))$",
      values_to = "mean_MSE"
    ) %>%
    dplyr::mutate(
      method = dplyr::recode(method,
                             "Sp_s" = "Data splitting (stratified)",
                             "th"   = "Data thinning"
      ),
      train_frac = factor(train_frac, levels = c("0.5", "0.7", "0.8", "0.9")),
      method = factor(method, levels = unique(method))
    ) %>%
    dplyr::filter(!method %in% c("sp_r", "thin_II"))
  
  # Create bar plot
  ggplot(df_long, aes(x = train_frac, y = mean_MSE, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(
      x = expression(plain("Training proportion")~"(" * epsilon * ")"),
      y = "Mean MSE",
      fill = "Method",
      title = title_expr
    ) +
    theme_minimal(base_size = 18)
}



plot_zGPS.AO_drug_eps <- function(drug_level_mse_wide) {
  
  mse_cols_pattern <- "mean_MSE_.+_(0\\.5|0\\.7|0\\.8|0\\.9)$"
  
  df_long_thinning_stratified <- drug_level_mse_wide %>%
    dplyr::select(drug, dplyr::matches(mse_cols_pattern)) %>%
    tidyr::pivot_longer(
      cols = -drug,
      names_to = c("method", "train_frac"),
      names_pattern = "^mean_MSE_(.+)_((?:0\\.5|0\\.7|0\\.8|0\\.9))$",
      values_to = "mean_MSE"
    ) %>%
    dplyr::mutate(
      method = dplyr::recode(method,
                             "Sp_s" = "Data splitting (stratified)",
                             "th"   = "Data thinning"
      ),
      train_frac = factor(train_frac, levels = c("0.5", "0.7", "0.8", "0.9")),
      method = factor(method, levels = unique(method))
    ) %>%
    dplyr::filter(!method %in% c("sp_r", "thin_II"))  # exclude both mean_MSE_sp_r and mean_MSE_thin_II
  
  ggplot(df_long_thinning_stratified, aes(x = train_frac, y = mean_MSE, fill = method)) +
    geom_boxplot(position = position_dodge(width = 0.7), width = 0.7, outlier.size = 0.8) +
    labs(
      title = expression("" * "MSE across drugs by  " * "" * epsilon),
      x = expression(plain("Training proportion")~"(" * epsilon * ")"),
      y = "MSE",
      fill = "Method"
    ) +
    theme_bw(base_size = 20) +
    theme(legend.position = "bottom", panel.border = element_blank())
}






plot_zGPS.AO_random <- function(drug_level_mse_wide, title_expr = "Drug-level MSE (Random splitting)", exclude_pattern = NULL) {
  
  # Optionally exclude columns
  if (!is.null(exclude_pattern)) {
    cols_to_keep <- !grepl(exclude_pattern, colnames(drug_level_mse_wide))
    drug_level_mse_wide <- drug_level_mse_wide[, cols_to_keep, drop = FALSE]
  }
  
  # Compute means across drugs (exclude first column 'drug')
  col_means <- colMeans(drug_level_mse_wide[, -1], na.rm = TRUE)
  
  tidy_df <- enframe(col_means, name = "metric_frac", value = "mean_MSE") %>%
    separate(metric_frac, into = c("metric", "train_frac"), sep = "_(?=\\d)", convert = TRUE) %>%
    mutate(metric = recode(metric,
                           mean_MSE_sp_r = "Data splitting (randomly)",
                           mean_MSE_Sp_s = "Data splitting (stratified)",
                           mean_MSE_th = "Data thinning")) %>%
    # Keep only Data splitting (randomly)
    filter(metric == "Data splitting (randomly)")
  
  # Plot
  ggplot(tidy_df, aes(x = factor(train_frac), y = mean_MSE, fill = metric)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(
      x = expression(plain("Training proportion")~"(" * epsilon * ")"),
      y = "MSE",
      title = title_expr
    ) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 18) +
    theme(
      legend.position = "none",
      panel.border = element_blank()
    )
}


plot_zGPS.AO_random_group <- function(df_wide, method_to_plot = NULL, title_expr = NULL) {
  # Regex pattern to detect mean MSE columns with train fractions 0.5,0.7,0.8,0.9
  mse_cols_pattern <- "mean_MSE_.+_(0\\.5|0\\.7|0\\.8|0\\.9)$"
  
  # Select and pivot longer into method and train_frac columns
  df_long <- df_wide %>%
    select(matches(mse_cols_pattern)) %>%
    pivot_longer(
      cols = everything(),
      names_to = c("method", "train_frac"),
      names_pattern = "^mean_MSE_(.+)_((?:0\\.5|0\\.7|0\\.8|0\\.9))$",
      values_to = "mean_MSE"
    ) %>%
    mutate(
      train_frac = factor(train_frac, levels = c("0.5", "0.7", "0.8", "0.9")),
      method = recode(method,
                      "SP_s_bth" = "Data splitting (stratified)",
                      "th_both" = "Data thinning",
                      "sp_r_bth" = "Data splitting (randomly)",
                      .default = method
      )
    )
  
  # Keep only the specified method if provided
  if (!is.null(method_to_plot)) {
    df_long <- df_long %>% filter(method == method_to_plot)
  }
  
  # Plot boxplot
  ggplot(df_long, aes(x = train_frac, y = mean_MSE, fill = method)) +
    geom_boxplot(fill = "#66C2A5",position = position_dodge(width = 0.7), width = 0.7, outlier.size = 0.8) +
    labs(
      title = title_expr,
      x = expression(plain("Training proportion")~"(" * epsilon * ")"),
      y = "MSE"
    ) +
    theme_bw(base_size = 18) +
    theme(legend.position = "none",panel.border = element_blank()  )  # Remove legend
}


plot_zGPS.AO_drug_eps_random <- function(drug_level_mse_wide) {
  
  mse_cols_pattern <- "mean_MSE_.+_(0\\.5|0\\.7|0\\.8|0\\.9)$"
  
  df_long_random_splitting <- drug_level_mse_wide %>%
    dplyr::select(drug, dplyr::matches(mse_cols_pattern)) %>%
    tidyr::pivot_longer(
      cols = -drug,
      names_to = c("method", "train_frac"),
      names_pattern = "^mean_MSE_(.+)_((?:0\\.5|0\\.7|0\\.8|0\\.9))$",
      values_to = "mean_MSE"
    ) %>%
    dplyr::mutate(
      method = dplyr::recode(method,
                             "sp_r" = "Data splitting (Random)"),
      train_frac = factor(train_frac, levels = c("0.5", "0.7", "0.8", "0.9")),
      method = factor(method, levels = unique(method))
    ) %>%
    dplyr::filter(!method %in% c("Sp_s", "th"))  # exclude both mean_MSE_sp_r and mean_MSE_thin_II
  
  ggplot(df_long_random_splitting, aes(x = train_frac, y = mean_MSE, fill = method)) +
    geom_boxplot(fill = "#66C2A5", position = position_dodge(width = 0.7), width = 0.7, outlier.size = 0.8) +
    labs(
      title = expression("" * "MSE across drugs by  " * "" * epsilon),
      x = expression(plain("Training proportion")~"(" * epsilon * ")"),
      y = "MSE",
      fill = "Method"
    ) +
    theme_bw(base_size = 18) +
    theme(legend.position = "none", panel.border = element_blank())
}











# -------------------------------
# Generic performance evaluation function for thinning results
# -------------------------------

evaluate_performance_thin <- function(res_obj, threshold = 0.5) {
  y_true <- res_obj$big_data$y_val
  y_pred <- res_obj$big_data$predicted_values
  
  y_binary <- ifelse(y_true > 0, 1, 0)
  pred_binary <- ifelse(y_pred > threshold, 1, 0)
  
  cm <- suppressMessages(
    suppressWarnings(confusionMatrix(
    factor(pred_binary, levels = c(0, 1)),
    factor(y_binary, levels = c(0, 1)),
    positive = "1"
  )))
  
  roc_obj <- roc(y_binary, y_pred)
  auc_val <- as.numeric(auc(roc_obj))
  
  precision <- cm$byClass["Pos Pred Value"]
  recall <- cm$byClass["Sensitivity"]
  
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA else 2 * (precision * recall) / (precision + recall)
  
  list(
    auc = auc_val,
    accuracy = cm$overall["Accuracy"],
    sensitivity = recall,
    specificity = cm$byClass["Specificity"],
    precision = precision,
    f1 = f1
  )
}

# --------------------------------
# Evaluate performance across thinning versions and seeds
# --------------------------------


perf_list_thin <- list()
thinning_versions <- c("50.0_50.0", "70.0_30.0", "80.0_20.0", "90.0_10.0")
seeds <- 123:132

for (tv in thinning_versions) {
  for (s in seeds) {
    obj_name <- paste0("res_epsilon", tv, "_seed", s)
    res_obj <- all_results_thin[[obj_name]]
    if (is.null(res_obj)) next
    
    perf <- evaluate_performance_thin(res_obj)
    
    perf_list_thin[[obj_name]] <- data.frame(
      thinning = tv,
      seed = s,
      auc = perf$auc,
      accuracy = perf$accuracy,
      sensitivity = perf$sensitivity,
      specificity = perf$specificity,
      precision = perf$precision,
      f1 = perf$f1
    )
  }
}

perf_all_thinning <- bind_rows(perf_list_thin)



# -------------------------------
# Compute performance summaries by thinning ratio
# ------------------------------- 

perf_summary_thin <- perf_all_thinning %>%
  group_by(thinning) %>%
  summarise(
    mean_auc = mean(auc, na.rm = TRUE),
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    mean_specificity = mean(specificity, na.rm = TRUE),
    mean_precision = mean(precision, na.rm = TRUE),
    mean_f1 = mean(f1, na.rm = TRUE))






# -------------------------------
# Generic performance evaluation function for stratified splitting approach
# -------------------------------



evaluate_performance_stratified <- function(res_obj, threshold = 0.5) {
  y_true <- res_obj$big_data$y_val
  y_pred <- res_obj$big_data$predicted_values
  
  
  y_binary <- ifelse(y_true > 0, 1, 0)
  pred_binary <- ifelse(y_pred > threshold, 1, 0)
  
  cm <- suppressMessages(
    suppressWarnings(confusionMatrix(
    factor(pred_binary, levels = c(0, 1)),
    factor(y_binary, levels = c(0, 1)),
    positive = "1"
  )))
  
  roc_obj <- roc(y_binary, y_pred)
  auc_val <- as.numeric(auc(roc_obj))
  
  precision <- cm$byClass["Pos Pred Value"]
  recall <- cm$byClass["Sensitivity"]
  
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA else 2 * (precision * recall) / (precision + recall)
  
  list(
    auc = auc_val,
    accuracy = cm$overall["Accuracy"],
    sensitivity = recall,
    specificity = cm$byClass["Specificity"],
    precision = precision,
    f1 = f1
  )
}

# --------------------------------
# Evaluate performance across training fractions and seeds
# --------------------------------

perf_list_stratified <- list()
train_fracs <- c(0.5, 0.7, 0.8, 0.9)
seeds <- 123:132

for (frac in train_fracs) {
  for (s in seeds) {
    obj_name <- sprintf("frac%.1f_seed%d", frac, s)
    res_obj <- all_results_stratified[[obj_name]]
    if (is.null(res_obj)) next
    
    perf <- evaluate_performance_stratified(res_obj)
    
    perf_list_stratified[[obj_name]] <- data.frame(
      train_frac = frac,
      seed = s,
      auc = perf$auc,
      accuracy = perf$accuracy,
      sensitivity = perf$sensitivity,
      specificity = perf$specificity,
      precision = perf$precision,
      f1 = perf$f1
    )
  }
}

perf_list_stratified <- bind_rows(perf_list_stratified)

# -------------------------------
# Compute performance summaries by training ratio
# ------------------------------- 

perf_summary_stratified<- perf_list_stratified %>%
  group_by(train_frac) %>%
  summarise(
    mean_auc = mean(auc, na.rm = TRUE),
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    mean_specificity = mean(specificity, na.rm = TRUE),
    mean_precision = mean(precision, na.rm = TRUE),
    mean_f1 = mean(f1, na.rm = TRUE))






# --------------------------------
# Random Data Splitting and Performance Evaluation
# --------------------------------


evaluate_performance_random <- function(res_obj, threshold = 0.5) {
  y_true <- res_obj$big_data$validation$y
  y_pred <- res_obj$big_data$validation$predicted_values
  
  y_binary <- ifelse(y_true > 0, 1, 0)
  pred_binary <- ifelse(y_pred > threshold, 1, 0)
  
  cm <- suppressMessages(
    suppressWarnings(confusionMatrix(
    factor(pred_binary, levels = c(0, 1)),
    factor(y_binary, levels = c(0, 1)),
    positive = "1"
  )))
  
  roc_obj <- roc(y_binary, y_pred)
  auc_val <- as.numeric(auc(roc_obj))
  
  precision <- cm$byClass["Pos Pred Value"]
  recall <- cm$byClass["Sensitivity"]
  
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA else 2 * (precision * recall) / (precision + recall)
  
  list(
    auc = auc_val,
    accuracy = cm$overall["Accuracy"],
    sensitivity = recall,
    specificity = cm$byClass["Specificity"],
    precision = precision,
    f1 = f1
  )
}


# --------------------------------
# Evaluate performance across training fractions and seeds
# --------------------------------
perf_list_random <- list()

for (frac in train_fracs) {
  for (s in seeds) {
    obj_name <- sprintf("frac%.1f_seed%d", frac, s)
    res_obj <- all_results_random_split[[obj_name]]
    if (is.null(res_obj)) next
    
    perf <- evaluate_performance_random(res_obj)
    
    perf_list_random[[obj_name]] <- data.frame(
      train_frac = frac,
      seed = s,
      auc = perf$auc,
      accuracy = perf$accuracy,
      sensitivity = perf$sensitivity,
      specificity = perf$specificity,
      precision = perf$precision,
      f1 = perf$f1
    )
  }
}



perf_list_random <- bind_rows(perf_list_random)

perf_summary_split_rand <- perf_list_random %>%
  group_by(train_frac) %>%
  summarise(
    mean_auc = mean(auc, na.rm = TRUE),
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    mean_specificity = mean(specificity, na.rm = TRUE),
    mean_precision = mean(precision, na.rm = TRUE),
    mean_f1 = mean(f1, na.rm = TRUE))







############################################################
seeds <- 123:132
thinning_versions <- c("50.0_50.0", "70.0_30.0", "80.0_20.0", "90.0_10.0")
train_frac <- c(0.5, 0.7, 0.8, 0.9)
version_labels <- c("50_50", "70_30", "80_20", "90_10")

roc_data_comb <- data.frame()

# ----- Thinning -----
for (i in 1:4) {
  for (s in seeds) {
    obj_name <- paste0("res_epsilon", thinning_versions[i], "_seed", s)
    res_obj <- all_results_thin[[obj_name]]
    if (is.null(res_obj)) next
    
    y_true <- res_obj$big_data$y_val
    y_pred <- res_obj$big_data$predicted_values
    y_binary <- ifelse(y_true > 0, 1, 0)
    
    roc_obj <- roc(y_binary, y_pred)
    
    roc_df <- data.frame(
      version = version_labels[i],
      seed = as.factor(s),
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Method = "Thinning",
      Approach = "Data thinning"
    )
    roc_data_comb <- bind_rows(roc_data_comb, roc_df)
  }
}

# ----- Stratified Splitting -----
for (i in 1:4) {
  for (s in seeds) {
    obj_name <- sprintf("frac%.1f_seed%d", train_frac[i], s)
    res_obj <- all_results_stratified[[obj_name]]
    if (is.null(res_obj)) next
    
    y_true <- res_obj$big_data$y_val
    y_pred <- res_obj$big_data$predicted_values
    y_binary <- ifelse(y_true > 0, 1, 0)
    
    roc_obj <- roc(y_binary, y_pred)
    
    roc_df <- data.frame(
      version = version_labels[i],
      seed = as.factor(s),
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Method = "Splitting_s",
      Approach = "Data splitting (stratified)"
    )
    roc_data_comb <- bind_rows(roc_data_comb, roc_df)
  }
}

# ----- Random Splitting -----
for (i in 1:4) {
  for (s in seeds) {
    obj_name <- sprintf("frac%.1f_seed%d", train_frac[i], s)
    res_obj <- all_results_random_split[[obj_name]]
    if (is.null(res_obj)) next
    
    y_true <- res_obj$big_data$validation$y
    y_pred <- res_obj$big_data$validation$predicted_values
    y_binary <- ifelse(y_true > 0, 1, 0)
    
    roc_obj <- roc(y_binary, y_pred)
    
    roc_df <- data.frame(
      version = version_labels[i],
      seed = as.factor(s),
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Method = "Splitting_r",
      Approach = "Data splitting (random)"
    )
    roc_data_comb <- bind_rows(roc_data_comb, roc_df)
  }
}

# ----- Overlay -----
for (i in 1:4) {
  for (s in seeds) {
    # Thinning
    obj_t <- paste0("res_epsilon", thinning_versions[i], "_seed", s)
    res_t <- all_results_thin[[obj_t]]
    if (!is.null(res_t)) {
      y_true_t <- res_t$big_data$y_val
      y_pred_t <- res_t$big_data$predicted_values
      roc_t <- roc(ifelse(y_true_t > 0, 1, 0), y_pred_t)
      
      roc_df_t <- data.frame(
        version = version_labels[i],
        seed = as.factor(s),
        FPR = 1 - roc_t$specificities,
        TPR = roc_t$sensitivities,
        Method = "Thinning",
        Approach = "Overlay"
      )
    } else roc_df_t <- NULL
    
    # Stratified
    obj_s <- sprintf("frac%.1f_seed%d", train_frac[i], s)
    res_s <- all_results_stratified[[obj_s]]
    if (!is.null(res_s)) {
      y_true_s <- res_s$big_data$y_val
      y_pred_s <- res_s$big_data$predicted_values
      roc_s <- roc(ifelse(y_true_s > 0, 1, 0), y_pred_s)
      
      roc_df_s <- data.frame(
        version = version_labels[i],
        seed = as.factor(s),
        FPR = 1 - roc_s$specificities,
        TPR = roc_s$sensitivities,
        Method = "Splitting_s",
        Approach = "Overlay"
      )
    } else roc_df_s <- NULL
    
    # Random
    obj_r <- sprintf("frac%.1f_seed%d", train_frac[i], s)
    res_r <- all_results_random_split[[obj_r]]
    if (!is.null(res_r)) {
      y_true_r <- res_r$big_data$validation$y
      y_pred_r <- res_r$big_data$validation$predicted_values
      roc_r <- roc(ifelse(y_true_r > 0, 1, 0), y_pred_r)
      
      roc_df_r <- data.frame(
        version = version_labels[i],
        seed = as.factor(s),
        FPR = 1 - roc_r$specificities,
        TPR = roc_r$sensitivities,
        Method = "Splitting_r",
        Approach = "Overlay"
      )
    } else roc_df_r <- NULL
    
    roc_data_comb <- bind_rows(roc_data_comb, roc_df_t, roc_df_s, roc_df_r)
  }
}





plot_roc_curves <- function(roc_data_comb) {
  
  # Rename methods
  roc_data_comb$Method <- recode(roc_data_comb$Method,
                                 "Thinning" = "Data thinning",
                                 "Splitting_s" = "Splitting (stratified)",
                                 "Splitting_r" = "Splitting (random)")
  
  # Rename final ROC dataset for clarity
  roc_combined_three_splits <- roc_data_comb
  
  # Plot
  ggplot(roc_combined_three_splits, aes(x = FPR, y = TPR, color = Method, group = interaction(seed, Method))) +
    geom_line(size = 0.8, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    facet_grid(Approach ~ version, labeller = labeller(
      version = c(
        "50_50" = "ε = 0.5",
        "70_30" = "ε = 0.7",
        "80_20" = "ε = 0.8",
        "90_10" = "ε = 0.9"
      )
    )) +
    labs(
      title = "ROC Curves: Comparison of Thinning, Random and Stratified Splitting, and Overlay",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Method"
    ) +
    scale_color_manual(values = c(
      "Data thinning" = "#E69F00",
      "Splitting (stratified)" = "#56B4E9",
      "Splitting (random)" = "#009E73"
    )) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      panel.spacing = unit(1.5, "lines")
    )
}








# Performance Evaluation from Three Data Splitting Methods Across Epsilon Values

results_summary <- data.frame(
  epsilon = c(0.5, 0.7, 0.8, 0.9),
  AUC_Th = c(0.879, 0.915, 0.915, 0.917),
  AUC_Sp = c(0.908, 0.932, 0.943, 0.960),
  AUC_Sr = c(0.814, 0.850, 0.865, 0.868),
  
  Accuracy_Th = c(0.919, 0.945, 0.957, 0.974),
  Accuracy_Sp = c(0.926, 0.945, 0.957, 0.973),
  Accuracy_Sr = c(0.879, 0.883, 0.880, 0.883),
  
  Sensitivity_Th = c(0.460, 0.439, 0.333, 0.118),
  Sensitivity_Sp = c(0.521, 0.447, 0.385, 0.239),
  Sensitivity_Sr = c(0.480, 0.510, 0.511, 0.485),
  
  Specificity_Th = c(0.961, 0.976, 0.984, 0.996),
  Specificity_Sp = c(0.962, 0.977, 0.984, 0.994),
  Specificity_Sr = c(0.937, 0.937, 0.935, 0.940),
  
  Precision_Th = c(0.518, 0.525, 0.488, 0.450),
  Precision_Sp = c(0.555, 0.548, 0.528, 0.509),
  Precision_Sr = c(0.526, 0.533, 0.534, 0.536),
  
  F1_Th = c(0.487, 0.477, 0.395, 0.185),
  F1_Sp = c(0.537, 0.492, 0.444, 0.323),
  F1_Sr = c(0.499, 0.520, 0.520, 0.503)
)

plot_performance_metrics <- function(results_summary) {
  
  # Reshape to long format for ggplot
  performance_long <- results_summary %>%
    pivot_longer(-epsilon, names_to = c("Metric","Method"), names_sep = "_")
  
  # Plot
  ggplot(performance_long, aes(x = epsilon, y = value, color = Method, shape = Method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ Metric, scales = "free_y") +
    labs(
      x = "Training Proportion (ε)",
      y = "Metric Value",
      color = "Method",
      shape = "Method",
      title = "Comparison of Performance Metrics Across Methods and Training Proportions"
    ) +
    scale_color_manual(
      values = c("Th" = "blue", "Sp" = "red", "Sr" = "green"),
      labels = c("Th" = "Data thinning ", 
                 "Sp" = "Data splitting (stratified)", 
                 "Sr" = "Data splitting (random)")
    ) +
    scale_shape_manual(
      values = c("Th" = 16, "Sp" = 17, "Sr" = 15),
      labels = c("Th" = "Data thinning ", 
                 "Sp" = "Data splitting (stratified)", 
                 "Sr" = "Data splitting (random)")
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}




