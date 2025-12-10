##########################################
# ===== Plotting Functions for zGPS.AO =====
# These functions generate heatmaps or scatter plots 
# of drug-AE relative risks (lambda_hat) and group-level RR (s) 
# from the zero–inflated Gamma–Poisson shrinker with AE ontology (zGPS.AO) results.
##########################################

# ----------------------------
# Function: plot.zGPS.AO
# Purpose: Main plotting method for zGPS.AO objects
# Arguments:
#   x: A zGPS.AO object returned by get_all_params()
#   drug_name: Optional, a specific drug to visualize
#   AE_grp_name: Optional, a specific AE group to visualize
#   interactive_plot: TRUE/FALSE for interactive plot (default FALSE)
#   ... : Additional plotting arguments passed to internal ggplot calls
# Returns: ggplot or plotly object
# ----------------------------


plot.zGPS.AO = function(x,
                     drug_name = NULL,
                     AE_grp_name = NULL,
                     interactive_plot = FALSE,
                     ...) {
  big_data = x$big_data
  cond = is.null(drug_name) & is.null(AE_grp_name)
  

  if (cond) {
    p = plot_heatmap(big_data)
  } else {
    p = plot_lambdas(drug_name, AE_grp_name, big_data)
  }
  
  if (interactive_plot) {
    p = ggplotly(p, tooltip = "text") %>%
      layout(
        yaxis = list(tickfont = list(size = 28)), 
        xaxis = list(tickfont = list(size = 26))   
      )
  }
  
  return(p)
}



# ----------------------------
# Function: plot_heatmap
# Purpose: Generate a heatmap of group-level RR (s) across drugs and AE groups
# Arguments:
#   big_data: Data.frame from zGPS.AO containing s and s_pval
# Returns: ggplot heatmap
# ----------------------------


plot_heatmap = function(big_data) {
  big_data %>%
    group_by(drug, AE_grp) %>%
    summarize(s = mean(s),
              q_val = mean(s_pval)) %>%
    ungroup() %>%
    mutate(
      text = paste0(
        'AE_grp: ',
        AE_grp,
        '\n',
        'drug: ',
        drug,
        '\n',
        's: ',
        round(s, 4),
        '\n'
        ,
        'q: ',
        round(q_val, 8)
      )
    ) %>%
    ggplot(aes(drug, AE_grp, fill = s, text = text)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme_ipsum() +
    labs(
      x = "Drugs",
      y = "AE groups",  # lowercase 'g'
      fill = "S"
    ) +
    theme(
      axis.title.y = element_text(size = 16, angle = 0, vjust = 0.5),
      axis.title.x = element_text(size = 16, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) ->
    p
  
  return(p)
}




# ----------------------------
# Function: plot_lambdas
# Purpose: Generate scatter plot of lambda_hat for a specific drug × AE group
# Arguments:
#   drug_name: Character, the drug to visualize
#   AE_grp_name: Character, the AE group to visualize
#   big_data: Data.frame from zGPS.AO containing lambda_hat, y, E, s, lambda_pval
# Returns: ggplot object
# ----------------------------



plot_lambdas = function(drug_name, AE_grp_name, big_data) {
  # browser()
  big_data %>%
    filter(drug == drug_name &
             AE_grp == AE_grp_name, ) %>%
    mutate(text = paste0(
      'y: ',
      round(y, 4) ,
      '\n',
      'E: ',
      round(E, 4),
      '\n',
      'lambda_hat: ',
      round(lambda_hat, 4),
      '\n',
      'q_val: ',
      round(lambda_pval,4)
    )) ->
    data
  
  data %>%
    ggplot(aes(AE, lambda_hat, text = text)) +
    geom_point(color = 'red', size = 2) +
    coord_flip() +
    geom_hline(aes(yintercept = mean(s), text =)) +
    theme_bw() +
    xlab(NULL) + 
    ggtitle(paste0(drug_name, ', ', AE_grp_name, ' (s=', round(mean(data$s), 2), ')')) ->
    p
  
  return(p)
}

# ----------------------------
# End of plotting functions for zGPS.AO
# ----------------------------
