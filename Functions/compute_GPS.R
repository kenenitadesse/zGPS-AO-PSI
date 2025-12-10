# =============================================================================
# Function: compute_GPS
# Description: Fit GPS (Gamma–Poisson Shrinkage) model on zGPS.AO results
# =============================================================================

compute_GPS <- function(zGPS.AO_res) {
  
  # Required libraries
  if(!require(openEBGM)) install.packages("openEBGM"); library(openEBGM)
  if(!require(dplyr)) install.packages("dplyr"); library(dplyr)
  

# ==========================================
# Define Custom Negative Binomial Mixture Likelihood
# ==========================================

# The GPS model is estimated using a two-component Negative Binomial mixture
# prior. The function below evaluates the mixture density for counts.


dbinbinom <- function(x, size1, prob1, size2, prob2, w) {
  w * dnbinom(x, size = size1, prob = prob1) +
    (1 - w) * dnbinom(x, size = size2, prob = prob2)
}

# Log-likelihood function for Maximum Likelihood Estimation.
# The parameters represent (alpha1, beta1, alpha2, beta2, mixture weight).
loglikelihood2NegativeBinomial <- function(p, a, E) {
  alpha1 <- p[1]
  beta1  <- p[2]
  alpha2 <- p[3]
  beta2  <- p[4]
  w      <- p[5]
  
  -sum(log(
    dbinbinom(a,
              size1 = alpha1,
              prob1 = beta1 / (beta1 + E),
              size2 = alpha2,
              prob2 = beta2 / (beta2 + E),
              w = w)
  ))
}

# -------------------------------
# Fit GPS Model
# -------------------------------
gps_input_data <- zGPS.AO_res$big_data

  a <- gps_input_data$y
  E <- gps_input_data$E
  
  # Initial values for parameters: (alpha1, beta1, alpha2, beta2, weight)
  init_params <- c(0.2, 0.1, 2, 4, 1/3)
  
  # Maximum Likelihood Estimation of the mixture prior using nlminb
  opt <- suppressWarnings(
    nlminb(
      start   = init_params,
      objective = loglikelihood2NegativeBinomial,
      a = a, E = E,
      lower = c(1e-4, 1e-4, 1e-4, 1e-4, 1e-4),
      upper = c(Inf, Inf, Inf, Inf, 1 - 1e-10)
    )
  )
  
  # Extract estimated parameters
  theta_hat <- opt$par
  
  # Posterior quantities: Qn, EBGM, and the 1% and 99% credibility bounds
  qn_vals     <- Qn(theta_hat, N = a, E = E)
  ebgm_vals   <- ebgm(theta_hat, N = a, E = E, qn = qn_vals)
  quan01 <- quantBisect(1, theta_hat, N = a, E = E, qn = qn_vals)
  quan99 <- quantBisect(99, theta_hat, N = a, E = E, qn = qn_vals)
  
  # -------------------------------
  # Return Results
  # -------------------------------
  
  result_df <-data.frame(
    AE_grp = gps_input_data$AE_grp,
    drug = gps_input_data$drug,
    AE = gps_input_data$AE,
    y = a,
    E = E,
    EBGM = ebgm_vals,
    quan01 = quan01,
    quan99 = quan99,
    stringsAsFactors = FALSE
  )
  return(result_df)
}

