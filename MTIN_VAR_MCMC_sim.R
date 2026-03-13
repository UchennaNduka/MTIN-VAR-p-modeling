# ============================================
# Bayesian MCMC for MTIN-VAR(p) Estimation
# ============================================

# Function to generate MTIN-distributed errors
rmtin <- function(n, mu, Omega, theta) {
  library(mvtnorm)
  d <- length(mu)
  w <- runif(n, 1 - theta, 1)
  t(sapply(w, function(wi) rmvnorm(1, mu, Omega/wi)))  
}

# Function to generate multiple samples of MTIN-VAR(p) processes
simulate_mtin_varp_data <- function(d, p, n, n_sim, theta_true, seed = 123) {
  # --- True Parameters ---
  set.seed(seed)
  intercept <- matrix(runif(d, -0.5, 0.5), nrow = d)
  AR_coefs <- array(runif(d * d * p, -0.5, 0.5), dim = c(d, d * p))
  true_Psi <- cbind(intercept, AR_coefs)
  
  true_Omega <- diag(d)
  true_Omega[lower.tri(true_Omega)] <- 0.3
  true_Omega <- true_Omega + t(true_Omega) - diag(diag(true_Omega))
  
  # --- Storage ---
  sim_results <- array(NA, dim = c(n, d, n_sim))
  
  # --- Simulation Loop ---
  for (i in 1:n_sim) {
    set.seed(seed + i)
    Y <- matrix(0, n, d)
    Y[1:p, ] <- matrix(rnorm(d * p), nrow = p)
    for (t in (p + 1):n) {
      lags <- as.vector(Y[(t - p):(t - 1), ])
      x_prev <- c(1, lags)
      if (length(x_prev) != ncol(true_Psi)) stop("Mismatch: x_prev and true_Psi dimensions")
      Y[t, ] <- as.vector(true_Psi %*% x_prev) + t(rmtin(1, rep(0, d), true_Omega, theta_true))
    }
    sim_results[,,i] <- Y
  }
  
  return(list(
    data = sim_results,
    true_Psi = true_Psi,
    true_Omega = true_Omega,
    true_theta = theta_true
  ))
}

# Log-likelihood function for VAR(p) with tail-inflated normal innovations
log_likelihood <- function(Y, p, Psi, Omega, theta) {
  n <- nrow(Y)
  d <- ncol(Y)
  inv_Omega <- solve(Omega)
  log_det_Omega <- log(det(Omega))
  quad_forms <- numeric(n - p)
  gamma_diffs <- numeric(n - p)
  
  for (t in (p + 1):n) {
    # Construct x_{t-1} = [1, y_{t-1}, ..., y_{t-p}]
    lag_indices <- (t - 1):(t - p)
    x_prev <- c(1, as.vector(t(Y[lag_indices, ])))  # Flatten row-wise
    x_prev <- matrix(x_prev, ncol = 1)
    
    # Compute residual
    y_t <- matrix(Y[t, ], ncol = 1)
    residual <- y_t - Psi %*% x_prev
    
    # Quadratic form
    quad_form <- as.numeric(t(residual) %*% inv_Omega %*% residual)
    quad_forms[t - p] <- quad_form
    
    # Gamma terms
    s <- d / 2 + 1
    x1 <- (1 - theta) * quad_form / 2
    x2 <- quad_form / 2
    gamma_term1 <- gamma(s) * pgamma(x1, shape = s, lower.tail = FALSE)
    gamma_term2 <- gamma(s) * pgamma(x2, shape = s, lower.tail = FALSE)
    gamma_diff <- gamma_term1 - gamma_term2
    
    gamma_diffs[t - p] <- max(gamma_diff, 1e-300)  # Avoid log(0)
  }
  
  # Log-likelihood components
  term1 <- -n/2 * log_det_Omega - log(theta)
  term2 <- -(d/2 + 1) * sum(log(quad_forms))
  term3 <- sum(log(gamma_diffs))
  
  term1 + term2 + term3
}

# ============================================
# Bayesian MCMC Estimation Function
# ============================================

bayesian_mcmc_estimator <- function(Y, p, true_Psi, true_Omega, true_theta,
                                    n_iter = 5000, burn_in = 1000,
                                    prior_params = list(
                                      Psi_sd = 1.0,
                                      Omega_nu = NULL,
                                      Omega_S = NULL,
                                      theta_a = 2,
                                      theta_b = 2
                                    ),
                                    prop_sd = list(Psi = 0.05, Omega = 0.1, theta = 0.05),
                                    verbose = FALSE, adapt = TRUE) {
  
  library(mvtnorm)
  
  # Start timing
  start_time <- proc.time()[1]
  
  d <- ncol(Y)
  n <- nrow(Y)
  T <- n - p
  
  # Dimensions of Psi
  k_psi <- d * (d * p + 1)
  
  # Default prior parameters
  if (is.null(prior_params$Omega_nu)) prior_params$Omega_nu <- d + 2
  if (is.null(prior_params$Omega_S)) prior_params$Omega_S <- diag(d) * (d + 2)
  
  # ============================================
  # Helper functions for Bayesian estimation
  # ============================================
  
  # Log-prior for Psi (matrix normal with independence)
  log_prior_Psi <- function(Psi) {
    sum(dnorm(as.vector(Psi), 0, prior_params$Psi_sd, log = TRUE))
  }
  
  # Log-prior for Omega (Inverse-Wishart)
  log_prior_Omega <- function(Omega) {
    nu <- prior_params$Omega_nu
    S <- prior_params$Omega_S
    
    # Check positive definiteness
    if (any(eigen(Omega, symmetric = TRUE)$values <= 0)) {
      return(-Inf)
    }
    
    p <- nrow(Omega)
    log_det_Omega <- log(det(Omega))
    inv_Omega <- solve(Omega)
    
    # Inverse-Wishart log-density
    constant <- (nu * p / 2) * log(2) + (p * (p - 1) / 4) * log(pi) +
      sum(lgamma((nu + 1 - 1:p)/2))
    log_density <- constant - (nu/2) * log(det(S)) - 
      ((nu + p + 1)/2) * log_det_Omega - 
      0.5 * sum(diag(S %*% inv_Omega))
    
    return(log_density)
  }
  
  # Log-prior for theta (Beta)
  log_prior_theta <- function(theta) {
    dbeta(theta, prior_params$theta_a, prior_params$theta_b, log = TRUE)
  }
  
  # Log-posterior (un-normalized)
  log_posterior <- function(Psi, Omega, theta) {
    log_lik <- log_likelihood(Y, p, Psi, Omega, theta)
    log_prior <- log_prior_Psi(Psi) + log_prior_Omega(Omega) + log_prior_theta(theta)
    
    if (!is.finite(log_lik) || !is.finite(log_prior)) {
      return(-Inf)
    }
    
    return(log_lik + log_prior)
  }
  
  # ============================================
  # Initialize parameters
  # ============================================
  
  # Initialize Psi using OLS
  X_ols <- matrix(1, nrow = T, ncol = 1)
  for (i in 1:p) {
    X_ols <- cbind(X_ols, Y[(p-i+1):(n-i), ])
  }
  Y_ols <- Y[(p+1):n, ]
  
  # Regularized OLS
  lambda <- 1e-6
  XtX <- t(X_ols) %*% X_ols
  Psi_current <- solve(XtX + lambda * diag(ncol(XtX))) %*% t(X_ols) %*% Y_ols
  Psi_current <- t(Psi_current)
  
  # Initialize Omega from OLS residuals
  residuals_ols <- Y_ols - X_ols %*% t(Psi_current)
  Omega_current <- (t(residuals_ols) %*% residuals_ols) / T
  Omega_current <- 0.5 * (Omega_current + t(Omega_current))
  Omega_current <- Omega_current + 1e-6 * diag(d)
  
  # Initialize theta
  theta_current <- 0.5
  
  # Storage
  samples <- list(
    Psi = array(NA, dim = c(d, d * p + 1, n_iter)),
    Omega = array(NA, dim = c(d, d, n_iter)),
    theta = numeric(n_iter),
    log_post = numeric(n_iter),
    accept = list(Psi = 0, Omega = 0, theta = 0)
  )
  
  # Initial log-posterior
  current_log_post <- log_posterior(Psi_current, Omega_current, theta_current)
  
  # Adaptive proposal tuning
  adapt_every <- 100
  n_accept_Psi <- 0
  n_accept_Omega <- 0
  n_accept_theta <- 0
  
  # ============================================
  # MCMC Sampling
  # ============================================
  
  if (verbose) cat("Starting MCMC sampling...\n")
  
  for (iter in 1:n_iter) {
    
    if (verbose && iter %% 500 == 0) {
      cat(sprintf("Iteration %d/%d, theta = %.4f, log-post = %.2f\n", 
                  iter, n_iter, theta_current, current_log_post))
    }
    
    # 1. Update Psi (Metropolis-Hastings)
    Psi_proposal <- Psi_current + 
      matrix(rnorm(k_psi, 0, prop_sd$Psi), nrow = d, ncol = d * p + 1)
    
    log_post_proposal <- log_posterior(Psi_proposal, Omega_current, theta_current)
    log_accept <- log_post_proposal - current_log_post
    
    if (is.finite(log_accept) && log(runif(1)) < min(log_accept, 0)) {
      Psi_current <- Psi_proposal
      current_log_post <- log_post_proposal
      samples$accept$Psi <- samples$accept$Psi + 1
      n_accept_Psi <- n_accept_Psi + 1
    }
    
    # 2. Update Omega (Metropolis-Hastings with symmetric proposal)
    # Generate symmetric perturbation
    perturb <- matrix(rnorm(d^2, 0, prop_sd$Omega), d, d)
    perturb <- (perturb + t(perturb)) / 2
    Omega_proposal <- Omega_current + perturb
    
    # Ensure positive definiteness
    if (all(eigen(Omega_proposal, symmetric = TRUE)$values > 0)) {
      log_post_proposal <- log_posterior(Psi_current, Omega_proposal, theta_current)
      log_accept <- log_post_proposal - current_log_post
      
      if (is.finite(log_accept) && log(runif(1)) < min(log_accept, 0)) {
        Omega_current <- Omega_proposal
        current_log_post <- log_post_proposal
        samples$accept$Omega <- samples$accept$Omega + 1
        n_accept_Omega <- n_accept_Omega + 1
      }
    }
    
    # 3. Update theta (Metropolis-Hastings)
    theta_proposal <- theta_current + rnorm(1, 0, prop_sd$theta)
    
    # Reflect at boundaries
    if (theta_proposal < 0) theta_proposal <- -theta_proposal
    if (theta_proposal > 1) theta_proposal <- 2 - theta_proposal
    theta_proposal <- max(min(theta_proposal, 0.999), 0.001)
    
    log_post_proposal <- log_posterior(Psi_current, Omega_current, theta_proposal)
    log_accept <- log_post_proposal - current_log_post
    
    if (is.finite(log_accept) && log(runif(1)) < min(log_accept, 0)) {
      theta_current <- theta_proposal
      current_log_post <- log_post_proposal
      samples$accept$theta <- samples$accept$theta + 1
      n_accept_theta <- n_accept_theta + 1
    }
    
    # Store samples
    samples$Psi[, , iter] <- Psi_current
    samples$Omega[, , iter] <- Omega_current
    samples$theta[iter] <- theta_current
    samples$log_post[iter] <- current_log_post
    
    # Adaptive tuning (during burn-in only)
    if (adapt && iter < burn_in && iter %% adapt_every == 0) {
      # Adjust proposal scales based on acceptance rates
      accept_rate_Psi <- n_accept_Psi / adapt_every
      accept_rate_Omega <- n_accept_Omega / adapt_every
      accept_rate_theta <- n_accept_theta / adapt_every
      
      # Target acceptance rate: 0.234 for random walk
      if (accept_rate_Psi < 0.15) prop_sd$Psi <- prop_sd$Psi * 0.9
      if (accept_rate_Psi > 0.35) prop_sd$Psi <- prop_sd$Psi * 1.1
      prop_sd$Psi <- max(prop_sd$Psi, 1e-6)
      
      if (accept_rate_Omega < 0.15) prop_sd$Omega <- prop_sd$Omega * 0.9
      if (accept_rate_Omega > 0.35) prop_sd$Omega <- prop_sd$Omega * 1.1
      prop_sd$Omega <- max(prop_sd$Omega, 1e-6)
      
      if (accept_rate_theta < 0.15) prop_sd$theta <- prop_sd$theta * 0.9
      if (accept_rate_theta > 0.35) prop_sd$theta <- prop_sd$theta * 1.1
      prop_sd$theta <- max(prop_sd$theta, 1e-6)
      
      if (verbose) {
        cat(sprintf("Adaptation (iter %d): Psi=%.3f, Omega=%.3f, theta=%.3f\n",
                    iter, accept_rate_Psi, accept_rate_Omega, accept_rate_theta))
      }
      
      # Reset counters
      n_accept_Psi <- 0
      n_accept_Omega <- 0
      n_accept_theta <- 0
    }
  }
  
  # End timing
  end_time <- proc.time()[1]
  runtime <- end_time - start_time
  
  # ============================================
  # Post-processing
  # ============================================
  
  # Discard burn-in
  keep <- (burn_in + 1):n_iter
  samples$Psi_post <- samples$Psi[, , keep]
  samples$Omega_post <- samples$Omega[, , keep]
  samples$theta_post <- samples$theta[keep]
  samples$log_post_post <- samples$log_post[keep]
  
  # Calculate posterior means
  post_mean_Psi <- apply(samples$Psi_post, 1:2, mean)
  post_mean_Omega <- apply(samples$Omega_post, 1:2, mean)
  post_mean_theta <- mean(samples$theta_post)
  
  # Calculate posterior standard deviations
  post_sd_Psi <- apply(samples$Psi_post, 1:2, sd)
  post_sd_Omega <- apply(samples$Omega_post, 1:2, sd)
  post_sd_theta <- sd(samples$theta_post)
  
  # Calculate 95% credible intervals
  post_ci_Psi <- apply(samples$Psi_post, 1:2, quantile, probs = c(0.025, 0.975))
  post_ci_Omega <- apply(samples$Omega_post, 1:2, quantile, probs = c(0.025, 0.975))
  post_ci_theta <- quantile(samples$theta_post, probs = c(0.025, 0.975))
  
  # Calculate estimation errors
  Error_Psi <- post_mean_Psi - true_Psi
  Error_Omega <- post_mean_Omega - true_Omega
  
  frob_Psi_errors <- tryCatch({
    norm(Error_Psi, type = "F")
  }, error = function(e) NA)
  
  frob_Omega_errors <- tryCatch({
    norm(Error_Omega, type = "F")
  }, error = function(e) NA)
  
  theta_errors <- post_mean_theta - true_theta
  
  # Calculate acceptance rates
  accept_rate_Psi <- samples$accept$Psi / n_iter
  accept_rate_Omega <- samples$accept$Omega / n_iter
  accept_rate_theta <- samples$accept$theta / n_iter
  
  # Effective sample size for theta (simplified)
  ess_theta <- tryCatch({
    acf_theta <- acf(samples$theta_post, plot = FALSE, lag.max = 100)
    mcmc_size <- length(samples$theta_post)
    mcmc_size / (1 + 2 * sum(acf_theta$acf[-1]))
  }, error = function(e) NA)
  
  # ============================================
  # Return results
  # ============================================
  
  results <- list(
    # Point estimates
    Psi_hat = post_mean_Psi,
    Omega_hat = post_mean_Omega,
    theta_hat = post_mean_theta,
    
    # Uncertainty estimates
    Psi_sd = post_sd_Psi,
    Omega_sd = post_sd_Omega,
    theta_sd = post_sd_theta,
    
    # Credible intervals
    Psi_ci = post_ci_Psi,
    Omega_ci = post_ci_Omega,
    theta_ci = post_ci_theta,
    
    # Estimation errors
    Psi_norm = frob_Psi_errors,
    Omega_norm = frob_Omega_errors,
    theta_error = theta_errors,
    
    # MCMC diagnostics
    runtime = runtime,
    n_iter = n_iter,
    burn_in = burn_in,
    n_post = length(keep),
    accept_rate_Psi = accept_rate_Psi,
    accept_rate_Omega = accept_rate_Omega,
    accept_rate_theta = accept_rate_theta,
    ess_theta = ess_theta,
    
    # Raw samples (optional, can be large)
    samples = if (verbose) samples else NULL,
    
    # Convergence diagnostics
    log_post_trace = samples$log_post,
    theta_trace = samples$theta
  )
  
  return(results)
}

# ============================================
# Main Simulation Study for Bayesian MCMC
# ============================================

library(doParallel)
library(foreach)

# Define simulation settings
sim_settings <- data.frame(
  id = c(
    "bayesian_s1", "bayesian_s2", "bayesian_s3", "bayesian_s4", "bayesian_s5", "bayesian_s6",
    "bayesian_s7", "bayesian_s8", "bayesian_s9", "bayesian_s10", "bayesian_s11", "bayesian_s12",
    "bayesian_s13", "bayesian_s14", "bayesian_s15", "bayesian_s16", "bayesian_s17", "bayesian_s18",
    "bayesian_s19", "bayesian_s20", "bayesian_s21", "bayesian_s22", "bayesian_s23", "bayesian_s24",
    "bayesian_s25", "bayesian_s26", "bayesian_s27"
  ),
  p = rep(2, 27),
  d = c(
    2, 3, 5, 2, 3, 5, 2, 3, 5,  # n=200
    2, 3, 5, 2, 3, 5, 2, 3, 5,  # n=500
    2, 3, 5, 2, 3, 5, 2, 3, 5   # n=1000
  ),
  theta = c(
    0.7, 0.7, 0.7, 0.8, 0.8, 0.8, 0.9, 0.9, 0.9,  # n=200
    0.7, 0.7, 0.7, 0.8, 0.8, 0.8, 0.9, 0.9, 0.9,  # n=500
    0.7, 0.7, 0.7, 0.8, 0.8, 0.8, 0.9, 0.9, 0.9   # n=1000
  ),
  n = c(
    rep(200, 9), rep(500, 9), rep(1000, 9)
  ),
  stringsAsFactors = FALSE
)

# Add n_sim column (number of simulations per setting)
sim_settings$n_sim <- 100

cat("============================================\n")
cat("Bayesian MCMC for MTIN-VAR(p) Simulation Study\n")
cat("============================================\n")
cat("Total settings to run:", nrow(sim_settings), "\n")
cat("Total simulations:", sum(sim_settings$n_sim), "\n")
cat("============================================\n")

# Create directory for results
results_dir <- "new_bayesian_mcmc_results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Setup parallel backend
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel MCMC estimation\n")

# Master results list
all_results <- list()
summary_results <- data.frame()

# Loop through each setting
for (setting_idx in 1:nrow(sim_settings)) {
  setting <- sim_settings[setting_idx, ]
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("Running Setting", setting_idx, "/", nrow(sim_settings), ":", setting$id, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("Parameters: p =", setting$p, ", d =", setting$d, 
      ", theta =", setting$theta, ", n =", setting$n, "\n")
  cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  # === Simulate Data for this setting ===
  cat("Generating simulated data...\n")
  set.seed(123 + setting_idx)
  simulated_data <- simulate_mtin_varp_data(
    d = setting$d, 
    p = setting$p, 
    n = setting$n, 
    n_sim = setting$n_sim, 
    theta_true = setting$theta
  )
  cat("Data generation completed.\n")
  
  # === Parallel Bayesian MCMC Estimation ===
  cat("Starting parallel Bayesian MCMC estimation...\n")
  start_time_setting <- proc.time()[1]
  
  estimation_results <- foreach(i = 1:setting$n_sim,
                                .packages = c("mvtnorm"),
                                .export = c("bayesian_mcmc_estimator", "rmtin", "log_likelihood"),
                                .combine = 'rbind') %dopar% {
                                  Y_i <- simulated_data$data[,,i]
                                  result <- tryCatch({
                                    # Run Bayesian MCMC
                                    mcmc_result <- bayesian_mcmc_estimator(
                                      Y = Y_i,
                                      p = setting$p,
                                      true_Psi = simulated_data$true_Psi,
                                      true_Omega = simulated_data$true_Omega,
                                      true_theta = setting$theta,
                                      n_iter = 5000,
                                      burn_in = 1000,
                                      verbose = FALSE,
                                      adapt = TRUE
                                    )
                                    
                                    # Extract relevant information for results data frame
                                    data.frame(
                                      id = setting$id,
                                      simulation = i,
                                      Psi_norm = mcmc_result$Psi_norm,
                                      Omega_norm = mcmc_result$Omega_norm,
                                      theta_error = mcmc_result$theta_error,
                                      theta_hat = mcmc_result$theta_hat,
                                      theta_sd = mcmc_result$theta_sd,
                                      runtime = mcmc_result$runtime,
                                      n_iter = mcmc_result$n_iter,
                                      burn_in = mcmc_result$burn_in,
                                      n_post = mcmc_result$n_post,
                                      accept_rate_Psi = mcmc_result$accept_rate_Psi,
                                      accept_rate_Omega = mcmc_result$accept_rate_Omega,
                                      accept_rate_theta = mcmc_result$accept_rate_theta,
                                      ess_theta = mcmc_result$ess_theta,
                                      theta_ci_lower = mcmc_result$theta_ci[1],
                                      theta_ci_upper = mcmc_result$theta_ci[2],
                                      stringsAsFactors = FALSE
                                    )
                                  }, error = function(e) {
                                    cat(sprintf("Error in simulation %d: %s\n", i, e$message))
                                    NULL
                                  })
                                  result
                                }
  
  end_time_setting <- proc.time()[1]
  setting_runtime <- end_time_setting - start_time_setting
  
  cat("Parallel MCMC estimation completed for setting", setting$id, "\n")
  cat("Setting runtime:", round(setting_runtime, 2), "seconds\n")
  
  # === Process results for this setting ===
  valid_results <- estimation_results[complete.cases(estimation_results), ]
  n_valid <- nrow(valid_results)
  
  cat("Successful MCMC estimations:", n_valid, "/", setting$n_sim, "\n")
  
  if (n_valid > 0) {
    # Add setting parameters to results
    valid_results$p <- setting$p
    valid_results$d <- setting$d
    valid_results$true_theta <- setting$theta
    valid_results$n <- setting$n
    
    # Store detailed results
    all_results[[setting$id]] <- valid_results
    
    # Calculate summary statistics for this setting
    setting_summary <- data.frame(
      id = setting$id,
      p = setting$p,
      d = setting$d,
      true_theta = setting$theta,
      n = setting$n,
      n_sim = setting$n_sim,
      n_valid = n_valid,
      success_rate = n_valid / setting$n_sim * 100,
      
      # Parameter estimation errors
      mean_Psi_norm = mean(valid_results$Psi_norm, na.rm = TRUE),
      sd_Psi_norm = sd(valid_results$Psi_norm, na.rm = TRUE),
      mean_Omega_norm = mean(valid_results$Omega_norm, na.rm = TRUE),
      sd_Omega_norm = sd(valid_results$Omega_norm, na.rm = TRUE),
      
      # Theta estimation
      mean_theta_hat = mean(valid_results$theta_hat, na.rm = TRUE),
      sd_theta_hat = sd(valid_results$theta_hat, na.rm = TRUE),
      mean_theta_sd = mean(valid_results$theta_sd, na.rm = TRUE),
      theta_bias = mean(valid_results$theta_hat, na.rm = TRUE) - setting$theta,
      theta_rmse = sqrt(mean(valid_results$theta_error^2, na.rm = TRUE)),
      coverage = mean(valid_results$theta_ci_lower <= setting$theta & 
                        setting$theta <= valid_results$theta_ci_upper, na.rm = TRUE) * 100,
      
      # Algorithm performance
      mean_runtime = mean(valid_results$runtime, na.rm = TRUE),
      sd_runtime = sd(valid_results$runtime, na.rm = TRUE),
      mean_ess_theta = mean(valid_results$ess_theta, na.rm = TRUE),
      sd_ess_theta = sd(valid_results$ess_theta, na.rm = TRUE),
      
      # MCMC diagnostics
      mean_accept_Psi = mean(valid_results$accept_rate_Psi, na.rm = TRUE),
      sd_accept_Psi = sd(valid_results$accept_rate_Psi, na.rm = TRUE),
      mean_accept_Omega = mean(valid_results$accept_rate_Omega, na.rm = TRUE),
      sd_accept_Omega = sd(valid_results$accept_rate_Omega, na.rm = TRUE),
      mean_accept_theta = mean(valid_results$accept_rate_theta, na.rm = TRUE),
      sd_accept_theta = sd(valid_results$accept_rate_theta, na.rm = TRUE),
      
      setting_runtime = setting_runtime,
      stringsAsFactors = FALSE
    )
    
    summary_results <- rbind(summary_results, setting_summary)
    
    # Save individual setting results
    filename <- paste0(results_dir, "/", setting$id, "_detailed.csv")
    write.csv(valid_results, filename, row.names = FALSE)
    cat("Detailed results saved to:", filename, "\n")
    
    # Display setting summary
    cat("\n--- Setting Summary ---\n")
    cat("Mean theta estimate:", round(mean(valid_results$theta_hat, na.rm = TRUE), 6), 
        "(True:", setting$theta, ")\n")
    cat("Bias:", round(mean(valid_results$theta_hat, na.rm = TRUE) - setting$theta, 6), "\n")
    cat("RMSE:", round(sqrt(mean(valid_results$theta_error^2, na.rm = TRUE)), 6), "\n")
    cat("Coverage:", round(setting_summary$coverage[1], 2), "%\n")
    cat("Mean runtime:", round(mean(valid_results$runtime, na.rm = TRUE), 4), "seconds\n")
    cat("Mean ESS (theta):", round(mean(valid_results$ess_theta, na.rm = TRUE), 2), "\n")
    cat("Mean acceptance rates: Psi=", round(setting_summary$mean_accept_Psi[1], 3),
        ", Omega=", round(setting_summary$mean_accept_Omega[1], 3),
        ", theta=", round(setting_summary$mean_accept_theta[1], 3), "\n", sep = "")
    
  } else {
    cat("WARNING: No valid MCMC estimations for setting", setting$id, "\n")
  }
  
  cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Setting", setting$id, "completed.\n")
}

# === Stop Cluster ===
stopCluster(cl)
cat("\nAll Bayesian MCMC simulations completed.\n")

# === Save Combined Results ===
if (length(all_results) > 0) {
  # Combine all detailed results
  all_detailed_results <- do.call(rbind, all_results)
  write.csv(all_detailed_results, 
            paste0(results_dir, "/all_detailed_results.csv"), 
            row.names = FALSE)
  cat("\nAll detailed results saved to: ", results_dir, "/all_detailed_results.csv\n", sep = "")
  
  # Save summary results
  write.csv(summary_results, 
            paste0(results_dir, "/summary_results.csv"), 
            row.names = FALSE)
  cat("Summary results saved to: ", results_dir, "/summary_results.csv\n", sep = "")
  
  # Create a formatted summary table
  formatted_summary <- summary_results[, c(
    "id", "p", "d", "true_theta", "n", "n_valid", "success_rate",
    "mean_theta_hat", "theta_bias", "theta_rmse", "coverage",
    "mean_theta_sd", "mean_ess_theta",
    "mean_runtime", "mean_accept_Psi", "mean_accept_Omega", "mean_accept_theta"
  )]
  
  # Round numeric columns
  numeric_cols <- sapply(formatted_summary, is.numeric)
  formatted_summary[, numeric_cols] <- round(formatted_summary[, numeric_cols], 6)
  formatted_summary$coverage <- round(formatted_summary$coverage, 2)
  formatted_summary$success_rate <- round(formatted_summary$success_rate, 2)
  
  write.csv(formatted_summary, 
            paste0(results_dir, "/formatted_summary.csv"), 
            row.names = FALSE)
  cat("Formatted summary saved to: ", results_dir, "/formatted_summary.csv\n", sep = "")
  
  # Display overall summary
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("OVERALL BAYESIAN MCMC SIMULATION STUDY SUMMARY\n")
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("Total settings:", nrow(sim_settings), "\n")
  cat("Total simulations attempted:", sum(sim_settings$n_sim), "\n")
  cat("Total successful simulations:", nrow(all_detailed_results), "\n")
  cat("Overall success rate:", 
      round(nrow(all_detailed_results) / sum(sim_settings$n_sim) * 100, 2), "%\n")
  
  # Overall performance metrics
  overall_theta_bias <- mean(all_detailed_results$theta_error, na.rm = TRUE)
  overall_theta_rmse <- sqrt(mean(all_detailed_results$theta_error^2, na.rm = TRUE))
  overall_coverage <- mean(all_detailed_results$theta_ci_lower <= all_detailed_results$true_theta & 
                             all_detailed_results$true_theta <= all_detailed_results$theta_ci_upper, na.rm = TRUE) * 100
  
  cat("\nOverall Performance Metrics:\n")
  cat("  Mean theta bias:", round(overall_theta_bias, 6), "\n")
  cat("  Mean theta RMSE:", round(overall_theta_rmse, 6), "\n")
  cat("  Overall coverage:", round(overall_coverage, 2), "%\n")
  cat("  Mean runtime:", round(mean(all_detailed_results$runtime, na.rm = TRUE), 4), "seconds\n")
  cat("  Mean ESS (theta):", round(mean(all_detailed_results$ess_theta, na.rm = TRUE), 2), "\n")
  
  # Display summary table
  cat("\nSummary Table (first 10 settings):\n")
  print(head(formatted_summary, 10))
  
  # Create a comprehensive report
  report_file <- paste0(results_dir, "/bayesian_mcmc_report.txt")
  sink(report_file)
  
  cat("BAYESIAN MCMC FOR MTIN-VAR(p) SIMULATION STUDY REPORT\n")
  cat("=======================================================\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("\nTotal settings:", nrow(sim_settings), "\n")
  cat("Total simulations:", sum(sim_settings$n_sim), "\n")
  cat("Successful simulations:", nrow(all_detailed_results), "\n")
  cat("Success rate:", round(nrow(all_detailed_results) / sum(sim_settings$n_sim) * 100, 2), "%\n")
  cat("\n")
  
  # MCMC diagnostics summary
  cat("MCMC Diagnostics Summary:\n")
  cat("  Mean acceptance rate (Psi):", round(mean(all_detailed_results$accept_rate_Psi, na.rm = TRUE), 3), "\n")
  cat("  Mean acceptance rate (Omega):", round(mean(all_detailed_results$accept_rate_Omega, na.rm = TRUE), 3), "\n")
  cat("  Mean acceptance rate (theta):", round(mean(all_detailed_results$accept_rate_theta, na.rm = TRUE), 3), "\n")
  cat("  Mean ESS (theta):", round(mean(all_detailed_results$ess_theta, na.rm = TRUE), 2), "\n")
  cat("\n")
  
  # By sample size
  cat("Performance by sample size (n):\n")
  for (n_val in unique(sim_settings$n)) {
    n_results <- all_detailed_results[all_detailed_results$n == n_val, ]
    if (nrow(n_results) > 0) {
      cat("  n =", n_val, ": Mean theta RMSE =", 
          round(sqrt(mean(n_results$theta_error^2, na.rm = TRUE)), 6),
          ", Coverage =", round(mean(n_results$theta_ci_lower <= n_results$true_theta & 
                                       n_results$true_theta <= n_results$theta_ci_upper, na.rm = TRUE) * 100, 2),
          "%, Mean runtime =", round(mean(n_results$runtime, na.rm = TRUE), 4), "s\n")
    }
  }
  
  # By dimension
  cat("\nPerformance by dimension (d):\n")
  for (d_val in unique(sim_settings$d)) {
    d_results <- all_detailed_results[all_detailed_results$d == d_val, ]
    if (nrow(d_results) > 0) {
      cat("  d =", d_val, ": Mean theta RMSE =", 
          round(sqrt(mean(d_results$theta_error^2, na.rm = TRUE)), 6),
          ", Coverage =", round(mean(d_results$theta_ci_lower <= d_results$true_theta & 
                                       d_results$true_theta <= d_results$theta_ci_upper, na.rm = TRUE) * 100, 2),
          "%, Mean runtime =", round(mean(d_results$runtime, na.rm = TRUE), 4), "s\n")
    }
  }
  
  # By theta
  cat("\nPerformance by true theta:\n")
  for (theta_val in unique(sim_settings$theta)) {
    theta_results <- all_detailed_results[all_detailed_results$true_theta == theta_val, ]
    if (nrow(theta_results) > 0) {
      cat("  theta =", theta_val, ": Mean theta RMSE =", 
          round(sqrt(mean(theta_results$theta_error^2, na.rm = TRUE)), 6),
          ", Coverage =", round(mean(theta_results$theta_ci_lower <= theta_results$true_theta & 
                                       theta_results$true_theta <= theta_results$theta_ci_upper, na.rm = TRUE) * 100, 2),
          "%, Mean runtime =", round(mean(theta_results$runtime, na.rm = TRUE), 4), "s\n")
    }
  }
  
  sink()
  cat("\nDetailed report saved to:", report_file, "\n")
  
} else {
  cat("\nWARNING: No successful MCMC simulations across all settings.\n")
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("BAYESIAN MCMC SIMULATION STUDY COMPLETED\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("Results saved in directory:", results_dir, "\n")
cat("Files created:\n")
cat("  - Individual setting results (bayesian_s1_detailed.csv, etc.)\n")
cat("  - all_detailed_results.csv (combined all settings)\n")
cat("  - summary_results.csv (statistical summary)\n")
cat("  - formatted_summary.csv (readable summary)\n")
cat("  - bayesian_mcmc_report.txt (comprehensive report)\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
