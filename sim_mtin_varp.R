
simulate_mtin_varp_data <- function(d = 3, p = 2, n = 500, n_sim = 100, theta_true = 0.7, seed = 123) {
  # --- Required MTIN Sampling Function ---
  if (!exists("rmtin")) stop("Function 'rmtin' is required but not found. Source 'rtmin.R' first.")

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


# First, source the rmtin generator
source("rtmin.R")
# Then run the simulation
#simulated_data <- simulate_mtin_varp_data(d = 3, p = 2, n = 500, n_sim = 100, theta_true = 0.7)

