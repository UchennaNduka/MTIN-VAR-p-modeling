# ============================================
# MTIN-VAR(p) Simulation and Estimation Script
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

# MTIN-VAR(p) Estimation Function with ECME Algorithm using matrix expressions
mtin.varp <- function(Data, p, h, max.iter, theta_init = 0.5, 
                      tol = 1e-8, param_tol = 1e-6, 
                      true_Psi = NULL, true_Omega = NULL, true_theta = NULL) {
  library(mvtnorm)  
  library(pracma)
  library(rootSolve)
  
  # Start timing
  start_time <- proc.time()[1]
  
  N <- nrow(Data)
  K <- ncol(Data)
  n <- N
  d <- K
  
  # Extract data
  Y <- Data
  T <- n - p
  
  # Initialize
  B.t <- matrix(rnorm(K*(K*p+1)), nrow = K, ncol = K*p+1)  # Ψ
  C.t <- diag(K)  # Ω
  theta <- theta_init
  
  iter <- 0
  converged <- FALSE
  lltrace <- numeric(max.iter)
  
  while (iter < max.iter && !converged) {
    iter <- iter + 1
    
    # ============================================
    # E-step: Calculate weights w_t
    # ============================================
    w <- numeric(T)  # w_t^{(k+1)}
    
    for (t in (p+1):n) {
      # Create x_{t-1} = [1, y_{t-1}, ..., y_{t-p}]
      x_prev <- matrix(1, ncol = 1)
      for (i in 1:p) {
        x_prev <- rbind(x_prev, matrix(Y[t-i, ], ncol = 1))
      }
      
      # Compute residual
      y_t <- matrix(Y[t, ], ncol = 1)
      residual <- y_t - B.t %*% x_prev
      
      # Quadratic form
      quad_form <- as.numeric(t(residual) %*% solve(C.t) %*% residual)
      
      # Compute weight using equation from MTIN distribution
      s1 <- d/2 + 2
      s2 <- d/2 + 1
      x1 <- (1 - theta) * quad_form / 2
      x2 <- quad_form / 2
      
      gamma_num <- gamma(s1) * (pgamma(x1, s1, lower=FALSE) - pgamma(x2, s1, lower=FALSE))
      gamma_den <- gamma(s2) * (pgamma(x1, s2, lower=FALSE) - pgamma(x2, s2, lower=FALSE))
      
      w[t-p] <- 2/quad_form * gamma_num / gamma_den
    }
    
    # ============================================
    # CM-step 1: Update Ψ and Ω using matrix expressions
    # ============================================
    
    # Create lagged data matrices
    # Create X matrix: [1, y_{t-1}, ..., y_{t-p}] for t = p+1,...,n
    X <- matrix(1, nrow = T, ncol = 1)  # Intercept column
    for (i in 1:p) {
      lag_data <- Y[(p-i+1):(n-i), ]  # y_{t-i}
      X <- cbind(X, lag_data)
    }
    
    # Current response matrix
    Y_curr <- Y[(p+1):n, ]
    
    # Calculate weighted sums using your notation
    # Initialize matrices
    R_T <- sum(w)  # R^T = sum(w_t)
    
    # R_i^{wy} = sum(w_t * y_{t-i}) for i = 0,...,p
    R_wy <- matrix(0, nrow = K, ncol = p+1)
    for (i in 0:p) {
      if (i == 0) {
        # i=0: y_t
        y_lagged <- Y_curr
      } else {
        # i>0: y_{t-i}
        y_lagged <- Y[(p-i+1):(n-i), ]
      }
      for (t in 1:T) {
        R_wy[, i+1] <- R_wy[, i+1] + w[t] * y_lagged[t, ]
      }
    }
    
    # R_{i,j}^{wyy} = sum(w_t * y_{t-i} * y_{t-j}^T) for i,j = 0,...,p
    R_wyy <- array(0, dim = c(K, K, p+1, p+1))
    for (i in 0:p) {
      for (j in 0:p) {
        # Get y_{t-i} and y_{t-j}
        if (i == 0) {
          y_i <- Y_curr
        } else {
          y_i <- Y[(p-i+1):(n-i), ]
        }
        
        if (j == 0) {
          y_j <- Y_curr
        } else {
          y_j <- Y[(p-j+1):(n-j), ]
        }
        
        # Compute weighted sum of outer products
        for (t in 1:T) {
          R_wyy[,, i+1, j+1] <- R_wyy[,, i+1, j+1] + w[t] * 
            (y_i[t, ] %*% t(y_j[t, ]))
        }
      }
    }
    
    # ============================================
    # Construct H_0 and H_1 matrices
    # ============================================
    
    # H_0 = [R_0^{wy}  R_0,1^{wyy}  ...  R_0,p^{wyy}]
    H_0 <- matrix(0, nrow = K, ncol = 1 + K*p)
    
    # First column: R_0^{wy}
    H_0[, 1] <- R_wy[, 1]  # i=0
    
    # Remaining columns: R_0,j^{wyy} for j=1,...,p
    for (j in 1:p) {
      H_0[, (1 + (j-1)*K + 1):(1 + j*K)] <- R_wyy[,, 1, j+1]  # i=0, j=j
    }
    
    # H_1 matrix
    H_1 <- matrix(0, nrow = 1 + K*p, ncol = 1 + K*p)
    
    # Top-left element: R^T
    H_1[1, 1] <- R_T
    
    # First row (except top-left): (R_j^{wy})^T for j=1,...,p
    for (j in 1:p) {
      H_1[1, (1 + (j-1)*K + 1):(1 + j*K)] <- t(R_wy[, j+1])
    }
    
    # First column (except top-left): R_j^{wy} for j=1,...,p
    for (j in 1:p) {
      H_1[(1 + (j-1)*K + 1):(1 + j*K), 1] <- R_wy[, j+1]
    }
    
    # Remaining blocks: R_{i,j}^{wyy} for i,j=1,...,p
    for (i in 1:p) {
      for (j in 1:p) {
        row_start <- 1 + (i-1)*K + 1
        row_end <- 1 + i*K
        col_start <- 1 + (j-1)*K + 1
        col_end <- 1 + j*K
        
        H_1[row_start:row_end, col_start:col_end] <- R_wyy[,, i+1, j+1]
      }
    }
    
    # ============================================
    # Update Ψ using Eq. 20: Ψ^{(k+1)} = H_0^{(k+1)} (H_1^{(k+1)})^{-1}
    # ============================================
    B.t1 <- H_0 %*% solve(H_1)
    
    # ============================================
    # Update Ω using Eq. 21
    # ============================================
    # Compute R_{0,0}^{wyy}
    R_00_wyy <- R_wyy[,, 1, 1]
    
    # Compute Ω^{(k+1)} using matrix expression
    term1 <- R_00_wyy
    term2 <- H_0 %*% t(B.t1)
    term3 <- B.t1 %*% t(H_0)
    term4 <- B.t1 %*% H_1 %*% t(B.t1)
    
    C.t1 <- (1/T) * (term1 - term2 - term3 + term4)
    
    # Ensure symmetry
    C.t1 <- 0.5 * (C.t1 + t(C.t1))
    
    # ============================================
    # CM-step 2: Update theta (same as before)
    # ============================================
    theta_eq <- function(theta) {
      sum_term <- 0
      for (t in 1:T) {
        # Create x_{t-1}
        x_prev <- matrix(1, ncol = 1)
        for (i in 1:p) {
          x_prev <- rbind(x_prev, matrix(Y[(p+t-i), ], ncol = 1))
        }
        
        # Compute residual with new B.t1
        y_t <- matrix(Y[p+t, ], ncol = 1)
        residual <- y_t - B.t1 %*% x_prev
        
        # Quadratic form with new C.t1
        quad_form <- as.numeric(t(residual) %*% solve(C.t1) %*% residual)
        x1 <- (1 - theta) * quad_form / 2
        x2 <- quad_form / 2
        
        gamma_num <- gamma(d/2 + 1) * (pgamma(x1, d/2 + 1, lower=FALSE) - 
                                         pgamma(x2, d/2 + 1, lower=FALSE))
        term <- (quad_form/2)^(d/2 + 1) / gamma_num
        sum_term <- sum_term + term
      }
      -T/theta + (1 - theta)^(d/2) * sum_term
    }
    
    # Find theta that makes theta_eq = 0
    theta_new <- uniroot(theta_eq, interval = c(1e-6, 1 - 1e-6))$root
    
    # ============================================
    # Check convergence
    # ============================================
    B_diff <- max(abs(B.t1 - B.t))
    C_diff <- max(abs(C.t1 - C.t))
    theta_diff <- abs(theta_new - theta)
    
    # Update parameters
    B.t <- B.t1
    C.t <- C.t1
    theta <- theta_new
    
    # Store log-likelihood
    lltrace[iter] <- log_likelihood(Y, p, B.t, C.t, theta)
    
    # Check convergence criteria
    if (iter > 1) {
      loglik_diff <- abs(lltrace[iter] - lltrace[iter-1])
      if (loglik_diff < tol && B_diff < param_tol && 
          C_diff < param_tol && theta_diff < param_tol) {
        converged <- TRUE
      }
    }
  }
  
  # End timing
  end_time <- proc.time()[1]
  runtime <- end_time - start_time
  
  # Calculate errors if true parameters are provided
  if (!is.null(true_Psi) && !is.null(true_Omega) && !is.null(true_theta)) {
    Error_Psi <- B.t - true_Psi
    Error_Omega <- C.t - true_Omega
    
    frob_Psi_errors   <- norm(Error_Psi, type = "F")
    frob_Omega_errors <- norm(Error_Omega, type = "F")
    theta_errors      <- theta - true_theta
    
    results <- data.frame(
      Psi_norm = frob_Psi_errors,
      Omega_norm = frob_Omega_errors,
      theta_error = theta_errors,
      theta_hat = theta,
      runtime = runtime,
      iterations = iter,
      converged = converged,
      final_loglik = lltrace[iter]
    )
  } else {
    results <- list(
      Psi = B.t,
      Omega = C.t,
      theta = theta,
      loglik = lltrace[iter],
      iterations = iter,
      runtime = runtime,
      converged = converged
    )
  }
  
  return(results)
}


# ============================================
# Main Simulation Study
# ============================================

library(doParallel)
library(foreach)

# === Parameters ===
d <- 2
p <- 2
n <- 200
n_sim <- 100
theta_true <- 0.7

cat("============================================\n")
cat("MTIN-VAR(p) Simulation Study\n")
cat("============================================\n")
cat("Parameters:\n")
cat("- Dimension (d):", d, "\n")
cat("- VAR order (p):", p, "\n")
cat("- Time series length (n):", n, "\n")
cat("- Number of simulations (n_sim):", n_sim, "\n")
cat("- True theta:", theta_true, "\n")
cat("============================================\n")

# === Simulate Data ===
cat("\nGenerating simulated data...\n")
set.seed(123)
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)
cat("Data generation completed.\n")

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("\nUsing", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
cat("\nStarting parallel estimation...\n")
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("mvtnorm", "pracma", "rootSolve"),
                              .export = c("mtin.varp", "rmtin", "log_likelihood",
                                          "simulated_data", "p", "d", "theta_true"),
                              .combine = 'rbind') %dopar% {
                                Y_i <- simulated_data$data[,,i]
                                result <- tryCatch({
                                  res <- mtin.varp(Data = Y_i, p = p, h = 28, max.iter = 500, 
                                                   theta_init = 0.5, 
                                                   true_Psi = simulated_data$true_Psi,
                                                   true_Omega = simulated_data$true_Omega,
                                                   true_theta = theta_true)
                                  # Return as a one-row data frame
                                  data.frame(
                                    simulation = i,
                                    Psi_norm = res$Psi_norm,
                                    Omega_norm = res$Omega_norm,
                                    theta_error = res$theta_error,
                                    theta_hat = res$theta_hat,
                                    runtime = res$runtime,
                                    iterations = res$iterations,
                                    converged = res$converged,
                                    final_loglik = res$final_loglik
                                  )
                                }, error = function(e) {
                                  cat("Error in simulation", i, ":", e$message, "\n")
                                  NULL
                                })
                                result
                              }

# === Stop Cluster ===
stopCluster(cl)
cat("Parallel estimation completed.\n")

# === Filter Valid Results ===
valid_results <- estimation_results[complete.cases(estimation_results), ]
cat("\nSuccessful estimations:", nrow(valid_results), "/", n_sim, "\n")

# === Summary Statistics ===
if (nrow(valid_results) > 0) {
  Output_ecme <- valid_results
  
  cat("\n============================================\n")
  cat("ESTIMATION RESULTS SUMMARY\n")
  cat("============================================\n")
  
  # Parameter estimation errors
  cat("\n--- Parameter Estimation Errors ---\n")
  cat("Mean Psi Frobenius norm:", round(mean(Output_ecme$Psi_norm), 6), "\n")
  cat("SD Psi Frobenius norm:", round(sd(Output_ecme$Psi_norm), 6), "\n")
  cat("Mean Omega Frobenius norm:", round(mean(Output_ecme$Omega_norm), 6), "\n")
  cat("SD Omega Frobenius norm:", round(sd(Output_ecme$Omega_norm), 6), "\n")
  
  # Theta estimation
  cat("\n--- Theta Estimation ---\n")
  cat("True theta:", theta_true, "\n")
  cat("Mean estimated theta:", round(mean(Output_ecme$theta_hat), 6), "\n")
  cat("SD estimated theta:", round(sd(Output_ecme$theta_hat), 6), "\n")
  cat("Mean theta error:", round(mean(Output_ecme$theta_error), 6), "\n")
  cat("Bias (theta):", round(mean(Output_ecme$theta_hat) - theta_true, 6), "\n")
  cat("RMSE (theta):", round(sqrt(mean(Output_ecme$theta_error^2)), 6), "\n")
  
  # Runtime and iterations
  cat("\n--- Algorithm Performance ---\n")
  cat("Mean runtime (seconds):", round(mean(Output_ecme$runtime), 4), "\n")
  cat("SD runtime (seconds):", round(sd(Output_ecme$runtime), 4), "\n")
  cat("Min runtime (seconds):", round(min(Output_ecme$runtime), 4), "\n")
  cat("Max runtime (seconds):", round(max(Output_ecme$runtime), 4), "\n")
  
  cat("\nMean number of iterations:", round(mean(Output_ecme$iterations), 2), "\n")
  cat("SD iterations:", round(sd(Output_ecme$iterations), 2), "\n")
  cat("Min iterations:", min(Output_ecme$iterations), "\n")
  cat("Max iterations:", max(Output_ecme$iterations), "\n")
  
  # Convergence
  cat("\n--- Convergence ---\n")
  convergence_rate <- sum(Output_ecme$converged) / nrow(Output_ecme) * 100
  cat("Convergence rate:", round(convergence_rate, 2), "%\n")
  
  # Log-likelihood
  cat("\n--- Log-Likelihood ---\n")
  cat("Mean final log-likelihood:", round(mean(Output_ecme$final_loglik), 4), "\n")
  cat("SD final log-likelihood:", round(sd(Output_ecme$final_loglik), 4), "\n")
  
  # Display first few results
  cat("\n============================================\n")
  cat("FIRST 5 SIMULATION RESULTS\n")
  cat("============================================\n")
  print(head(Output_ecme, 5))
  
  # Save results to file
  write.csv(Output_ecme, "mtin_varp_simulation_results1.csv", row.names = FALSE)
  cat("\nResults saved to 'mtin_varp_simulation_results.csv'\n")
  
    
  # Create a comprehensive summary table
  summary_table <- data.frame(
    Metric = c(
      "Mean Psi Norm", "SD Psi Norm",
      "Mean Omega Norm", "SD Omega Norm",
      "Mean Theta Estimate", "SD Theta Estimate",
      "Theta Bias", "Theta RMSE",
      "Mean Runtime (s)", "SD Runtime (s)",
      "Mean Iterations", "SD Iterations",
      "Convergence Rate (%)"
    ),
    Value = c(
      round(mean(Output_ecme$Psi_norm), 6),
      round(sd(Output_ecme$Psi_norm), 6),
      round(mean(Output_ecme$Omega_norm), 6),
      round(sd(Output_ecme$Omega_norm), 6),
      round(mean(Output_ecme$theta_hat), 6),
      round(sd(Output_ecme$theta_hat), 6),
      round(mean(Output_ecme$theta_hat) - theta_true, 6),
      round(sqrt(mean(Output_ecme$theta_error^2)), 6),
      round(mean(Output_ecme$runtime), 4),
      round(sd(Output_ecme$runtime), 4),
      round(mean(Output_ecme$iterations), 2),
      round(sd(Output_ecme$iterations), 2),
      round(convergence_rate, 2)
    )
  )
  
  cat("\n============================================\n")
  cat("COMPREHENSIVE SUMMARY TABLE\n")
  cat("============================================\n")
  print(summary_table)
  
  # Save summary table
  write.csv(summary_table, "mtin_varp_summary_table1.csv", row.names = FALSE)
  cat("\nSummary table saved to 'mtin_varp_summary_table.csv'\n")
  
} else {
  cat("\nNo successful estimations to analyze.\n")
}

cat("\n============================================\n")
cat("SIMULATION STUDY COMPLETED\n")
cat("============================================\n")
