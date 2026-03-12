mtin_varp_using_VAR <- function(Data, p) {
  # Load required package
  if (!requireNamespace("vars", quietly = TRUE)) {
    stop("Package 'vars' is required but not installed.")
  }
  library(vars)

  # Prepare data
  Data_df <- as.data.frame(Data)
  colnames(Data_df) <- paste0("Y", 1:ncol(Data_df))

  # Fit VAR model with specified lag order p
  var_model <- vars::VAR(Data_df, p = p, type = "const")

  # Extract coefficients
  Psi_hat_raw <- Bcoef(var_model)  # d x (d*p + 1), const is last column
  d <- ncol(Data)
  
  # Reorder: move const column to the first position
  Psi_hat <- cbind(Psi_hat_raw[, ncol(Psi_hat_raw)],  # const
                   Psi_hat_raw[, -ncol(Psi_hat_raw)])  # rest (lag matrices)

  # Extract estimated residual covariance matrix
  Omega_hat <- summary(var_model)$covres
  
  set.seed(123)
  intercept <- matrix(runif(d, -0.5, 0.5), nrow = d)
  AR_coefs <- array(runif(d * d * p, -0.5, 0.5), dim = c(d, d * p))
  true_Psi <- cbind(intercept, AR_coefs)

  true_Omega <- diag(d)
  true_Omega[lower.tri(true_Omega)] <- 0.3
  true_Omega <- true_Omega + t(true_Omega) - diag(diag(true_Omega))


  # Compute errors
  frob_Psi_errors   <- norm(Psi_hat - true_Psi, type = "F")
  frob_Omega_errors <- norm(Omega_hat - true_Omega, type = "F")
  #theta_errors      <- NA  # No tail inflation in VAR

  results <- cbind(Psi_norm = frob_Psi_errors,
                   Omega_norm = frob_Omega_errors)
                   
  return(results)
}

#Try <- mtin_varp_using_VAR(simulated_data$data[,,1], p=2)
#Try

