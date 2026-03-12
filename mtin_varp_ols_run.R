# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 2
p <- 2
n <- 200
n_sim <- 100
theta_true <- 0.7

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=2_n=200_theta=07.csv", row.names = TRUE)


# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 2
p <- 2
n <- 200
n_sim <- 100
theta_true <- 0.8

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=2_n=200_theta=08.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 2
p <- 2
n <- 200
n_sim <- 100
theta_true <- 0.9

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=2_n=200_theta=09.csv", row.names = TRUE)


# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 2
p <- 2
n <- 500
n_sim <- 100
theta_true <- 0.7

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=2_n=500_theta=07.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 2
p <- 2
n <- 500
n_sim <- 100
theta_true <- 0.8

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=2_n=500_theta=08.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 2
p <- 2
n <- 500
n_sim <- 100
theta_true <- 0.9

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=2_n=500_theta=09.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 2
p <- 2
n <- 1000
n_sim <- 100
theta_true <- 0.7

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=2_n=1000_theta=07.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 2
p <- 2
n <- 1000
n_sim <- 100
theta_true <- 0.8

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=2_n=1000_theta=08.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 2
p <- 2
n <- 1000
n_sim <- 100
theta_true <- 0.9

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=2_n=1000_theta=09.csv", row.names = TRUE)




# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 3
p <- 2
n <- 200
n_sim <- 100
theta_true <- 0.7

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=3_n=200_theta=07.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 3
p <- 2
n <- 200
n_sim <- 100
theta_true <- 0.8

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=3_n=200_theta=08.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 3
p <- 2
n <- 200
n_sim <- 100
theta_true <- 0.9

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=3_n=200_theta=09.csv", row.names = TRUE)


# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 3
p <- 2
n <- 500
n_sim <- 100
theta_true <- 0.7

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=3_n=500_theta=07.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 3
p <- 2
n <- 500
n_sim <- 100
theta_true <- 0.8

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=3_n=500_theta=08.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 3
p <- 2
n <- 500
n_sim <- 100
theta_true <- 0.9

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=3_n=500_theta=09.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 3
p <- 2
n <- 1000
n_sim <- 100
theta_true <- 0.7

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=3_n=1000_theta=07.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 3
p <- 2
n <- 1000
n_sim <- 100
theta_true <- 0.8

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=3_n=1000_theta=08.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 3
p <- 2
n <- 1000
n_sim <- 100
theta_true <- 0.9

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=3_n=1000_theta=09.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 5
p <- 2
n <- 200
n_sim <- 100
theta_true <- 0.7

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=5_n=200_theta=07.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 5
p <- 2
n <- 200
n_sim <- 100
theta_true <- 0.8

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=5_n=200_theta=08.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 5
p <- 2
n <- 200
n_sim <- 100
theta_true <- 0.9

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=5_n=200_theta=09.csv", row.names = TRUE)


# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 5
p <- 2
n <- 500
n_sim <- 100
theta_true <- 0.7

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=5_n=500_theta=07.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 5
p <- 2
n <- 500
n_sim <- 100
theta_true <- 0.8

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=5_n=500_theta=08.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 5
p <- 2
n <- 500
n_sim <- 100
theta_true <- 0.9

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=5_n=500_theta=09.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 5
p <- 2
n <- 1000
n_sim <- 100
theta_true <- 0.7

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=5_n=1000_theta=07.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 5
p <- 2
n <- 1000
n_sim <- 100
theta_true <- 0.8

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=5_n=1000_theta=08.csv", row.names = TRUE)



# === Load Dependencies ===
source("rtmin.R")
source("mtin_varp_est_ols.R")
source("sim_mtin_varp.R")

library(doParallel)
library(foreach)

# === Parameters ===
d <- 5
p <- 2
n <- 1000
n_sim <- 100
theta_true <- 0.9

# === Simulate Data ===
simulated_data <- simulate_mtin_varp_data(d = d, p = p, n = n, n_sim = n_sim, theta_true = theta_true)

# === Setup Parallel Backend ===
n_cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel estimation\n")

# === Parallel Estimation ===
estimation_results <- foreach(i = 1:n_sim,
                              .packages = c("vars"),
                              .export = c("mtin_varp_using_VAR", "simulated_data", "p")) %dopar% {
  Y_i <- simulated_data$data[,,i]
  result <- tryCatch({
    mtin_varp_using_VAR(Y_i, p = p)
  }, error = function(e) NULL)
  result
}

# === Stop Cluster ===
stopCluster(cl)

# === Filter Valid Results ===
valid_results <- Filter(Negate(is.null), estimation_results)
cat("Successful estimations:", length(valid_results), "/", n_sim, "\n")

# === Example: Access results ===
if (length(valid_results) > 0) {
  print(valid_results)  # Print the first successful estimation
}
Output_ols <- do.call(rbind,valid_results)
Output_ols

write.csv(Output_ols, "C:\\Users\\November First\\Documents\\ols_mtin_var_p=2_d=5_n=1000_theta=09.csv", row.names = TRUE)


