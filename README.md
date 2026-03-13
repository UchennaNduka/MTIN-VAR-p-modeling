This is a collection of R scripts that provides a framework for parameter estimation in vector autoregressive models driven by Multivariate Tail Inflated Normal 
(MTIN) innovations. These scripts are used for experiments involving the MTIN-VAR(p) model. The ECME algorithm and Bayesian MCMC procedure are proposed for parameter
estimation in the model. These two methods were compared with standard parameter estimation procedures for VAR(p) that can be implemented in R. The standard procedures are captured in R
packages called VAR and BVAR.

**Overview**

The scripts performs the following functions:

1. **Data simulation:** Generates datasets based on a specified MTIN-VAR(p) structure.

2. **Parameter estimation:** Estimates the parameters of the MTIN-VAR(p) model based on each of the simulated samples using ECME algorithm, Bayesian MCMC, VAR, and BVAR.

3. **Boxplot visualization:** Plots the box plot of the Frobenius norm of estimation errors for each method across experimental settings.

4. **Data export:** Filters successful estimations and saves the aggregated results to a CSV file.

**Usage**
   
 1. Set your working directory to the folder containing the scripts.
 2. Update the file path in the write.csv function at the end of the script to your desired local directory.
 3. Run the relevant scripts.
    
(i) For ECME algorithm run **MTIN_VAR_ECME_sim.R**,

(ii) For Bayesian MCMC run **MTIN_VAR_MCMC_sim.R**,

(iii) For OLS using VAR package run **mtin_varp_ols_run.R**,

(iv) For Bayesian VAR using BVAR package run **mtin_varp_bvar_run.R**, 

(v) **Psi_boxplots_rcodes.R** gives the boxplots for $\mathbf{\Psi}-\mathbf{\hat{\Psi}}$, 

(vi) **Omega_boxplots_rcodes.R** gives the boxplots for $\mathbf{\Omega}-\mathbf{\hat{\Omega}}$, 

(vii) **Theta_boxplots_rcodes.R** gives the boxplots for $\theta-\hat{\theta}$. 

