This is a collection of R scripts that provides a framework for parameter estimation in vector autoregressive models driven by Multivariate Tail Inflated Normal 
(MTIN) innovations. These scripts are used for experiments involving the MTIN-VAR(p) model. The ECME algorithm and Bayesian MCMC procedure are proposed for parameter
estimation in the model. These two methods were compared with standard parameter estimation procedures for VAR(p) that can be implemented in R. The standard procedures are captured in R
packages called VAR and BVAR.

**Overview**

The scripts performs the following functions:

1. **Data simulation:** Generates datasets based on a specified MTIN-VAR(p) structure.
