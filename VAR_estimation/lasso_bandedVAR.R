# In papers DGP notation: (F0),(V1),(S1)
# standard and truncated results for:
# Figures 3 + 4.
# Table 9

# Generates results for VAR coefficient estimation error of lasso with truncation and without truncation. 
# data generated from banded VAR model with innovation errors with diagonal covariance structure. 

# Here the autocovariance is included in the CV measure, for the choice of tau.
# And no standardisation is carried out.


# In order to compute errors for errors from different distributions requires uncommenting 
# And to generate results for against dimension or sample size, and (n,p) settings used in the tables.

setwd("/user/work/zl22291/simulations/truncation")
source(file = "functions/data_generation.R")
source(file = "functions/estimation.R")

#### sourcing functions ####
# source(file = "truncation/functions/data_generation.R")
# source(file = "truncation/functions/estimation.R")
#### creating list to store objects ####
list_setup = list(
  "truncated" = list(
    "t_2.1" = list(),
    "t_3" = list(),
    "t_4" = list(),
    "log_normal" = list(),
    "gaussian" = list()
  ),
  "not_truncated" = list(
    "t_2.1" = list(),
    "t_3" = list(),
    "t_4" = list(),
    "log_normal" = list(),
    "gaussian" = list()
  )
)
list_setup_2 <- list(
  "M" = list_setup,
  "F" = list_setup,
  "l_2_inf" = list_setup # maximal euclidean norm of a column
)
lasso_table = list_setup_2


############################# uncomment whether need to run against sample size n, or against dimension, or (n,p) settings used for tables #####################
#### lasso_vs_n ####
n_p = cbind(n = seq(100, 300, 50), p = 50)
#### lasso_vs_p ####
# n_p = cbind(n = 200, p = seq(50, 300, 50))
#### table values ####
# n_p = cbind(n = c(100,100,200,200,500,500), p = c(50,100,50,100,100,200))

nsim = 200
#data = VAR_1_data_ind_t(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_df = 2.1)

A_coeff = A_coeff_banded(n_p)

library(Rcpp)
maxEuclideanNormR <- cppFunction('double maxEuclideanNormR(NumericMatrix m) {
  int n = m.nrow();
  double max_norm = 0.0;
  
  for (int i = 0; i < n; i++) {
    NumericVector row = m(i, _);
    double sum_of_squares = 0.0;
    for (int j = 0; j < row.size(); j++) {
      sum_of_squares += row[j] * row[j];
    }
    double norm = sqrt(sum_of_squares);
    if (i == 0 || norm > max_norm) {
      max_norm = norm;
    }
  }
  
  return max_norm;
}')

# For each dimension type of data need to loop through to get vector of errors for that norm

# Given a data matrix this function truncates and then fits LASSO and returns matrix with difference to true coefficients.
library(sparsevar)
sparsevar_trunc_fit_diff = function(data, n_tau = 60, lag = 0, trim = 1, trim_d=1, trim_all=F, A_coeff, max = T, m_norm, cv_lag=T){
  data = cross_val_and_trunc(data, n_tau, lag, trim, trim_d, trim_all, max = max, cv_lag=cv_lag, standardise = F)
  fit_diff = fitVAR(data, p = 1, parallel = T, ncores = 4)$A[[1]] - A_coeff
  return(fit_diff)
}

error_func = function(m_norm, diff_mat){
  if(m_norm == "l_2_inf"){
    error = maxEuclideanNormR(diff_mat)
  } else {
    error = norm(diff_mat, paste0(m_norm))
  }
  return(error)
}

compute_errors = function(data, table, distr){
  # truncated
  for(i in names(data)){
    fit_diff = lapply(data[[i]], sparsevar_trunc_fit_diff, A_coeff = A_coeff[[i]])
    
    for(j in c("M", "F", "l_2_inf")){
      trunc_lasso_error = sapply(fit_diff, error_func, m_norm = j)
      table[[j]]$truncated[[distr]][[i]] = trunc_lasso_error
    }
  }
  
  # not-truncated
  for(i in names(data)){
    fit_diff = lapply(data[[i]], function(x){fitVAR(x, p = 1,parallel = T, ncores = 4)$A[[1]] - A_coeff[[i]]})
    
    for(j in c("M", "F", "l_2_inf")){
      trunc_lasso_error = sapply(fit_diff, error_func, m_norm = j)
      table[[j]]$not_truncated[[distr]][[i]] = trunc_lasso_error
    }
  }
  return(table)
}

# gaussian
set.seed(200)
lasso_table = compute_errors(
  data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "gauss"),
  table = lasso_table, distr = "gaussian"
)
#### save ####
saveRDS(lasso_table, file = "lasso_cv_max_norm_auto_V2/gaussian_vs_n.rds")
rm(lasso_table)
gc()









