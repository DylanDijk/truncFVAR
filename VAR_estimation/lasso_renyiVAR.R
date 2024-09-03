# In papers DGP notation: (F0),(V2),(S1)
# standard and truncated results for:
# Table 10

# Generates results for VAR coefficient estimation error of lasso with truncation and without truncation. 
# data generated from renyi VAR model with innovation errors with diagonal covariance structure. 

# Here the autocovariance is included in the CV measure, for the choice of tau.
# And no standardisation is carried out.


# In order to compute errors for errors from different distributions requires uncommenting 
# And to generate results for against dimension or sample size, and (n,p) settings used in the tables.

setwd("/user/work/zl22291/simulations/truncation")
source(file = "functions/data_generation.R")
source(file = "functions/estimation.R")

##### sourcing functions ####
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
# list to store objects
list_setup_2 <- list(
  "M" = list_setup,
  "F" = list_setup,
  "l_2_inf" = list_setup # maximal euclidean norm of a column
)

lasso_ke_table = list_setup_2


############################# uncomment whether need to run against sample size n, or against dimension, or (n,p) settings used for tables #####################
#### table values ####
n_p = cbind(n = c(100,100,200,200,500,500), p = c(50,100,50,100,100,200))
nsim = 200

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
sparsevar_trunc_fit_diff = function(data, n_tau = 60, lag = 0, trim = 1, trim_d=1, trim_all=F, A_coeff, max = T, cv_lag=T){
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



# i = names(data)[1]
# j = "M"
# distr = "t_2.1"

compute_errors = function(data, table, distr){
  print(paste(distr))
  # truncated
  for(i in names(data)){
    fit_diff = Map(sparsevar_trunc_fit_diff, data = data[[i]],  A_coeff = attributes(data[[i]])$A)
    
    for(j in c("M", "F", "l_2_inf")){
      trunc_lasso_error = sapply(fit_diff, error_func, m_norm = j)
      table[[j]]$truncated[[distr]][[i]] = trunc_lasso_error
    }
  }
  
  # not-truncated
  for(i in names(data)){
    fit_diff = Map( function(x,y){fitVAR(x, p = 1,parallel = T, ncores = 4)$A[[1]] - y}, data[[i]],  attributes(data[[i]])$A)
    
    for(j in c("M", "F", "l_2_inf")){
      trunc_lasso_error = sapply(fit_diff, error_func, m_norm = j)
      table[[j]]$not_truncated[[distr]][[i]] = trunc_lasso_error
    }
  }
  return(table)
}

# t_2.1

# set.seed(200)
# lasso_ke_table = compute_errors(
#   data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = renyi_coeff, innov_df = 4, innov_dist = "t"),
#   table = lasso_ke_table, distr = "t_4"
# )

#log_normal
set.seed(200)
lasso_ke_table = compute_errors(
data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = renyi_coeff, innov_dist = "log"),
  table = lasso_ke_table, distr = "log_normal"
)


#gaussian
# set.seed(200)
# lasso_ke_table = compute_errors(
#   data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = renyi_coeff, innov_dist = "gauss"),
#   table = lasso_ke_table, distr = "gaussian"
# )

#### save ####
saveRDS(lasso_ke_table, file = "lasso_cv_max_norm_auto_V2/gauss_ke_table_renyi.rds")
# save(lasso_vs_n, file = "heavytailed_sparsevar_lasso_cv_with_max_norm/lasso_vs_n_cv_with_max_norm.RData")
rm(lasso_ke_table)
gc()









