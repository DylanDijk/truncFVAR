##########################################################################################################################################
# This script gets the VAR coefficient estimates with and without truncation, when the VAR order used for fitting is 2, but the true VAR order is 1.
# Results are uses to generate results in Table - looking at performance under VAR order misspecification
##########################################################################################################################################

# Need to first get functions to generate data
# This is now sourcing function that can generate stable banded VAR process of lags greater than 1
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/3c4339ea8bf00a735408be0659d0c88bcf5dac0e/functions/data_generation.R")

# source(file = "functions/estimation.R")
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/84f22825f88656a6b6efdb8cf628d1e7610397b0/functions/estimation.R")


########################################################################################################################################################
# Below I am modifying the script https://github.com/DylanDijk/truncFVAR/blob/84f22825f88656a6b6efdb8cf628d1e7610397b0/VAR_estimation/lasso_bandedVAR.R
########################################################################################################################################################

# what lags to generate data for
d_0 = 1

# what lags to fit 
d = 2

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

#### table values ####
n_p = cbind(n = c(100,100,200,200,500,500), p = c(50,100,50,100,100,200))

nsim = 200
#data = VAR_1_data_ind_t(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_df = 2.1)

A_coeff = A_coeff_banded_d(n_p, d = d_0)
# now adding zero matrix to side of var coefficient matrices.
if(length(A_coeff[[1]]) < d){
  A_coeff_ex <- lapply(A_coeff, function(A_list) {
    A_combined <- do.call(cbind, A_list)
    p <- nrow(A_combined)
    cbind(A_combined, matrix(0, nrow = p, ncol = p))
  })
}



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
sparsevar_trunc_fit_diff = function(data, n_tau = 60, lag = d, trim = 1, trim_d=1, trim_all=F, A_coeff, max = T, m_norm, cv_lag=T){
  data = cross_val_and_trunc(data, n_tau, lag, trim, trim_d, trim_all, max = max, cv_lag=cv_lag, standardise = F)
  fit_diff = fitVAR(data, p = d, parallel = T, ncores = 4)
  
  fit_diff = do.call(cbind, fit_diff$A[1:2]) - A_coeff
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

norms <- c("M", "F", "l_2_inf")

compute_errors = function(data, table, distr){
  # truncated
  for(i in names(data)){
    print(i)
    fit_diff = parallel::mclapply(
      data[[i]],
      sparsevar_trunc_fit_diff,
      A_coeff = A_coeff_ex[[i]], mc.cores = 4
    )
    
    for(j in norms){
      trunc_lasso_error = vapply(fit_diff, error_func, numeric(1), m_norm = j)
      table[[j]]$truncated[[distr]][[i]] = trunc_lasso_error
    }
  }
  
  # not-truncated
  for(i in names(data)){
    fit_diff = parallel::mclapply(
      data[[i]],
      function(x){
        fit = fitVAR(x, p = d, parallel = T, ncores = 4)
        do.call(cbind, fit$A[1:2]) - A_coeff_ex[[i]]
      },
      mc.cores =4)
    
    for(j in norms){
      trunc_lasso_error = vapply(fit_diff, error_func, numeric(1), m_norm = j)
      table[[j]]$not_truncated[[distr]][[i]] = trunc_lasso_error
    }
  }
  return(table)
}

# gaussian
set.seed(200)
system.time({
  lasso_table = compute_errors(
    data = VAR_d_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "gauss"),
    table = lasso_table, distr = "gaussian"
  )
})

# t-2.1
# set.seed(200)
# lasso_table = compute_errors(
#   data = VAR_d_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "t", innov_df = 2.1),
#   table = lasso_table, distr = "t_2.1"
# )

# t-3
# set.seed(200)
# lasso_table = compute_errors(
#   data = VAR_d_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "t", innov_df = 3),
#   table = lasso_table, distr = "t_3"
# )

# t-4
# set.seed(200)
# lasso_table = compute_errors(
#   data = VAR_d_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "t", innov_df = 4),
#   table = lasso_table, distr = "t_4"
# )

# log-normal
# set.seed(200)
# system.time({
# lasso_table = compute_errors(
#   data = VAR_d_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "log"),
#   table = lasso_table, distr = "log_normal"
# )
# })


#### save ####
setwd("/user/work/zl22291/simulations/truncation/proj1_revisions/comment_4")
saveRDS(lasso_table, file = "VARorder_misspecific_gaussian_d0_1_d_2.rds")
rm(lasso_table)
gc()

