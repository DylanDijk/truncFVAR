##########################################################################################################################################
# This script gets the VAR coefficient estimates with and without truncation, when the VAR order used for fitting is 3, but the true VAR order is 2.
# Results are used to create the heatmap plot of the average truncated estimate.
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
d_0 = 1:3
d_0 = 2

# what lags to fit 
d = 3

#### creating list to store objects ####
list_setup = list(
  "truncated" = list(
    "t_2.1" = list()
  ),
  "not_truncated" = list(
    "t_2.1" = list()
  )
)

lasso_table = list_setup


############################# uncomment whether need to run against sample size n, or against dimension, or (n,p) settings used for tables #####################
#### lasso_vs_n ####
# n_p = cbind(n = seq(100, 300, 50), p = 50)
#### lasso_vs_p ####
# n_p = cbind(n = 200, p = seq(50, 300, 50))
#### table values ####
# n_p = cbind(n = c(100,100,200,200,500,500), p = c(50,100,50,100,100,200))


n_p = cbind(n = c(200), p = c(50))


# n_p = cbind(n = c(100,100), p = c(50,100))
# n_p = cbind(n = c(100,100), p = c(5,7))

nsim = 200
#data = VAR_1_data_ind_t(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_df = 2.1)

A_coeff = A_coeff_banded_d(n_p, d = d_0)
# now adding zero matrix to side of var coefficient matrices.
if(length(A_coeff[[1]]) < d){
  A_coeff_ex <- lapply(A_coeff, function(A_list) {
    
    # Bind all matrices inside the list horizontally
    A_combined <- do.call(cbind, A_list)
    
    p <- nrow(A_combined)
    
    # Append zero block
    cbind(A_combined, matrix(0, nrow = p, ncol = p))
  })
}



# For each dimension type of data need to loop through to get vector of errors for that norm



# Given a data matrix this function truncates and then fits LASSO and returns matrix with difference to true coefficients.
library(sparsevar)
sparsevar_trunc_fit_diff = function(data, n_tau = 60, lag = d, trim = 1, trim_d=1, trim_all=F, max = T, m_norm, cv_lag=T){
  data = cross_val_and_trunc(data, n_tau, lag, trim, trim_d, trim_all, max = max, cv_lag=cv_lag, standardise = F)
  fit_diff = fitVAR(data, p = d, parallel = T, ncores = 4)
  
  fit_diff = cbind(fit_diff$A[[1]], fit_diff$A[[2]], fit_diff$A[[3]])
  return(fit_diff)
}



compute_errors = function(data, table, distr){
  # truncated
  for(i in names(data)){
    fit_diff = lapply(data[[i]], sparsevar_trunc_fit_diff)
  }
  
  table$truncated$t_2.1[[i]] = fit_diff
  
  # not-truncated
  for(i in names(data)){
    fit_diff = lapply(data[[i]], function(x){
      fit = fitVAR(x, p = d, parallel = T, ncores = 4)
      cbind(fit$A[[1]], fit$A[[2]], fit$A[[3]])
    })
  }
  
  table$not_truncated$t_2.1[[i]] = fit_diff
  
  return(table)
}


# t-2.1
set.seed(200)
system.time({
  lasso_table = compute_errors(
    data = VAR_d_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "t", innov_df = 2.1),
    table = lasso_table, distr = "t_2.1"
  )
})


#### save ####
setwd("/user/work/zl22291/simulations/truncation/proj1_revisions/comment_4")
saveRDS(lasso_table, file = "heatmap_results__d0_2_d_3.rds")
rm(lasso_table)
gc()


