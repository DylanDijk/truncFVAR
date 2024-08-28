# In papers DGP notation: (F0),(V2),(S1)
# Table 6

# Generates results comparing truncated and sample covariance and autocovariance matrix estimates for
# dependent samples, generated from VAR(1) model with renyi random graph coefficient matrix, innovation errors with diagonal covariance structure. 

# Here the autocovariance is included in the CV measure and in the error.

# These simulations standardised the data, by scaling up the truncation parameter tau by the MAD of the corresponding variable.
# This makes little difference, as the data has been simulated so that each variable has variance 1.

#### sourcing functions ####
source(file = "truncation/functions/data_generation.R")
source(file = "truncation/functions/estimation.R")
#### creating list to store objects ####
list_setup = list(
  "truncated" = list(
    "t_2.1" = list(),
    "t_3" = list(),
    "t_4" = list(),
    "gaussian" = list(),
    "lognormal" = list()
  ),
  "not_truncated" = list(
    "t_2.1" = list(),
    "t_3" = list(),
    "t_4" = list(),
    "gaussian" = list(),
    "lognormal" = list()
  )
)

# list to store objects
list_setup_2 <- list(
  "M" = list_setup,
  "F" = list_setup
)


n_p = cbind(n = c(50, 50, 100, 200, 200), p = c(100, 200, 200, 50, 100))
nsim = 200

# cv_lag set to TRUE, so including autocov in CV measure
# cov_error_comb set to TRUE, so computing error using both cov and autocov
compute_errors = function(data, distr, table, m_norm, cv_lag = T, cov_error_comb = T){
  
  # true covariance matrices
  cov_var_l = vector(mode = "list", length = nrow(n_p))
  names(cov_var_l) = names(A_coeff_banded(n_p))
  
  for(i in 1:nrow(n_p)){
    cov_var_l[[i]] = vector(mode = "list", length = nsim)
    for(k in 1:nsim){
      cov_var_l[[i]][[k]] = as.matrix(cov_of_var(A = attributes(data[[i]])$A[[k]]))
    }  
  } 
  # true autocovariance matrices
  autocov_var_l = vector(mode = "list", length = nrow(n_p)); names(autocov_var_l) = names(cov_var_l)
  for(i in names(autocov_var_l)){
    autocov_var_l[[i]] = vector(mode = "list", length = nsim)
    for(k in 1:nsim){
      autocov_var_l[[i]][[k]] = attributes(data[[i]])$A[[k]] %*% cov_var_l[[i]][[k]]
    }
  }
  
  combcov_var_l = vector(mode = "list", length = nrow(n_p)); names(combcov_var_l) = names(cov_var_l)
  # combcov_var_l = mapply(function(m1, m2) cbind(m1, m2), cov_var_l, autocov_var_l, SIMPLIFY = FALSE)
  combcov_var_l = mapply(function(l1,l2) mapply(function(m1, m2) cbind(m1, m2), l1, l2, SIMPLIFY = FALSE), cov_var_l, autocov_var_l, SIMPLIFY = FALSE)
  
  trunc_cov_error = vector(mode = "list", length = length(names(data)))
  names(trunc_cov_error) = names(data)
  for(i in names(data)){
    for(j in c("F", "M")){
      table[[j]]$truncated[[distr]][[i]] = 
        cross_val_and_error(data = data[[i]], cv_lag = cv_lag, cov_error_comb = cov_error_comb,
                            n_tau = 60, lag = 0, trim = 1, trim_d = 1, trim_all = F, max = T,
                            error = "list", true_cov = combcov_var_l[[i]], error_norm = paste(j))
    }
  }  
  
  cov_list = lapply(data, function(x){lapply(x,cov)})
  autocov_list = lapply(data, function(x){lapply(x,acf_no_center, lag = 1)})
  combcov_list = mapply(function(l1,l2) mapply(function(m1, m2) cbind(m1, m2), l1, l2, SIMPLIFY = FALSE), cov_list, autocov_list, SIMPLIFY = FALSE)
  norm_cov_diff = vector(mode = "list", length = length(names(data)))
  names(norm_cov_diff) = names(data)
  for(i in names(data)){
    for(j in c("F", "M")){
      table[[j]]$not_truncated[[distr]][[i]] = Map(function(x,y){norm(x-y, paste(j))}, combcov_list[[i]], combcov_var_l[[paste0(i)]])
    }
  }
  
  return(table)
}

# gaussian

ke_table = list_setup_2
set.seed(200)

ke_table = compute_errors(data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = renyi_coeff, innov_dist = "gauss"),
                          table = ke_table, distr = "gaussian")


# t distribution

library(stringr)
for(i in c("t_2.1", "t_3", "t_4")){
  set.seed(200)
  ke_table = compute_errors(data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = renyi_coeff, innov_dist = "t",
                                                  innov_df = as.numeric(str_extract_all(i, "\\d+\\.\\d+|\\d+"))),
                            table = ke_table, distr = i)
}

# lognormal
set.seed(200)
ke_table = compute_errors(data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = renyi_coeff, innov_dist = "log"),
                          table = ke_table, distr = "lognormal")


###############################################################################################################################

