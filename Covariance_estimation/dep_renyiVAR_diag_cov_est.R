# In papers DGP notation: (F0),(V2),(S1)
# Table 4

# Generates results comparing truncated and sample covariance matrix estimates for
# dependent samples, generated from VAR(1) model with Renyi coefficient matrix, innovation errors with diagonal covariance structure. 

# These simulations standardised the data, by scaling up the truncation parameter tau by the MAD of the corresponding variable.
# This makes little difference, as the data has been generated from a VAR(1) model.

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

compute_errors = function(data, distr, table, m_norm, cv_lag = F){
  
  # covariances of VAR
  cov_var = vector(mode = "list", length = nrow(n_p))
  names(cov_var) = names(A_coeff_banded(n_p))
  
  for(i in 1:nrow(n_p)){
    cov_var[[i]] = vector(mode = "list", length = nsim)
    for(k in 1:nsim){
      cov_var[[i]][[k]] = cov_of_var(A = attributes(data[[i]])$A[[k]])   
    }  
  } 
  
  
  trunc_cov_error = vector(mode = "list", length = length(names(data)))
  names(trunc_cov_error) = names(data)
  for(i in names(data)){
    for(j in c("F", "M")){
      table[[j]]$truncated[[distr]][[i]] = 
        cross_val_and_error(data = data[[i]], cv_lag = cv_lag,
                            n_tau = 60, lag = 0, trim = 1, trim_d = 1, trim_all = F, max = T, error = "list", true_cov = cov_var[[i]], error_norm = paste(j))
    }
  }  
  
  # error = function(y){norm(y - combined_cov[[paste0(i)]], paste0(j))}
  
  cov_list = lapply(data, function(x){lapply(x,acf_no_center)})
  norm_cov_diff = vector(mode = "list", length = length(names(data)))
  names(norm_cov_diff) = names(data)
  for(i in names(data)){
    for(j in c("F", "M")){
      table[[j]]$not_truncated[[distr]][[i]] =  unlist(Map(function(x,y){norm(x-y, paste(j))}, cov_list[[i]], cov_var[[i]]))
    }
  }
  
  return(table)
}


# gaussian

ke_table = list_setup_2
set.seed(200)
# ke_table = compute_errors(data = VAR_1_data_ind_gauss(nsim = nsim, n_p = n_p, A_coeff = A_coeff),
#                             table = ke_table, distr = "gaussian")
ke_table = compute_errors(data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = renyi_coeff, innov_dist = "gauss"),
                          table = ke_table, distr = "gaussian")


# t distribution

library(stringr)
for(i in c("t_2.1", "t_3", "t_4")){
  set.seed(200)
  # ke_table = compute_errors(data = VAR_1_data_ind_t(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_df = as.numeric(str_extract_all(i, "\\d+\\.\\d+|\\d+"))),
  #                           table = ke_table, distr = i)
  ke_table = compute_errors(data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = renyi_coeff, innov_dist = "t",
                                                  innov_df = as.numeric(str_extract_all(i, "\\d+\\.\\d+|\\d+"))),
                            table = ke_table, distr = i)
}

# lognormal
set.seed(200)
# ke_table = compute_errors(data = VAR_1_data_ind_lognorm(nsim = nsim, n_p = n_p, A_coeff = A_coeff),
#                           table = ke_table, distr = "lognormal")
ke_table = compute_errors(data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = renyi_coeff, innov_dist = "log"),
                          table = ke_table, distr = "lognormal")
