# In papers DGP notation: (F0),(V1),(S1)
# Table 3

# Generates results comparing truncated and sample covariance matrix estimates for
# dependent samples, generated from VAR(1) model with banded coefficient matrix, innovation errors with diagonal covariance structure. 

# These simulations standardised the data, by scaling up the truncation parameter tau by the MAD of the corresponding variable.
# This makes little difference, as the data has been generated from a VAR(1) model.

#### sourcing functions ####
source(file = "functions/data_generation.R")
source(file = "functions/estimation.R")
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
list_setup_2 <- list(
  "M" = list_setup,
  "F" = list_setup
)


n_p = cbind(n = c(50, 50, 100, 200, 200), p = c(100, 200, 200, 50, 100))
nsim = 200

# Banded VAR coefficient matrix
A_coeff = A_coeff_banded(n_p)

cov_var_l = vector(mode = "list", length = nrow(n_p))
for(k in 1:nrow(n_p)){
  names(cov_var_l)[k] <- paste("(", n_p[k, "n"], ",", n_p[k, "p"], ")", sep = "")
}
for(i in names(cov_var_l)){
  cov_var_l[[i]] = round(as.matrix(cov_of_var(A = A_coeff[[i]])),3)
}

compute_errors = function(data, distr, table, m_norm, cv_lag = F){
  
  trunc_cov_error = vector(mode = "list", length = length(names(data)))
  names(trunc_cov_error) = names(data)
  for(i in names(data)){
    for(j in c("F", "M")){
      table[[j]]$truncated[[distr]][[i]] = 
        cross_val_and_error(data = data[[i]], cv_lag = cv_lag,
                            n_tau = 60, lag = 0, trim = 1, trim_d = 1, trim_all = F, max = T, error = function(y){norm(y - cov_var_l[[paste0(i)]], paste0(j))})
    }
  }  
  
  cov_list = lapply(data, function(x){lapply(x,cov)})
  norm_cov_diff = vector(mode = "list", length = length(names(data)))
  names(norm_cov_diff) = names(data)
  for(i in names(data)){
    for(j in c("F", "M")){
      table[[j]]$not_truncated[[distr]][[i]] = sapply(cov_list[[i]], function(y){norm(y - cov_var_l[[paste0(i)]], paste0(j))})
    }
  }
  
  return(table)
}

ke_table = list_setup_2
# gaussian
set.seed(200)
ke_table = compute_errors(data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "gauss"),
                          table = ke_table, distr = "gaussian")


# t distribution
library(stringr)
for(i in c("t_2.1", "t_3", "t_4")){
  set.seed(200)
  ke_table = compute_errors(data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "t",
                                                  innov_df = as.numeric(str_extract_all(i, "\\d+\\.\\d+|\\d+"))),
                            table = ke_table, distr = i)
}

# lognormal
set.seed(200)

ke_table = compute_errors(data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "log"),
                          table = ke_table, distr = "lognormal")
