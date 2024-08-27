# In papers DGP notation: (F0),(VO),(S1)
# Table 1

# Generates results comparing truncated and sample covariance matrix estimates for
# independent samples with diagonal covariance structure. 

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
list_setup_2 <- list(
  "M" = list_setup,
  "F" = list_setup
)
ke_table = list_setup_2


# Setting max = T, so that max norm will be used in cross validation as opposed to frobenius norm.
compute_errors = function(data, distr, table, m_norm, max = T){
  error = function(x){norm(x - diag(nrow(x)), paste0(m_norm))}
  trunc_cov_error = lapply(data, cross_val_and_error, n_tau = 60, lag = 0, trim = 1, trim_d = 1, trim_all = F, max, error = error)
  table[[m_norm]]$truncated[[distr]] = trunc_cov_error
  
  cov_list = lapply(data, function(x){lapply(x,cov)})
  norm_cov_diff = lapply(cov_list, function(x){sapply(x, error)})
  table[[m_norm]]$not_truncated[[distr]] = norm_cov_diff
  
  return(table)
}



n_p_values = cbind(n = c(50, 50, 100, 200, 200), p = c(100, 200, 200, 50, 100))
nsim = 200

set.seed(250)
library(stringr)
for(j in c("M", "F")){
  ke_table = compute_errors(data = ind_data_ind_gauss(nsim = nsim, n_p = n_p_values),
                            table = ke_table, distr = "gaussian", max = T, m_norm = j)
  
  # t distribution
  for(i in c("t_2.1", "t_3", "t_4")){
    ke_table = compute_errors(data = ind_data_ind_t(nsim = nsim, n_p = n_p_values, df = as.numeric(str_extract_all(i, "\\d+\\.\\d+|\\d+"))),
                              table = ke_table, distr = i, max = T, m_norm = j)
  }
  #log-normal
  ke_table = compute_errors(data = ind_data_ind_lognorm(nsim = nsim, n_p = n_p_values),
                            table = ke_table, distr = "lognormal", max = T, m_norm = j)
}










