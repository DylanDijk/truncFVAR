# In papers DGP notation: (F0),(VO),(S2)
# Table 2

# Generates results comparing truncated and sample covariance matrix estimates for
# independent samples with power decay covariance structure, 0.9^{|i-j|}. 

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

ke_table_power_decay = list_setup_2


n_p_values = cbind(n = c(50, 50, 100, 200, 200), p = c(100, 200, 200, 50, 100))
nsim = 200

##### power decay #####
covstr = vector(mode = "list", length = nrow(n_p_values))
power_dec = vector(mode = "list", length = nrow(n_p_values)) 
for(k in 1:nrow(n_p_values)){
  p = n_p_values[k, "p"]
  power_dec[[k]] = matrix(nrow = p, ncol = p)
  for(i in 1:p){
    for(j in 1:p){
      power_dec[[k]][i,j] = 0.9^(abs(i-j))
    }
  }
  covstr[[k]] = with(
    eigen(power_dec[[k]]),
    vectors %*% sqrt(diag(values)) %*% t(vectors)
  )
}
for (i in 1:nrow(n_p_values)){
  names(covstr)[i] <- paste("(", n_p_values[i, "n"], ",", n_p_values[i, "p"], ")", sep = "")
  names(power_dec)[i] <- paste("(", n_p_values[i, "n"], ",", n_p_values[i, "p"], ")", sep = "")
}
# heatmap(covstr$`(200,50)`, Rowv = NA, Colv = NA)
##########

compute_errors = function(data, distr, table, m_norm, max = T){
  
  trunc_cov_error = vector(mode = "list", length = length(names(data)))
  names(trunc_cov_error) = names(data)
  for(i in names(data)){
    trunc_cov_error[[i]] = cross_val_and_error(data = data[[i]], 
                                               n_tau = 60, lag = 0, trim = 1, trim_d = 1, trim_all = F,max = max, error = function(y){norm(y - power_dec[[paste0(i)]], paste0(m_norm))})
  }
  table[[m_norm]]$truncated[[distr]] = trunc_cov_error
  
  cov_list = lapply(data, function(x){lapply(x,cov)})
  norm_cov_diff = vector(mode = "list", length = length(names(data)))
  names(norm_cov_diff) = names(data)
  for(i in names(data)){
    norm_cov_diff[[i]] = sapply(cov_list[[i]], function(y){norm(y - power_dec[[paste0(i)]], paste0(m_norm))})
  }
  
  table[[m_norm]]$not_truncated[[distr]] = norm_cov_diff
  
  return(table)
}


set.seed(350)
library(stringr)
for(j in c("M", "F")){
  print(j)
  ke_table_power_decay = compute_errors(data = ind_data_ind_gauss(nsim = nsim, n_p = n_p_values, covstr = covstr),
                                        table = ke_table_power_decay, distr = "gaussian", m_norm = j)
  
  # t distribution
  for(i in c("t_2.1", "t_3", "t_4")){
    print(i)
    ke_table_power_decay = compute_errors(
      data = ind_data_ind_t(nsim = nsim, n_p = n_p_values, df = as.numeric(str_extract_all(i, "\\d+\\.\\d+|\\d+")), covstr = covstr),
      table = ke_table_power_decay, distr = i, m_norm = j)
  }
  #log-normal
  ke_table_power_decay = compute_errors(data = ind_data_ind_lognorm(nsim = nsim, n_p = n_p_values, covstr = covstr),
                                        table = ke_table_power_decay, distr = "lognormal", m_norm = j)
}



