# In papers DGP notation: (F1),(V1),(S1)
# Table 7

# Generates results comparing truncated and sample covariance and autocovariance matrix estimates for
# dependent samples, generated from Factor plus VAR(1) model with banded coefficient matrix, innovation errors with diagonal covariance structure. 
# The factor model has the number of factors fixed at 3.

# Here the autocovariance is included in the CV measure and in the error.

# I now do not standardise the truncation with MAD. And I am using the normalising constant with the loading matrix,
# as is done in fnets simulations (C2).

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


# lasso np values
n = c(100,100,200,200,500,500)
p = c(50,100,50,100,100,200)
n_p = cbind(n = n, p = p)
r = 3
n_r = cbind(n = n, p = r)
nsim = 200


# Banded VAR coefficient matrix
A_coeff = A_coeff_banded(n_p)
library(Matrix)

# setting seed before randomly generating other fixed quantities
set.seed(100)
# Lambda
lambda_list = vector(mode = "list", length = nrow(n_p))
names(lambda_list) <- names(A_coeff_banded(n_p))

for(i in 1:nrow(n_p)){
  lambda_list[[i]] = vector(length = nsim, mode = "list")
  p_l = n_p[i,"p"]
  
  for(k in 1:nsim){
    lambda_list[[i]][[k]] = matrix(rnorm(r*p_l), nrow = r,  ncol = p_l)
  }
}



#### Computing true covariances ####
# covariances of idio
cov_var_i = vector(mode = "list", length = nrow(n_p))
names(cov_var_i) = names(A_coeff)

for(i in names(cov_var_i)){
  cov_var_i[[i]] = round(as.matrix(cov_of_var(A = A_coeff[[i]])),3)
}



#### Computing errors ####
compute_errors = function(data, distr, table, m_norm, norm_loading = T){
  
  
  # covariances of factor VAR
  cov_var_fac = vector(mode = "list", length = nrow(n_p))
  names(cov_var_fac) = names(cov_var_i)
  
  for(i in 1:nrow(n_p)){
    cov_var_fac[[i]] = vector(mode = "list", length = nsim)
    for(k in 1:nsim){
      cov_var_fac[[i]][[k]] = cov_of_var_manual(attr(data, "A_fac")[[i]][[k]])   
    }  
  }
  
  # autocovariance of factor var
  autocov_var_fac = vector(mode = "list", length = nrow(n_p)); names(autocov_var_fac) = names(cov_var_fac)
  for(i in names(autocov_var_fac)){
    autocov_var_fac[[i]] = vector(mode = "list", length = nsim)
    for(k in 1:nsim){
      autocov_var_fac[[i]][[k]] = ( attr(data, "A_fac")[[i]][[k]] ) %*% cov_var_fac[[i]][[k]]
    }
  }
  
  # covariances of combined
  cov_var_comb = vector(mode = "list", length = nrow(n_p))
  names(cov_var_comb) = names(cov_var_i)
  
  for(i in 1:nrow(n_p)){
    cov_var_comb[[i]] = vector(mode = "list", length = nsim)
    for(k in 1:nsim){
      if(norm_loading){
        lambda = t(sweep(lambda_list[[i]][[k]], 2, attr(data, "A_norm")[[i]][k,], `*`))
      }else{
        lambda = t(lambda_list[[i]][[k]])
      }  
      cov_var_comb[[i]][[k]] = (lambda %*% cov_var_fac[[i]][[k]]) %*% t(lambda) + cov_var_i[[i]]
    }
  }
  
  
  #autocovariances of combined
  autocov_var_comb = vector(mode = "list", length = nrow(n_p))
  names(autocov_var_comb) = names(cov_var_i)
  
  for(i in 1:nrow(n_p)){
    autocov_var_comb[[i]] = vector(mode = "list", length = nsim)
    for(k in 1:nsim){
      if(norm_loading){
        lambda = t(sweep(lambda_list[[i]][[k]], 2, attr(data, "A_norm")[[i]][k,], `*`))
      }else{
        lambda = t(lambda_list[[i]][[k]])
      }  
      autocov_var_comb[[i]][[k]] = (lambda %*% autocov_var_fac[[i]][[k]]) %*% t(lambda)  +  
        (A_coeff[[i]]) %*% cov_var_i[[i]]
    }
  }
  
  combcov_var = vector(mode = "list", length = nrow(n_p)); names(combcov_var) = names(cov_var_i)
  # combcov_var_l = mapply(function(m1, m2) cbind(m1, m2), cov_var_l, autocov_var_l, SIMPLIFY = FALSE)
  combcov_var = mapply(function(l1,l2) mapply(function(m1, m2) cbind(m1, m2), l1, l2, SIMPLIFY = FALSE), cov_var_comb, autocov_var_comb, SIMPLIFY = FALSE)
  
  
  
  trunc_cov_error = vector(mode = "list", length = length(names(data)))
  names(trunc_cov_error) = names(data)
  for(i in names(data)){
    for(j in c("F", "M")){
      table[[j]]$truncated[[distr]][[i]] = 
        cross_val_and_error(data = data[[i]], cv_lag = T, cov_error_comb = T, standardise = F,
                            n_tau = 60, lag = 0, trim = 1, trim_d = 1, trim_all = F, max = T,
                            error = "list", true_cov = combcov_var[[i]], error_norm = paste(j))
    }
  }  
  
  # error = function(y){norm(y - combined_cov[[paste0(i)]], paste0(j))}
  
  cov_list = lapply(data, function(x){lapply(x,cov)})
  autocov_list = lapply(data, function(x){lapply(x,acf_no_center, lag = 1)})
  combcov_list = mapply(function(l1,l2) mapply(function(m1, m2) cbind(m1, m2), l1, l2, SIMPLIFY = FALSE), cov_list, autocov_list, SIMPLIFY = FALSE)
  norm_cov_diff = vector(mode = "list", length = length(names(data)))
  names(norm_cov_diff) = names(data)
  for(i in names(data)){
    for(j in c("F", "M")){
      table[[j]]$not_truncated[[distr]][[i]] = Map(function(x,y){norm(x-y, paste(j))}, combcov_list[[i]], combcov_var[[paste0(i)]])
    }
  }
  
  return(table)
}





# gaussian
ke_table = list_setup_2
set.seed(200)
ke_table = compute_errors(data = factor_mod_data_VAR_idio(dist = "gauss", n_p = n_p, norm_loading = T),
                          table = ke_table, distr = "gaussian")


# t distribution
for(i in c(2.1, 3, 4)){
  set.seed(200)
  ke_table = compute_errors(data = factor_mod_data_VAR_idio(dist = "t", innov_df = i, n_p = n_p, norm_loading = T),
                            table = ke_table, distr = paste("t",sep = "_", i))
}

# log normal
set.seed(200)
ke_table = compute_errors(data = factor_mod_data_VAR_idio(dist = "log", n_p = n_p, norm_loading = T),
                          table = ke_table, distr = "lognormal")


