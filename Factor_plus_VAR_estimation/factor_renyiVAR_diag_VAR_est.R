# In papers DGP notation: (F1),(V2),(S1)
# Figure 6

# Generates results comparing VAR coefficient estimation error when using truncated and not truncated data for
# data generated from Factor plus Renyi random graph VAR model with innovation errors with diagonal covariance structure. 

# Here the autocovariance is included in the CV measure and in the error.

# No standardisation is carried out here.

# In order to compute errors for errors from different distributions requires uncommenting 

#### sourcing functions ####
source(file = "truncation/functions/data_generation.R")
source(file = "truncation/functions/estimation.R")
source(file = "truncation/functions/fac_adj_funcs/fac_adj_funcs.R")
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

list_setup_3 = list(
  "comm" = list_setup_2,
  "A_est" = list_setup_2
)

lasso_fac_ke_table = list_setup_3


# lasso_vs_n
# n_p = cbind(n = rep(seq(100, 300, 50),2), p = c(rep(30,5),rep(50,5)))
n = c(100,100,200,200,500,500)
p = c(50,100,50,100,100,200)
# n = c(100,100)
# p = c(50,100)
n_p = cbind(n = n, p = p)
r = 3
n_r = cbind(n = n, p = r)
nsim = 200
# lasso_vs_p
# n_p = cbind(n = 200, p = seq(50, 300, 50))
# nsim = 1
#data = VAR_1_data_ind_t(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_df = 2.1)

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



glm_pca_fun = function(x, q){
  if(q > 0){
    mpca = m.static.pca(x, q = q)
    xx_c = mpca$f %*% t(mpca$lam)
    xx_i = x - xx_c
    glmn = sparsevar::fitVAR(xx_i, p = 1, method = "cv", parallel = TRUE, ncores = 4)
  } else {
    glmn = sparsevar::fitVAR(x, p = 1)
  }  
  glmn_A = glmn$A[[1]]
  return(list(A = glmn_A, glmn = glmn, common = xx_c))
}





# Given a data matrix this function either truncates or does not truncate data and then
# computes the fitted coefficient error and the common component error
library(sparsevar)
trunc_fit_diff_comm = function(data, common_comp, A_coeff, trunc = T){
  if(trunc){
    data = cross_val_and_trunc(data, n_tau = 60, lag = 0, trim = 1, trim_d=1, trim_all=F, max = T, cv_lag=T, standardise = F)
  }
  glmnf = glm_pca_fun(x = data, q = 3)
  
  fit_diff = glmnf$A - A_coeff
  common_est = norm(glmnf$common - common_comp, "2") / norm(common_comp, "2")
  return(list(fit_diff = fit_diff, common_est = common_est))
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
  print(paste(distr))
  ## compute common components of a set of simulations for (n,p)
  # getting normalised loadig matrices
  lambda_a_norm = vector(mode = "list", length = nrow(n_p))
  for(i in 1:nrow(n_p)){
    for(k in 1:nsim){
      lambda_a_norm[[i]][[k]] = sweep(lambda_list[[i]][[k]], 2, attr(data, "A_norm")[[i]][k,], `*`) 
    }
  }
  common_comps = vector(mode = "list", length = nrow(n_p))
  common_comps = mapply(function(l1,l2) mapply(function(m1,m2) m2 %*% m1 ,l1, l2, SIMPLIFY = F),
                        lambda_a_norm, attributes(data)$fac_var, SIMPLIFY = F)
  names(common_comps) = names(A_coeff_banded(n_p))
  
  rm(lambda_a_norm)
  gc()
  
  # truncated
  for(truncs in c(T,F)){
    trunc_name = ifelse(truncs, "truncated", "not_truncated")
    for(i in names(data)){
      pca_lasso = mapply(trunc_fit_diff_comm, data[[i]], common_comps[[i]], A_coeff = attr(data, "A_var")[[i]], trunc = truncs, SIMPLIFY = F)
      for(j in c("M", "F", "l_2_inf")){
        trunc_lasso_error = sapply(pca_lasso, function(x){error_func(m_norm = j, diff_mat = x$fit_diff)})
        table$A_est[[j]][[trunc_name]][[distr]][[i]] = trunc_lasso_error
      }
      table$comm[["F"]][[trunc_name]][[distr]][[i]] = sapply(pca_lasso, function(x){x$common_est})
    }
  }
  return(table)
}



#################### uncomment corresponding to which distribution you want to run #############################
# for t dist, you also need to modify the innov_df argument.

# t
# set.seed(200)
# lasso_fac_ke_table = compute_errors(
#   data = factor_mod_data_VAR_idio(dist = "t", innov_df = 4, n_p = n_p, fac_var_ret = T, norm_loading = T,
#                                   indo_A = renyi_coeff),
#   table = lasso_fac_ke_table, distr = "t_4"
# )


# gaussian
# set.seed(200)
# lasso_fac_ke_table = compute_errors(
#   data = factor_mod_data_VAR_idio(dist = "gauss", n_p = n_p, fac_var_ret = T, norm_loading = T,
#                                   indo_A = renyi_coeff),
#   table = lasso_fac_ke_table, distr = "gaussian"
# )

# lognormal
set.seed(200)
lasso_fac_ke_table = compute_errors(
  data = factor_mod_data_VAR_idio(dist = "log", n_p = n_p, fac_var_ret = T, norm_loading = T,
                                  indo_A = renyi_coeff),
  table = lasso_fac_ke_table, distr = "log_normal"
)


# #### save ####
# saveRDS(lasso_fac_ke_table, file = "factor_plus_VAR/fac_var_lasso_cv_with_auto_renyi_lognormal.rds")
# # save(lasso_vs_n, file = "heavytailed_sparsevar_lasso_cv_with_max_norm/lasso_vs_n_cv_with_max_norm.RData")
# rm(lasso_fac_ke_table)
# gc()

