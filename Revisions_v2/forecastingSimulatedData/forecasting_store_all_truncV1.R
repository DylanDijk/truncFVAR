##########################################################################################################################################
# This script gets the next step forecast for simulated factor plus VAR data. And saves the "oracle" best linear predictor as well. 
  # Need to uncomment the code corresponding to the distribution of interest for generating the data.
  # These results are used to generate table containing forecasting performance results.
##########################################################################################################################################

# In this script I want to compute the forecast errors for simulated data.

# Need to first get functions to generate data
# This is now sourcing function that can generate stable banded VAR process of lags greater than 1
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/eb4b1a4dfbb040c82822cd6047fcb1484e674ebb/functions/data_generation.R")

# source(file = "functions/estimation.R")
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/84f22825f88656a6b6efdb8cf628d1e7610397b0/functions/estimation.R")
###################################################################################################################################################


library(fnets)

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

lasso_table = list_setup


##### Generate data #####
set.seed(500)
n = c(100,100,200,200,500,500)
p = c(50,100,50,100,100,200)
# n = c(100)
# p = c(50)
n_p = cbind(n = n, p = p)
r = 3
n_r = cbind(n = n, p = r)
nsim = 200

# n_p = cbind(n = c(100,100), p = c(100,200))
#data = VAR_1_data_ind_t(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_df = 2.1)

A_coeff = A_coeff_banded(n_p)

# setting seed before randomly generating other fixed quantities
set.seed(100)
# Lambda
lambda_list = vector(mode = "list", length = nrow(n_p))
names(lambda_list) <- names(A_coeff)

for(i in 1:nrow(n_p)){
  lambda_list[[i]] = vector(length = nsim, mode = "list")
  p_l = n_p[i,"p"]
  
  for(k in 1:nsim){
    lambda_list[[i]][[k]] = matrix(rnorm(r*p_l), nrow = r,  ncol = p_l)
  }
}
#####################################################################################################################################################
# I want to now look at forecasting where the forecasting for the common component projects the original data point X_n
# And the estimate of the idiosyncratic component is based on subtracting the insample estimate of the common component (with truncation) from the original data 

# modified version of fnets:::predict.fnets
predict_fnets_mod  = function(object, newdata = NULL, n.ahead = 1, fc.restricted = TRUE, 
                              r = c("ic", "er"), data) 
{
  if (is.null(newdata)) {
    newdata <- attr(object, "args")$x
  } else if (object$q >= 1) {
    stop("To produce forecasts when a common component is present, estimate a model on the new data. \n")
  }
  newdata <- as.ts(newdata)
  n.ahead <- fnets:::posint(n.ahead)
  if (nrow(newdata) < n.ahead) {
    n.ahead <- nrow(newdata)
    warning("Forecast horizon restricted by number of observations")
  }
  cpre <- common.predict_mod(object, t(newdata), n.ahead, fc.restricted,r, data_n = data[nrow(data),])
  ipre <- idio.predict_mod(object, t(newdata), cpre, n.ahead, og_data = data)
  out <- list(forecast = cpre$fc + ipre$fc, common.pred = cpre, 
              idio.pred = ipre, mean.x = object$mean.x, r = cpre$r)
  return(out)
}

# modified version of common.predict
common.predict_mod = function(object, x, n.ahead = 1, fc.restricted = TRUE, r = c("ic","er"), data_n) 
{
  xx <- x - object$mean.x
  p <- dim(x)[1]
  if (!is.numeric(r)) {
    r.method <- match.arg(r, c("ic", "er"))
    r <- NULL
  }
  else r.method <- NULL
  pre <- list(is = 0 * t(x), fc = matrix(0, nrow = n.ahead, 
                                         ncol = p))
  if (attr(object, "factor") == "unrestricted") {
    if (object$q < 1) {
      warning("There should be at least one factor for common component estimation!")
    }
    else {
      if (fc.restricted) 
        pre <- common.restricted.predict_mod(xx = xx, Gamma_x = object$acv$Gamma_x, 
                                             Gamma_c = object$acv$Gamma_c, q = object$q, 
                                             r = r, r.method = r.method, n.ahead = n.ahead, data_n = data_n)
      else pre <- common.unrestricted.predict(xx = xx, 
                                              cve = object, n.ahead = n.ahead)
    }
  }
  if (attr(object, "factor") == "restricted") {
    if (object$q < 1) {
      warning("There should be at least one factor for common component estimation!")
    }
    else if (object$q >= 1) {
      if (!fc.restricted) 
        warning("fc.restricted is being set to TRUE, as fnets object is generated with fm.restricted = TRUE")
      pre <- common.restricted.predict_mod(xx = xx, Gamma_x = object$acv$Gamma_x, 
                                           Gamma_c = object$acv$Gamma_c, q = object$q, r = object$q, 
                                           r.method = r.method, n.ahead = n.ahead, data_n = data_n)
    }
  }
  return(pre)
}

common.restricted.predict_mod = function(xx, Gamma_x, Gamma_c, q, r = NULL, max.r = NULL, r.method = NULL, n.ahead = 1, data_n = data_n) 
{
  p <- dim(xx)[1]
  n <- dim(xx)[2]
  if (is.null(max.r)) 
    max.r <- max(q, min(50, round(sqrt(min(n, p)))))
  if (n.ahead >= dim(Gamma_c)[3]) {
    warning("At most ", (dim(Gamma_c)[3] - 1)/2, "-step ahead forecast is available!")
    n.ahead <- (dim(Gamma_c)[3] - 1)/2
  }
  if (is.null(r)) {
    if (r.method == "ic") {
      abc <- abc.factor.number(xx, covx = Gamma_x[, , 1], 
                               q.max = max.r)
      r <- max(q, abc$q.hat[5])
      sv <- svd(Gamma_x[, , 1], nu = max.r, nv = 0)
    }
    else if (r.method == "er") {
      sv <- svd(Gamma_x[, , 1], nu = max.r, nv = 0)
      r <- which.max(sv$d[q:max.r]/sv$d[1 + q:max.r]) + 
        q - 1
    }
  }
  else sv <- svd(Gamma_x[, , 1], nu = max.r, nv = 0)
  is <- sv$u[, 1:r, drop = FALSE] %*% t(sv$u[, 1:r, drop = FALSE]) %*% 
    xx
  if (n.ahead >= 1) {
    fc <- matrix(0, nrow = p, ncol = n.ahead)
    # proj.x <- t(t(sv$u[, 1:r, drop = FALSE])/sv$d[1:r]) %*% 
    #   t(sv$u[, 1:r, drop = FALSE]) %*% xx[, n]
    proj.x <- t(t(sv$u[, 1:r, drop = FALSE])/sv$d[1:r]) %*% 
      t(sv$u[, 1:r, drop = FALSE]) %*% data_n
    for (hh in 1:n.ahead) fc[, hh] <- t(Gamma_c[, , hh + 
                                                  1]) %*% proj.x
  }
  else {
    fc <- NA
  }
  out <- list(is = as.ts(t(is)), fc = as.ts(t(fc)), r = r, 
              n.ahead = n.ahead)
  return(out)
}

# modified version of idio.predict that gets estimate of idiosyncratic component by subtracting estimated common component (still with truncation) from original data (without truncation)
idio.predict_mod = function (object, x, cpre, n.ahead = 1, og_data) 
{
  p <- dim(x)[1]
  n <- dim(x)[2]
  xx <- x - object$mean.x
  beta <- object$idio.var$beta
  if (is.null(beta)) 
    beta <- object$beta
  d <- dim(beta)[1]/p
  A <- t(beta)
  og_data <- as.ts(og_data)
  og_data = t(og_data)
  is <- og_data - t(cpre$is)
  if (n.ahead >= 1) {
    fc <- matrix(0, nrow = p, ncol = n.ahead)
    for (ii in 1:n.ahead) {
      for (ll in 1:d) fc[, ii] <- fc[, ii] + A[, p * (ll - 
                                                        1) + 1:p] %*% is[, n + ii - ll]
      is <- cbind(is, fc[, ii])
    }
  }
  else {
    fc <- NA
  }
  out <- list(is = as.ts(t(is[, 1:n])), fc = as.ts(t(fc)), 
              n.ahead = n.ahead)
  return(out)
}

#####################################################################################################################################################

# function that returns forecast after truncating and not truncating
# takes single matrix from simulation as input
library(Matrix)
forcast_tr_ntr = function(tr = FALSE, data, q = r, h = 1, A, D, fac_t, xi_n, lambda, A_norm){
  if(tr){
    fnet = fnets(data, fm.restricted = TRUE, robust = TRUE, robust.standardise = TRUE, var.args = list(n.cores = 4), center = FALSE, q = q)
  } else {
    fnet = fnets(data, fm.restricted = TRUE, robust = FALSE, var.args = list(n.cores = 4), center = FALSE, q = q)
  }
  forcast_list = vector(mode = "list", length = 3)
  forcast_list <- lapply(forcast_list, function(x) NA)
  names(forcast_list) = c("combined_pred", "best_linear_forecast", "factor_r")
  
  pred = predict_fnets_mod(fnet, n.ahead = h, fc.restricted = TRUE, data = data)
  
  # scale lambda if needed
  lambda_scaled = sweep(lambda, 2, A_norm, "*")
  # factor part
  factor_part = as.numeric(fac_t[nrow(fac_t), ] %*% t(D) %*% lambda_scaled)
  # oracle forecast
  best_linear = as.numeric((A %*% xi_n)) + factor_part
  
  forcast_list[[1]] = as.numeric(pred$forecast[1, ])
  forcast_list[[2]] = best_linear
  forcast_list[[3]] = pred$r
  
  return(forcast_list)
}


# forecast_error = function(pred, oracle){
#   sum((pred - oracle)^2) / sum(oracle^2)
# }


# takes the full data object, then for each n_p element in list, computes the forecast for each simulation.
compute_errors = function(data, table, distr){
  
  for(i in names(data)){
    
    print(i)
    
    F_list  = attr(data, "fac_var")[[i]]
    A_norm  = attr(data, "A_norm")[[i]]
    A_norm = split(A_norm, row(A_norm))
    lambda  = lambda_list[[i]]
    xi_list = attr(data, "xi_n")[[i]]
    D_list = attr(data, "A_fac")[[i]]
    
    for(truncs in c(T,F)){
      trunc_name = ifelse(truncs, "truncated", "not_truncated")
      
      table[[trunc_name]][[distr]][[i]] = mapply(function(X, D, fac_t, lam, norm_vec, xi_n){
        
        fc = forcast_tr_ntr(
          tr = truncs,
          data = X,
          A = A_coeff[[i]],
          D = D,
          fac_t = fac_t,
          xi_n = xi_n,
          lambda = lam,
          A_norm = norm_vec
        )
        
        pred   = as.numeric(fc$combined_pred)
        oracle = as.numeric(fc$best_linear_forecast)
        
        list(pred = pred, oracle = oracle, factor_r = fc$factor_r)
        
      }, data[[i]], D_list, F_list, lambda, A_norm, xi_list, SIMPLIFY = FALSE)
    }
    
  }
  
  return(table)
}



# gaussian
# system.time({
# set.seed(200)
# data = factor_mod_data_VAR_idio(dist = "gauss", n_p = n_p, fac_var_ret = T, norm_loading = T, xi_n = TRUE)
# lasso_fac_ke_table = compute_errors(data = data,
#                                     table = lasso_table, distr = "gaussian"
# )
# })

# t21
# set.seed(200)
# data = factor_mod_data_VAR_idio(dist = "t", innov_df = 2.1, n_p = n_p, fac_var_ret = T, norm_loading = T, xi_n = TRUE)
# lasso_fac_ke_table = compute_errors(data = data,
#                                     table = lasso_table, distr = "t_2.1")


# t3
# set.seed(200)
# data = factor_mod_data_VAR_idio(dist = "t", innov_df = 3, n_p = n_p, fac_var_ret = T, norm_loading = T, xi_n = TRUE)
# lasso_fac_ke_table = compute_errors(data = data,
#                                     table = lasso_table, distr = "t_3")

# t4
# set.seed(200)
# data = factor_mod_data_VAR_idio(dist = "t", innov_df = 4, n_p = n_p, fac_var_ret = T, norm_loading = T, xi_n = TRUE)
# lasso_fac_ke_table = compute_errors(data = data,
#                                     table = lasso_table, distr = "t_4")

# log_normal
set.seed(200)
data = factor_mod_data_VAR_idio(dist = "log", n_p = n_p, fac_var_ret = T, norm_loading = T, xi_n = TRUE)
lasso_fac_ke_table = compute_errors(data = data,
                                    table = lasso_table, distr = "log_normal")

setwd("/user/work/zl22291/simulations/truncation/proj1_revisions/comment_5")
# saveRDS(lasso_fac_ke_table, file = "fac_var_forecasting_t21_truncV1.Rds")
# saveRDS(lasso_fac_ke_table, file = "fac_var_forecasting_t3_truncV1.Rds")
# saveRDS(lasso_fac_ke_table, file = "fac_var_forecasting_t4_truncV1.Rds")
saveRDS(lasso_fac_ke_table, file = "fac_var_forecasting_lognormal_truncV1.Rds")
# saveRDS(lasso_fac_ke_table, file = "fac_var_forecasting_Gaussian_truncV1.Rds")



