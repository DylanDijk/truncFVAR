##########################################################################################################################################
# This script gets the time for our proposed method for VAR(1) data generated with a Banded coefficient matrix
##########################################################################################################################################

# setwd("/user/work/zl22291/simulations/truncation")
# source(file = "functions/data_generation.R")
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/84f22825f88656a6b6efdb8cf628d1e7610397b0/functions/data_generation.R")

# source(file = "functions/estimation.R")
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/84f22825f88656a6b6efdb8cf628d1e7610397b0/functions/estimation.R")


#### lasso_vs_p ####
n_p = cbind(n = 200, p = seq(50, 300, 50))

nsim = 200
A_coeff = A_coeff_banded(n_p)

gc()

#######################################################################################
# Helper: build (yy, xxs)
build_var_design <- function(xx, var.order = 1) {
  n <- nrow(xx)
  p <- ncol(xx)
  
  yy <- as.vector(xx[n:(var.order + 1), ])         # keeping your exact indexing
  xxs <- xx[(n - 1):(1), , drop = FALSE]           # keeping your exact indexing
  
  xxs <- methods::as(
    Matrix::Matrix(
      kronecker(Matrix::Diagonal(p), xxs),
      sparse = TRUE
    ),
    "dgCMatrix"
  )
  list(yy = yy, xxs = xxs)
}
# function that uses cv.glmnet directly, to avoid using fitVAR
bench_one_p <- function(data) {
  
  xx <- cross_val_and_trunc(
    data,
    n_tau = 60, lag = 0, trim = 1, trim_d = 1,
    trim_all = F, max = T, cv_lag = T,
    standardise = FALSE
  )
  
  des <- build_var_design(xx, var.order = 1)
  
  
  # fit = glmnet::cv.glmnet(
  #   y = des$yy, x = des$xxs,
  #   intercept = FALSE, standardize = FALSE, parallel = TRUE, lambda = lambda_p200
  # )
  fit = glmnet::cv.glmnet(
    y = des$yy, x = des$xxs,
    intercept = FALSE, standardize = FALSE, parallel = TRUE
  )
}



# Given a data matrix this function truncates and then fits LASSO and returns matrix with difference to true coefficients.
# with parallelisation:
sparsevar_trunc_fit_diff = function(data, n_tau = 60, lag = 0, trim = 1, trim_d=1, trim_all=F, max = T, cv_lag=T){
  
  el_time = system.time({
    data = cross_val_and_trunc(data, n_tau, lag, trim, trim_d, trim_all, max = max, cv_lag=cv_lag, standardise = F)
    fit_diff = bench_one_p(data)
  })["elapsed"]
  
  return(as.numeric(el_time))
}


# gaussian
set.seed(200)

# banded VAR
data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "gauss")

time_list = setNames(
  lapply(names(data), function(x) vector("numeric", nsim)),
  names(data)
)

library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
system.time(
  for(i in names(data)){
    print(i)
    for(k in 1:nsim){
      print(k)  
      time_list[[i]][k] = sparsevar_trunc_fit_diff(data = data[[i]][[k]])
    }
  }  
)

#### save ####
setwd("/user/work/zl22291/simulations/truncation/proj1_revisions")
saveRDS(time_list, file = "trunc_time_list_banded_cvglmnet.rds")


