# In papers DGP notation: (F0),(V1),(S1)
# adaHuber results for:
# Figures 3 + 4.

# Generates results for VAR coefficient estimation error of adaHuber::adaHuber.cv.lasso method. 
# From Wang, Lili, Chao Zheng, et al. (2018). A New Principle for Tuning-Free Huber Regression.

# data generated from banded VAR model with innovation errors with diagonal covariance structure. 

# In order to compute errors for errors from different distributions requires uncommenting 
# And to generate results for against dimension or sample size.

#### sourcing functions ####
source(file = "truncation/functions/data_generation.R")
source(file = "truncation/functions/estimation.R")
#### creating list to store objects ####
list_setup = list(
  "adaHuber" = list(
    "t_2.1" = list(),
    "log_normal" = list(),
    "gaussian" = list()
  )
) 
list_setup_2 <- list(
  "M" = list_setup,
  "F" = list_setup,
  "l_2_inf" = list_setup # maximal euclidean norm of a column
)

############################# uncomment whether need to run against sample size n, or against dimension. #####################
#### lasso_vs_n ####
n_p = cbind(n = seq(100, 300, 50), p = 50)
#### lasso_vs_p ####
# n_p = cbind(n = 200, p = seq(50, 300, 50))
nsim = 200



A_coeff = A_coeff_banded(n_p)


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

# For each dimension type of data need to loop through to get vector of errors for that norm
error_func = function(m_norm, diff_mat){
  if(m_norm == "l_2_inf"){
    error = maxEuclideanNormR(diff_mat)
  } else {
    error = norm(diff_mat, paste0(m_norm))
  }
  return(error)
}


# Given a data matrix this function fits adaHuber.cv.lasso and returns matrix with difference to true coefficients.
library(doParallel)
num_cores <- detectCores()
ncores = 4
cl = makeCluster(ncores, type = "FORK")
registerDoParallel(cl)
# stopCluster(cl)


adahuber_var_parallel_fit_diff = function(data, A_coeff, ncores = 4){
  n = nrow(data)
  p = ncol(data)
  
  res = 
    foreach(i=1:p, .packages = c("adaHuber"))%dopar%{
      adaHuber::adaHuber.cv.lasso(X = data[(n-1):1,], Y = c(data[n:2,i]))$coef[-1]
      
    }
  
  A_adahuber_parallel = do.call(rbind, res)
  rm(res); gc()
  fit_diff = A_adahuber_parallel - A_coeff
  rm(A_adahuber_parallel); gc()
  return(fit_diff)
}

compute_errors = function(data, table, distr){
  for(i in names(data)){
    for(k in 1:nsim){
      fit_diff = adahuber_var_parallel_fit_diff(data = data[[i]][[k]], A_coeff = A_coeff[[i]])
      
      for(j in c("M", "F", "l_2_inf")){
        table[[j]]$adaHuber[[distr]][[i]][k] = error_func(fit_diff, m_norm = j)
      }
      rm(fit_diff)
      gc()
    }
  }  
  return(table)
}

############################ uncommonent lines related to the distribution of interest ################################################
# In order to replicate results the seed needs to be set before each time data is generated.
set.seed(200)
# gaussian
data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "gauss")
# t_21
# data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_df = 2.1, innov_dist = "t")
# lognormal
# data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "log")

for(j in c("M", "F", "l_2_inf")){
  for(i in names(data)){
    list_setup_2[[j]]$adaHuber$gaussian[[i]]= vector(length = nsim)
    # list_setup_2[[j]]$adaHuber$t_2.1[[i]]= vector(length = nsim)
    # list_setup_2[[j]]$adaHuber$log_normal[[i]]= vector(length = nsim)
  }
}
adaHuber_lasso_vs_n_gaussian = compute_errors(data = data, table = list_setup_2, distr = "gaussian")
# adaHuber_lasso_vs_n_t_21 = compute_errors(data = data, table = adaHuber_lasso_vs_n_t_21, distr = "t_2.1")
# adaHuber_lasso_vs_n_lognorm = compute_errors(data = data, table = adaHuber_lasso_vs_n_lognorm, distr = "log_normal")


