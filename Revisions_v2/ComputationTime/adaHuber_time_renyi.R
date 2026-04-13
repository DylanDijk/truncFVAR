##########################################################################################################################################
# This script gets the time for adaHuber for VAR(1) data generated with a Renyi coefficient matrix
##########################################################################################################################################

source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/84f22825f88656a6b6efdb8cf628d1e7610397b0/functions/data_generation.R")
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/84f22825f88656a6b6efdb8cf628d1e7610397b0/functions/estimation.R")

#### lasso_vs_p ####
n_p = cbind(n = 200, p = seq(50, 300, 50))
nsim = 200


# Given a data matrix this function fits adaHuber.cv.lasso and returns matrix with difference to true coefficients.
library(doParallel)
num_cores <- detectCores()
ncores = 4
cl = makeCluster(ncores, type = "FORK")
registerDoParallel(cl)
# stopCluster(cl)


adahuber_var_parallel_fit_diff = function(data, ncores = 4){
  n = nrow(data)
  p = ncol(data)
  
  el_time = system.time({
    res = 
      foreach(i=1:p, .packages = c("adaHuber"))%dopar%{
        adaHuber::adaHuber.cv.lasso(X = data[(n-1):1,], Y = c(data[n:2,i]))
      }})["elapsed"]    
  
  return(as.numeric(el_time))
  
}


set.seed(200)
# gaussian
# renyi VAR
data =VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = renyi_coeff, innov_dist = "gauss")

time_list = setNames(
  lapply(names(data), function(x) vector("numeric", nsim)),
  names(data)
)

system.time(
  for(i in names(data)){
    print(i)
    for(k in 1:nsim){
      print(k)  
      time_list[[i]][k] = adahuber_var_parallel_fit_diff(data = data[[i]][[k]])
    }
  }  
)

### save ####
setwd("/user/work/zl22291/simulations/truncation/proj1_revisions")
saveRDS(time_list, file = "adaHuber_time_list_renyi.rds")



