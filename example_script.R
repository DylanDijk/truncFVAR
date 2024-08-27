# example I could put on github
source(file = "truncation/functions/data_generation.R")
source(file = "truncation/functions/estimation.R")
n_p = cbind(n = c(200), p = c(50))
nsim = 1
A_coeff = A_coeff_banded(n_p)
data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "t", innov_df = 4)
data = data$`(200,50)`[[1]]
data_tr = cross_val_and_trunc(data, n_tau = 60, cv_lag = T, lag = 0, trim = 1, trim_d=1, trim_all=F, max = T, standardise = F)
for(i in 1:50){
  plot(data_tr[,i], type = "l", ylim = c(min(data), max(data)), main = paste(i))
  plot(data[,i], type = "l", ylim = c(min(data), max(data)), main = paste(i))
}
