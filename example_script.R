# example I could put on github
source(file = "truncation/functions/data_generation.R")
source(file = "truncation/functions/estimation.R")
n = 200; p = 50
n = 200; p = 20
n_p = cbind(n = c(n), p = c(p))
nsim = 1
A_coeff = A_coeff_banded(n_p)
data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "t", innov_df = 4)
A = A_coeff$`(200,20)`

data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = renyi_coeff, innov_dist = "t", innov_df = 4)
A = attributes(data$`(200,20)`)$A[[1]]
data = data$`(200,50)`[[1]]
data = data$`(200,20)`[[1]]
data_tr = cross_val_and_trunc(data, n_tau = 60, cv_lag = T, lag = 0, trim = 1, trim_d=1, trim_all=F, max = T, standardise = F)
for(i in 1:50){
  plot(data_tr[,i], type = "l", ylim = c(min(data), max(data)), main = paste(i))
  plot(data[,i], type = "l", ylim = c(min(data), max(data)), main = paste(i))
}

set.seed(200)
spvar = sparsevar::fitVAR(data = data, p =1, method = "cv")
set.seed(200)
glmn = glmnet::cv.glmnet(y = as.vector(data[n:2,]), x = kronecker(diag(1, p), data[(n-1):1,]), intercept = T, standardize = FALSE)
glmn_A = matrix(coef(glmn, s = glmn$lambda.min)[-1], nrow = p, ncol = p, byrow = T)

spvar$A[[1]]
spvar$A[[1]][1,1]
glmn_A
glmn_A[1,1]

norm(glmn_A - A, "F")
norm(spvar$A[[1]] - A, "F")
heatmap(glmn_A, Rowv = NA, Colv = NA)
heatmap(spvar$A[[1]], Rowv = NA, Colv = NA)
heatmap(attr(data), Rowv = NA, Colv = NA)

as.matrix(as.vector(m), p, p)
matrix(as.vector(m), nrow = 3, ncol = 3, byrow = T)
as.vector(m)
