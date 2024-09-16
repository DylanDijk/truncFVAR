# truncFVAR

Code used for "Heavy-tailed robust estimation of factor-adjusted vector autoregressive models for high-dimensional time series"

This repository contains all of the code used to generate the figures and tables for the simulation results.

## Example
Below we give code to truncate simulated data, and then comparing the errors of covarianc estimation and Lasso estimation.

#### sourcing functions
```r
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/master/functions/data_generation.R")
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/master/functions/estimation.R")
```
#### generate VAR data
```r
n_p = cbind(n = 200, p = 50)
A = A_coeff_banded(n_p)[[1]]
VAR_data = VAR_1_data_ind(nsim = 1, n_p = n_p, innov_dist = "t", innov_df = 2.1, A_coeff = A)$`(200,50)`[[1]]
```
#### truncate
```r
trunc_VAR_data = cross_val_and_trunc(data = VAR_data, cv_lag = T, standardise = F)
```
#### plot of original series and truncated series
```r
for(i in 1:n_p[,"p"]){
  plot(trunc_VAR_data[,i], type = "l", ylim = c(min(VAR_data), max(VAR_data)), main = paste(i))
  plot(VAR_data[,i], type = "l", ylim = c(min(VAR_data), max(VAR_data)), main = paste(i))
}
```
#### cov estimation
```r
true_cov = cov_of_var(A = A)
norm(acf_no_center(VAR_data)-true_cov, "M") 
norm(acf_no_center(trunc_VAR_data)-true_cov, "M")
```
#### lasso
```r
norm(sparsevar::fitVAR(VAR_data, p = 1, parallel = T, ncores = 4)$A[[1]] - A, "M")
norm(sparsevar::fitVAR(trunc_VAR_data, p = 1, parallel = T, ncores = 4)$A[[1]] - A, "M")
```


## Corresponding code to simulation results




