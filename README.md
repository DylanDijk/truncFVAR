# truncFVAR

This repository contains the code used to generate the figures and tables for the simulation results, for the paper:  
"Heavy-tailed robust estimation of factor-adjusted vector autoregressive models for high-dimensional time series"


## Example
<details>
  <summary>Example of using functions from this repository, to truncate simulated data, and then compare errors from covariance estimation and Lasso estimation, to estimation without truncation.
</summary>

#### sourcing functions
```r
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/master/functions/data_generation.R")
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/master/functions/estimation.R")
```
#### generate VAR data
```r
n_p = cbind(n = 200, p = 50)
A = A_coeff_banded(n_p)
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
A = A[[1]]
true_cov = cov_of_var(A = A)
norm(acf_no_center(VAR_data)-true_cov, "M") 
norm(acf_no_center(trunc_VAR_data)-true_cov, "M")
```
#### Lasso
```r
norm(sparsevar::fitVAR(VAR_data, p = 1, parallel = T, ncores = 4)$A[[1]] - A, "M")
norm(sparsevar::fitVAR(trunc_VAR_data, p = 1, parallel = T, ncores = 4)$A[[1]] - A, "M")
```
</details>

## fnets R package

The [fnets](https://github.com/haeran-cho/fnets) R package has been updated on
GitHub to include truncation, implemented as described
in this paper. The package was used to carry out the real data application in the paper.

To install the version used for the paper:
```r
devtools::install_github("haeran-cho/fnets@89da3c3")
```

<details>
  <summary>Covariance estimation example</summary>
  
  
This example looks at covariance estimation, with truncation, for data generated from a VAR model.

I source the `estimation.r` script from this repo, as it contains a function to calculate the true covariance matrix of a VAR process given the coefficient matrix A
```r
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/master/functions/estimation.R")
```
#### Generating data
```r
VAR_data = fnets::sim.var(n = 200, p = 50, heavy = TRUE, df = 2.1)
```
#### Truncating the data
```r
VAR_data_tr = fnets::cv_trunc(data = VAR_data$data, cv_lag = 1, standardise = FALSE)$data
```
#### Computing the true covariance of the VAR(1) process
```r
true_cov = cov_of_var(A = VAR_data$A)
```
#### Computing sample covariance for the truncated and original data
```r
cov_est = fnets:::acf_no_center(data = VAR_data$data, lag = 0)
cov_rob_est = fnets:::acf_no_center(data = VAR_data_tr, lag = 0)
```

#### Covariance max norm estimation error
```r
norm(cov_est - true_cov, "M")
norm(cov_rob_est - true_cov, "M")
```
</details>

<details>
  <summary>VAR estimation example</summary>

#### Generate heavy-tailed VAR data and estimating A
A is estimated with fnets with `robust = TRUE`

```r
VAR_data = fnets::sim.var(n = 200, p = 50, heavy = TRUE, df = 2.1)
fnet_fit = fnets::fnets(x = VAR_data$data, center = FALSE, q = 0, robust = TRUE, fm.restricted = TRUE)
```
#### Looking at model fit
Comparing robust estimate of A to the true A
```r
par(mfrow = c(1,2), mar = c(1,1,2,1))
zlim <- range(c(fnet_fit$idio.var$beta, VAR_data$A))
image(t(fnet_fit$idio.var$beta), col = heat.colors(10), axes = FALSE, zlim = zlim, main = "fnets estimate")
image((VAR_data$A), col = heat.colors(10), axes = FALSE, zlim = zlim, main = "Ground truth")
```

</details>


<details>
  <summary>Factor-adjusted VAR estimation example</summary>

#### Generate heavy-tailed factor plus VAR data 

This process is generated as described in (F2) in the simulations section of the paper
```r
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/master/functions/data_generation.R")
facvar_dat = fac_var_dat(n = 200, p = 50, r = 3, dist = "t", innov_df = 2.1)
```

#### Estimating A

```r
fnet_fit = fnets::fnets(x = facvar_dat, center = FALSE, q = 3, robust = TRUE, fm.restricted = TRUE)
```

#### Looking at model fit
Comparing robust estimate of A to the true A
```r
par(mfrow = c(1,2), mar = c(1,1,2,1))
A = attributes(facvar_dat)$A[[1]]
zlim <- range(c(A, fnet_fit$idio.var$beta))
image(t(fnet_fit$idio.var$beta), col = heat.colors(10), axes = FALSE, zlim = zlim, main = "fnets estimate")
image(A, col = heat.colors(10), axes = FALSE, zlim = zlim, main = "Ground truth")
```

</details>

## Corresponding code to simulation results

To add here, the scripts and corresponding figure and table numbers


