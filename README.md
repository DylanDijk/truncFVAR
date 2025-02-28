# truncFVAR

Code used for the paper: "Heavy-tailed robust estimation of factor-adjusted vector autoregressive models for high-dimensional time series"

This repository contains the code used to generate the figures and tables for the simulation results.

## Example
<details>
  <summary>(Click to expand) Here we give an example of using the functions used for the simulations results, to truncate simulated data, and then comparing the errors of covariance estimation and Lasso estimation.</summary>

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


## Corresponding code to simulation results




