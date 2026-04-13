# This script uses the fnets package, of the version on github with commit 89da3c3 made on Nov 20, 2024
# To install:
# devtools::install_github("haeran-cho/fnets@89da3c3")
library(fnets)

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


#### Load data. ####
# description of variables: https://users.ssc.wisc.edu/~bhansen/econometrics/FRED-MD_description.pdf
set.seed(500)
library(fbi) # version 0.7.0
# filepath <- "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2023-12.csv"
setwd("/user/work/zl22291/simulations/truncation/proj1_revisions/redo_of_Real_data_forecasting")
# setwd("~/Documents/Projects/PhD/Journal_submission_project_1/Journal_econometrics_and_statistics/Revisions/revision_scripts/redo_of_Real_data_forecasting")
filepath <- "FREDMD-2023-12.csv"
data <- fredmd(filepath, date_start = NULL, date_end = NULL, transform = TRUE)
data = data[13:nrow(data),] # starting dataset from 1960-01-01
dates = data$date
years = as.numeric(substr(dates, 1, 4))
data <- data[, colSums(is.na(data)) == 0]
ncol(data)
data = data[,-1]
#####

#### Stationarity transformations. ####
# After fbi transformations there are still series that look like there is trend

## based on adf.test
library(tseries)
ind_diff = vector(length = ncol(data))
for(i in 1:ncol(data)){
  if(adf.test(data[,i])$p.value > 0.01){
    ind_diff[i] = i
  }
}
ind_diff = ind_diff[ind_diff != 0]
diff_data = data[-1,]
diff_data[,ind_diff] = diff(as.matrix(data[,ind_diff]))
data = diff_data

#####

#### Center and scale. ####
center_d = function(x){
  mean.x <- apply(x, 2, mean) 
  x <- t(x) - mean.x
  x = t(x)
  return(x)
}
scale_d = function(x){
  sd.x <- apply(x, 2, sd) 
  x <- t(x)/sd.x
  x = t(x)
  return(x)
}

data = center_d(data)
data = scale_d(data)

#####

#### Get number of factors, based on truncated data. ####
data_tr = fnets::cv_trunc(data = data, cv_lag = 1)

# https://cran.r-project.org/web/packages/dfms/vignettes/introduction.html
ic = dfms::ICr(data_tr$data, max.r =20)
q = ic$r.star[2]
# plot(ic)
# screeplot(ic)
#####
# Then apply rolling forecasting procedure, re-centering each window, and allowing truncation parameter to be scaled.

# function that returns forecast after truncating and not truncating
library(Matrix)
forcast_tr_ntr = function(tr = FALSE, data, q = c("ic", "er"), h = h){
  if(tr){
    fnet = fnets(data, fm.restricted = TRUE, robust = TRUE, robust.standardise = TRUE, var.args = list(n.cores = 4), center = FALSE, q = q)
  } else {
    fnet = fnets(data, fm.restricted = TRUE, robust = FALSE, var.args = list(n.cores = 4), center = FALSE, q = q)
  }
  forcast_list = vector(mode = "list", length = 8)
  forcast_list <- lapply(forcast_list, function(x) NA)
  names(forcast_list) = c("combined_pred", "factor_pred", "idio_pred", "number_of_factors", "loadings", "idio_beta","lambda", "truncation")
  
  pred = predict_fnets_mod(fnet, n.ahead = h, fc.restricted = TRUE, data = data)
  
  forcast_list[[1]] = pred$forecast
  forcast_list[[2]] = pred$common.pred$fc
  forcast_list[[3]] = pred$idio.pred$fc
  forcast_list[[4]] = fnet$q
  if(is.null(fnet$loadings)){
    forcast_list[[5]] = NA
  } else {
    forcast_list[[5]] = fnet$loadings
  }
  # forcast_list[[5]] = ifelse(is.null(fnet$loadings), NA, fnet$loadings)
  forcast_list[[6]] = Matrix(fnet$idio.var$beta, sparse = TRUE)
  forcast_list[[7]] = fnet$idio.var$lambda
  if(tr){
    forcast_list[[8]] = attributes(fnet)$truncation
  }
  return(forcast_list)
}

# rolling forecast 
ws = 120
h = 3
forec_results = vector(mode = "list", length = 3)
names(forec_results) = c("tr", "ntr", "true")

forec_results$tr = vector(mode = "list", length = nrow(data)-ws - (h-1))
forec_results$ntr = vector(mode = "list", length = nrow(data)-ws - (h-1))
forec_results$true = vector(mode = "list", length = nrow(data)-ws - (h-1))

# system.time({
for(i in 1:(nrow(data)- ws - (h-1))){
  # for(i in c(1)){
  dat = data[i:((ws)+i+h-1),]
  dat = center_d(dat)
  dat = as.matrix(dat)
  
  # debug(forcast_tr_ntr)
  forec_results$ntr[[i]] = forcast_tr_ntr(tr = FALSE, data = dat[1:(ws),], h = h, q = q)
  forec_results$tr[[i]] = forcast_tr_ntr(tr = TRUE, data = dat[1:(ws),], h = h, q = q)
  forec_results$true[[i]] = dat[(ws+1):(ws+3),]
}
# })

saveRDS(forec_results, file = "realdata_forecasting_truncV1.Rds")
