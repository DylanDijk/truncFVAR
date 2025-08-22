# This script uses the fnets package, of the version on github with commit 89da3c3 made on Nov 20, 2024
# To install:
# devtools::install_github("haeran-cho/fnets@89da3c3")
library(fnets)


#### Load data. ####
# description of variables: https://users.ssc.wisc.edu/~bhansen/econometrics/FRED-MD_description.pdf
set.seed(500)
library(fbi) # version 0.7.0
# filepath <- "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2023-12.csv"
filepath <- "Real_data_forecasting/FREDMD-2023-12.csv"
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
  
  pred = predict(fnet, n.ahead = h, fc.restricted = TRUE)
  
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

for(i in 1:(nrow(data)- ws - (h-1))){
  # for(i in 1:2){
  dat = data[i:((ws)+i+h-1),]
  dat = center_d(dat)
  dat = as.matrix(dat)

  forec_results$ntr[[i]] = forcast_tr_ntr(tr = FALSE, data = dat[1:(ws),], h = h, q = q)
  forec_results$tr[[i]] = forcast_tr_ntr(tr = TRUE, data = dat[1:(ws),], h = h, q = q)
  forec_results$true[[i]] = dat[(ws+1):(ws+3),]
}

# saveRDS(forec_results, file = "forecast_real_data/FRED_MD_CLEAN_forcast_errors_diff_ADF.Rds")









