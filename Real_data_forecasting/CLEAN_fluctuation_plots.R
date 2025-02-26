set.seed(500)
library(fbi) # version 0.7.0
# description of variables: https://users.ssc.wisc.edu/~bhansen/econometrics/FRED-MD_description.pdf
filepath <- "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2023-12.csv"
data <- fredmd(filepath, date_start = NULL, date_end = NULL, transform = TRUE)
data = data[13:nrow(data),] # starting dataset from 1960-01-01
dates = data$date
years = as.numeric(substr(dates, 1, 4))
data <- data[, colSums(is.na(data)) == 0]
ncol(data)
data = data[,-1]
years = years[-1]
###################################################################################

forec_errors = readRDS(file = "Real_data_forecasting/FRED_MD_CLEAN_forcast_errors_diff_ADF.Rds")

nsim = length(forec_errors$tr)
h = 1

library(murphydiagram)
for(j in c(66,67,21,57,64,73)){

  forec_errors_tr = vector(length = nsim)
  for(i in 1:nsim){
    forec_errors_tr[i] = abs((forec_errors$tr[[i]])[["combined_pred"]][h,j] - (forec_errors$true[[i]])[h,j])
  }
  forec_errors_ntr = vector(length = nsim)
  for(i in 1:nsim){
    forec_errors_ntr[i] = abs((forec_errors$ntr[[i]])[["combined_pred"]][h,j] - (forec_errors$true[[i]])[h,j])
  }
  
  fluctuation_test(forec_errors_tr, forec_errors_ntr, mu = 0.3, dmv_fullsample = T,
                   lag_truncate =0,
                   conf_level = 0.1, time_labels = years[121:764])
  title(paste(colnames(data)[j]))
}

