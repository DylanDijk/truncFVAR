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

data_tr = fnets::cv_trunc(data = data, cv_lag = 1)


#################### Now want to get eigenvalues of covariance matrix
library(Matrix)   
library(ggplot2)


# Parameters
n_rep <- 100                  # number of resamples per dimension
p_seq <- 5:108                 # dimensions to try
n <- nrow(data)                # time length

get_eigs <- function(X) {
  Sigma_hat <- cov(X)
  ev <- eigen(Sigma_hat, symmetric = TRUE, only.values = TRUE)$values
  # sort(ev, decreasing = TRUE)[1:2]
  sort(ev, decreasing = TRUE)[1]
}

# Storage
results <- data.frame()

# Loop over cross-sectional dimensions
set.seed(3)
for (p in p_seq) {
  print(p)
  for (rep in 1:n_rep) {
    cols <- sample(ncol(data), p)         
    X_sub <- data[, cols]
    eigs <- get_eigs(X_sub)
    # results <- rbind(results,
    #                  data.frame(p = p, eigenvalue = eigs, rank = 1:2))
    results <- rbind(results,
                     data.frame(p = p, eigenvalue = eigs, rank = 1))
  }
}

# Plot (similar to your figure)

ggplot(results, aes(x = factor(p), y = eigenvalue, fill = factor(rank))) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, size = 0.3, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("white", "grey70"), name = "Eigenvalue rank",
                    labels = c("1st largest", "2nd largest")) +
  scale_x_discrete(
    breaks = levels(factor(results$p))[seq(1, length(unique(results$p)), by = 5)]
  ) +
  labs(x = "Dimension", y = "Eigenvalue") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")   # remove legend if you want same style as paper


# gplot = 
# ggplot(results, aes(x = factor(p), y = eigenvalue, fill = factor(rank))) +
#   geom_boxplot(outlier.size = 0.5, width = 0.6, size = 0.3, position = position_dodge(width = 0.7)) +
#   scale_fill_manual(values = c("white", "grey70"), name = "Eigenvalue rank",
#                     labels = c("1st largest", "2nd largest")) +
#   scale_x_discrete(
#     breaks = levels(factor(results$p))[seq(1, length(unique(results$p)), by = 5)]
#   ) +
#   labs(x = "Dimension", y = "Eigenvalue") +
#   theme_bw(base_size = 14) +
#   theme(legend.position = "none")   # remove legend if you want same style as paper
# 
# gdat = ggplot_build(gplot)$data[[1]]
# 
# gplot + geom_segment(data=gdat, aes(x=xmin, xend=xmax, 
#                                y=middle, yend=middle), colour="red", size=1.1, inherit.aes = FALSE)



