# In papers DGP notation: (F0),(V1),(S1)
# Figure 1 + 2

# Creates plots that display how choice of truncation parameter tau varies under different settings.
# The plots are created using dependent samples, generated from VAR(1) model with banded coefficient matrix, innovation errors with diagonal covariance structure. 

# These simulations standardised the data, by scaling up the truncation parameter tau by the MAD of the corresponding variable.
# This makes little difference, as the data has been generated from a VAR(1) model.


# In order to generate all the plots in Figures 1 + 2, requires modifying the script.
# To generate the plots for both against dimension or sample size, requires uncommenting the # vs n block code or # vs p.
# To generate the plots for the choice of tau with the autocovariance in the CV measure (Figure 2) requires setting cv_lag = T,
# in the compute_errors() function.


#### sourcing functions ####
source(file = "truncation/functions/data_generation.R")
source(file = "truncation/functions/estimation.R")


list_setup = list(
  "t_2.1" = list(),
  "t_3" = list(),
  "t_4" = list(),
  "gaussian" = list()
)

# vs p
n_p = cbind(n = c(200), p = seq(50, 300, 50))
nsim = 200

# vs n
# n_p = cbind(n = seq(50, 300, 50), p = 100)
# nsim = 200

# Banded VAR coefficient matrix
A_coeff = A_coeff_banded(n_p)

compute_errors = function(data, distr, table, max){
  
  for(i in names(data)){
    table[[distr]][[i]] = 
      cross_val_and_error(data = data[[i]], 
                          n_tau = 60, lag = 0, trim = 1, trim_d = 1, trim_all = F, max = max, cv_lag = T)
  }  
  
  return(table)
}



tau_vs_p = list_setup
# tau_vs_n = list_setup

# gaussian
set.seed(200)
tau_vs_p = compute_errors(data = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "gauss"),
                          table = tau_vs_p, distr = "gaussian", max = T)
# t
library(stringr)
for(i in c("t_2.1", "t_3", "t_4")){
  set.seed(200)
  tau_vs_p = compute_errors(data =  VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_dist = "t",
                                                   innov_df = as.numeric(str_extract_all(i, "\\d+\\.\\d+|\\d+"))),
                            table = tau_vs_p, distr = i, max = T)
  
}

################## line graphs ################## 

dat = tau_vs_n
x = log(seq(from = 50, to = 300, by = 50))
xlab = "log(n)"
margin = 0.1
ylab=""

min_taus_t_2.1 = lapply(dat$t_2.1, function(x){log(mean(as.numeric(x)))})
min_quart_t_2.1 = lapply(dat$t_2.1, function(x){round(mean(as.numeric(names(x))),3)})
min_taus_gauss = lapply(dat$gaussian, function(x){log(mean(as.numeric(x)))})
min_quart_gauss = lapply(dat$gaussian, function(x){round(mean(as.numeric(names(x))),3)})

# ylim = c(min(unlist(min_taus_gauss)) - 0.1,max(unlist(min_taus_t_2.1)) + 0.1),

par(mar=c(5,5,4,1)+.1)
plot(y = unlist(min_taus_t_2.1), x = x, xaxt='n',
     ylab = ylab, xlab = xlab, type = "b",
     main = "Sample size", pch = 16,  col = "purple", 
     ylim = c(0,2),
     cex = 1.5, cex.lab = 2, cex.axis = 1.45, cex.main =2.3,
     xlim = c(min(x)-margin, max(x)+margin))
text(y = unlist(min_taus_t_2.1) + 0.05, x = x, as.character(min_quart_t_2.1), cex = 1.5)
for(i in 1:3){
  dist_i = c("t_3", "t_4", "gaussian")[i]
  colr = c("blue", "red", "green")[i]
  min_taus = unlist(lapply(dat[[dist_i]], function(x){log(mean(as.numeric(x)))}))
  min_quart = lapply(dat[[dist_i]], function(x){round(mean(as.numeric(names(x))),3)})
  points(x = x, y = min_taus, col = colr, type = "b", pch = 16, cex = 1.5)
  text(y = min_taus + 0.1, x = x, as.character(min_quart), cex = 1.5)
}
axis(1,at=x, labels = as.character(seq(from = 50, to = 300, by = 50)), cex.axis = 1.45)

margin = 3
legend("top", legend=c("T-2.1", "T-3", "T-4", "Gaussian"),
       col=c("purple", "blue", "red", "green"), pch=rep(16,4), cex=2,  bty = "n", horiz = T, text.width = 0, x.intersp = 0.2)


################## Box-plots ##################    
dat = tau_vs_n
dat = tau_vs_p

dat = lapply(dat, function(x){lapply(x, log)})
library(reshape2)
long_data <- melt(dat, id.vars = c("group", "type"), variable.name = "element", value.name = "value")
long_data$L1 <- factor(long_data$L1, levels = c("gaussian", "t_4", "t_3", "t_2.1"))

long_data$L2 <- factor(long_data$L2, levels = c("(50,100)", "(100,100)", "(150,100)", "(200,100)", "(250,100)", "(300,100)"))
long_data$L2 <- factor(long_data$L2, levels = c("(200,50)", "(200,100)", "(200,150)", "(200,200)", "(200,250)", "(200,300)"))

# Create the boxplot using ggplot2
library(ggplot2)
ggplot(long_data, aes(x = L2, y = value, fill = L1)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  scale_fill_manual(values =  c("green", "red", "blue", "purple")) + 
  scale_x_discrete(
    labels = c("(50,100)" = "50", "(100,100)" = "100", "(150,100)" = "150",
               "(200,100)" = "200", "(250,100)" = "250", "(300,100)" = "300")
  ) +
  # scale_x_discrete(
  #   labels = c("(200,50)" = "50", "(200,100)" = "100", "(200,150)"  = "150", "(200,200)"  = "200", "(200,250)"= "250", "(200,300)"= "300")
  # ) + 
  labs(x = "log(n)", y = NULL) +
  # labs(x = "log(2log(p))", y = NULL) +
  coord_cartesian(ylim = c(NA, 2.8)) +  # Set y-axis limit to 15
  theme_minimal() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)
  )

################################################################################################################################################ 






