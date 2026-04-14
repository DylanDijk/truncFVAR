##########################################################################################################################################
# This script generates the table containing forecasting performance results 
  # - using the results generated from the forecasting_store_all_truncV1.R script
##########################################################################################################################################


####### trunc V1 results #############

res_t21 = readRDS("revision_scripts/comment_5/no_trunc_Xn//fac_var_forecasting_t21_truncV1.Rds")
res_t3 = readRDS("revision_scripts/comment_5/no_trunc_Xn/fac_var_forecasting_t3_truncV1.Rds")
res_t4 = readRDS("revision_scripts/comment_5/no_trunc_Xn/fac_var_forecasting_t4_truncV1.Rds")
res_log = readRDS("revision_scripts/comment_5/no_trunc_Xn/fac_var_forecasting_lognormal_truncV1.Rds")
res_gauss = readRDS("revision_scripts/comment_5/no_trunc_Xn/fac_var_forecasting_Gaussian_truncV1.Rds")

res = res_t21

###################################################################################################################################


layers1 <- c("truncated", "not_truncated")

for (l1 in layers1) {
  res[[l1]][["t_3"]]       <- res_t3[[l1]][["t_3"]]
  res[[l1]][["t_4"]]       <- res_t4[[l1]][["t_4"]]
  res[[l1]][["log_normal"]]     <- res_log[[l1]][["log_normal"]]
  res[[l1]][["gaussian"]]       <- res_gauss[[l1]][["gaussian"]]
}

res_error <- lapply(res, function(trunc_level) {
  lapply(trunc_level, function(dist_level) {
    lapply(dist_level, function(np_level) {
      sapply(np_level, function(x) {
        max(abs(x$pred - x$oracle)) / max(abs(x$oracle))
        # max(abs(x$pred - x$oracle))
        # sum((x$pred - x$oracle)^2)/ sum((x$oracle)^2)
      })
    })
  })
})


n_p <- c("(100,50)", "(100,100)", "(200,50)", "(200,100)", "(500,100)", "(500,200)")
dist <- c("log_normal","t_2.1","t_3","t_4","gaussian")
# dist <- c("t_2.1","t_3","t_4")
methods <- c("truncated", "not_truncated")


##################################################################################################################################################
# not including t4 distribution column
dist <- dist[dist != "t_4"]
matrix_tab <- matrix(
  "",
  nrow = length(n_p) * length(methods),
  ncol = length(dist) * 3
)

row_labels <- c()
row_idx <- 1

for (k in seq_along(n_p)) {
  
  k_n_p <- n_p[k]
  
  # store values temporarily for comparison
  temp_means <- matrix(NA, nrow = 2, ncol = length(dist))
  temp_max   <- matrix(NA, nrow = 2, ncol = length(dist))
  temp_sds   <- matrix(NA, nrow = 2, ncol = length(dist))
  
  # first pass: compute values
  for (m_idx in seq_along(methods)) {
    m <- methods[m_idx]
    
    for (i in seq_along(dist)) {
      i_dist <- dist[i]
      
      vals <- res_error[[m]][[i_dist]][[k_n_p]]
      
      temp_means[m_idx, i] <- mean(vals)
      temp_max[m_idx, i]   <- max(vals)
      temp_sds[m_idx, i]   <- sd(vals)
    }
  }
  
  # second pass: format + bold
  for (m_idx in seq_along(methods)) {
    
    m <- methods[m_idx]
    
    row_labels <- c(
      row_labels,
      paste0(k_n_p, " ", ifelse(m == "truncated", "Truncated", "Untruncated"))
    )
    
    for (i in seq_along(dist)) {
      
      mean_val <- temp_means[m_idx, i]
      max_val  <- temp_max[m_idx, i]
      sd_val   <- temp_sds[m_idx, i]
      
      other_idx <- ifelse(m_idx == 1, 2, 1)
      
      # ---- MEAN ----
      if (!is.na(mean_val) && mean_val < temp_means[other_idx, i]) {
        mean_str <- paste0("\\textbf{", sprintf("%.2f", mean_val), "}")
      } else {
        mean_str <- sprintf("%.2f", mean_val)
      }
      
      # ---- MAX ----
      if (!is.na(max_val) && max_val < temp_max[other_idx, i]) {
        max_str <- paste0("\\textbf{", sprintf("%.2f", max_val), "}")
      } else {
        max_str <- sprintf("%.2f", max_val)
      }
      
      # ---- SD (in brackets) ----
      sd_str <- paste0("(", sprintf("%.2f", sd_val), ")")
      
      matrix_tab[row_idx, 3*(i-1) + 1] <- mean_str
      matrix_tab[row_idx, 3*(i-1) + 2] <- max_str
      matrix_tab[row_idx, 3*(i-1) + 3] <- sd_str
    }
    
    row_idx <- row_idx + 1
  }
}

rownames(matrix_tab) <- row_labels

colnames(matrix_tab) <- as.vector(rbind(
  paste0(dist, "_Mean"),
  paste0(dist, "_Max"),
  paste0(dist, "_SD")
))

library(xtable)
print(
  xtable(matrix_tab),
  sanitize.text.function = identity,
  sanitize.rownames.function = identity,
  sanitize.colnames.function = identity
)
