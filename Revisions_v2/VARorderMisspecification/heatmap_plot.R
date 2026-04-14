##########################################################################################################################################
# Generates heatmap of average estimate truncated mean estimate when the VAR order used for fitting is 3, but the true VAR order is 2.
  # The true VAR coefficient matrix is banded, so is the same across all the simulations.
##########################################################################################################################################


########################### FINAL clean plots for the response letter #####################################################

res = readRDS("revision_scripts/comment_4/heatmap_results__d0_2_d_3.rds")
# res is a list of 200 matrices

#### First need to subtract from the fitted values in res, the true VAR coefficient matrix. ####
# True VAR coefficient matrix. A_coeff_ex.
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/3c4339ea8bf00a735408be0659d0c88bcf5dac0e/functions/data_generation.R")
source("https://raw.githubusercontent.com/DylanDijk/truncFVAR/84f22825f88656a6b6efdb8cf628d1e7610397b0/functions/estimation.R")
d_0 = 2
# what lags to fit 
d = 3

n_p = cbind(n = c(200), p = c(50))

A_coeff = A_coeff_banded_d(n_p, d = d_0)
# now adding zero matrix to side of var coefficient matrices.
if(length(A_coeff[[1]]) < d){
  A_coeff_ex <- lapply(A_coeff, function(A_list) {
    
    # Bind all matrices inside the list horizontally
    A_combined <- do.call(cbind, A_list)
    
    p <- nrow(A_combined)
    
    # Append zero block
    cbind(A_combined, matrix(0, nrow = p, ncol = p))
  })
}

dim(A_coeff_ex$`(200,50`)





library(ggplot2)
library(reshape2)

# Extract full true matrix
true_mat <- A_coeff_ex$`(200,50`

# Symmetric limits (recommended for comparison with estimates)
max_abs <- max(abs(true_mat))
lims <- c(-max_abs, max_abs)

# Convert to long format
df_true <- melt(true_mat)
colnames(df_true) <- c("row", "col", "val")


common_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    
    # ↓ Slightly smaller fonts
    plot.title = element_text(size = 35, hjust = 0.5),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    
    # ↓ Softer legend (not bold)
    legend.title = element_text(size = 24, face = "plain"),
    legend.text  = element_text(size = 20),
    
    # ↓ THIS IS KEY: reduce margins so panel is bigger
    plot.margin = margin(5, 5, 5, 5)
  )



# Plot heatmap
p_true <- ggplot(df_true, aes(x = col, y = row, fill = val)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    midpoint = 0, limits = lims
  ) +
  labs(
    title = "True Coefficient Matrix",
    x = "Columns",
    y = "Rows",
    fill = "Value"
  ) +
  coord_fixed() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  common_theme +
  guides(
    fill = guide_colorbar(
      barheight = 12,
      barwidth = 1,
      frame.colour = "black",
      ticks.colour = "black"
    )
  )
p_true



library(ggplot2)
library(reshape2)

# Function to compute mean estimate
compute_mean <- function(mat_list) {
  arr <- simplify2array(mat_list)
  apply(arr, c(1, 2), mean)
}


# Compute mean matrices
mean_trunc <- compute_mean(res$truncated$t_2.1$`(200,50`)
mean_not   <- compute_mean(res$not_truncated$t_2.1$`(200,50`)


# Convert to long format
df_trunc <- melt(mean_trunc)
df_not   <- melt(mean_not)

colnames(df_trunc) <- c("row", "col", "val")
colnames(df_not)   <- c("row", "col", "val")

# Plot: Truncated
p1 <- ggplot(df_trunc, aes(x = col, y = row, fill = val)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    midpoint = 0, limits = lims
  ) +
  labs(
    title = "Truncated (Mean Estimate)",
    x = "Columns",
    y = "Rows",
    fill = "Value"
  ) +
  coord_fixed() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  common_theme +
  guides(
    fill = guide_colorbar(
      barheight = 12,
      barwidth = 1,
      frame.colour = "black",
      ticks.colour = "black"
    )
  )

# Plot: Not Truncated
p2 <- ggplot(df_not, aes(x = col, y = row, fill = val)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    midpoint = 0, limits = lims
  ) +
  labs(
    title = "Not Truncated (Mean Estimate)",
    x = "Columns",
    y = "Rows",
    fill = "Value"
  ) +
  coord_fixed() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    plot.title = element_text(
      hjust = 0.5,
      size = 25,
      margin = margin(b = 10)
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 14)
  ) +
  guides(
    fill = guide_colorbar(
      barheight = 12,
      barwidth = 1,
      frame.colour = "black",
      ticks.colour = "black"
    )
  )

# Display plots
p1