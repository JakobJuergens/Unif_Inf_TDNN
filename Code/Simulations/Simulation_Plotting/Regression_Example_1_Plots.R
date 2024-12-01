##### Setup #####
# Load packages from CRAN
library(tidyverse)
library(patchwork)
library(metR)

# Set up regression function and heteroskedasticity function
reg_f <- function(covariates) {
  # return(1)
  return(sum(cos(5 * covariates)))
}

htscty_f <- function(covariates) {
  # return(0)
  return(0.25 * abs(sum(covariates^2)))
}

##### Generate Plots for Regression Function and Heteroskedasticity #####
# Create data-frame for plotting purposes
plot_df <- as_tibble(expand.grid(seq(0, 1, 0.01), seq(0, 1, 0.01)))
plot_df$reg_vals <- unlist(purrr::map(.f = ~ reg_f(plot_df[.x, 1:2]), .x = 1:nrow(plot_df)))
plot_df$var_vals <- unlist(purrr::map(.f = ~ htscty_f(plot_df[.x, 1:2]), .x = 1:nrow(plot_df)))

# Plot response surface
reg_plot <- ggplot(data = plot_df) +
  geom_contour_fill(aes(x = Var1, y = Var2, z = reg_vals),
    bins = 100
  ) +
  scale_fill_distiller(palette = "Spectral") +
  xlab("Covariate 1") +
  ylab("Covariate 2") +
  labs(fill = "Value of Regression Function") +
  theme_light() +
  theme(
    legend.position = "bottom", legend.key.size = unit(1.5, "cm"),
    text = element_text(size = 16)
  )

# Plot Heteroskedasticity
var_plot <- ggplot(data = plot_df) +
  geom_contour_fill(aes(x = Var1, y = Var2, z = var_vals),
    bins = 100
  ) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  xlab("Covariate 1") +
  ylab(element_blank()) +
  labs(fill = "Error Term Variance") +
  theme_light() +
  theme(
    legend.position = "bottom", legend.key.size = unit(1.5, "cm"),
    text = element_text(size = 16)
  )

# Combine Plots for Paper
comb_plot <- reg_plot + var_plot
ggsave(
  filename = "../../Graphics/Reg_Exmp1.pdf", plot = comb_plot,
  width = 16, height = 6, units = "in"
)

##### Generate Plots for Estimator #####
# Load Estimation Results
results <- list(
  readRDS(file = "Simulation_Results/reg_exp_1_n10000s_10_25.RDS"),
  readRDS(file = "Simulation_Results/reg_exp_1_n10000s_100_250.RDS"),
  readRDS(file = "Simulation_Results/reg_exp_1_n10000s_500_5000.RDS")
)

# Create Tibble and matrix versions of Estimation results
indices <- 1:length(results)
n_reps <- results[[1]]$n_reps
points <- results[[1]]$points
values <- results[[1]]$values

estimates_mat <- purrr::map(.x = results, .f = ~ do.call(cbind, .x$estimates))
estimates_tibbles <- purrr::map(.x = indices, .f = ~ cbind(tibble(
  cov1 = points[, 1], cov2 = points[, 2],
  vals = values
), estimates_mat[[.x]]) %>%
  pivot_longer(cols = !c("cov1", "cov2", "vals"), names_to = "Simulation Run", values_to = "Estimate"))

# Calculate Statistics for each point
estimates_avgs <- purrr::map(.x = estimates_mat, .f = ~ rowMeans(.x))
estimates_bias <- purrr::map(.x = indices, .f = ~ estimates_avgs[[.x]] - values)
estimates_vars <- purrr::map(
  .x = indices,
  .f = function(x) unlist(purrr::map(1:nrow(points), .f = function(y) var(estimates_mat[[x]][y, ])))
)

# Create Tibble with Statistics
stat_tibble <- purrr::map(
  .x = indices,
  .f = ~ tibble(
    cov1 = points[, 1], cov2 = points[, 2],
    avgs = estimates_avgs[[.x]], bias = estimates_bias[[.x]], vars = estimates_vars[[.x]]
  )
)

# Calculate Scales for Fill-Legends
limits_bias_u <- unlist(purrr::map(.x = indices, .f = ~ max(stat_tibble[[.x]]$bias)))
limits_bias_l <- unlist(purrr::map(.x = indices, .f = ~ min(stat_tibble[[.x]]$bias)))
limits_bias <- unlist(purrr::map(.x = indices, .f = ~ max(abs(c(limits_bias_l[.x], limits_bias_u[.x])))))

limits_vars <- unlist(purrr::map(.x = indices, .f = ~ max(stat_tibble[[.x]]$vars)))

bias_colors <- list(c("green", "red"), c("blue", "orange"), c("purple", "yellow"))

# Create Plots of Statistics
# Plot Bias Surface
est_bias_plots <- purrr::map(
  .x = indices,
  .f = ~ ggplot(data = stat_tibble[[.x]], aes(x = cov1, y = cov2, z = bias)) +
    geom_contour_fill(bins = 100) +
    geom_contour(color = "black", size = 0.1, bins = 20) +
    # scale_fill_distiller(palette = "Spectral", limit = c(-1,1)*limits_bias[.x]) +
    scale_fill_gradient2(low = bias_colors[[.x]][1], high = bias_colors[[.x]][2], midpoint = 0, limit = c(-1, 1) * limits_bias[.x]) +
    xlab("Covariate 1") +
    ylab("Covariate 2") +
    labs(fill = "Estimator Bias") +
    theme_light() +
    theme(
      legend.position = "bottom", legend.key.size = unit(1.5, "cm"),
      text = element_text(size = 16)
    )
)

# Plot Estimator Variance Surface
est_var_plots <- purrr::map(
  .x = indices,
  .f = ~ ggplot(data = stat_tibble[[.x]], aes(x = cov1, y = cov2, z = vars)) +
    geom_contour_fill(bins = 100) +
    geom_contour(color = "black", size = 0.1, bins = 20) +
    # scale_fill_distiller(palette = "RdBu", direction = -1) +
    scale_fill_gradient2(low = "white", high = bias_colors[[.x]][2], limit = c(0, limits_vars[.x])) +
    xlab("Covariate 1") +
    ylab(element_blank()) +
    labs(fill = "Estimator Variance") +
    theme_light() +
    theme(
      legend.position = "bottom", legend.key.size = unit(1.5, "cm"),
      text = element_text(size = 16)
    )
)

# Combine Plots for Paper
comb_plots <- purrr::map(
  .x = indices,
  .f = ~ est_bias_plots[[.x]] + est_var_plots[[.x]]
)
ggsave(
  filename = paste0("../../Graphics/Reg_Exmp1_", results$kernel_orders[1], "_", results$kernel_orders[2], "_Est.pdf"),
  plot = comb_plot2, width = 16, height = 6, units = "in"
)
