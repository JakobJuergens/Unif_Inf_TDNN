##### Setup #####
# Load packages from CRAN
library(tidyverse)
library(patchwork)
library(metR)

# List simulation output files
sim_files <- list.files(path = 'Simulation_Results/Reg_Exp1/')

# Load Estimation Results
Estimation_Results <- purrr::map(.x = 1:length(sim_files),
                                  .f = ~readRDS(file = paste0('Simulation_Results/Reg_Exp1/', sim_files[[.x]])))

DNN1_estimates <- purrr::map(.x = 1:length(sim_files),
                      .f = ~ Estimation_Results[[.x]]$DNN_estimates[[1]])

DNN2_estimates <- purrr::map(.x = 1:length(sim_files),
                             .f = ~ Estimation_Results[[.x]]$DNN_estimates[[2]])

TDNN_estimates <- purrr::map(.x = 1:length(sim_files),
                             .f = ~ Estimation_Results[[.x]]$TDNN_estimates)

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
  labs(fill = "Error Term Standard Deviation") +
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
# Create Tibble and matrix versions of Estimation results
indices <- 1:length(sim_files)
n_reps <- Estimation_Results[[1]]$n_reps
points <- Estimation_Results[[1]]$points
values <- Estimation_Results[[1]]$values

# Create Matrices of Estimates
DNN1_mat <- purrr::map(.x = indices,
                       .f = ~ matrix(data = do.call(cbind, DNN1_estimates[[.x]]), nrow = nrow(points)))
DNN2_mat <- purrr::map(.x = indices,
                       .f = ~ matrix(do.call(cbind, DNN2_estimates[[.x]]), nrow = nrow(points)))
TDNN_mat <- purrr::map(.x = indices,
                       .f = ~ matrix(do.call(cbind, TDNN_estimates[[.x]]), nrow = nrow(points)))

# Remove List of Data to save memory
rm(Estimation_Results, DNN1_estimates, DNN2_estimates, TDNN_estimates)
gc()

# Calculate Statistics for each point
DNN1_avgs <- purrr::map(.x = DNN1_mat, .f = ~ rowMeans(.x))
DNN1_bias <- purrr::map(.x = indices, .f = ~ DNN1_avgs[[.x]] - values)
DNN1_vars <- purrr::map(
  .x = indices, .f = ~unlist(purrr::map(1:nrow(points), .f = function(y) var(DNN1_mat[[.x]][y, ])))
)

DNN2_avgs <- purrr::map(.x = DNN2_mat, .f = ~ rowMeans(.x))
DNN2_bias <- purrr::map(.x = indices, .f = ~ DNN2_avgs[[.x]] - values)
DNN2_vars <- purrr::map(
  indices, .f = function(x) unlist(purrr::map(1:nrow(points), .f = function(y) var(DNN2_mat[[x]][y, ])))
)

TDNN_avgs <- purrr::map(.x = TDNN_mat, .f = ~ rowMeans(.x))
TDNN_bias <- purrr::map(.x = indices, .f = ~ TDNN_avgs[[.x]] - values)
TDNN_vars <- purrr::map(
  indices, .f = function(x) unlist(purrr::map(1:nrow(points), .f = function(y) var(TDNN_mat[[x]][y, ])))
)

# Clean Memory
rm(DNN1_mat, DNN2_mat, TDNN_mat)
gc()

# Create Tibble with Statistics
stat_tibble <- purrr::map(
  .x = indices,
  .f = ~ tibble(
    cov1 = points[, 1], cov2 = points[, 2],
    DNN1_avgs = DNN1_avgs[[.x]], DNN1_bias = DNN1_bias[[.x]], DNN1_vars = DNN1_vars[[.x]],
    DNN2_avgs = DNN2_avgs[[.x]], DNN2_bias = DNN2_bias[[.x]], DNN2_vars = DNN2_vars[[.x]],
    TDNN_avgs = TDNN_avgs[[.x]], TDNN_bias = TDNN_bias[[.x]], TDNN_vars = TDNN_vars[[.x]]
  )
)

# Create Plots of Statistics
# Plot Bias Surface
DNN1_bias_plots <- purrr::map(
  .x = indices,
  .f = ~ ggplot(data = stat_tibble[[.x]], aes(x = cov1, y = cov2, z = DNN1_bias)) +
    geom_contour_fill(bins = 100) +
    geom_contour(color = "black", linewidth = 0.1, bins = 20) +
    scale_fill_distiller(palette = "Spectral") +
    xlab("Covariate 1") +
    ylab("Covariate 2") +
    labs(fill = "Estimator Bias") +
    theme_light() +
    theme(
      legend.position = "bottom", legend.key.size = unit(1.5, "cm"),
      text = element_text(size = 16)
    )
)

DNN2_bias_plots <- purrr::map(
  .x = indices,
  .f = ~ ggplot(data = stat_tibble[[.x]], aes(x = cov1, y = cov2, z = DNN2_bias)) +
    geom_contour_fill(bins = 100) +
    geom_contour(color = "black", linewidth = 0.1, bins = 20) +
    scale_fill_distiller(palette = "Spectral") +
    xlab("Covariate 1") +
    ylab("Covariate 2") +
    labs(fill = "Estimator Bias") +
    theme_light() +
    theme(
      legend.position = "bottom", legend.key.size = unit(1.5, "cm"),
      text = element_text(size = 16)
    )
)

TDNN_bias_plots <- purrr::map(
  .x = indices,
  .f = ~ ggplot(data = stat_tibble[[.x]], aes(x = cov1, y = cov2, z = TDNN_bias)) +
    geom_contour_fill(bins = 100) +
    geom_contour(color = "black", linewidth = 0.1, bins = 20) +
    scale_fill_distiller(palette = "Spectral") +
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
DNN1_var_plots <- purrr::map(
  .x = indices,
  .f = ~ ggplot(data = stat_tibble[[.x]], aes(x = cov1, y = cov2, z = DNN1_vars)) +
    geom_contour_fill(bins = 100) +
    geom_contour(color = "black", linewidth = 0.1, bins = 20) +
    scale_fill_distiller(palette = "RdBu", direction = -1) +
    xlab("Covariate 1") +
    ylab(element_blank()) +
    labs(fill = "Estimator Variance") +
    theme_light() +
    theme(
      legend.position = "bottom", legend.key.size = unit(1.5, "cm"),
      text = element_text(size = 16)
    )
)

DNN2_var_plots <- purrr::map(
  .x = indices,
  .f = ~ ggplot(data = stat_tibble[[.x]], aes(x = cov1, y = cov2, z = DNN2_vars)) +
    geom_contour_fill(bins = 100) +
    geom_contour(color = "black", linewidth = 0.1, bins = 20) +
    scale_fill_distiller(palette = "RdBu", direction = -1) +
    xlab("Covariate 1") +
    ylab(element_blank()) +
    labs(fill = "Estimator Variance") +
    theme_light() +
    theme(
      legend.position = "bottom", legend.key.size = unit(1.5, "cm"),
      text = element_text(size = 16)
    )
)

TDNN_var_plots <- purrr::map(
  .x = indices,
  .f = ~ ggplot(data = stat_tibble[[.x]], aes(x = cov1, y = cov2, z = TDNN_vars)) +
    geom_contour_fill(bins = 100) +
    geom_contour(color = "black", linewidth = 0.1, bins = 20) +
    scale_fill_distiller(palette = "RdBu", direction = -1) +
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
DNN1_comb_plots <- purrr::map(
  .x = indices,
  .f = ~ DNN1_bias_plots[[.x]] + DNN1_var_plots[[.x]]
)

DNN2_comb_plots <- purrr::map(
  .x = indices,
  .f = ~ DNN2_bias_plots[[.x]] + DNN2_var_plots[[.x]]
)

TDNN_comb_plots <- purrr::map(
  .x = indices,
  .f = ~ TDNN_bias_plots[[.x]] + TDNN_var_plots[[.x]]
)

# Generate names for files
DNN1_file_names <- purrr::map(.x = sim_files,
                              .f = ~ paste0('NPR_Plot_DNN1_', str_split_i(string = .x, pattern = 'reg_exp_1_', i = 2), '.pdf'))
DNN2_file_names <- purrr::map(.x = sim_files,
                              .f = ~ paste0('NPR_Plot_DNN2_', str_split_i(string = .x, pattern = 'reg_exp_1_', i = 2), '.pdf'))
TDNN_file_names <- purrr::map(.x = sim_files,
                              .f = ~ paste0('NPR_Plot_TDNN_', str_split_i(string = .x, pattern = 'reg_exp_1_', i = 2), '.pdf'))

# Save plots
save <- purrr::map(
  .x = indices,
  .f = ~ ggsave(filename = paste0('Graphics/Reg_Exp1/DNN/', DNN1_file_names[[.x]]),
                plot = DNN1_comb_plots[[.x]], width = 16, height = 6, units = "in")
)
save <- purrr::map(
  .x = indices,
  .f = ~ ggsave(filename = paste0('Graphics/Reg_Exp1/DNN/', DNN2_file_names[[.x]]),
                plot = DNN2_comb_plots[[.x]], width = 16, height = 6, units = "in")
)
save <- purrr::map(
  .x = indices,
  .f = ~ ggsave(filename = paste0('Graphics/Reg_Exp1/TDNN/', TDNN_file_names[[.x]]),
                plot = TDNN_comb_plots[[.x]], width = 16, height = 6, units = "in")
)

