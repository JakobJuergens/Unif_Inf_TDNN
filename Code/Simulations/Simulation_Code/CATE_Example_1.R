# Load packages from CRAN
library(tidyverse)
library(patchwork)
library(metR)

# Load my DNN package
install.packages("D:/Jakob_Clouds/dropbox/Research/Individual_Projects/U_Statistics/Unif_Inf_TDNN/Code/tdnnR_0.1.0.tar.gz")
library(tdnnR)

# General Setup
n_obs <- 10000

# Generate Covariates
cov_dim <- 2
cov_mat <- matrix(data = runif(n = n_obs*cov_dim, 0, 1), nrow = n_obs, ncol = cov_dim)

# Calculate Response Values
responses <- rep(NA, times = n_obs)

reg_f0 <- function(covariates){
  return(sum(cos(5*covariates)))
}

reg_f1 <- function(covariates){
  return(sum(cos(5*covariates)) + exp(-sum(covariates)))
}

htscty_f <- function(covariates){
  return(0.5*abs(sum(covariates)))
}

reg_vals <- cbind(unlist(purrr::map(.f = ~ reg_f0(cov_mat[.x, ]), .x = 1:n_obs)),
                  unlist(purrr::map(.f = ~ reg_f1(cov_mat[.x, ]), .x = 1:n_obs)))

errors <-  unlist(
  purrr::map(.f = ~ rnorm(n = 1, mean = 0, sd = htscty_f(cov_mat[.x, ])), .x = 1:n_obs))


# Create data-frame for plotting purposes
plot_df <- as_tibble(expand.grid(seq(0,1,0.01), seq(0,1,0.01)))
plot_df$reg_vals0 <- unlist(purrr::map(.f = ~ reg_f0(plot_df[.x,]), .x = 1:nrow(plot_df)))
plot_df$reg_vals1 <- unlist(purrr::map(.f = ~ reg_f1(plot_df[.x,]), .x = 1:nrow(plot_df)))
plot_df$var_vals <- unlist(purrr::map(.f = ~ htscty_f(plot_df[.x,]), .x = 1:nrow(plot_df)))

# Plot response surfaces
reg0_plot <- ggplot(data = plot_df) +
  geom_contour_fill(aes(x = Var1, y = Var2, z = reg_vals0),
                    bins = 100) +
  scale_fill_distiller(palette = "Spectral") +
  scale_colour_continuous(limits = c(-3,3)) +
  xlab("Covariate 1") + ylab("Covariate 2") + labs(fill = 'Value of Regression Function') +
  theme_light()

reg1_plot <- ggplot(data = plot_df) +
  geom_contour_fill(aes(x = Var1, y = Var2, z = reg_vals1),
                    bins = 100) +
  scale_fill_distiller(palette = "Spectral") +
  xlab("Covariate 1") + ylab("Covariate 2") + labs(fill = 'Value of Regression Function') +
  theme_light()


# Plot Heteroskedasticity
var_plot <- ggplot(data = plot_df) +
  geom_contour_fill(aes(x = Var1, y = Var2, z = var_vals),
                    bins = 100) +
  scale_fill_distiller(palette='RdBu', direction=-1) +
  xlab("Covariate 1") + ylab(element_blank()) + labs(fill = 'Error Term Variance') +
  theme_light() +
  theme(legend.position="bottom", legend.key.size = unit(1.5, 'cm'),
        text = element_text(size = 16))

# Combine Plots for Paper
comb_plot <- reg0_plot / reg1_plot + plot_layout(guides = "collect", axes = "collect") &
  scale_fill_distiller(palette = "Spectral", limits = c(-2,3)) &
  theme(legend.position="bottom", legend.key.size = unit(1.5, 'cm'),
        text = element_text(size = 16))

ggsave(filename = '../../Graphics/CATE_Exmp1.pdf', plot = comb_plot,
       width = 16, height = 12, units = 'in')
