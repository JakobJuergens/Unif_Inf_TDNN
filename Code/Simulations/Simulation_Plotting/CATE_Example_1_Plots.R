# Load packages from CRAN
library(tidyverse)
library(patchwork)
library(metR)

# Load my DNN package
install.packages("D:/Jakob_Clouds/dropbox/Research/Individual_Projects/U_Statistics/Unif_Inf_TDNN/Code/tdnnR_0.1.0.tar.gz")
library(tdnnR)

reg_f0 <- function(covariates){
  return(sum(cos(5*covariates)))
}

reg_f1 <- function(covariates){
  return(sum(cos(5*covariates)) + exp(-sum(covariates)))
}

htscty_f <- function(covariates){
  return(0.25 * sum(covariates^2))
}

propsc_f <- function(covariates){
  return(0.1 + 0.8*(1/(1 + exp(-10*(sum(covariates^0.8) - 1)))))
}

# Create data-frame for plotting purposes
plot_df <- as_tibble(expand.grid(seq(0,1,0.01), seq(0,1,0.01)))
plot_df$reg_vals0 <- unlist(purrr::map(.f = ~ reg_f0(plot_df[.x, 1:2]), .x = 1:nrow(plot_df)))
plot_df$reg_vals1 <- unlist(purrr::map(.f = ~ reg_f1(plot_df[.x,1:2]), .x = 1:nrow(plot_df)))
plot_df$var_vals <- unlist(purrr::map(.f = ~ htscty_f(plot_df[.x,1:2]), .x = 1:nrow(plot_df)))
plot_df$propsc_vals <- unlist(purrr::map(.f = ~ propsc_f(plot_df[.x,1:2]), .x = 1:nrow(plot_df)))

# Plot response surfaces
reg0_plot <- ggplot(data = plot_df) +
  ggtitle(label = 'Without Treatment') +
  geom_contour_fill(aes(x = Var1, y = Var2, z = reg_vals0),
                    bins = 100) +
  xlab("Covariate 1") + ylab("Covariate 2") + labs(fill = 'Value of Regression Function') +
  theme_light()

reg1_plot <- ggplot(data = plot_df) +
  ggtitle(label = 'With Treatment') +
  geom_contour_fill(aes(x = Var1, y = Var2, z = reg_vals1),
                    bins = 100) +
  xlab("Covariate 1") + ylab("Covariate 2") + labs(fill = 'Value of Regression Function') +
  theme_light()

# Plot Propensity Score
propsc_plot <- ggplot(data = plot_df) +
  ggtitle(label = 'Propensity Score') +
  geom_contour_fill(aes(x = Var1, y = Var2, z = propsc_vals),
                    bins = 100) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(0,1)) +
  xlab("Covariate 1") + ylab("Covariate 2") + labs(fill = 'Propensity Score') +
  theme_light() +
  theme(legend.position="bottom", legend.key.size = unit(1, 'cm'),
        text = element_text(size = 16))


# Plot Heteroskedasticity
var_plot <- ggplot(data = plot_df) +
  ggtitle(label = 'Error Variance') +
  geom_contour_fill(aes(x = Var1, y = Var2, z = var_vals),
                    bins = 100) +
  scale_fill_distiller(palette='RdBu', direction=-1) +
  xlab("Covariate 1") + ylab(element_blank()) + labs(fill = 'Error Term Variance') +
  theme_light() +
  theme(legend.position="bottom", legend.key.size = unit(1, 'cm'),
        text = element_text(size = 16))

# Combine Plots for Paper
comb_plot <- ((reg0_plot + reg1_plot) +
                plot_layout(guides = "collect", axes = "collect") &
                scale_fill_distiller(palette = "Spectral", limits = c(-2,3)) &
                theme(legend.position="bottom", legend.key.size = unit(1, 'cm'),
                      text = element_text(size = 16))) /
  ((propsc_plot + var_plot) +
     plot_layout(guides = "collect", axes = "collect") &
     theme(legend.position="bottom", legend.key.size = unit(1, 'cm'),
           text = element_text(size = 16)))

ggsave(filename = '../../Graphics/CATE_Exmp1.pdf', plot = comb_plot,
       width = 16, height = 12, units = 'in')
