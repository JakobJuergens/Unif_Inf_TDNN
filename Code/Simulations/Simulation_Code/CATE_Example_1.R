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
  return(0.25 * sum(covariates^2))
}

propsc_f <- function(covariates){
  return(0.1 + 0.8*(1/(1 + exp(-10*(sum(covariates^0.8) - 1)))))
}

reg_vals <- cbind(unlist(purrr::map(.f = ~ reg_f0(cov_mat[.x, ]), .x = 1:n_obs)),
                  unlist(purrr::map(.f = ~ reg_f1(cov_mat[.x, ]), .x = 1:n_obs)))

errors <-  unlist(
  purrr::map(.f = ~ rnorm(n = 1, mean = 0, sd = htscty_f(cov_mat[.x, ])), .x = 1:n_obs))


