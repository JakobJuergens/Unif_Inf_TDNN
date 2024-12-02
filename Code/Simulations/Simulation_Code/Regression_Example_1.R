##### Setup #####
# Load packages from CRAN
library(tidyverse)
library(patchwork)
library(parallel)
library(metR)
RNGkind("L'Ecuyer-CMRG")

# Load my DNN package
install.packages("D:/Jakob_Clouds/dropbox/Research/Individual_Projects/U_Statistics/Unif_Inf_TDNN/Code/tdnnR_0.1.0.tar.gz")
library(tdnnR)

# General Setup
seed <- 394856
set.seed(seed)
n_obs <- 10000
n_reps <- 10000
cov_dim <- 2
s_1 <- 2500
s_2 <- 5000
n_workers <- 20

# Set up regression function and heteroskedasticity function
reg_f <- function(covariates) {
  return(sum(cos(5 * covariates)))
}

htscty_f <- function(covariates) {
  return(0.25 * sum(covariates^2))
}

# Define Simulation function
sim_function <- function(i) {
  # Generate Data
  cov_mat <- matrix(data = runif(n = n_obs * cov_dim, 0, 1), nrow = n_obs, ncol = cov_dim)
  responses <- unlist(purrr::map(
    .f = ~ reg_f(cov_mat[.x, 1:cov_dim]) + rnorm(n = 1, mean = 0, sd = htscty_f(cov_mat[.x, 1:cov_dim])),
    .x = 1:n_obs
  ))

  # Combine data into suitable matrix
  data_mat <- cbind(responses, cov_mat)


  # Estimate regression function at points of interest
  # using TDNN
  estimates <- purrr::map(
    .x = 1:nrow(points),
    .f = ~ tdnnR::TDNN(
      x = points[.x, ], data = data_mat, s1 = s_1, s2 = s_2,
      presorted = FALSE, standardize = FALSE, asymp_approx_weights = TRUE
    )
  )

  # return estimates
  return(estimates)
}

##### Estimation Procedure #####
# Create grid of points of interest
points <- data.matrix(expand.grid(seq(0, 1, 0.05), seq(0, 1, 0.05)))
values <- unlist(purrr::map(
  .x = 1:nrow(points),
  .f = ~ reg_f(covariates = points[.x, ])
))

# Create Cluster and export relevant objects
cl <- makePSOCKcluster(names = n_workers)
junk <- clusterEvalQ(cl, library(metR))
junk <- clusterEvalQ(cl, library(purrr))
clusterExport(cl = cl, varlist = c("points", "TDNN", "reg_f", "htscty_f", "n_obs", "cov_dim", "sim_function", "s_1", "s_2"))
clusterSetRNGStream(cl = cl, iseed = seed)

# Run simulation
results <- clusterApplyLB(cl = cl, x = 1:n_reps, fun = sim_function)

# Stop Cluster
stopCluster(cl)

# Recombine estimates in a useful fashion
TDNN_Estimates <- do.call(purrr::map(.x = 1:n_reps, .f = ~ unlist(purrr::map(1:nrow(points), function(y) results[[.x]][[y]]$TDNN_res))),
  what = cbind
)

DNN1_Estimates <- do.call(purrr::map(.x = 1:n_reps, .f = ~ unlist(purrr::map(1:nrow(points), function(y) results[[.x]][[y]]$DNN1_res))),
  what = cbind
)

DNN2_Estimates <- do.call(purrr::map(.x = 1:n_reps, .f = ~ unlist(purrr::map(1:nrow(points), function(y) results[[.x]][[y]]$DNN2_res))),
  what = cbind
)

# Save Simulation
saveRDS(
  object = list(
    "points" = points, "values" = values,
    "TDNN_estimates" = TDNN_Estimates, "DNN_estimates" = list(DNN1_Estimates, DNN2_Estimates),
    "cov_dim" = cov_dim, "kernel_orders" = c(s_1, s_2), "n_obs" = n_obs, "n_reps" = n_reps
  ),
  file = paste0("Simulation_Results/reg_exp_1_n", n_obs, "s_", s_1, "_", s_2, ".RDS")
)
