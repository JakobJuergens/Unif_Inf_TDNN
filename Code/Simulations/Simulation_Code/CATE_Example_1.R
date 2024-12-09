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
n_reps <- 5000
cov_dim <- 2
s_1 <- 120
s_2 <- 250
n_folds <- 50
n_workers <- 20

# Define Data generating functions
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

CATE_f <- function(covariates){
  return(exp(-sum(covariates)))
}

# Define Simulation function
sim_function <- function(i) {
  cov_mat <- matrix(
    data = runif(n = n_obs * cov_dim, 0, 1),
    nrow = n_obs, ncol = cov_dim
  )

  prop_scores <- unlist(purrr::map(
    .x = 1:n_obs,
    .f = ~ propsc_f(cov_mat[.x, 1:cov_dim])
  ))

  treatments <- unlist(purrr::map(
    .x = 1:n_obs,
    .f = ~ rbinom(n = 1, size = 1, prob = prop_scores[.x])
  ))

  responses <- unlist(purrr::map(
    .x = 1:n_obs,
    .f = ~ ifelse(treatments[.x] == 1, yes = reg_f1(cov_mat[.x, 1:cov_dim]), no = reg_f0(cov_mat[.x, 1:cov_dim]))
    + rnorm(n = 1, mean = 0, sd = htscty_f(cov_mat[.x, 1:cov_dim]))
  ))

  # Combine data into suitable matrix
  data_mat <- cbind(responses, treatments, cov_mat)

  # Estimate regression function at points of interest
  # using TDNN
  estimates <- tdnnR::TDNN_DML2_mult(points = points, data = data_mat, s1 = s_1, s2 = s_2, n_folds = n_folds,
       asymp_approx_weights = TRUE)

  # return estimates
  return(estimates)
}

##### Estimation Procedure #####
# Create grid of points of interest
points <- data.matrix(expand.grid(seq(0, 1, 0.1), seq(0, 1, 0.1)))
values <- unlist(purrr::map(
  .x = 1:nrow(points),
  .f = ~ CATE_f(covariates = points[.x, ])
))

# Create Cluster and export relevant objects
cl <- makePSOCKcluster(names = n_workers)
junk <- clusterEvalQ(cl, library(metR))
junk <- clusterEvalQ(cl, library(purrr))
clusterExport(cl = cl, varlist = c("points", "TDNN", "reg_f0", "reg_f1", "propsc_f",
                                   "htscty_f", "n_obs", "cov_dim", "n_folds", "sim_function", "s_1", "s_2"))
clusterSetRNGStream(cl = cl, iseed = seed)

# Run simulation
results <- clusterApplyLB(cl = cl, x = 1:n_reps, fun = sim_function)

# Stop Cluster
stopCluster(cl)

# Recombine estimates in a useful fashion
TDNN_Estimates <- do.call(purrr::map(.x = 1:n_reps, .f = ~ unlist(purrr::map(1:nrow(points), function(y) results[[.x]][y,1]))),
                          what = cbind
)

DNN1_Estimates <- do.call(purrr::map(.x = 1:n_reps, .f = ~ unlist(purrr::map(1:nrow(points), function(y) results[[.x]][y,2]))),
                          what = cbind
)

DNN2_Estimates <- do.call(purrr::map(.x = 1:n_reps, .f = ~ unlist(purrr::map(1:nrow(points), function(y) results[[.x]][y,3]))),
                          what = cbind
)

# Save Simulation
saveRDS(
  object = list(
    "points" = points, "values" = values,
    "TDNN_estimates" = TDNN_Estimates, "DNN_estimates" = list(DNN1_Estimates, DNN2_Estimates),
    "cov_dim" = cov_dim, "kernel_orders" = c(s_1, s_2), "n_obs" = n_obs, "n_reps" = n_reps, "n_folds" = n_folds
  ),
  file = paste0("Simulation_Results/CATE_Exp1/CATE_exp_1_n", n_obs, "s_", s_1, "_", s_2, ".RDS")
)
