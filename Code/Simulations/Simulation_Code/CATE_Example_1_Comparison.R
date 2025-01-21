##### Setup #####
# Load packages from CRAN
library(grf)
library(tidyverse)
library(patchwork)
library(parallel)
library(metR)
RNGkind("L'Ecuyer-CMRG")

# General Setup
seed <- 394856
set.seed(seed)
n_obs <- 10000
n_reps <- 5000
cov_dim <- 2
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
  cforest <- causal_forest(X = cov_mat, Y = responses, W = treatments)
  estimates <- predict(object = cforest, points)
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
junk <- clusterEvalQ(cl, library(grf))
clusterExport(cl = cl, varlist = c("points", "reg_f0", "reg_f1", "propsc_f",
                                   "htscty_f", "n_obs", "cov_dim", "sim_function"))
clusterSetRNGStream(cl = cl, iseed = seed)

# Run simulation
results <- clusterApplyLB(cl = cl, x = 1:n_reps, fun = sim_function)

# Stop Cluster
stopCluster(cl)

# Recombine estimates in a useful fashion
CF_Estimates <- do.call(purrr::map(.x = 1:n_reps, .f = ~ unlist(purrr::map(1:nrow(points), function(y) results[[.x]]))),
                          what = cbind
)

# Save Simulation
saveRDS(
  object = list(
    "points" = points, "values" = values,
    "CF_estimates" = CF_Estimates,
    "cov_dim" = cov_dim, "n_obs" = n_obs, "n_reps" = n_reps
  ),
  file = paste0("Simulation_Results/CATE_Exp1/CATE_exp_1_n", n_obs, "_CF", ".RDS")
)
