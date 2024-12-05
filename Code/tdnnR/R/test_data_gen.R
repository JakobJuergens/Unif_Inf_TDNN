reg_f0 <- function(covariates) {
  return(sum(cos(5 * covariates)))
}

reg_f1 <- function(covariates) {
  return(sum(cos(5 * covariates)) + exp(-sum(covariates)))
}

htscty_f <- function(covariates) {
  return(0.25 * sum(covariates^2))
}

propsc_f <- function(covariates) {
  return(0.1 + 0.8 * (1 / (1 + exp(-10 * (sum(covariates^0.8) - 1)))))
}

CATE_f <- function(covariates){
  return(exp(-sum(covariates)))
}

#' Generate a data set for the nonparametric regression setup.
#'
#' @param n_obs Number of desired observations
#' @param cov_dim Dimensionality of covariates
#' @return A suitable matrix of data.
#'
reg_test_data <- function(n_obs, cov_dim) {
  cov_mat <- matrix(data = runif(n = n_obs * cov_dim, 0, 1), nrow = n_obs, ncol = cov_dim)

  responses <- unlist(purrr::map(
    .x = 1:n_obs,
    .f = ~ reg_f0(cov_mat[.x, 1:cov_dim]) + rnorm(n = 1, mean = 0, sd = htscty_f(cov_mat[.x, 1:cov_dim]))
  ))

  # Combine data into suitable matrix
  data_mat <- cbind(responses, cov_mat)

  return(data_mat)
}

#' Generate a data set for the nonparametric regression setup.
#'
#' @param n_obs Number of desired observations
#' @param cov_dim Dimensionality of covariates
#' @return A suitable matrix of data.
#'
CATE_test_data <- function(n_obs, cov_dim) {
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

  return(data_mat)
}
