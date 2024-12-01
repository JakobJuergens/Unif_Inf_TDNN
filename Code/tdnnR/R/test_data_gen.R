reg_f <- function(covariates) {
  return(sum(cos(5 * covariates)))
}

htscty_f <- function(covariates) {
  return(0.25 * abs(sum(covariates^2)))
}

#' Generate a data set for the nonparametric regression setup.
#'
#' @param n_obs Number of desired observations
#' @param cov_dim Dimensionality of covariates
#' @return A number.
#'
reg_test_data <- function(n_obs, cov_dim) {
  cov_mat <- matrix(data = runif(n = n_obs * cov_dim, 0, 1), nrow = n_obs, ncol = cov_dim)
  responses <- unlist(purrr::map(
    .f = ~ reg_f(cov_mat[.x, 1:cov_dim]) + rnorm(n = 1, mean = 0, sd = htscty_f(cov_mat[.x, 1:cov_dim])),
    .x = 1:n_obs
  ))

  # Combine data into suitable matrix
  data_mat <- cbind(responses, cov_mat)

  return(data_mat)
}
