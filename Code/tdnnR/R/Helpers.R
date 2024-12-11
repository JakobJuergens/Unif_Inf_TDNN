#' Calculate the distance between two points of interest
#'
#' @param x Point A of interest
#' @param y Point B of interest
#' @return The Euclidean distance between point A and B
dist <- function(x, y) {
  sqrt(sum((x - y)^2))
}

#' Sorts a sample according to its distance relative
#' to a point of interest (for NPR Setup)
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
NPR_pre_sort <- function(x, data) {
  n <- nrow(data)
  dists <- sqrt(rowSums((sweep(x = data[ , -1], MARGIN = 2, FUN = '-', x)^2)))
  dist_ranks <- order(dists)
  return(data[dist_ranks, ])
}

#' Sorts a sample according to its distance relative
#' to a point of interest (for CATE Setup)
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
CATE_pre_sort <- function(x, data) {
  n <- nrow(data)
  dists <- sqrt(rowSums((sweep(x = data[ , -c(1,2)], MARGIN = 2, FUN = '-', x)^2)))
  dist_ranks <- order(dists)
  return(list('data' = data[dist_ranks, ], 'order' = dist_ranks))
}

#' Sorts a sample according to its distance relative
#' to a point of interest (for CATE Setup)
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
pre_sort <- function(x, cov_mat) {
  dists <- rowSums((sweep(x = cov_mat, MARGIN = 2, FUN = '-', x)^2))
  dist_ranks <- order(dists)
  return(dist_ranks)
}

#' Standardizes a point of interest and a sample
#' to improve comparability between dimensions (for NPR Setup)
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
NPR_standardize <- function(x, data) {
  d <- ncol(data) - 1
  for (i in 1:d) {
    sd <- sd(data[, i + 1])
    data[, i + 1] <- data[, i + 1] / sd
    x[i] <- x[i] / sd
  }
  return(list(x = x, data = data))
}

#' Standardizes a point of interest and a sample
#' to improve comparability between dimensions (for CATE Setup)
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
CATE_standardize <- function(x, data) {
  d <- ncol(data) - 2
  for (i in 2:d) {
    sd <- sd(data[, i + 1])
    data[, i + 1] <- data[, i + 1] / sd
    x[i] <- x[i] / sd
  }
  return(list(x = x, data = data))
}

#' Given a number of observations and a number of desired folds,
#' returns the a list of vectors containing the indices of each fold
#'
#' @param n_obs Number of Observations
#' @param n_folds Number of desired folds
generate_folds <- function(n_obs, n_folds){

  # Randomize Order of Indices
  indices <- sample(x = 1:n_obs, size = n_obs, replace = FALSE)
  # Generate Folds
  obs_p_fold <- ceiling(n_obs / n_folds)
  folds <- purrr::map(
    .x = 1:n_folds,
    .f = ~ sort(indices[((.x - 1)*obs_p_fold + 1):(min(c(.x * obs_p_fold, n_obs)))])
  )
  # Return folds list
  return(folds)
}

#' Evaluates the Neyman-Orthogonal Score Function for the CATE
#' at a given point
#'
#' @param X Covariates
#' @param mu_0 (Estimated) value of the untreated regression function at X
#' @param mu_1 (Estimated) value of the treated regression function at X
#' @param Y Response
#' @param W Treatment Status (coded as 0 for untreated and 1 for treated)
#' @param pi (Estimated) Propensity Score
CATE_score <- function(X, mu_0, mu_1, Y, W, pi){
  return((mu_1 - mu_0) + (W/pi - (1-W)/(1-pi))*(Y - ifelse(W == 1, yes = mu_1, no = mu_0)))
}
