#' Calculate the DNN-DML2 CATE-Estimator
#'
#' @param x The point of interest.
#' @param data The data set containing the observations.
#' As a matrix containing the response value in the first column
#' and each covariate in a subsequent columns
#' @param n_folds The number of folds used for crossfitting
#' @param s The subsampling scale.
#' @param presorted True or False whether the data is sorted according to its
#' distance to the point of interest (default value = FALSE)
#' @param standardize True or False whether to standardize (default value = FALSE)
#' @param asymp_approx_weights True or False whether to use an asymptotic
#' approximation for the weights (default value = FALSE)
#' @return A number.
#'
#' @export
DNN_DML2 <- function(x, data, s, n_folds,
                     presorted = FALSE, standardize = FALSE,
                     asymp_approx_weights = TRUE) {

  # Standardize and Sort if chosen
  if(standardize == TRUE){
    stand <- CATE_standardize(x = x, data = data)
    x <- stand$x
    data <- stand$data
  }
  d <- ncol(data) - 2
  if (presorted == FALSE) {
    order <- pre_sort(x = x, cov_mat = data[, 3:(d+2)])
    data <- data[order, ]
  }

  n <- nrow(data)

  # Create List of Fold indices to Estimate Propensity Scores and Regression Functions
  folds <- generate_folds(n_obs = nrow(data), n_folds = n_folds)
  nuisance_par_ests <- matrix(data = NA, nrow = nrow(data), ncol = 3)

  # Identify treated and untreated units
  untreated <- which(data[,2] == 0)
  treated <- which(data[,2] == 1)

  # Estimate Nuisance Parameters
  for (i in 1:n_folds) {
    tmp <- purrr::map(
        .x = folds[[i]],
        .f = ~ DNN_Nuisance(x = data[.x, -(1:2)], data = data[-folds[[i]], ], s = s,
                            asymp_approx_weights = asymp_approx_weights)
      )
    for(j in 1:length(folds[[i]])){
      nuisance_par_ests[folds[[i]][j], ] <- tmp[[j]]
    }
  }

  # Estimate Parameter of Interest using DNN estimator weights
  res <- 0
  factor <- 1
  prefactor <- 0

  # using exact weights
  if (asymp_approx_weights == FALSE) {
    for (i in 1:(n - s + 1)) {
      index <- n - s + 2 - i
      res <- res + factor * CATE_score(X = data[i, -c(1,2)], mu_0 = nuisance_par_ests[i, 2], mu_1 = nuisance_par_ests[i, 3],
                                       Y = data[i, 1], W = data[i, 2], pi = nuisance_par_ests[i, 1])
      prefactor <- prefactor + factor
      factor <- factor * ((n - index + 1) / i)
    }
    res <- res / prefactor
  }

  # using approximate weights
  if (asymp_approx_weights == TRUE) {
    alpha <- s / n

    for (i in 1:(n - s + 1)) {
      if (alpha * (1 - alpha)^(i - 1) == 0) {
        break
      }
      res <- res + alpha * (1 - alpha)^(i - 1) * CATE_score(X = data[i, -c(1,2)], mu_0 = nuisance_par_ests[i, 2], mu_1 = nuisance_par_ests[i, 3],
                                                            Y = data[i, 1], W = data[i, 2], pi = nuisance_par_ests[i, 1])
      prefactor <- prefactor + alpha * (1 - alpha)^(i - 1)
    }

    res <- res / prefactor
  }

  # return calculated result
  names(res) <- NULL
  return(res)
}
