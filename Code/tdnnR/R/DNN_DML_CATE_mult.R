#' Calculate the DNN-DML2 CATE-Estimator for a list of points of interest
#'
#' @param points A matrix containing points of interest as its rows.
#' @param data The data set containing the observations.
#' As a matrix containing the response value in the first column
#' and each covariate in a subsequent columns
#' @param n_folds The number of folds used for crossfitting
#' @param s The subsampling scale.
#' @param asymp_approx_weights True or False whether to use an asymptotic
#' approximation for the weights (default value = FALSE)
#' @return A number.
#'
#' @export
DNN_DML2_mult <- function(points, data, s, n_folds, #standardize = FALSE,
                     asymp_approx_weights = TRUE) {

  # Standardize and Sort if chosen
  #if(standardize == TRUE){
  #  stand <- CATE_standardize(x = x, data = data)
  #  x <- stand$x
  #  data <- stand$data
  #}

  n <- nrow(data)
  n_points <- nrow(points)

  # Create List of Fold indices to Estimate Propensity Scores and Regression Functions
  folds <- generate_folds(n_obs = n, n_folds = n_folds)
  nuisance_par_ests <- matrix(data = NA, nrow = n, ncol = 3)

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
  estimates <- vector(mode = 'numeric', length = n_points)

  for(j in 1:n_points){
    x <- points[j,]
    ordering <- CATE_pre_sort(x = x, data = data)$order
    data_j <- data[ordering, ]
    nuisance_par_ests_j <- nuisance_par_ests[ordering, ]

    res <- 0
    factor <- 1
    prefactor <- 0

    # using exact weights
    if (asymp_approx_weights == FALSE) {
      for (i in 1:(n - s + 1)) {
        index <- n - s + 2 - i
        res <- res + factor * CATE_score(X = data_j[i, -c(1,2)], mu_0 = nuisance_par_ests_j[i, 2], mu_1 = nuisance_par_ests_j[i, 3],
                                         Y = data_j[i, 1], W = data_j[i, 2], pi = nuisance_par_ests_j[i, 1])
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
        res <- res + alpha * (1 - alpha)^(i - 1) * CATE_score(X = data_j[i, -c(1,2)], mu_0 = nuisance_par_ests_j[i, 2], mu_1 = nuisance_par_ests_j[i, 3],
                                                              Y = data_j[i, 1], W = data_j[i, 2], pi = nuisance_par_ests_j[i, 1])
        prefactor <- prefactor + alpha * (1 - alpha)^(i - 1)
      }
      res <- res / prefactor
    }

    estimates[j] <-res
  }


  # return calculated result
  names(estimates) <- NULL
  return(estimates)
}
