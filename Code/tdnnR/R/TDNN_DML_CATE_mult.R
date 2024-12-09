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
TDNN_DML2_mult <- function(points, data, s1, s2, n_folds, #standardize = FALSE,
                          asymp_approx_weights = TRUE) {

  # Standardize and Sort if chosen
  #if(standardize == TRUE){
  #  stand <- CATE_standardize(x = x, data = data)
  #  x <- stand$x
  #  data <- stand$data
  #}

  # Check that s2 > s1 and reverse if not
  if (s1 > s2) {
    tmp <- s1
    s1 <- s2
    s2 <- tmp
  }

  d <- ncol(data) - 2
  n <- nrow(data)
  n_points <- nrow(points)

  # Calculate weights for two-scale DNN
  w1 <- (1 - (s1 / s2)^(-2 / d))^(-1)
  w2 <- 1 - w1

  # Create List of Fold indices to Estimate Propensity Scores and Regression Functions
  folds <- generate_folds(n_obs = n, n_folds = n_folds)
  nuisance_par_ests <- matrix(data = NA, nrow = n, ncol = 3)

  # Estimate Nuisance Parameters
  for (i in 1:n_folds) {
    tmp <- purrr::map(
      .x = folds[[i]],
      .f = ~ TDNN_Nuisance(x = data[.x, -(1:2)], data = data[-folds[[i]], ], s1 = s1, s2 = s2,
                          asymp_approx_weights = asymp_approx_weights)
    )
    for(j in 1:length(folds[[i]])){
      nuisance_par_ests[folds[[i]][j], ] <- tmp[[j]]
    }
  }

  # Calculate Neyman-Orthogonal Score at each Observation
  NO_mom <- unlist(purrr::map(.x = 1:n,
                       .f = ~ CATE_score(X = data[.x, 3:(d+2)], mu_0 = nuisance_par_ests[.x, 2], mu_1 = nuisance_par_ests[.x, 3],
                                         Y = data[.x, 1], W = data[.x, 2], pi = nuisance_par_ests[.x, 1])))

  # Estimate Parameter of Interest using DNN estimator weights
  estimates <- matrix(NA, nrow = n_points, ncol = 3)

  for(j in 1:n_points){
    x <- points[j,]
    ordering <- CATE_pre_sort(x = x, data = data)$order
    NO_mom_j <- NO_mom[ordering]

    res_1 <- 0
    res_2 <- 0
    factor_1 <- 1
    factor_2 <- 1
    prefactor_1 <- 0
    prefactor_2 <- 0

    # using approximate weights
    if (asymp_approx_weights == TRUE) {
      alpha_1 <- s1/n
      alpha_2 <- s2/n

      for (i in 1:(n - s1 + 1)) {
        if(alpha_1*(1-alpha_1)^(i-1) == 0){break}
        res_1 <- res_1 + alpha_1*(1-alpha_1)^(i-1) * NO_mom_j[i]
        prefactor_1 <- prefactor_1 + alpha_1*(1-alpha_1)^(i-1)
      }

      for (i in 1:(n - s2 + 1)) {
        if(alpha_2*(1-alpha_2)^(i-1) == 0){break}
        res_2 <- res_2 + alpha_2*(1-alpha_2)^(i-1) * NO_mom_j[i]
        prefactor_2 <- prefactor_2 + alpha_2*(1-alpha_2)^(i-1)
      }
    }

    res_1 <- res_1/prefactor_1
    res_2 <- res_2/prefactor_2

    names(res_1) <- NULL
    names(res_2) <- NULL

    # Combine results
    res <- w1 * res_1 + w2 * res_2

    estimates[j,] <- c(res, res_1, res_2)
  }


  # return calculated result
  names(estimates) <- NULL
  return(estimates)
}
