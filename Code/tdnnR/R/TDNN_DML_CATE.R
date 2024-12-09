#' Calculate the TDNN-DML2 CATE-Estimator
#'
#' @param x The point of interest.
#' @param data The data set containing the observations.
#' As a matrix containing the response value in the first column
#' and each covariate in a subsequent columns
#' @param s1 The first subsampling scale.
#' @param s2 The second subsampling scale.
#' @param n_folds The number of folds used for crossfitting
#' @param presorted True or False whether the data is sorted according to its
#' distance to the point of interest (default value = FALSE)
#' @param standardize True or False whether to standardize (default value = FALSE)
#' @param asymp_approx_weights True or False whether to use an asymptotic
#' approximation for the weights (default value = FALSE)
#' @return A number.
#'
#' @export
TDNN_DML2 <- function(x, data, s1, s2, n_folds,
                      presorted = FALSE, standardize = FALSE,
                      asymp_approx_weights = TRUE) {

    # Standardize and Sort if chosen
  if(standardize == TRUE){
    stand <- CATE_standardize(x = x, data = data)
    x <- stand$x
    data <- stand$data
  }
  if (presorted == FALSE) {
    data <- CATE_pre_sort(x = x, data = data)
  }

  # Check that s2 > s1 and reverse if not
  if (s1 > s2) {
    tmp <- s1
    s1 <- s2
    s2 <- tmp
  }

  # Calculate weights for two-scale DNN
  d <- ncol(data) - 2
  w1 <- (1 - (s1 / s2)^(-2 / d))^(-1)
  w2 <- 1 - w1

  n <- nrow(data)

  # Create List of Fold indices to Estimate Propensity Scores and Regression Functions
  folds <- generate_folds(n_obs = nrow(data), n_folds = n_folds)
  nuisance_par_ests <- matrix(data = NA, nrow = n, ncol = 3)

  # Identify treated and untreated units
  untreated <- which(data[,2] == 0)
  treated <- which(data[,2] == 1)

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

  # Estimate Parameter of Interest using TDNN estimator weights
  res_1 <- 0
  res_2 <- 0
  factor_1 <- 1
  factor_2 <- 1
  prefactor_1 <- 0
  prefactor_2 <- 0

  # using exact weights
  if (asymp_approx_weights == FALSE) {
    for (i in 1:(s2 - s1)) {
      index <- n - s1 + 2 - i
      res_1 <- res_1 + factor_1 * CATE_score(X = data[index, -c(1,2)], mu_0 = nuisance_par_ests[index, 2], mu_1 = nuisance_par_ests[index, 3],
                                             Y = data[index, 1], W = data[index, 2], pi = nuisance_par_ests[index, 1])
      prefactor_1 <- prefactor_1 + factor_1
      factor_1 <- factor_1 * ((n - index + 1) / i)

    }

    for (i in 1:(n - s2 + 1)) {
      index <- n - s2 + 2 - i

      res_1 <- res_1 + factor_1 * CATE_score(X = data[index, -c(1,2)], mu_0 = nuisance_par_ests[index, 2], mu_1 = nuisance_par_ests[index, 3],
                                             Y = data[index, 1], W = data[index, 2], pi = nuisance_par_ests[index, 1])
      res_2 <- res_2 + factor_2 * CATE_score(X = data[index, -c(1,2)], mu_0 = nuisance_par_ests[i, 2], mu_1 = nuisance_par_ests[i, 3],
                                             Y = data[index, 1], W = data[index, 2], pi = nuisance_par_ests[index, 1])
      prefactor_1 <- prefactor_1 + factor_1
      prefactor_2 <- prefactor_2 + factor_2

      factor_1 <- factor_1 * ((n - index + 1) / (i + s2 - s1))
      factor_2 <- factor_2 * ((n - index + 1) / i)
    }
  }

  # using approximate weights
  if (asymp_approx_weights == TRUE) {
    alpha_1 <- s1/n
    alpha_2 <- s2/n

    for (i in 1:(n - s1 + 1)) {
      if(alpha_1*(1-alpha_1)^(i-1) == 0){break}
      res_1 <- res_1 + alpha_1*(1-alpha_1)^(i-1) * CATE_score(X = data[i, -c(1,2)], mu_0 = nuisance_par_ests[i, 2], mu_1 = nuisance_par_ests[i, 3],
                                                              Y = data[i, 1], W = data[i, 2], pi = nuisance_par_ests[i, 1])
      prefactor_1 <- prefactor_1 + alpha_1*(1-alpha_1)^(i-1)
    }

    for (i in 1:(n - s2 + 1)) {
      if(alpha_2*(1-alpha_2)^(i-1) == 0){break}
      res_2 <- res_2 + alpha_2*(1-alpha_2)^(i-1) * CATE_score(X = data[i, -c(1,2)], mu_0 = nuisance_par_ests[i, 2], mu_1 = nuisance_par_ests[i, 3],
                                                              Y = data[i, 1], W = data[i, 2], pi = nuisance_par_ests[i, 1])
      prefactor_2 <- prefactor_2 + alpha_2*(1-alpha_2)^(i-1)
    }
  }

  res_1 <- res_1/prefactor_1
  res_2 <- res_2/prefactor_2

  names(res_1) <- NULL
  names(res_2) <- NULL

  # Combine results
  res <- w1 * res_1 + w2 * res_2

  # return calculated result
  return(list("TDNN_res" = res, "DNN1_res" = res_1, "DNN2_res" = res_2))
}
