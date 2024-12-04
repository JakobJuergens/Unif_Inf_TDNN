#' Calculate the distributional nearest neighbor estimate for a point x
#' based on a sample data with subsampling scale s.
#' The implementation is based on the equivalent representation as an L-statistic
#' described in Steele (2009) and the TDNN estimator from Demirkaya et al. (2024)
#'
#' @param x The point of interest.
#' @param data The data set containing the observations.
#' As a matrix containing the response value in the first column
#' and each covariate in a subsequent columns
#' @param s1 The first subsampling scale.
#' @param s2 The second subsampling scale.
#' @param presorted True or False whether the data is sorted according to its
#' distance to the point of interest (default value = FALSE)
#' @param standardize True or False whether to standardize (default value = FALSE)
#' @param asymp_approx_weights True or False whether to use an asymptotic
#' approximation for the weights (default value = FALSE)
#' @return A number.
#'
#' @export
TDNN <- function(x, data, s1, s2,
                 presorted = FALSE, standardize = FALSE,
                 asymp_approx_weights = TRUE, verbose = FALSE) {
  # Check inputs for executing function
  point_check(x)
  data_check(data)
  data_point_compatibility_check(x, data)
  samplingscale_check2(s1, s2, data)

  # Standardize Data
  if(standardize == TRUE){
    stand <- NPR_standardize(x = x, data = data)
    x <- stand$x
    data <- stand$data
  }

  # Pre-sort data relative to point of interest
  if (presorted == FALSE) {
    data <- NPR_pre_sort(x = x, data = data)
  }

  Y <- as.vector(data[, 1])
  d <- ncol(data) - 1
  n <- length(Y)

  if (verbose) {
    print(paste0("d = ", d, ", n = ", n))
  }

  # Check that s2 > s1 and reverse if not
  if (s1 > s2) {
    tmp <- s1
    s1 <- s2
    s2 <- tmp
  }

  if (verbose) {
    print(paste0("s1 = ", s1, ", s2 = ", s2))
  }
  # Calculate weights for two-scale DNN
  w1 <- (1 - (s1 / s2)^(-2 / d))^(-1)
  w2 <- 1 - w1

  if (verbose) {
    print(paste0("w1 = ", w1, ", w2 = ", w2))
  }

  # Calculate the two DNN estimators
  res_1 <- 0
  res_2 <- 0
  prefactor_1 <- 0
  prefactor_2 <- 0

  # using exact weights
  if(asymp_approx_weights == FALSE){
    factor_1 <- 1
    factor_2 <- 1

    for (i in 1:(s2 - s1)) {
      index <- n - s1 + 2 - i
      res_1 <- res_1 + factor_1 * Y[index]
      prefactor_1 <- prefactor_1 + factor_1
      factor_1 <- factor_1 * ((n - index + 1) / i)

    }

    for (i in 1:(n - s2 + 1)) {
      index <- n - s2 + 2 - i

      res_1 <- res_1 + factor_1 * Y[index]
      res_2 <- res_2 + factor_2 * Y[index]
      prefactor_1 <- prefactor_1 + factor_1
      prefactor_2 <- prefactor_2 + factor_2

      factor_1 <- factor_1 * ((n - index + 1) / (i + s2 - s1))
      factor_2 <- factor_2 * ((n - index + 1) / i)
    }
  }

  # using asymptotically approximated weights
  if(asymp_approx_weights == TRUE){
    alpha_1 <- s_1/n
    alpha_2 <- s_2/n

    for (i in 1:(n - s_1 + 1)) {
      if(alpha_1*(1-alpha_1)^(i-1) == 0){break}
      res_1 <- res_1 + alpha_1*(1-alpha_1)^(i-1) * Y[i]
      prefactor_1 <- prefactor_1 + alpha_1*(1-alpha_1)^(i-1)
    }

    for (i in 1:(n - s_2 + 1)) {
      if(alpha_2*(1-alpha_2)^(i-1) == 0){break}
      res_2 <- res_2 + alpha_2*(1-alpha_2)^(i-1) * Y[i]
      prefactor_2 <- prefactor_2 + alpha_2*(1-alpha_2)^(i-1)
    }
  }

  if (verbose) {
    print(paste0("prefactor_1 = ", prefactor_1, ", prefactor_2 = ", prefactor_2))
  }

  res_1 <- res_1/prefactor_1
  res_2 <- res_2/prefactor_2

  if (verbose) {
    print(paste0("res_1 = ", res_1, ", res_2 = ", res_2))
  }

  # Combine results
  res <- w1 * res_1 + w2 * res_2

  # return results
  return(list("TDNN_res" = res, "DNN1_res" = res_1, "DNN2_res" = res_2))
}
