#' Calculate the distributional nearest neighbor estimate for a point x
#' based on a sample data with subsampling scale s.
#' The implementation is based on the equivalent representation as an L-statistic
#' described in Steele (2009)
#'
#' @param x The point of interest.
#' @param data The data set containing the observations.
#' As a matrix containing the response value in the first column
#' and each covariate in a subsequent columns
#' @param s The subsampling scale.
#' @param presorted True or False whether the data is sorted according to its
#' distance to the point of interest (default value = FALSE)
#' @param standardize True or False whether to standardize (default value = FALSE)
#' @param asymp_approx_weights True or False whether to use an asymptotic
#' approximation for the weights (default value = FALSE)
#' @return A number.
#'
#' @export
DNN <- function(x, data, s,
                presorted = FALSE, standardize = FALSE,
                asymp_approx_weights = TRUE, verbose = FALSE) {

  # Check inputs for executing function
  point_check(x)
  data_check(data)
  data_point_compatibility_check(x, data)
  samplingscale_check(s, data)

  # Standardize Data
  if (standardize == TRUE) {
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

  # Calculate DNN estimator
  res <- 0
  factor <- 1
  prefactor <- 0

  # using exact weights
  if (asymp_approx_weights == FALSE) {
    for (i in 1:(n - s + 1)) {
      index <- n - s + 2 - i
      res <- res + factor * Y[index]
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
      res <- res + alpha * (1 - alpha)^(i - 1) * Y[i]
      prefactor <- prefactor + alpha * (1 - alpha)^(i - 1)
    }

    res <- res / prefactor
  }

  # return calculated result
  return(res)
}
