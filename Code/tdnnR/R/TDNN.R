#' Calculate the distributional nearest neighbor estimate for a point x
#' based on a sample data with subsampling scale s.
#' The implementation is based on the equivalent representation as an L-statistic
#' described in Steele (2009)
#'
#' @param x The point of interest.
#' @param data The data set containing the observations.
#' @param s The subsampling scale.
#' @return A number.
#'
#' @export
DNN <- function(x, data, s, presorted = FALSE) {
  # Check inputs for executing function
  point_check(x)
  data_check(data)
  data_point_compatibility_check(x, data)
  samplingscale_check(s, data)

  # Pre-sort data relative to point of interest
  if (presorted == FALSE) {
    data <- pre_sort(x = x, data = data)
  }

  n <- nrow(data)
  d <- ncol(data) - 1

  # Calculate DNN estimator
  res <- 0
  factor <- 1
  prefactor <- 1 / choose(n, s)
  for (i in 1:(n - s + 1)) {
    index <- n - s + 2 - i
    res <- res + factor * data[index, 1]
    factor <- factor * ((n - index + 1) / (n - index - s + 2))
  }
  res <- prefactor * res

  # return calculated result
  return(res)
}

#' Calculate the distributional nearest neighbor estimate for a point x
#' based on a sample data with subsampling scale s.
#' The implementation is based on the equivalent representation as an L-statistic
#' described in Steele (2009) and the TDNN estimator from Demirkaya et al. (2024)
#'
#' @param x The point of interest.
#' @param data The data set containing the observations.
#' @param s1 The first subsampling scale.
#' @param s2 The second subsampling scale.
#' @return A number.
#'
#' @export
TDNN <- function(x, data, s1, s2, presorted = FALSE) {
  # Check inputs for executing function
  point_check(x)
  data_check(data)
  data_point_compatibility_check(x, data)
  samplingscale_check2(s1, s2, data)

  # Pre-sort data relative to point of interest
  if (presorted == FALSE) {
    data <- pre_sort(x = x, data = data)
  }

  # Check that s2 > s1 and reverse if not
  if (s1 > s2) {
    tmp <- s1
    s_1 <- s2
    s_2 <- tmp
  }

  # Calculate weights for two-scale DNN
  w1 <- 0.5
  w2 <- 1 - w_1

  # Calculate the two DNN estimators
  res1 <- 0
  res2 <- 0
  factor1 <- 1
  factor2 <- 1
  prefactor1 <- 1 / choose(n, s1)
  prefactor2 <- 1 / choose(n, s2)

  for (i in 1:(s2 - s1)) {
    index <- n - s2 + 2 - i
    res2 <- res2 + factor2 * data[index, 1]
    factor2 <- factor2 * ((n - index + 1) / (n - index - s2 + 2))
  }

  for (i in 1:(n - s1 + 1)) {
    index <- n - s1 + 2 - i

    res1 <- res1 + factor1 * data[index, 1]
    res2 <- res2 + factor2 * data[index, 1]

    factor1 <- factor1 * ((n - index + 1) / (n - index - s1 + 2))
    factor2 <- factor2 * ((n - index + 1) / (n - index - s2 + 2))
  }

  res1 <- prefactor1 * res1
  res2 <- prefactor2 * res2

  # Combine results
  res <- w1 * res1 + w2 * res2

  # return results
  return(res)
}
