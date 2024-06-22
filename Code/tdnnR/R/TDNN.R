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
TDNN <- function(x, data, s1, s2, presorted = FALSE, verbose = FALSE) {
  # Check inputs for executing function
  point_check(x)
  data_check(data)
  data_point_compatibility_check(x, data)
  samplingscale_check2(s1, s2, data)

  # Pre-sort data relative to point of interest
  if (presorted == FALSE) {
    data <- pre_sort(x = x, data = data)
  }

  Y <- as.vector(data[, 1])
  d <- ncol(data) - 1
  n <- length(Y)

  if(verbose){
    print(paste0("d = ", d, ", n = ", n))
  }

  # Check that s2 > s1 and reverse if not
  if (s1 > s2) {
    tmp <- s1
    s1 <- s2
    s2 <- tmp
  }

  if(verbose)
  {
    print(paste0("s1 = ", s1, ", s2 = ", s2))
  }
  # Calculate weights for two-scale DNN
  w1 <- (1 - (s1 / s2)^(-2 / d))^(-1)
  w2 <- 1 - w1

  if (verbose) {
    print(paste0("w1 = ", w1, ", w2 = ", w2))
  }

  # Calculate the two DNN estimators
  res1 <- 0
  res2 <- 0
  factor1 <- 1
  factor2 <- 1

  for (i in 1:(s2 - s1)) {
    index <- n - s1 + 2 - i
    res1 <- res1 + factor1 * Y[index]
    factor1 <- factor1 * ((n - index + 1) / i)
  }

  for (i in 1:(n - s2 + 1)) {
    index <- n - s2 + 2 - i

    res1 <- res1 + factor1 * Y[index]
    res2 <- res2 + factor2 * Y[index]

    factor1 <- factor1 * ((n - index + 1) / (i + s2 - s1))
    factor2 <- factor2 * ((n - index + 1) / i)
  }

  if(verbose){
    print(paste0("res1 = ", res1, ", res2 = ", res2))
  }

  # At this point the factors are (n-1 choose s-1)
  # so we ca simplify for the prefactors
  prefactor1 <- (factor1*(1 + n/s1))^(-1)
  prefactor2 <- (factor2*(1 + n/s2))^(-1)

  if (verbose) {
    print(paste0("prefactor1 = ", prefactor1, ", prefactor2 = ", prefactor2))
  }

  res1 <- prefactor1 * res1
  res2 <- prefactor2 * res2

  # Combine results
  res <- w1 * res1 + w2 * res2

  # return results
  return(res)
}
