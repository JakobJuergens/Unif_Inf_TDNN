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
DNN <- function(x, data, s){
  # Check inputs for executing function
  point_check(x)
  data_check(data)
  data_point_compatibility_check(x, data)

  # Pre-sort data relative to point of interest
  data_sorted <- pre_sort(x = x, data = data)
  n <- nrow(data_sorted)
  d <- ncol(data_sorted) - 1

  # Calculate DNN estimator
  res <- 0
  factor <- 1
  prefactor <- 1/choose(n,s)
  for(i in 1:(n-s+1)){
    index <- n-s+2-i
    res <- res + factor*data_sorted[index,1]
    factor <- factor*((n-index+1)/(n - index - s + 2))
  }
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
TDNN <- function(x, data, s1, s2){
  # Check inputs for executing function
  point_check(x)
  data_check(data)
  data_point_compatibility_check(x, data)

  # Pre-sort data relative to point of interest
  data_sorted <- pre_sort(x = x, data = data)

  # Calculate weights for two-scale DNN
  w1 <- 0.5
  w2 <- 1 - w_1
}


