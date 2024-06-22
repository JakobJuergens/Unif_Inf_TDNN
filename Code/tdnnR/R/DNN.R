#' Calculate the distributional nearest neighbor estimate for a point x
#' based on a sample data with subsampling scale s.
#' The implementation is based on the equivalent representation as an L-statistic
#' described in Steele (2009)
#'
#' @param x The point of interest.
#' @param data The data set containing the observations.
#' As a matrix containing the response value in the first column.
#' @param s The subsampling scale.
#' @return A number.
#'
#' @export
DNN <- function(x, data, s,
                presorted = FALSE, standardize = FALSE,
                verbose = FALSE) {
  # Check inputs for executing function
  point_check(x)
  data_check(data)
  data_point_compatibility_check(x, data)
  samplingscale_check(s, data)

  # Standardize Data
  if(standardize == TRUE){
    stand <- standardize(x = x, data = data)
    x <- stand$x
    data <- stand$data
  }

  # Pre-sort data relative to point of interest
  if (presorted == FALSE) {
    data <- pre_sort(x = x, data = data)
  }

  Y <- as.vector(data[,1])
  d <- ncol(data) - 1
  n <- length(Y)

  if (verbose) {
    print(paste0("d = ", d, ", n = ", n))
  }

  # Calculate DNN estimator
  res <- 0
  factor <- 1
  prefactor <- 1 / choose(n, s)
  for (i in 1:(n - s + 1)) {
    index <- n - s + 2 - i
    res <- res + factor * Y[index]
    factor <- factor * ((n - index + 1) / i)
  }
  res <- prefactor * res

  # return calculated result
  return(res)
}
