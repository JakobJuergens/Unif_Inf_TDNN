#' Calculate the Jackknife Variance Estimator for the
#' distributional nearest neighbor estimate for a point x
#' based on a sample data with subsampling scale s.
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
#' approximation for the weights (default value = TRUE)
#' @return A number.
#'
#' @export

DNN_JK_var <- function(x, data, s,
                       presorted = FALSE, standardized = FALSE,
                       asymp_approx_weights = TRUE) {
  # Standardize Data
  if (standardize == TRUE) {
    stand <- NPR_standardize(x = x, data = data)
    x <- stand$x
    data <- stand$data
  }

  d <- ncol(data) - 1
  # Pre-sort data relative to point of interest
  if (presorted == FALSE) {
    order <- pre_sort(x = x, cov_mat = data[, -1])
    data <- data[order, ]
  }

  Y <- as.vector(data[, 1])
  n <- length(Y)

  # Calculate JK Variance estimator
  res <- 0
  # using exact weights
  if (asymp_approx_weights == FALSE) {
    res <- NA
  }
  # using approximate weights
  if (asymp_approx_weights == TRUE) {
    alpha <- s / n
    alpha_tilde <- s / (n - 1)

    for(i in 1:(n - s + 1)) {

    }

    for(i in (n-s+2):n){

    }

    res <- (n - 1) / n * res
  }
  return(NA)
}
