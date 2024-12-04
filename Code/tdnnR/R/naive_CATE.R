#' LOREM IPSUM
#'
#' @param x The point of interest.
#' @param data The data set containing the observations.
#' As a matrix containing the response value in the first column, the
#' treatment status in the second and each covariate in subsequent columns
#' @param s The subsampling scale.
#' @param presorted True or False whether the data is sorted according to its
#' distance to the point of interest (default value = FALSE)
#' @param standardize True or False whether to standardize (default value = FALSE)
#' @param asymp_approx_weights True or False whether to use an asymptotic
#' approximation for the weights (default value = FALSE)
#' @return A number.
#'
#' @export
naive_CATE_DNN <- function(x, data, s,
                           presorted = FALSE, standardize = FALSE,
                           asymp_approx_weights = TRUE) {

  # Standardize and presort if desired
  if (standardize == TRUE) {
    stand <- standardize(x = x, data = data)
    x <- stand$x
    data <- stand$data
  }
  if (presorted == FALSE) {
    data <- pre_sort(x = x, data = data)
  }

  # Split data matrix according to treatment status
  untreated_data <- data[data[, 2] == 1, ]
  treated_data <- data[data[, 2] == 0, ]

  # Run DNN Estimation on each sub-data set separately
  untreated_reg_est <- DNN(
    x = x, data = untreated_data, s = s,
    presorted = TRUE, standardize = FALSE,
    asymp_approx_weights = asymp_approx_weights,
    verbose = FALSE
  )

  treated_reg_est <- DNN(
    x = x, data = treated_data, s = s,
    presorted = TRUE, standardize = FALSE,
    asymp_approx_weights = asymp_approx_weights,
    verbose = FALSE
  )

  # Return difference
  return(treated_reg_est - untreated_reg_est)
}

#' LOREM IPSUM
#'
#' @param x The point of interest.
#' @param data The data set containing the observations.
#' As a matrix containing the response value in the first column, the
#' treatment status in the second and each covariate in subsequent columns
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
naive_CATE_TDNN <- function(x, data, s1, s2,
                            presorted = FALSE, standardize = FALSE,
                            asymp_approx_weights = TRUE) {

  # Standardize and presort if desired
  if (standardize == TRUE) {
    stand <- standardize(x = x, data = data)
    x <- stand$x
    data <- stand$data
  }
  if (presorted == FALSE) {
    data <- pre_sort(x = x, data = data)
  }

  # Split data matrix according to treatment status
  untreated_data <- data[data[, 2] == 1, ]
  treated_data <- data[data[, 2] == 0, ]

  # Run DNN Estimation on each sub-data set separately
  untreated_reg_est <- TDNN(
    x = x, data = untreated_data, s_1 = s_1, s_2 = s_2,
    presorted = TRUE, standardize = FALSE,
    asymp_approx_weights = asymp_approx_weights,
    verbose = FALSE
  )

  treated_reg_est <- TDNN(
    x = x, data = treated_data, s_1 = s_1, s_2 = s_2,
    presorted = TRUE, standardize = FALSE,
    asymp_approx_weights = asymp_approx_weights,
    verbose = FALSE
  )

  # Return difference
  return(treated_reg_est - untreated_reg_est)
}
