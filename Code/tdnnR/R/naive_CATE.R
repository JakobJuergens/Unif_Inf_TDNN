#' Naive estimator for conditional average treatment effects based on the DNN
#' estimator.
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
    stand <- CATE_standardize(x = x, data = data)
    x <- stand$x
    data <- stand$data
  }
  if (presorted == FALSE) {
    data <- CATE_pre_sort(x = x, data = data)$data
  }

  # Split data matrix according to treatment status
  untreated_data <- data[data[, 2] == 0, ]
  treated_data <- data[data[, 2] == 1, ]

  # Run DNN Estimation on each sub-data set separately
  untreated_reg_est <- DNN(
    x = x, data = untreated_data[,-2], s = s,
    presorted = TRUE, standardize = FALSE,
    asymp_approx_weights = asymp_approx_weights
  )

  treated_reg_est <- DNN(
    x = x, data = treated_data[,-2], s = s,
    presorted = TRUE, standardize = FALSE,
    asymp_approx_weights = asymp_approx_weights
  )

  # Return difference
  return(treated_reg_est - untreated_reg_est)
}

#' Naive estimator for conditional average treatment effects based on the TDNN
#' estimator.
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
    stand <- CATE_standardize(x = x, data = data)
    x <- stand$x
    data <- stand$data
  }
  if (presorted == FALSE) {
    data <- CATE_pre_sort(x = x, data = data)$data
  }

  # Split data matrix according to treatment status
  untreated_data <- data[data[, 2] == 0, ]
  treated_data <- data[data[, 2] == 1, ]

  # Run DNN Estimation on each sub-data set separately
  untreated_reg_est <- TDNN(
    x = x, data = untreated_data[,-2], s1 = s1, s2 = s2,
    presorted = TRUE, standardize = FALSE,
    asymp_approx_weights = asymp_approx_weights
  )$TDNN_res

  treated_reg_est <- TDNN(
    x = x, data = treated_data[,-2], s1 = s1, s2 = s2,
    presorted = TRUE, standardize = FALSE,
    asymp_approx_weights = asymp_approx_weights
  )$TDNN_res

  # Return difference
  return(treated_reg_est - untreated_reg_est)
}
