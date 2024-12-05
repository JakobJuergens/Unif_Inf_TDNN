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

  n <- nrow(data)

  # Create List of Fold indices to Estimate Propensity Scores and Regression Functions
  folds <- generate_folds(n_obs = nrow(data), n_folds = n_folds)
  nuisance_par_ests <- matrix(data = NA, nrow = nrow(data), ncol = 3)

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

  # Estimate Parameter of Interest using DNN estimator weights
  res <- 0
  factor <- 1
  prefactor <- 0

  # using exact weights
  if (asymp_approx_weights == FALSE) {
    res <- NA
  }

  # using approximate weights
  if (asymp_approx_weights == TRUE) {
    res <- NA
  }

  # return calculated result
  names(res) <- NULL
  return(res)
}
