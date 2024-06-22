#' Check the provided data for correct format etc
#'
#' @param data The data set containing the observations.
#' @return TRUE or FALSE
data_check <- function(data) {
  if (!is.matrix(data)) {
    stop("Data needs to be supplied as a matrix.")
  }
  if (any(is.na(data))) {
    stop("Data matrix cannot contain NA.")
  }
  if (!all(is.numeric(data))) {
    stop("Data matrix can only contain real numbers.")
  }
  return(TRUE)
}

#' Check the provided point of interest for correct format etc
#'
#' @param x The point of interest
#' @return TRUE or FALSE
point_check <- function(x) {
  if (!all(is.numeric(x))) {
    stop("Point of interest must be given as coordinate vector.")
  }
  return(TRUE)
}

#' Check the provided data and point of interest for compatibility
#'
#' @param x The point of interest
#' @param data The data set containing the observations.
#' @return TRUE or FALSE
data_point_compatibility_check <- function(x, data) {
  if (length(x) != (ncol(data) - 1)) {
    stop("Point of interest and provided data are of different dimensionality.")
  }
  return(TRUE)
}

#' Check the provided data and sampling scale of interest for compatibility
#'
#' @param s The subsampling scale
#' @param data The data set containing the observations.
#' @return TRUE or FALSE
samplingscale_check <- function(s, data) {
  if (1 == 0) {
    stop("Error 1 occured.")
  }
  return(TRUE)
}

#' Check the provided data and subsampling scales for compatibility
#'
#' @param s1 The first subsampling scale
#' @param s2 The second subsampling scale
#' @param data The data set containing the observations.
#' @return TRUE or FALSE
samplingscale_check2 <- function(s1, s2, data) {
  n <- nrow(data)
  if (!(all(c(s1,s2)%%1==0))) {
    stop("s1 and s2 need to be integers.")
  }
  if (s1 > n) {
    stop("s1 exceeds sample size.")
  }
  if (s2 > n) {
    stop("s2 exceeds sample size.")
  }
  return(TRUE)
}
