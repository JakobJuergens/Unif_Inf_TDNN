#' Calculate the distance between two points of interest
#'
#' @param x Point A of interest
#' @param y Point B of interest
#' @return The Euclidean distance between point A and B
dist <- function(x, y) {
  sqrt(sum((x - y)^2))
}

#' Sorts a sample according to its distance relative
#' to a point of interest (for NPR Setup)
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
NPR_pre_sort <- function(x, data) {
  n <- nrow(data)
  dists <- sqrt(rowSums((sweep(x = data[ , -1], MARGIN = 2, FUN = '-', x)^2)))
  dist_ranks <- order(dists)
  return(data[dist_ranks, ])
}

#' Sorts a sample according to its distance relative
#' to a point of interest (for CATE Setup)
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
CATE_pre_sort <- function(x, data) {
  n <- nrow(data)
  dists <- sqrt(rowSums((sweep(x = data[ , -c(1,2)], MARGIN = 2, FUN = '-', x)^2)))
  dist_ranks <- order(dists)
  return(data[dist_ranks, ])
}

#' Standardizes a point of interest and a sample
#' to improve comparability between dimensions (for NPR Setup)
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
NPR_standardize <- function(x, data) {
  d <- ncol(data) - 1
  for (i in 1:d) {
    sd <- sd(data[, i + 1])
    data[, i + 1] <- data[, i + 1] / sd
    x[i] <- x[i] / sd
  }
  return(list(x = x, data = data))
}

#' Standardizes a point of interest and a sample
#' to improve comparability between dimensions (for CATE Setup)
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
CATE_standardize <- function(x, data) {
  d <- ncol(data) - 2
  for (i in 2:d) {
    sd <- sd(data[, i + 1])
    data[, i + 1] <- data[, i + 1] / sd
    x[i] <- x[i] / sd
  }
  return(list(x = x, data = data))
}
