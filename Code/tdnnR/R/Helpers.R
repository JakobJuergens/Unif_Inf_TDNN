#' Calculate the distance between two points of interest
#'
#' @param x Point A of interest
#' @param y Point B of interest
#' @return The Euclidean distance between point A and B
dist <- function(x, y) {
  sqrt(sum((x - y)^2))
}

#' Sorts a sample according to its distance relative
#' to a point of interest
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
pre_sort <- function(x, data) {
  n <- nrow(data)
  dist_ranks <- order(unlist(purrr::map(
    .x = 1:n,
    .f = ~ philentropy::euclidean(x, data[.x, -1], testNA = FALSE)
  )))
  return(data[dist_ranks, ])
}

#' Standardizes a point of interest and a sample
#' to improve comparability between dimensions
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
standardize <- function(x, data) {
  d <- ncol(data) - 1
  for (i in 1:d) {
    sd <- sd(data[, i + 1])
    data[, i + 1] <- data[, i + 1] / sd
    x[i] <- x[i] / sd
  }
  return(list(x = x, data = data))
}
