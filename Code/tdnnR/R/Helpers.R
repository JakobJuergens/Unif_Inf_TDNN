#' Calculate the distance between two points of interest
#'
#' @param x Point A of interest
#' @param y Point B of interest
#' @return The Euclidean distance between point A and B
dist <- function(x, y){
  sqrt(sum((x-y)^2))
}

#' Sorts a sample according to its distance relative
#' to a point of interest
#'
#' @param x Point of interest
#' @param data Sample under Consideration
#' @return The sample sorted relative to x
pre_sort <- function(x, data){

}
