#' Vector of reciprocal distances
#'
#' If there is a custom definicion of distances, it uses it. Otherwise, it
#' computes the reciprocal distances of points in x.
#'
#' If there exists a matrix object named
#' \code{.personal.definition.of.distances} in the global environment, the
#' function returns the lower triangular part of it as a vector, regardless of
#' \code{x}. Otherwise, the Euclidean distances between points in \code{x} are
#' returned as a vector.
#'
#' @param x a \code{data.frame} with planar coordinates
#' @export
vecdist <- function(x) {
  if(exists(".personal.definition.of.distances"))
  {
    M <- .personal.definition.of.distances
    return(M[lower.tri(M)])
  }
  else
  {
      return(as.vector(dist(x)))
  }
}
