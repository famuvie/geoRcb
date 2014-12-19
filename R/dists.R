#' Vector of reciprocal distances
#'
#' If there is a custom definicion of distances, it uses it.
#' Otherwise, it computes the reciprocal distances of points in x.
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
