#' Default transition function in geoR
#'
#' Default function to use with \code{\link[gdistance]{transition}}.
#'
#' To use another function, override this one by defining it using the same
#' name.
#'
#' @param x numeric vector of length 2 with cost/conductivity non-negative
#'   numbers.
#'
#' @return cost/conductivity value of the transition between neighbouring cells.
#'
#' @examples
#'    geoR_trf(c(0, 1))
#'
#'    # redefine function
#'    geoR_trf <- mean
#'
#'    geoR_trf(c(0, 1))
geoR_trf <- function(x) {
  if (any(x == 0)) return(0)
  else return(mean(x))
}
