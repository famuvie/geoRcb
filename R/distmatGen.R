#' Compute cost-based distances
#'
#' Generate a cost-based distance matrix among observations and/or between
#' observations and prediction locations.
#'
#' Use the centroids of the conductivity or conductivity surface raster cells as
#' prediction locations.
#'
#' @export
#' @import gdistance
#' @param pts A SpatialPoints[DataFrame] with a defined projection or a
#'   two-column data.frame/matrix of coordinates
#' @param condsurf raster. Conductivity surface.
#' @param ret specify whether to return distances between obs-obs ("obs"),
#'   obs-loc ("loc"), or both ("both") (default)
#' @param directions See \code{\link[gdistance]{transition}}.
#' @param silent logical. If TRUE avoids any warnings or messages.
#' @examples
#' data(noise)
#' if (require(raster)) {
#'   r <- raster(extent(c(0,3,0,3)), nrows = 3, ncols = 3)
#'   wall.idx <- 4:5
#'   values(r)[-wall.idx] <- 1
#'   obs <- coordinates(r)[-wall.idx, ]
#'
#'   ddm <- distmatGen(obs, r, ret = "obs", directions = 8)
#'
#'   par(mfrow = c(2, 1))
#'   plot(r)
#'   text(obs[, 1], obs[, 2],
#'        col = c('red', rep('black', 6)))
#'
#'   plot(dist(obs)[1:6], ddm[-1, 1],
#'        type = 'n',
#'        main = 'Distances to point 1',
#'        xlab = 'Euclidean distance',
#'        ylab = 'Cost-based distance')
#'   text(dist(obs)[1:6], ddm[-1, 1],
#'        labels = 2:7)
#'   abline(0, 1)
#' }
distmatGen <- function(pts,
                       condsurf,
                       ret = c('both', 'obs', 'loc'),
                       directions = 16,
                       silent = FALSE) {

  ## Derive fromcoords as a two-column matrix
  if (inherits(pts, 'matrix')) {
    fromcoords <- pts
  } else if (inherits(pts, 'data.frame')) {
    fromcoords <- as.matrix(pts)
  } else if (inherits(pts, 'SpatialPoints')) {
    fromcoords <- coordinates(pts)
  } else {
    stop(paste('pts must be either a two column matrix/data.frame',
               'or a SpatialPoints object.'))
  }
  stopifnot(ncol(fromcoords) == 2)

  r <- condsurf

  ret <- match.arg(ret)

  ## Extend surface if needed
  if (!extent(r) >= extent(pts)) {
    xn <- floor(extent(pts)[1]-1)
    xx <- ceiling(extent(pts)[2]+1)
    yn <- floor(extent(pts)[3]-1)
    yx <- ceiling(extent(pts)[4]+1)
    ext <- unlist(lapply(list(xn, xx, yn, yx), function(x) round(x, 0)))
    r <- extend(r, extent(ext))
  }

  ## Transition matrix with diagonal-correction
  tr <- transition(r, geoR_trf, directions)
  tr <- geoCorrection(tr, type = 'c', multpl = FALSE)

  ## Obs-Obs distance matrix
  if (ret == "both" || ret == "obs") {
    oret <- diag(0, nrow = nrow(fromcoords))
    oret.c <- costDistance(tr, fromcoords)
    oret[lower.tri(oret)] <- oret.c
    oret <- oret+t(oret)

    # Warn about isolated points
    if (!silent && length(idx <- idx_isol(oret)) > 0) {
      msg <- paste("The",
                   ifelse(length(idx) > 1, 'points', 'point'),
                   paste(idx, collapse = ', '),
                   ifelse(length(idx) > 1, 'look', 'looks'),
                   "isolated wrt other observations.\n",
                   "Please check.")
      warning(msg)
    }
  }

  if(ret == "both" || ret == "loc"){
    # tocoords <- xyFromCell(tr,which(values(r)==res(r)))
    tocoords <- coordinates(r)
    lret <- t(costDistance(tr, fromcoords, tocoords))

    # Warn about isolated points
    if (!silent && length(idx <- idx_isol(lret)) > 0) {
      msg <- paste("The",
                   ifelse(length(idx) > 1, 'points', 'point'),
                   paste(idx, collapse = ', '),
                   ifelse(length(idx) > 1, 'look', 'looks'),
                   "isolated wrt all prediction locations.\n",
                   "Please check.")
      warning(msg)
    }
  }

  ans <- switch(ret,
                both = list(obs = oret,
                            loc = lret),
                obs  = oret,
                loc  = lret)

  return(ans)
}

#' Identify isolated locations
#'
#' Defined as points whith cost distances to every other
#' point of either 0 or Infinite.
#'
#' param x distance matrix.
idx_isol <- function(x) {
  is.zero_or_inf <- function(x) all(is.infinite(x) | x == 0)
  which(apply(x, 2, is.zero_or_inf))
}
