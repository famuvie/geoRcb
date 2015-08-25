#' Compute cost-based distances
#'
#' Generate a cost-based distance matrix among observations and/or between
#' observations and prediction locations.
#'
#' Use the centroids of the cost surface raster cells as prediction locations.
#'
#' If a cost surface raster is not supplied a uniform (cost 1) grid will be
#' generated to match to bounding box of pts
#'
#' @export
#' @import gdistance
#' @param pts A SpatialPoints[DataFrame] with a defined projection or a
#'   two-column data.frame/matrix of coordinates
#' @param costsurf cost surface Raster
#' @param ret specify whether to return distances between obs-obs ("obs"), obs-loc
#'   ("loc"), or both ("both") (default)
#' @param directions. See \code{\link[gdistance]{transition}}.
#' @param silent logical. If TRUE avoids any warnings or messages.
#' @examples
#'  data(noise)
#'  r <- raster(nrows=40,ncols=40,resolution=1)
#'  r <- setExtent(r,extent(obs),keepres=T)
#'  aggfac <- 5
#'  costsurf <- rasterize(malilla,r,background=aggfac,field=rep(10000,length(polygons(malilla))))
#'  costsurf <- reclassify(costsurf,c(aggfac+1,Inf,NA))
#'  costsurf <- aggregate(costsurf,aggfac,expand=T,fun=max,na.rm=T)
#'
#'  res <- distmatGen(obs,costsurf,ret="obs")
#'  plot(dd.distmat[,51],res[,51],ylim=c(0,600),xlim=c(0,600))
distmatGen <- function(pts,
                       costsurf,
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

  ret <- match.arg(ret)

  ## Rename cost surface (?)
  r <- costsurf

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
idx_isol <- function(x) {
  is.zero_or_inf <- function(x) all(is.infinite(x) | x == 0)
  which(apply(x, 2, is.zero_or_inf))
}
