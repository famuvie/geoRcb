#' Noise data
#'
#' Noise measurements, cartography and cost-based distance matrices from a
#' pilot-study in Valencia, Spain (see reference).
#'
#' @name noise
#' @aliases dd.distmat dl.distmat loc malilla obs
#' @references Antonio López-Quílez and Facundo Muñoz (2009). Geostatistical
#'   computing of acoustic maps in the presence of barriers. \emph{Mathematical
#'   and Computer Modelling} \strong{50}(5-6): 929–938.
#'   DOI:10.1016/j.mcm.2009.05.021
#' @docType data
#' @examples
#'    data(noise)
#'  if (require(sp)) {
#'    plot(malilla)
#'    points(obs, pch=19, col='red')
#'
#'    ## geodata structure with transformed covariates
#'    covarnames=sapply(1:3, function(x) paste("d2TV", x, sep=""))
#'    obs1.df <- as.data.frame(cbind(Leq=obs$Leq,
#'                             1/(1+(as.data.frame(obs)[covarnames]/20)^2)))
#'    obs.gd <- as.geodata(cbind(coordinates(obs), obs1.df),
#'                         data.col="Leq",
#'                         covar.col=c('d2TV1','d2TV2','d2TV3'))
#'    plot(obs.gd)
#'  }
NULL

#' @name dd.distmat
#' @docType data
#' @describeIn noise cost-based distance matrix from observations to observations
