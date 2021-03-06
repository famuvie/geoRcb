% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/likGRF.R
\name{likfit}
\alias{likfit}
\title{Fit a cost-based variogram model}
\usage{
likfit(geodata, coords = geodata$coords, data = geodata$data,
  trend = "cte", ini.cov.pars, fix.nugget = FALSE, nugget = 0,
  fix.kappa = TRUE, kappa = 0.5, fix.lambda = TRUE, lambda = 1,
  fix.psiA = TRUE, psiA = 0, fix.psiR = TRUE, psiR = 1, cov.model,
  realisations, lik.method = "ML", components = TRUE, nospatial = TRUE,
  limits = pars.limits(), dists.mat, print.pars = FALSE, messages, ...)
}
\arguments{
\item{geodata}{a list containing elements \code{coords} and
    \code{data} as described next. 
    Typically an object of the class \code{"geodata"}.
        If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }

\item{coords}{an \eqn{n \times 2}{n x 2} matrix where each row has the 2-D
    coordinates of the \eqn{n} data locations.
    By default it takes the
    component \code{coords} of the argument \code{geodata}, if provided.
  }

\item{data}{a vector with \emph{n} data values. By default it takes the
    component \code{data} of the argument \code{geodata}, if provided.  }

\item{trend}{specifies the mean part of the model. See documentation
    of \code{\link{trend.spatial}} for further details.
    Defaults to \code{"cte"}.  }

\item{ini.cov.pars}{initial values for the covariance parameters:
    \eqn{\sigma^2}{sigma^2} (partial sill) and \eqn{\phi}{phi} (range
    parameter).  Typically a vector with two components. However a
    matrix can be used to provide several initial values. See
    \code{DETAILS} below.}

\item{fix.nugget}{logical, indicating whether the parameter
    \eqn{\tau^2}{tau^2} (nugget variance) should be regarded as fixed
    (\code{fix.nugget = TRUE}) or should be
    estimated (\code{fix.nugget = FALSE}). Defaults to
    \code{FALSE}.  }

\item{nugget}{value of the nugget parameter.
    Regarded as a fixed value if \code{fix.nugget = TRUE} otherwise
    as the initial value for the minimisation algorithm. Defaults to zero.  }

\item{fix.kappa}{logical, indicating whether the extra parameter
    \eqn{\kappa}{kappa} should be regarded as fixed
    (\code{fix.kappa = TRUE}) or should be 
    estimated (\code{fix.kappa = FALSE}). Defaults to \code{TRUE}.  }

\item{kappa}{value of the extra parameter \eqn{\kappa}{kappa}.
    Regarded as a fixed value if \code{fix.kappa = TRUE}
    otherwise as the initial value for the minimisation algorithm. Defaults to
    \eqn{0.5}. This parameter is valid only if the covariance function is one
    of: \code{"matern"}, \code{"powered.exponential"}, \code{"cauchy"} or
    \code{"gneiting.matern"}.  For more details on covariance functions
    see documentation for \code{\link{cov.spatial}}.  }

\item{fix.lambda}{logical, indicating whether the Box-Cox transformation parameter
    \eqn{\lambda}{lambda} should be regarded as fixed
    (\code{fix.lambda = TRUE}) or should be
    be estimated (\code{fix.lambda = FALSE}). Defaults to \code{TRUE}.  }

\item{lambda}{value of the Box-Cox transformation parameter
    \eqn{\lambda}{lambda}.
    Regarded as a fixed value if \code{fix.lambda = TRUE} otherwise 
    as the initial value for the minimisation algorithm. Defaults to
    \eqn{1}. Two particular cases are \eqn{\lambda = 1}{lambda = 1}
    indicating no transformation and \eqn{\lambda = 0}{lambda
      = 0} indicating log-transformation.  }

\item{fix.psiA}{logical, indicating whether the anisotropy angle parameter
    \eqn{\psi_R}{psi_R} should be regarded as fixed
    (\code{fix.psiA = TRUE}) or should 
    be estimated (\code{fix.psiA = FALSE}). Defaults to
    \code{TRUE}.  }

\item{psiA}{value (in radians) for the anisotropy angle parameter
    \eqn{\psi_A}{psi_A}. 
    Regarded as a fixed value if \code{fix.psiA = TRUE}
    otherwise  as the initial value for the minimisation algorithm.
    Defaults to \eqn{0}. See \code{\link{coords.aniso}} for further
    details on anisotropy correction.  }

\item{fix.psiR}{logical, indicating whether the anisotropy ratio parameter
    \eqn{\psi_R}{psi_R} should be regarded as fixed
    (\code{fix.psiR = TRUE}) or should be estimated (\code{fix.psiR = FALSE}). Defaults to
    \code{TRUE}.  }

\item{psiR}{value, always greater than 1, for the anisotropy ratio  parameter
    \eqn{\psi_R}{psi_R}. 
    Regarded as a fixed value if \code{fix.psiR = TRUE}
    otherwise as the initial value for the minimisation algorithm.
    Defaults to \eqn{1}. See \code{\link{coords.aniso}} for further
    details on anisotropy correction.  }

\item{cov.model}{a string specifying the model for the correlation
    function. For further details see documentation for \code{\link{cov.spatial}}.
    Reads values from an \code{variomodel} object passed to \code{ini.cov.pars} if any, otherwise
    defaults to the \emph{exponential} model.  }

\item{realisations}{optional. Logical or a vector indicating the number of replication
    for each datum. For further information see \code{DETAILS} below and
    documentation for \code{\link{as.geodata}}. }

\item{lik.method}{(formely method.lik) options are \code{"ML"} for maximum likelihood and \code{"REML"} for
    restricted maximum likelihood. Defaults to \code{"ML"}.  }

\item{components}{an \eqn{n \times 3}{n x 3} data-frame with fitted
    values for the three model components: trend, spatial and residuals.
    See the section \code{DETAILS} below for the model specification.  }

\item{nospatial}{logical. If \code{TRUE} parameter estimates for the
    model without spatial component are included in the output.  }

\item{limits}{values defining lower and upper limits for the model
    parameters used in the numerical minimisation.
    The auxiliary function \code{\link{pars.limits}} is called to set
    the limits.
    See also \bold{Limits} in DETAILS below.}

\item{dists.mat}{n x n symmetric matrix with cost-based distances between
observations}

\item{print.pars}{logical. If \code{TRUE} the parameters and the value
  of the negative log-likelihood (up to a constant) are printed each
  time the function to be minimised is called.  }

\item{messages}{logical. Indicates whether
    status messages should be printed on the screen (or output device)
    while the function is running.  }

\item{...}{additional parameters to be passed to the minimisation
    function. Typically  arguments of the type  \code{control()} which controls the
    behavior of the minimisation algorithm. For further details see documentation 
    for the minimisation function \code{\link{optim}}.  }
}
\description{
All the arguments work as in \code{\link[geoR]{likfit}}, except the
additional argument \code{dists.mat}, which takes a symmetric matrix of
distances between observation locations
}
\examples{
## geodata structure with transformed covariates
data(noise)
if (require(sp)) {
  covarnames=sapply(1:3, function(x) paste("d2TV", x, sep=""))
  obs.df <- data.frame(Leq=obs$Leq,
                       1/(1+(as.data.frame(obs)[covarnames]/20)^2))
  obs.gd <- as.geodata(cbind(coordinates(obs), obs.df),
                       data.col="Leq",
                       covar.col=c('d2TV1','d2TV2','d2TV3'))
  trend=~d2TV1*(d2TV2+d2TV3)

  ## compute euclidean and cost-based empirical variograms
  vg.std <- variog(obs.gd, trend=trend)
  vg.dmat <- variog(obs.gd, trend=trend, dists.mat=dd.distmat)

  ## fitting variogram models
  vgmdl.std  <- likfit(geodata = obs.gd, trend=trend,
                       ini = c(8,300), cov.model = "matern")
  vgmdl.dmat <- likfit(geodata = obs.gd, trend=trend,
                       ini = c(8,300), cov.model = "matern",
                       dists.mat=dd.distmat)

  ## Fitted parameters
  data.frame(
    parameters=c("tausq","sigmasq","phi"),
    Euclidean=c(round(vgmdl.std$tausq,2),round(vgmdl.std$sigmasq,2),round(vgmdl.std$phi,0)),
    Cost_based=c(round(vgmdl.dmat$tausq,2),round(vgmdl.dmat$sigmasq,2),round(vgmdl.dmat$phi,0)))

  ## practical range
  ## defined as the value for which the correlation function
  ## decays to 5\% of its value at 0
  x=seq(0,800)
  y=cov.spatial(x,cov.pars=vgmdl.std$cov.pars)
  min(x[y<0.05*y[1]])    # 358
  y=cov.spatial(x,cov.pars=vgmdl.dmat$cov.pars)
  min(x[y<0.05*y[1]])    # 502
  # Note that the cost-based  analysis detects a
  # longer-ranged correlation structure

  ## plotting and comparing empirical variograms
  ## (with classical and robust estimation)
  ## and fitted variogram models
  op <- par(mfrow=c(2,2))
  u = 13  # binning
  for( est in c('classical','modulus') )
  {
    vg.std <- variog(obs.gd, trend=trend, estimator.type=est, uvec=u)
    vg.dmat <- variog(obs.gd, trend=trend, dists.mat=dd.distmat, estimator.type=est, uvec=u)
    plot(vg.std,pch=20, cex=1.2, max.dist=max(vg.dmat$u), ylim=c(0,20), col='gray')
    lines(vg.std,pch=20, col='gray')
    lines(vg.dmat,pch=20,cex=2,col='red')
    legend('topleft',c('Euclidean','Cost'),lty=1,lwd=2,col=c('gray','red'))
    title(paste('binning:',u,'   estimator:',est))

    plot(vg.std, pch=20, cex=1.2, max.dist=800, col='gray')   # empirical standard
    lines(vgmdl.std,lwd=2,col='gray')         # fitted model standard
    points(as.data.frame(vg.dmat[c(1,2)]),pch=20,cex=2,col='red')   # empirical cost-based
    lines(vgmdl.dmat,lwd=2,col='darkred')             # fitted model cost-based
    legend('topleft',c('Euclidean','Cost'),lty=1,lwd=2,col=c('gray','red'))
    title(paste('binning:',u,'   estimator:',est))
  }
  par(op)
}
}

