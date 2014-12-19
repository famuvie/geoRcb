## Only modified to prevent computing the reciprocal distances between
## coordinates, in case a custom distance matrix is in use.
#' Variance-covariance
"varcov.spatial" <-
  function(coords = NULL, dists.lowertri = NULL, cov.model = "matern",
           kappa = 0.5, nugget = 0, cov.pars = stop("no cov.pars argument"),
           inv = FALSE, det = FALSE,
           func.inv = c("cholesky", "eigen", "svd", "solve"),
           scaled = FALSE, only.decomposition = FALSE,
           sqrt.inv = FALSE, try.another.decomposition = TRUE,
           only.inv.lower.diag = FALSE, ...)
{
  func.inv <- match.arg(func.inv)
  cov.model <- sapply(cov.model, match.arg, choices = .geoR.cov.models)
  if(only.inv.lower.diag)  inv <- TRUE
  if(is.null(coords) & is.null(dists.lowertri))
    stop("one of the arguments, coords or dists.lowertri must be provided")
  if (!is.null(coords) & !is.null(dists.lowertri))
    stop("only ONE argument, either coords or dists.lowertri must be provided")
  if (!is.null(coords))  n <- nrow(coords)
  if (!is.null(dists.lowertri))
    n <- as.integer(round(0.5 * (1 + sqrt(1 + 8 * length(dists.lowertri)))))
  tausq <- nugget
  if (is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
  }
  else {
    sigmasq <- cov.pars[, 1]
    phi <- cov.pars[, 2]
  }
##  print(c(tausq=tausq, sigmasq=sigmasq, phi=phi, kappa=kappa))
  #if (!is.null(coords)) dists.lowertri <- as.vector(dist(coords))  # (GUENMAP - original uncommented)
  if (!is.null(coords)) dists.lowertri <- vecdist(coords)           # (GUENMAP)
  if (round(1e+12 * min(dists.lowertri)) == 0)
    warning("Two or more pairs of data at coincident (or very close) locations. \nThis may cause crashes in some matrices operations.\n")
  varcov <- matrix(0, n, n)
  if (scaled) {
    if (all(phi < 1e-12))
      varcov <- diag(x = (1 + (tausq/sum(sigmasq))), n)
    else {
      if (is.vector(cov.pars)) cov.pars.sc <- c(1, phi)
      else cov.pars.sc <- cbind(1, phi)
      covvec <- cov.spatial(obj = dists.lowertri, cov.model = cov.model,
                            kappa = kappa, cov.pars = cov.pars.sc)
      varcov[lower.tri(varcov)] <- covvec
      varcov <- t(varcov)
      varcov[lower.tri(varcov)] <- covvec
      remove("covvec")
      if(sum(sigmasq) < 1e-16) diag(varcov) <- 1
      else diag(varcov) <- 1 + (tausq/sum(sigmasq))
    }
  }
  else {
    if (all(sigmasq < 1e-10) | all(phi < 1e-10)) {
      varcov <- diag(x = (tausq + sum(sigmasq)), n)
     }
    else {
      covvec <- cov.spatial(obj = dists.lowertri, cov.model = cov.model,
                            kappa = kappa, cov.pars = cov.pars)
      varcov[lower.tri(varcov)] <- covvec
      varcov <- t(varcov)
      varcov[lower.tri(varcov)] <- covvec
      remove("covvec")
      diag(varcov) <- tausq + sum(sigmasq)
    }
  }
  if (inv | det | only.decomposition | sqrt.inv | only.inv.lower.diag) {
    if (func.inv == "cholesky") {
        varcov.sqrt <- try(chol(varcov), silent=TRUE)
      if (inherits(varcov.sqrt, "try-error")) {
        if (try.another.decomposition){
          cat("trying another decomposition (svd)\n")
          func.inv <- "svd"
        }
        else {
          print(varcov.sqrt[1])
          stop()
        }
      }
      else {
        if (only.decomposition | inv) remove("varcov")
        if (!only.decomposition) {
          if (det) cov.logdeth <- sum(log(diag(varcov.sqrt)))
          if (sqrt.inv) inverse.sqrt <- solve(varcov.sqrt)
          if (inv) {
            invcov <- chol2inv(varcov.sqrt)
            if (!sqrt.inv) remove("varcov.sqrt")
          }
        }
      }
    }
    if (func.inv == "svd") {
      varcov.svd <- svd(varcov, nv = 0)
      cov.logdeth <- try(sum(log(sqrt(varcov.svd$d))), silent=TRUE)
      if (inherits(cov.logdeth, "try-error")) {
        if (try.another.decomposition){
          cat("trying another decomposition (eigen)\n")
          func.inv <- "eigen"
        }
        else {
          print(cov.logdeth[1])
          stop()
        }
      }
      else {
        if (only.decomposition | inv) remove("varcov")
        if (only.decomposition)
          varcov.sqrt <- crossprod(t(varcov.svd$u) * sqrt(sqrt(varcov.svd$d)))
        if (inv) {
          invcov <- crossprod(t(varcov.svd$u)/sqrt(varcov.svd$d))
        }
        if (sqrt.inv)
          inverse.sqrt <- crossprod(t(varcov.svd$u)/sqrt(sqrt(varcov.svd$d)))
      }
    }
    if (func.inv == "solve") {
      if (det)
        stop("the option func.inv == \"solve\" does not allow computation of determinants. \nUse func.inv = \"chol\",\"svd\" or \"eigen\"\n")
      invcov <- try(solve(varcov), silent=TRUE)
      if (inherits(cov.logdeth, "try-error")) {
        if (try.another.decomposition)
          func.inv <- "eigen"
        else {
          print(invcov[1])
          stop()
        }
      }
      remove("varcov")
    }
    if (func.inv == "eigen") {
      varcov.eig <- try(eigen(varcov, symmetric = TRUE), silent=TRUE)
      cov.logdeth <- try(sum(log(sqrt(varcov.eig$val))), silent=TRUE)
      if (inherits(cov.logdeth, "try.error") | inherits(varcov.eig, "try-error")) {
        diag(varcov) <- 1.0001 * diag(varcov)
        varcov.eig <- try(eigen(varcov, symmetric = TRUE), silent=TRUE)
        cov.logdeth <- try(sum(log(sqrt(varcov.eig$val))), silent=TRUE)
        if (inherits(cov.logdeth, "try.error") | inherits(varcov.eig, "try-error")) {
          return(list(crash.parms = c(tausq=tausq, sigmasq=sigmasq, phi=phi, kappa=kappa)))
        }
      }
      else {
        if (only.decomposition | inv) remove("varcov")
        if (only.decomposition)
          varcov.sqrt <- crossprod(t(varcov.eig$vec)* sqrt(sqrt(varcov.eig$val)))
        if (inv) invcov <- crossprod(t(varcov.eig$vec)/sqrt(varcov.eig$val))
        if (sqrt.inv)
          inverse.sqrt <- crossprod(t(varcov.eig$vec)/sqrt(sqrt(varcov.eig$val)))
      }
    }
  }
  if (!only.decomposition) {
    if (det) {
      if (inv) {
        if (only.inv.lower.diag)
          result <- list(lower.inverse = invcov[lower.tri(invcov)],
                         diag.inverse = diag(invcov), log.det.to.half = cov.logdeth)
        else result <- list(inverse = invcov, log.det.to.half = cov.logdeth)
      }
      else {
        result <- list(varcov = varcov, log.det.to.half = cov.logdeth)
      }
      if (sqrt.inv)
        result$sqrt.inverse <- inverse.sqrt
    }
    else {
      if (inv) {
        if (only.inv.lower.diag)
          result <- list(lower.inverse = invcov[lower.tri(invcov)],
                         diag.inverse = diag(invcov))
        else {
          if (sqrt.inv)
            result <- list(inverse = invcov, sqrt.inverse = inverse.sqrt)
          else result <- list(inverse = invcov)
        }
      }
      else result <- list(varcov = varcov)
    }
  }
  else result <- list(sqrt.varcov = varcov.sqrt)
  result$crash.parms <- NULL
  return(result)
}


