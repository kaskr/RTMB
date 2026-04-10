## Helper: sdreport that works for all integration methods
sdreport_xtra <- function(obj,
                          sdr, ## Append to this object
                          ignore.parm.uncertainty = FALSE,
                          getReportCovariance = TRUE,
                          type = c("mean", "mode"),
                          what = c("reportvector", "raneffvector"),
                          ...) {
  applicable <- obj$env$intern || length(obj$env$integrate)
  if (!applicable) {
    return(list())
  }
  if (!is.null(obj$env$random)) {
    warning("sdreport: Cannot mix inner and outer integration methods (yet)")
    return(list())
  }
  type <- match.arg(type)
  what <- match.arg(what)
  if (what == "reportvector" && length(sdr$env$ADreportDims)==0) {
    ## Nothing to report - drop out
    return(list())
  }
  getReportCovariance <- getReportCovariance && (what == "reportvector")
  ans <- list()
  par.fixed <- sdr$par.fixed
  cov.fixed <- sdr$cov.fixed
  ## Get parameter list
  pl <- lapply(sdr$env$parameters, function(x) attr(x, "shape") %||% x)
  data <- obj$env$data
  func <- attr(data, "func")
  ## Get number vars to adreport
  ## Or the number of parameters
  if (what == "reportvector") {
    tmp <- MakeADFun(data, pl, ADreport=TRUE, type="ADFun")
    lgts <- unlist(lapply(tmp$env$ADreportDims, prod))
    n <- sum(lgts)
    statistic <- GetTape(tmp)
    nam <- rep(names(lgts), lgts)
  } else {
    random <- which(rep( names(pl), sapply(pl, length) ) %in% obj$env$.random)
    n <- length(random)
    statistic <- function(parvec) parvec[random]
    nam <- rep(names(pl), lengths(pl))[random]
  }
  ## For epsilon method
  eps <- rep(0, n)
  p <- c(eps, par.fixed)
  getObj <- function(moment=1, scale=1) {
    g <- function(par) {
      eps <- par[[1]]
      par <- par[-1]
      parvec <- do.call("c", par)
      func(par) * scale + sum(eps * statistic(parvec)^moment)
    }
    MakeADFun(g,
              c(list(eps), pl),
              random=obj$env$.random,
              integrate = obj$env$integrate,
              intern=obj$env$intern,
              map=obj$env$map,
              silent=TRUE)
  }
  ## First moment object
  obj1 <- getObj(moment=1, scale=1)
  ## Some helper objects
  all <- seq_len(length(obj1$par))
  lpar <- tail(all, length(par.fixed))
  leps <- head(all, n)
  ## Mean always needed
  a_mean <- head(as.vector(obj1$gr(p)), n)
  a_mode <- NULL
  ## Full cov requested ?
  cov <- matrix(NA)
  if (getReportCovariance) {
    obj1$env$retape_adgrad()
    ## V(phi(u,x)|x) = - D2_eps (f)
    B <- -obj1$env$f(theta=p, type="ADGrad", order=1, keepx=leps, keepy=leps)
  }
  ## If type=="mode" overwrite 'obj1'
  if (type == "mode") {
    obj1 <- getObj(moment=1, scale=1e6)
    a_mode <- head(as.vector(obj1$gr(p)), n)
  }
  ## Get (transposed) jacobian of *either* a_mean or a_mode:
  ## Note: interchanging 'keepx' and 'keepy' gives same result but usually slower
  if (!ignore.parm.uncertainty) {
    obj1$env$retape_adgrad()
    JT <- obj1$env$f(theta=p, type="ADGrad", order=1, keepx=leps, keepy=lpar)
  } else {
    JT <- matrix(0, length(lpar), length(leps))
  }
  ## Second moment object
  ## (FIXME: Could be avoided if we could calculate hessian diagonals wrt epsilon)
  obj2 <- getObj(moment=2, scale=1)
  b <- head(as.vector(obj2$gr(p)), n)
  ## NOTE: zero matrix if isTRUE(ignore.theta.uncertainty)
  Vtheta <- cov.fixed
  ## Full covariance
  if (getReportCovariance) {
    cov <- t(JT) %*% Vtheta %*% JT + B
  }
  ## Put pieces together
  diag.term1 <- colSums( JT * (Vtheta %*% JT) )
  diag.term2 <- b-a_mean^2
  ans$value <- if (type == "mean") a_mean else a_mode
  names(ans$value) <- nam
  ans$sd <- sqrt(diag.term1 + diag.term2)
  ans$cov <- cov
  ans
}
