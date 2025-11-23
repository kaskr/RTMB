## Translation of QUADPACK routine
## Gauss-Kronrod K21/G10
rdqk21 <- function(f, a, b) {
  ## Kronrod nodes
  xk <- c(.995657163025808080735527280689003,
          .973906528517171720077964012084452,
          .930157491355708226001207180059508,
          .865063366688984510732096688423493,
          .780817726586416897063717578345042,
          .679409568299024406234327365114874,
          .562757134668604683339000099272694,
          .433395394129247190799265943165784,
          .294392862701460198131126603103866,
          .14887433898163121088482600112972,0.)
  xk <- c(-xk , rev(xk)[-1])
  ## Kronrod weights
  wk <- c(.011694638867371874278064396062192,
          .03255816230796472747881897245939,
          .05475589657435199603138130024458,
          .07503967481091995276704314091619,
          .093125454583697605535065465083366,
          .109387158802297641899210590325805,
          .123491976262065851077958109831074,
          .134709217311473325928054001771707,
          .142775938577060080797094273138717,
          .147739104901338491374841515972068,
          .149445554002916905664936468389821)
  wk <- c(wk , rev(wk)[-1])
  ## Gauss nodes (xk[ig])
  ig <- seq(from=2, by=2, length=10)
  ## Gauss weights
  wg <- c(.066671344308688137593568809893332,
          .149451349150580593145776339657697,
          .219086362515982043995534934228163,
          .269266719309996355091226921569469,
          .295524224714752870173892994651338)
  wg <- c(wg, rev(wg))
  ## Integrate
  x <- 0.5 * ((b - a) * xk + b + a)
  fx <- f(x)
  K <- sum(wk * fx) * (b - a)/2
  G <- sum(wg * fx[ig]) * (b - a)/2
  c( result = K, abserr = abs(K-G) )
}

## Translation of QUADPACK routine
## Gauss-Kronrod K15/G7
rdqk15 <- function(f, a, b) {
  ## Kronrod nodes
  xk <- c(.991455371120812639206854697526329,
          .949107912342758524526189684047851,
          .864864423359769072789712788640926,
          .741531185599394439863864773280788,
          .58608723546769113029414483825873,
          .405845151377397166906606412076961,
          .207784955007898467600689403773245, 0.)
  xk <- c(-xk , rev(xk)[-1])
  ## Kronrod weights
  wk <- c(.02293532201052922496373200805897,
          .063092092629978553290700663189204,
          .104790010322250183839876322541518,
          .140653259715525918745189590510238,
          .16900472663926790282658342659855,
          .190350578064785409913256402421014,
          .204432940075298892414161999234649,
          .209482141084727828012999174891714)
  wk <- c(wk , rev(wk)[-1])
  ## Gauss nodes (xk[ig])
  ig <- seq(from=2, by=2, length=7)
  ## Gauss weights
  wg <- c(.129484966168869693270611432679082,
          .27970539148927666790146777142378,
          .381830050505118944950369775488975,
          .417959183673469387755102040816327)
  wg <- c(wg, rev(wg)[-1])
  ## Integrate
  x <- 0.5 * ((b - a) * xk + b + a)
  fx <- f(x)
  K <- sum(wk * fx) * (b - a)/2
  G <- sum(wg * fx[ig]) * (b - a)/2
  c( result = K, abserr = abs(K-G) )
}

isConstant <- function(x) !getVariables(x)
Value <- function(x) getValues(x)
## f <- function(x)sin(exp(x))*exp(x)
## MakeTape(function(x)RTMB:::adintegrate(f,1*x,7*x),1)(1)
## F <- MakeTape(function(mu)RTMB:::adintegrate(function(x)dnorm(x,mu),-2,2),1)
dotfun <- function(f,...) {
  f <- match.fun(f)
  function(x) f(x,...)
}
adintegrate <- function(f, a, b, cfg) {
  if (!ad_context()) stop("No active AD context")
  a <- advector(a)
  b <- advector(b)
  ainf <- isConstant(a) && is.infinite(a)
  binf <- isConstant(b) && is.infinite(b)
  if (ainf && binf) {
    return(adintegrate(f, a, 0, cfg) + adintegrate(f, 0, b, cfg))
  } else if (ainf || binf) {
    ## Transformation:
    if (isConstant(b) && Value(b) == Inf)
      g <- function(u) { t <- 1-u ;  f(a + u / t) / (t*t) } ## [a, Inf]
    if (isConstant(b) && Value(b) == -Inf)
      g <- function(u) { t <- 1-u ; -f(a - u / t) / (t*t) } ## [a, -Inf]
    if (isConstant(a) && Value(a) == Inf)
      g <- function(u) { t <- 1-u ; -f(b + u / t) / (t*t) } ## [Inf, b]
    if (isConstant(a) && Value(a) == -Inf)
      g <- function(u) { t <- 1-u ;  f(b - u / t) / (t*t) } ## [-Inf, b]
    return(adintegrate(g, 0, 1, cfg))
  }
  old <- TapeConfig()
  TapeConfig(vectorize="enable")
  F <- MakeTape(function(ab) rdqk21(f, ab[1], ab[2]), numeric(2))
  TapeConfig(old)
  F$simplify() ## FIXME: Remove many RefOp
  decompose_refs(.pointer(environment(F)$mod))
  F$reorder() ## Move c(a,b) parameters up front
  ## NOTE: F may depend on parameters that have to be 'resolved'
  vars <- resolve_refs(.pointer(environment(F)$mod))
  F <- F$atomic()
  ptr <- .pointer(environment(F)$mod)
  bisect_atom(ptr, c(a, b, vars), cfg)
}

##' AD adaptive integration.
##'
##' Univariate adaptive integration closely following \link[stats]{integrate}.
##' The R (stats) version is used in standard (non AD) evaluation mode,
##' while the AD version uses a simplified re-implementation.
##'
##' @param f Integrand.
##' @param lower Lower integration limit.
##' @param upper Upper integration limit.
##' @param ... Passed to `f`.
##' @param subdivisions Max number of subdivisions.
##' @param rel.tol Relative tolerance (not used by the AD version)
##' @param abs.tol Absolute tolerance.
##' @param stop.on.error Stop on error?
##' @param keep.xy Not used.
##' @param aux Not used.
##' @examples
##' ## Example with many sub-divisions
##' f <- function(x) sin(exp(x))
##' F <- MakeTape(function(x) integrate(f, 0, x)$value, 0)
##' F(7)
##' integrate(f, 0, 7)
##' ## Example with singularity
##' f <- function(x) dbeta(x, shape1=.1, shape2=.1)
##' F <- MakeTape(function(x) integrate(f, 0, x)$value, 0)
##' F(.1)
##' integrate(f, 0, .1)
setMethod("integrate",
          c(f="ANY"),
          function(f, lower, upper, ..., subdivisions, rel.tol, abs.tol, stop.on.error) {
  if (!ad_context()) {
    stats::integrate(f=f, lower=lower, upper=upper, ..., subdivisions=subdivisions, rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=stop.on.error)
  } else {
    f <- dotfun(f, ...)
    cfg <- list(subdivisions, rel.tol, abs.tol, stop.on.error)
    names(cfg) <- c("subdivisions", "rel.tol", "abs.tol", "stop.on.error")
    ans <- as.list(adintegrate(f=f, a=lower, b=upper, cfg=cfg))
    names(ans) <- c("value", "abs.error", "subdivisions")
    ans
  }
})
