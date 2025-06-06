##' Interpolation
##'
##' Some interpolation methods are available to be used as part of 'RTMB' objective functions.
##'
##' `interpol1Dfun` and `interpol2Dfun` are kernel smoothers useful in the case where you need a 3rd order \emph{smooth} representation of a \emph{data} vector or matrix.
##' A typical use case is when a high-resolution map needs to be accessed along a random effect trajectory.
##' Both 1D and 2D cases accept an 'interpolation radius' parameter (default R=2) controlling the degree of smoothness. Note, that only the value R=1 will match the data exactly, while higher radius trades accuracy for smoothness. Note also that these smoothers do not attempt to extrapolate: The returned value will be `NaN` outside the valid range (`xlim` / `ylim`).
##'
##' `splinefun` imitates the corresponding `stats` function. The AD implementation (in contrast to `interpol1Dfun`) works for parameter dependent y-coordinates.
##'
##' @rdname Interpolation
##' @name Interpolation
##' @param xlim Domain of x
##' @param ylim Domain of y
##' @param z Matrix to be interpolated
##' @param ... Configuration parameters
##' @examples
##' ## ======= interpol1D
##' ## R=1 => exact match of observations
##' f <- interpol1Dfun(sin(1:10), R=1)
##' layout(t(1:2))
##' plot(sin(1:10))
##' plot(f, 1, 10, add=TRUE)
##' title("R=1")
##' F <- MakeTape(f, 0)
##' F3 <- F$jacfun()$jacfun()$jacfun()
##' plot(Vectorize(F3), 1, 10)
##' title("3rd derivative")
##' ## ======= interpol2D
##' ## R=1 => exact match of observations
##' f <- interpol2Dfun(volcano, xlim=c(0,1), ylim=c(0,1), R=1)
##' f(0,0) == volcano[1,1]   ## Top-left corner
##' f(1,1) == volcano[87,61] ## Bottom-right corner
##' ## R=2 => trades accuracy for smoothness
##' f <- interpol2Dfun(volcano, xlim=c(0,1), ylim=c(0,1), R=2)
##' f(0,0) - volcano[1,1]    ## Error Top-left corner
##' F <- MakeTape(function(x) f(x[1],x[2]), c(.5,.5))
##' ## ======= splinefun
##' T <- MakeTape(function(x){
##'    S <- splinefun(sin(x))
##'    S(4:6)
##' }, 1:10)
NULL

##' @describeIn Interpolation Construct a kernel smoothed representation of a vector.
##' @return function of x.
interpol1Dfun <- function(z, xlim=c(1,length(z)), ...) {
    stopifnot(is.null(dim(z)))
    f2 <- interpol2Dfun(as.matrix(z), xlim, c(-.5,.5), ...)
    function(x) {
        f2(x, 0)
    }
}

##' @describeIn Interpolation Construct a kernel smoothed representation of a matrix.
##' @return function of x and y.
interpol2Dfun <- function(z, xlim=c(1,nrow(z)), ylim=c(1,ncol(z)), ...) {
    con <- list(...)
    if (is.null(con$R))
        con$R <- 2
    ptr <- ip2D(z, xlim, ylim, con)
    function(x, y) {
        if (inherits(x, "advector") || inherits(y, "advector")) {
            ip2D_eval_ad(ptr, advector(x), advector(y))
        } else {
            ip2D_eval_num(ptr, x, y)
        }
    }
}

## Modified stats::splinefun to handle the fixed spline case with possibly AD input.
splinefun_stats <- function(x., y., method) {
    s_stats <- stats::splinefun(x., y., method)
    s_tmb <- NULL
    function(x, deriv = 0L) {
        if (!inherits(x, "advector")) {
            s_stats(x, deriv)
        } else {
            if (is.null(s_tmb)) {
                s_tmb <- splinefun(advector(x.),
                                   advector(y.),
                                   method)
            }
            s_tmb(x, deriv)
        }
    }
}

setGeneric("splinefun")
##' @describeIn Interpolation Construct a spline function.
##' @param x spline x coordinates
##' @param y spline y coordinates
##' @param method Same as for the stats version, however only the three first are available.
setMethod("splinefun", signature(x="ad",
                                 y="ad",
                                 ties="missing"),
          function(x, y, method=c("fmm", "periodic", "natural")) {
              method <- match.arg(method)
              if (!inherits(x, "advector") &&
                  !inherits(y, "advector") &&
                  !ad_context()) {
                  return (splinefun_stats(x, y, method=method))
              }
              x <- advector(x)
              y <- advector(y)
              ## C code choices =>
              ##   case 1: periodic_spline
              ##   case 2: natural_spline
              ##   case 3: fmm_spline
              iMeth <- match(method, c("periodic", "natural", "fmm"))
              ptr <- splineptr(x, y, iMeth)
              S <- function(x) {
                  x <- advector(x)
                  splineptr_eval(ptr, x)
              }
              function(x, deriv=0L) {
                  if (deriv > 0) {
                      S <- MakeTape(S, numeric(length(x)))
                      for (i in seq_len(deriv)) {
                          S <- MakeTape(function(x) sum(S(x)), S$par())
                          S <- S$jacfun()
                          S$simplify()
                      }
                  }
                  S(x)
              }
          })
##' @describeIn Interpolation Construct a spline function.
setMethod("splinefun", signature(x="ad",
                                 y="missing",
                                 ties="missing"),
          function(x, method=c("fmm", "periodic", "natural")) {
              y <- x
              x <- seq_along(y)
              splinefun(x, y, method)
          })
