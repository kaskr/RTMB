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
##' f <- interpol2Dfun(volcano, xlim=c(0,1), ylim=c(0,1))
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

setGeneric("splinefun")
##' @describeIn Interpolation Construct a spline function.
##' @param x spline x coordinates
##' @param y spline y coordinates
##' @param method Same as for the stats version, however only the three first are available.
setMethod("splinefun", signature(y="advector",
                                 ties="missing"),
          function(x, y, method=c("fmm", "periodic", "natural")) {
              if (!is.numeric(x)) stop("'x' must be numeric")
              i <- order(x)
              x <- x[i]
              y <- y[i]
              method <- match.arg(method)
              ## C code choices =>
              ##   case 1: periodic_spline
              ##   case 2: natural_spline
              ##   case 3: fmm_spline
              iMeth <- match(method, c("periodic", "natural", "fmm"))
              ptr <- splineptr(x, y, iMeth)
              function(x) {
                  if (inherits(x, "advector"))
                      stop("AD spline must have constant (numeric) input")
                  splineptr_eval(ptr, x)
              }
          })
setMethod("splinefun", signature(x="advector",
                                 y="missing",
                                 ties="missing"),
          function(x, method=c("fmm", "periodic", "natural")) {
              y <- x
              x <- seq_along(y)
              splinefun(x, y, method)
          })
