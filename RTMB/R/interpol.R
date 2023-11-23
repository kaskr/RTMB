##' Interpolation
##'
##' Some interpolation methods are available to be used as part of 'RTMB' objective functions.
##' @rdname Interpolation
##' @name Interpolation
##' @param xlim Domain of x
##' @param ylim Domain of y
##' @param ... Configuration parameters
##' @examples
##' f <- interpol2Dfun(volcano, xlim=c(0,1), ylim=c(0,1))
##' F <- MakeTape(function(x) f(x[1],x[2]), c(.5,.5))
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
