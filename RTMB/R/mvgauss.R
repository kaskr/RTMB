##' @describeIn MVgauss Multivariate normal distribution. \link{OSA-residuals} can be used for argument \code{x}.
##' @details Multivariate normal density evaluation is done using `dmvnorm()`. This is meant for dense covariance matrices. If *many evaluations* are needed for the *same covariance matrix* please note that you can pass matrix arguments: When `x` is a matrix the density is applied to each row of `x` and the return value will be a vector (length = `nrow(x)`) of densities.
##' @param x Density evaluation point
##' @param Sigma Covariance matrix
##' @param mu Mean parameter vector
##' @param log Logical; Return log density?
##' @param scale Extra scale parameter - see section 'Scaling'.
dmvnorm <- function(x, mu=0, Sigma, log=FALSE, scale=1) {
    if (!unit(scale)) {
        return (dscale("dmvnorm", x, mu, Sigma,
                       log=log, scale=scale, vectorize=TRUE))
    }
    if (inherits(x, "simref")) {
        if (!log) stop("'simref' is for *log* density evaluation only")
        nr <- nrow(as.matrix(Sigma))
        n <- length(x) / nr
        if (length(mu)<nr) mu <- rep(mu, length.out=nr)
        x[] <- MASS::mvrnorm(n, mu, Sigma)
        return( rep(0, n) )
    }
    if (inherits(x, "osa")) {
        keep <- x@keep
        x <- x@x
        dim(keep) <- dim(x)
    } else {
        keep <- NULL
    }
    ## R convention is to have samples by row
    x <- t(as.matrix(x))
    if (!is.null(keep))
        keep <- t(as.matrix(keep))
    mu <- t(as.matrix(mu))
    Sigma <- as.matrix(Sigma)
    d <- nrow(Sigma)
    x0 <- as.vector(x) - as.vector(mu)
    dim(x0) <- c(d, length(x0) / d)
    anstype <- .anstype(x0, Sigma)
    anstype( dmvnorm0(advector(x0), advector(Sigma), log, keep) )
}

##' @describeIn MVgauss Multivariate normal distribution. OSA is \emph{not} implemented.
##' @details The function `dgmrf()` is essentially identical to `dmvnorm()` with the only difference that `dgmrf()` is specified via the *precision* matrix (inverse covariance) assuming that this matrix is *sparse*.
##' @param Q Sparse precision matrix
dgmrf <- function(x, mu=0, Q, log=FALSE, scale=1) {
    if (!unit(scale)) {
        return (dscale("dgmrf", x, mu, Q,
                       log=log, scale=scale, vectorize=TRUE))
    }

    if (!ad_context()) { ## Workaround: see C++ code 'gmrf0'
        F <- .MakeTape(function(...)advector(dgmrf(x,mu,Q,log)),numeric(0))
        return (F$eval(numeric(0)))
    }
    ## R convention is to have samples by row
    x <- t(as.matrix(x))
    mu <- t(as.matrix(mu))
    d <- attr(Q, "Dim")[1]
    x0 <- as.vector(x) - as.vector(mu)
    dim(x0) <- c(d, length(x0) / d)
    anstype <- .anstype(x0, Q)
    anstype( dgmrf0(advector(x0), as(Q, "adsparse"), log) )
}

##' @describeIn MVgauss Gaussian stationary mean zero AR(k) density
##' @details Autoregressive density evaluation is implemented for all orders via `dautoreg()` (including the simplest AR1).
##' We note that this variant is for a *stationary*, *mean zero* and *variance one* process.
##' FIXME: Provide parameterization via partial correlations.
##' @param phi Autoregressive parameters
dautoreg <- function(x, mu=0, phi, log=FALSE, scale=1) {
    "[<-" <- ADoverload("[<-")
    if (!zero(mu) || !unit(scale)) {
        return (dscale("dautoreg", x, 0, phi,
                       log=log, center=mu, scale=scale))
    }
    k <- length(phi)
    M <- matrix(0, k, k)
    for (i in 1:k) {
        for (j in 1:k) {
            d <- abs(i-j)
            if (i != j) {
                M[i, d] <- M[i, d] + phi[j]
            }
        }
    }
    I <- diag(k)
    gamma <- solve(I-M, phi)
    sigma <- sqrt(1-sum(phi*gamma))
    V0 <- diag(k)
    for (i in 1:k) {
        for (j in 1:k) {
            d <- abs(i-j)
            if (i != j){
                V0[i, j] <- gamma[d]
            }
        }
    }
    k <- min(length(x), k)
    V0 <- V0[1:k, 1:k]
    ## Speedup (code would work the same without)
    if (inherits(x, "simref")) {
        xref <- x
        x <- x$value
        x[1:k] <- MASS::mvrnorm(1, 0, V0)
        for (i in (tail(seq_along(x), -k))) {
            x[i] <- rnorm(1, sum(phi * x[i - (1:k)]), sigma)
        }
        xref[] <- x
        return (0)
    }
    ans <- dmvnorm(x[1:k], 0, V0, log=TRUE)
    for (i in (tail(seq_along(x), -k))) {
        ans <- ans + dnorm(x[i], sum(phi * x[i - (1:k)]), sigma, log=TRUE)
    }
    if (!log) ans <- exp(ans)
    ans
}

##' @describeIn MVgauss Separable extension of Gaussian log-densities
##' @details Separable extension can be constructed for an unlimited number of inputs. Each input must be a function returning a *gaussian* *mean zero* **log** density. The output of `dseparable` is another **log** density which can be evaluated for array arguments. For example `dseparable(f1,f2,f3)` takes as input a 3D array `x`. `f1` acts in 1st array dimension of `x`, `f2` in 2nd dimension and so on. In addition to `x`, parameters `mu` and `scale` can be supplied - see below.
##' @section Scaling:
##' All the densities accept a `scale` argument which replaces `SCALE` and `VECSCALE` functionality of TMB.
##' Scaling is applied elementwise on the residual `x-mu`. This works as expected when `scale` is a *scalar* or a *vector* object of the same length as `x`.
##' In addition, `dmvnorm` and `dgmrf` can be scaled by a vector of length equal to the covariance/precision dimension. In this case the `scale` parameter is recycled by row to meet the special row-wise vectorization of these densities.
##' @param ... Log densities
##' @examples
##' func <- function(x, sd, parm, phi) {
##'    ## IID N(0, sd^2)
##'    f1 <- function(x)sum(dnorm(x, sd=sd, log=TRUE))
##'    Sigma <- diag(2) + parm
##'    ## MVNORM(0, Sigma)
##'    f2 <- function(x)dmvnorm(x, Sigma=Sigma, log=TRUE)
##'    ## AR(2) process
##'    f3 <- function(x)dautoreg(x, phi=phi, log=TRUE)
##'    ## Separable extension (implicit log=TRUE)
##'    -dseparable(f1, f2, f3)(x)
##' }
##' parameters <- list(x = array(0, c(10, 2, 10)), sd=2, parm=1, phi=c(.9, -.2))
##' obj <- MakeADFun(function(p)do.call(func, p), parameters, random="x")
##' ## Check that density integrates to 1
##' obj$fn()
##' ## Check that integral is independent of the outer parameters
##' obj$gr()
##' ## Check that we can simulate from this density
##' s <- obj$simulate()
dseparable <- function(...) {
    f <- list(...)
    stopifnot(all(sapply(f, is.function)))
    dsep <- function(x, mu=0, scale=1, log=TRUE) {
        if (!log) stop("'dseparable' is for *log* density evaluation only")
        if (!zero(mu) || !unit(scale)) {
            return (dscale(dsep, x,
                           log=log, center=mu, scale=scale))
        }
        if (inherits(x, "osa")) {
            ok <-
                (length(x@x) == length(x@keep)) &&  ## no CDF method!
                all(!getVariables(x@keep)) &&       ## All 'keep' are constant
                all(getValues(x@keep) == 1)         ## and equal to one.
            if (!ok)
                stop("'osa' (marginalization) is not fully implemented for separable densities")
            x <- x@x
        }
        dosim <- inherits(x, "simref")
        d <- dim(x)
        stopifnot(length(d) == length(f))
        rot <- c(seq_along(d)[-1L], 1L)
        x0 <- x
        ans <- 0
        if (dosim)
            x <- array(rnorm(length(x)), d)
        for (i in seq_along(f)) {
            xmat <- x
            dim(xmat) <- c(nrow(x), length(x) / nrow(x))
            zero <- rep(0, nrow(x))
            J <- MakeTape(f[[i]], zero)$jacfun()
            ans <- ans + ( f[[i]](zero) + d[i] * log(sqrt(2 * pi)) ) * prod(d[-i])
            ## FIXME: Check J linear and J(0)=0
            if (dosim) {
                Q <- -J$jacobian(zero)
                Lt <- chol(Q)
                x[] <- solve(Lt, xmat)
            } else {
                x[] <- -apply(xmat, 2, J)
            }
            x <- aperm(x, rot)
        }
        if (dosim) {
            x0[] <- x
            return (0)
        }
        ans <- ans - .5 * sum(x * x0)
        ans <- ans - length(x) * log(sqrt(2 * pi))
        ans
    }
    dsep
}

## Utility to scale a density:
dscale <- function(f, x, ...,
                   log=FALSE, center=0, scale, vectorize=FALSE) {
    f <- match.fun(f)
    ## Handle the special 'byrow' vectorization offered by 'mvnorm' and 'gmrf'
    if (vectorize) {
        if (length(scale) > 1) {
            if (length(scale) != length(x)) {
                ## weird 'byrow' case
                if (!is.matrix(x))
                    stop("'x' must be a matrix")
                nc <- ncol(x)
                if (length(scale) != nc)
                    stop("Vector 'scale' must be compatible with *rows* of 'x'")
                scale <- matrix(scale, nrow(x), ncol(x), byrow=TRUE)
            } else {
                if (!is.null(dim(x))) {
                    if (!identical(dim(x), dim(scale))) {
                        stop("'dim(scale)' must equal 'dim(x)'")
                    }
                }
            }
        }
    }
    ## Aggregation 'byrow' in vectorized case
    if (vectorize)
        Sum <- rowSums
    else
        Sum <- sum
    ## Check 'center' and 'scale'
    if (length(scale) != 1L)
        if (length(scale) != length(x))
            stop("Vector 'scale' must have same length as 'x'")
    if (length(center) != 1L)
        if (length(center) != length(x))
            stop("Vector 'center' must have same length as 'x'")
    do.center <- !zero(center)
    if (inherits(x, "osa")) {
        if (!log) stop("'osa' is for *log* density evaluation only")
        if (do.center) x@x <- x@x - center
        x@x <- x@x / scale
        keep <- as.matrix(x@keep)[,1] ## Ignore CDF adjustment (it is zero)
        dim(keep) <- dim(x)
        ans <- f(x, ..., log=TRUE) - Sum(keep * log(scale))
        return (ans)
    }
    if (do.center) x <- x - center
    nrep <- length(x) / length(scale)
    ans <- f(x / scale, ..., log=TRUE) - nrep * Sum(log(scale))
    if (log) ans else exp(ans)
}

zero <- function(x) identical(x, 0)
unit <- function(x) identical(x, 1)

##' @describeIn MVgauss Helper to generate an unstructured correlation matrix to use with `dmvnorm`
##' @section Unstructured correlation:
##' Replacement of `UNSTRUCTURED_CORR` functionality of TMB. Constuct object using `us <- unstructured(k)`.
##' Now `us` has two methods: `x <- us$parms()` gives the parameter vector used as input to the objective function, and `us$corr(x)` turns the parameter vector into an unstructured correlation matrix.
##' @param k Dimension
unstructured <- function(k) {
    N <- (k * k - k) / 2
    list(
        parms = function() rep(0, N),
        corr = function(x) {
            "[<-" <- ADoverload("[<-")
            if (length(x) != N) stop("Expected ", N, " parameters")
            L <- diag(k)
            L[lower.tri(L)] <- x
            cov2cor( L %*% t(L) )
        })
}
