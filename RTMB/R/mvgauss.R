##' @describeIn MVgauss Multivariate normal distribution. \link{OSA-residuals} can be used for argument \code{x}.
##' @details Multivariate normal density evaluation is done using `dmvnorm()`. This is meant for dense covariance matrices. If *many evaluations* are needed for the *same covariance matrix* please note that you can pass matrix arguments: When `x` is a matrix the density is applied to each row of `x` and the return value will be a vector (length = `nrow(x)`) of densities.
##' @param x Density evaluation point
##' @param Sigma Covariance matrix
##' @param mu Mean parameter vector
##' @param log Logical; Return log density?
dmvnorm <- function(x, mu=0, Sigma, log=FALSE) {
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
dgmrf <- function(x, mu=0, Q, log=FALSE) {
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
dautoreg <- function(x, phi, log=FALSE) {
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
    ans <- dmvnorm(x[1:k], 0, V0, log=TRUE)
    for (i in (tail(seq_along(x), -k))) {
        ans <- ans + dnorm(x[i], sum(phi * x[i - (1:k)]), sigma, log=TRUE)
    }
    if (!log) ans <- exp(ans)
    ans
}

##' @describeIn MVgauss Separable extension of Gaussian log-densities
##' @details Separable extension can be constructed for an unlimited number of inputs. Each input must be a function returning a *gaussian* *mean zero* **log** density. The output of `dseparable` is another **log** density which can be evaluated for array arguments. For example `dseparable(f1,f2,f3)` takes as input a 3D array. `f1` acts in 1st array dimension, `f2' in 2nd dimension and so on.
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
    function(x) {
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
}
