##' Distributional assignment operator
##'
##' @details Provides a slightly simplified syntax *inspired by*, but *not* compatible with, other probabilistic programming languages (e.g. BUGS/JAGS):
##'
##' - \code{x %~% distribution(...)} is syntactic sugar for \code{.nll <- .nll - sum(distribution(x,...,log=TRUE))}
##' - The variable \code{.nll} is automatically initialized to \code{0} and returned on exit.
##'
##' @note If the shorter name `~` is preferred, it can be locally overloaded using \code{"~" <- RTMB::"%~%"}.
##' @examples
##' f <- function(parms) {
##'   getAll(parms)
##'   x %~% dnorm(mu, 1)
##'   y %~% dpois(exp(x))
##' }
##' p <- list(mu=0, x=numeric(10))
##' y <- 1:10
##' obj <- MakeADFun(f, p, random="x")
##' @param x LHS; Random effect or data for which distribution assignment applies
##' @param distr RHS; Distribution expression
"%~%" <- function(x, distr) {
    pf <- parent.frame()
    if (is.null(pf$.nll))
        delayedAssign(".nll", {on.exit(return(.nll)); 0}, pf, pf)
    pf$.nll
    distr <- substitute(distr)
    distr <- as.list(distr)
    x <- substitute(x)
    distr <- c(distr[1], x, distr[-1])
    distr$log <- TRUE
    ## Have sum-optimized version?
    sumdistr <- as.character(.sum_optimized_table[as.character(distr[[1]])])
    if (is.na(sumdistr)) {
        cl <- as.call(distr)
        cl <- substitute(.nll <- .nll - sum(x), list(x=cl))
    }
    else {
        distr[[1]] <- substitute(":::"("RTMB", x), list(x=sumdistr))
        cl <- as.call(distr)
        cl <- substitute(.nll <- .nll - x, list(x=cl))
    }
    ## print(cl)
    eval.parent(cl)
}

.sum_optimized_table <- c("dnorm" = "sumdnorm")
