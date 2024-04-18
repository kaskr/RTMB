##' AD adjoint code from R
##'
##' Writing custom AD adjoint derivatives from R
##'
##' @details Reverse mode derivatives (adjoint code) can be implemented from R using the function `ADjoint`. It takes as input a function of a single argument `f(x)` representing the function value, and another function of *three* arguments `df(x, y, dy)` representing the adjoint derivative wrt `x` defined as `d/dx sum( f(x) * dy )`. Both `y` and `dy` have the same length as `f(x)`. The argument `y` can be assumed equal to `f(x)` to avoid recalculation during the reverse pass. It should be assumed that all arguments `x`, `y`, `dy` are vectors without any attributes. In case of matrix functions, the argument dimensions therefore have to be recovered from the lengths (see `logdet` example).
##' Higher order derivatives automatically work provided that `df` is composed by functions that `RTMB` already knows how to differentiate.
##' @note `ADjoint` may be useful when you need a special atomic function which is not yet available in `RTMB`, or just to experiment with reverse mode derivatives.
##' However, the approach may cause a *significant overhead* compared to native `RTMB` derivatives. In addition, the approach is *not thread safe*, i.e. calling R functions cannot be done in parallel using OpenMP.
##' @param f R function representing the function value.
##' @param df R function representing the reverse mode derivative.
##' @param name Internal name of this atomic.
##' @examples
##' ############################################################################
##' ## Lambert W-function defined by W(y*exp(y))=y
##' W <- function(x) {
##'   logx <- log(x)
##'   y <- pmax(logx, 0)
##'   while (any(abs(logx - log(y) - y) > 1e-9, na.rm = TRUE)) {
##'       y <- y - (y - exp(logx - y)) / (1 + y)
##'   }
##'   y
##' }
##' ## Derivatives
##' dW <- function(x, y, dy) {
##'    dy / (exp(y) * (1. + y))
##' }
##' ## Define new derivative symbol
##' LamW <- ADjoint(W, dW)
##' ## Test derivatives
##' (F <- MakeTape(function(x)sum(LamW(x)), numeric(3)))
##' F(1:3)
##' F$print()                ## Note the 'name'
##' F$jacobian(1:3)          ## gradient
##' F$jacfun()$jacobian(1:3) ## hessian
##' ############################################################################
##' ## Log determinant
##' logdet <- ADjoint(
##'    function(x) {
##'        dim(x) <- rep(sqrt(length(x)), 2)
##'        determinant(x, log=TRUE)$modulus
##'    },
##'    function(x, y, dy) {
##'        dim(x) <- rep(sqrt(length(x)), 2)
##'        t(solve(x)) * dy
##'    },
##'    name = "logdet")
##' MakeTape(logdet, diag(2))
##' @rdname ADjoint
##' @name ADjoint
##' @return A function that allows for numeric and taped evaluation.
ADjoint <- function(f, df, name=NULL) {
    if (is.null(name)) {
        name <- substring(deparse(substitute(f))[1], 1, 7)
    }
    F <- function(x) {
        f(x)
    }
    attr(F, "reverse") <- df
    attr(F, "name") <- name
    function(x) {
        if (ad_context())
            TapedEval(F, x)
        else f(x)
    }
}
