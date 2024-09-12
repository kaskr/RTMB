##' AD adjoint code from R
##'
##' Writing custom AD adjoint derivatives from R
##'
##' @details Reverse mode derivatives (adjoint code) can be implemented from R using the function `ADjoint`. It takes as input a function of a single argument `f(x)` representing the function value, and another function of *three* arguments `df(x, y, dy)` representing the adjoint derivative wrt `x` defined as `d/dx sum( f(x) * dy )`. Both `y` and `dy` have the same length as `f(x)`. The argument `y` can be assumed equal to `f(x)` to avoid recalculation during the reverse pass. It should be assumed that all arguments `x`, `y`, `dy` are vectors without any attributes *except* for dimensions, which are stored on first evaluation. The latter is convenient when implementing matrix functions (see `logdet` example).
##' Higher order derivatives automatically work provided that `df` is composed by functions that `RTMB` already knows how to differentiate.
##' @section Complex case:
##' The argument \code{complex=TRUE} specifies that the functions `f` and `df` are complex differentiable (holomorphic) and that arguments `x`, `y` and `dy` should be assumed complex (or \link{adcomplex}). Recall that complex differentiability is a strong condition excluding many continuous functions e.g. `Re`, `Im`, `Conj` (see example).
##' @note `ADjoint` may be useful when you need a special atomic function which is not yet available in `RTMB`, or just to experiment with reverse mode derivatives.
##' However, the approach may cause a *significant overhead* compared to native `RTMB` derivatives. In addition, the approach is *not thread safe*, i.e. calling R functions cannot be done in parallel using OpenMP.
##' @param f R function representing the function value.
##' @param df R function representing the reverse mode derivative.
##' @param name Internal name of this atomic.
##' @param complex Logical; Assume complex and \link{adcomplex} types for all arguments?
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
##'    function(x) determinant(x, log=TRUE)$modulus,
##'    function(x, y, dy) t(solve(x)) * dy,
##'    name = "logdet")
##' (F <- MakeTape(logdet, diag(2)))
##' ## Test derivatives
##' ## Compare with numDeriv::hessian(F, matrix(1:4,2))
##' F$jacfun()$jacobian(matrix(1:4,2)) ## Hessian
##' ############################################################################
##' ## Holomorphic extension of 'solve'
##' matinv <- ADjoint(
##'    solve,
##'    function(x,y,dy) -t(y) %*% dy %*% t(y),
##'    complex=TRUE)
##' (F <- MakeTape(function(x) Im(matinv(x+AD(1i))), diag(2)))
##' ## Test derivatives
##' ## Compare with numDeriv::jacobian(F, matrix(1:4,2))
##' F$jacobian(matrix(1:4,2))
##' @rdname ADjoint
##' @name ADjoint
##' @return A function that allows for numeric and taped evaluation.
ADjoint <- function(f, df, name=NULL, complex=FALSE) {
    if (is.null(name)) {
        name <- substring(deparse(substitute(f))[1], 1, 7)
    }
    if (complex) {
        f <- cplxify1(f)
        df <- cplxify3(df)
    }
    F <- function(x) {
        f(x)
    }
    attr(F, "reverse") <- df
    attr(F, "name") <- name
    function(x) {
        if (complex)
            x <- unsplit(x)
        ans <-
            if (ad_context())
                TapedEval(F, x)
            else f(x)
        if (complex)
            ans <- resplit(ans)
        ans
    }
}
## Helpers
cplxify1 <- function(f) {
    f ## Lazy eval
    function(x) unsplit(f(resplit(x)))
}
cplxify3 <- function(df) {
    df ## Lazy eval
    function(x, y, dy) unsplit(Conj(df(resplit(x),
                                       resplit(y),
                                       Conj(resplit(dy)))))
}
