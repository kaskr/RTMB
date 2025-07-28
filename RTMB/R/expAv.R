## Helper: Bound eigen values of real matrix by Disc(C, R)
eigenDisc <- function(A) {
    B <- A
    diag(B) <- 0
    B <- abs(B)
    ## Circle theorem - row version
    center <- diag(A)
    radius <- rowSums(B)
    Min <- min(center - radius)
    Max <- max(center + radius)
    c( C = (Max+Min)/2 , R = (Max-Min)/2 )
}

##' Matrix exponential of sparse matrix multiplied by a vector.
##'
##' Calculates `expm(A) %*% v` using plain series summation. The number of terms is determined adaptively when `uniformization=TRUE`.
##' The uniformization method essentially pushes the spectrum of the operator inside a zero centered disc, within which a tight uniform error bound is available. This effectively reduces the number of terms needed to calculate the series to a given accuracy.
##' If `A` is a generator matrix (i.e. `expm(A)` is a probability matrix) and if `v` is a probability vector, then the relative error of the result is bounded by `tol`.
##' However, note that series summation may be unstable depending on the spectral radius of A. If you get `NaN` values, consider setting `rescale_freq=1` for better stability (see details).
##'
##' Additional supported arguments via `...` currently include:
##'
##' - `Nmax` Integer (`2e9` by default). Use no more than this number of terms even if the specified accuracy cannot be met. When using `expAv` as part of likelihood optimization, you can set a lower value to avoid long unnecessary computation when the optimizer tries extreme parameters. For example, if the spectral radius of `A` cannot realistically exceed some known value `rhomax` one can set `Nmax=qpois(tol, rhomax, lower.tail = FALSE)`.
##' - `warn` Logical (`TRUE` by default). Give warning if number of terms is truncated by `Nmax` (autodiff code only).
##' - `trace` Logical (`FALSE` by default). Trace the number of terms when it adaptively changes (autodiff code only).
##' - `rescale_freq` Integer (`50` by default) controlling how often to rescale for numerical stability. Set to a lower value for more frequent rescaling at the cost of longer computation time. The default value should suffice for a generator matrix of spectral radius up to at least `1e6` (`.Machine$double.xmax^(1/50)`).
##' - `rescale` Logical; Set to `FALSE` to disable rescaling for higher speed. All other values are ignored.
##'
##' @references
##' Grassmann, W. K. (1977). Transient solutions in Markovian queueing systems. \emph{Computers & Operations Research}, 4(1), 47--53.
##'
##' Sherlock, C. (2021). Direct statistical inference for finite Markov jump processes via the matrix exponential. \emph{Computational Statistics}, 36(4), 2863--2887.
##' @param A Sparse matrix (usually a generator)
##' @param v Vector (or matrix)
##' @param transpose Calculate `expm(t(A)) %*% v` ? (faster due to the way sparse matrices are stored)
##' @param uniformization Use uniformization method?
##' @param tol Accuracy if A is a generator matrix and v a probability vector.
##' @param ... Extra configuration parameters
##' @param cache Re-use internal AD calculations by setting an attribute on this object (`A` by default - use NULL to disable caching).
##' @return Vector (or matrix)
##' @examples
##' set.seed(1); A <- Matrix::rsparsematrix(5, 5, .5)
##' expAv(A, 1:5) ## Matrix::expm(A) %*% 1:5
##' F <- MakeTape(function(x) expAv(A*x, 1:5, trace=TRUE), 1)
##' F(1)
##' F(2) ## More terms needed => trigger retaping
expAv <- function(A, v, transpose=FALSE, uniformization=TRUE, tol=1e-8, ..., cache=A) {
    force(cache)
    cfg <- list(Nmax=2e9, warn=TRUE, trace=FALSE, rescale_freq=50)
    dotargs <- list(...)
    cfg[names(dotargs)] <- dotargs
    if (isFALSE(cfg[["rescale"]]))
      cfg[["rescale_freq"]] <- cfg[["Nmax"]] ## never rescale
    if (ad_context()) A <- as(A, "adsparse")
    N <- cfg$Nmax ## Default
    C <- 0 ## Default
    if (uniformization) {
        disc <- eigenDisc(A)
        diag(A) <- diag(A) - disc["C"]
        getN <- function(rho) qpois(tol, rho, lower.tail=FALSE)
        N <- DataEval( getN, disc["R"] )
        C <- disc["C"]
    }
    v <- as.matrix(v)
    if (ad_context()) {
        v <- advector(v)
        N <- advector(N)
        C <- advector(C)
        if (!transpose) A <- t(A)
        ans <- expATv(A, v, N, C, cfg, cache) ## A transposed internally
    } else {
        if (transpose) A <- Matrix::t(A)
        ans <- term <- v
        log_scale <- 0
        if (N > cfg[["Nmax"]]) {
          N <- cfg[["Nmax"]]
          warning(sprintf("Number of terms reduced to 'Nmax' (%d) - result might be inaccurate", cfg[["Nmax"]]))
        }
        for (n in seq_len(N)) {
            term <- A %*% term / n
            ans <- ans + term
            if (!(n %% cfg[["rescale_freq"]])) {
              s <- sum(abs(ans))
              term <- term / s
              ans <- ans / s
              log_scale <- log_scale + log(s)
            }
        }
        ans <- exp(C + log_scale) * ans
    }
    drop(as.matrix(ans))
}
