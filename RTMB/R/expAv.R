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
##' The uniformization method essentially pushes the spectrum of the operator inside a zero centered disc, within which a uniform error bound is available.
##' If `A` is a generator matrix (i.e. `expm(A)` is a probability matrix) and if `v` is a probability vector, then the relative error of the result is bounded by `tol`.
##'
##' Additional supported arguments via `...` currently include:
##'
##' - `Nmax` Use no more than this number of terms even if the spcified accuracy cannot be met.
##' - `warn` Give warning if number of terms is truncated by `Nmax`.
##' - `trace` Trace the number of terms when it adaptively changes.
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
##' F <- MakeTape(function(x) expAv(A*x, 1:5), 1)
##' F(1)
##' F(2) ## More terms needed => trigger retaping
expAv <- function(A, v, transpose=FALSE, uniformization=TRUE, tol=1e-8, ..., cache=A) {
    force(cache)
    cfg <- list(Nmax=100, warn=FALSE, trace=TRUE)
    dotargs <- list(...)
    cfg[names(dotargs)] <- dotargs
    if (ad_context()) A <- as(A, "adsparse")
    N <- cfg$Nmax ## Default
    if (uniformization) {
        ## Template of A
        .A <- new("dgCMatrix",
                  x=numeric(length(A@x)),
                  i=A@i,
                  p=A@p,
                  Dim=A@Dim
                  )
        disc <- DataEval(function(x){.A@x[] <- x; eigenDisc(.A)}, A@x)
        diag(A) <- diag(A) - disc["C"]
        getN <- function(rho) qpois(tol, rho, lower.tail=FALSE)
        N <- DataEval( getN, disc["R"] )
    }
    v <- as.matrix(v)
    if (ad_context()) {
        v <- advector(v)
        N <- advector(N)
        if (!transpose) A <- t(A)
        ans <- expATv(A, v, N, cfg, cache) ## A transposed internally
    } else {
        if (transpose) A <- Matrix::t(A)
        ans <- term <- v
        for (n in seq_len(N)) {
            term = A %*% term / n
            ans <- ans + term
        }
    }
    if (uniformization) ans <- exp(disc["C"]) * ans
    drop(as.matrix(ans))
}
