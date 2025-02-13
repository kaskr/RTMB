lchol <- ADjoint(
    function(x) {
        t(chol(t(x)))
    },
    function(x, y, dy) {
        if (!inherits(x, "advector")) {
            ## Fast lapack based kernel
            S <- sytrisol(y, dy)
            diag(S) <- .5 * diag(S)
            S
        } else {
            ## Slow rewritten version
            T <- t(y) %*% dy
            u <- upper.tri(T, FALSE)
            T[u] <- t(T)[u]
            Linv <- solve(y)
            S <- t(Linv) %*% T %*% Linv
            S[u] <- 0
            diag(S) <- .5 * diag(S)
            S
        }
    },
    name = "lchol")

##' @describeIn ADmatrix AD matrix cholesky
chol.advector <- function(x, ...) {
    t(lchol(t(x)))
}
