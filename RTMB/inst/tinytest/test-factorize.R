## Test matrix factorizations

library(RTMB)

################################################################################
## General eigen decomposition
################################################################################

matlog <- function(x) {
    e <- eigen(x, symmetric=FALSE)
    e$vec %*% (log(e$val) * solve(e$vec))
}
matexp <- function(x) {
    e <- eigen(x, symmetric=FALSE)
    e$vec %*% (exp(e$val) * solve(e$vec))
}
f <- function(x) Re(matlog(matexp(x)))

m <- structure(c(-0.626, 0.184, -0.836, 1.595, 0.33, -0.82, 0.487,  0.738, 0.576, -0.305, 1.512, 0.39, -0.621, -2.215, 1.125, -0.045,  -0.016, 0.944, 0.821, 0.594, 0.919, 0.782, 0.075, -1.989, 0.62 ), dim = c(5L, 5L))

F <- MakeTape(f, m)
expect_equal(f(m), F(m), info="eigen general")
expect_equal(F$jacobian(m), diag(25), info="eigen general derivatives")

################################################################################
## Cholesky decomposition
################################################################################

x <- m %*% t(m)
F <- MakeTape(chol, x)
expect_equal(chol(x), F(x), info="chol")
ndchol <- numDeriv::jacobian(chol, x)
expect_equal(F$jacobian(x), ndchol, info="chol derivatives", tol = 1e-6)
