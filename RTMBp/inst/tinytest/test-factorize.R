## Test matrix factorizations

library(RTMBp)

################################################################################
## General eigen decomposition
################################################################################

## NOTE: Remember that eigen decomposition is not unique. We test on
## some unique derived matrix quantities instead.

matlog <- function(x) {
    e <- eigen(x, symmetric=FALSE)
    e$vec %*% (log(e$val) * solve(e$vec))
}
matexp <- function(x) {
    e <- eigen(x, symmetric=FALSE)
    e$vec %*% (exp(e$val) * solve(e$vec))
}
f <- function(x) Re(matlog(matexp(x)))

## Test matrix (Not symmetric - eigen values and vectors are complex)
m <- structure(c(-0.626, 0.184, -0.836, 1.595, 0.33, -0.82, 0.487,
                 0.738, 0.576, -0.305, 1.512, 0.39, -0.621, -2.215,
                 1.125, -0.045, -0.016, 0.944, 0.821, 0.594, 0.919,
                 0.782, 0.075, -1.989, 0.62 ),
               dim = c(5L, 5L))

F <- MakeTape(f, m)
expect_equal(f(m), F(m), info="general real eigen")
expect_equal(F$jacobian(m), diag(25), info="general real eigen derivatives")

f <- function(x) Re(matlog(matexp(x+AD(2i)*x)))
F <- MakeTape(f, m)
expect_equal(f(m), F(m), info="general complex eigen")
expect_equal(F$jacobian(m), diag(25), info="general complex eigen derivatives")

################################################################################
## Symmetric eigen decomposition (lower triangle access)
################################################################################

matlog <- function(x) {
    e <- eigen(x, symmetric=TRUE)
    e$vec %*% (log(e$val) * t(e$vec))
}
f <- function(x) Re(matlog(x))
x <- m %*% t(m)
F <- MakeTape(f, x)
expect_equal(f(x), F(x), info="symmetric real eigen")
ndf <- numDeriv::jacobian(F, x)
J <- F$jacobian(x)
## NOTE: Observed to pass with default tol - but it's tight
expect_equal(J, ndf, info="symmetric real eigen derivatives", tol=1e-6)

## complex case (Hermitian)
matlog <- function(x) {
    e <- eigen(x, symmetric=TRUE)
    e$vec %*% (log(e$val) * Conj(t(e$vec)))
}
f <- function(x) Re(matlog(x+(1-diag(5))*AD(1i)+diag(5)*2))
F <- MakeTape(f, x)
expect_equal(f(x), F(x), info="symmetric complex eigen")
ndf <- numDeriv::jacobian(F, x)
J <- F$jacobian(x)
## NOTE: Observed to pass with default tol - but it's tight
expect_equal(J, ndf, info="symmetric complex eigen derivatives", tol=1e-6)

################################################################################
## Cholesky decomposition
################################################################################

x <- m %*% t(m)
F <- MakeTape(chol, x)
expect_equal(chol(x), F(x), info="chol")
ndchol <- numDeriv::jacobian(chol, x)
J <- F$jacobian(x)
## NOTE: Observed to pass with default tol - but it's tight
expect_equal(J, ndchol, info="chol derivatives", tol=1e-6)
Jfun <- F$jacfun()
expect_equal(J, Jfun(x), info="chol derivative replay")
