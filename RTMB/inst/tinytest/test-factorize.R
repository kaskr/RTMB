## Test matrix factorizations

library(RTMB)

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

## Defective and repeated-eigenvalue matrices: NaN rather than an R error
## or a silently wrong derivative.
## 3-state chain generator; defective (Jordan block) when the rates are equal
Qchain <- function(x) {
    "[<-" <- ADoverload("[<-")
    r <- exp(x)
    Q <- matrix(0, 3, 3)
    Q[1, 2] <- r[1]; Q[2, 3] <- r[2]
    Q[1, 1] <- -r[1]; Q[2, 2] <- -r[2]
    Q
}
f <- function(x) Re(matexp(Qchain(x))[1, 3])
F <- MakeTape(f, c(0, 1))
expect_equal(F(c(0, 1)), f(c(0, 1)), info="defective eigen: value away from tie")
expect_true(is.nan(F(c(0, 0))), info="defective eigen: NaN value, no error")
expect_true(all(is.nan(F$jacobian(c(0, 0)))), info="defective eigen: NaN derivatives")
## 3-state generator with all rates equal: eigenvalues (0, -3, -3) with a
## full set of eigenvectors
Q3 <- function(x) {
    "[<-" <- ADoverload("[<-")
    Q <- matrix(0, 3, 3)
    off <- which(row(Q) != col(Q))
    Q[off] <- exp(x)
    for (i in 1:3) Q[i, i] <- -sum(Q[i, ])
    Q
}
f <- function(x) Re(matexp(Q3(x))[1, 2])
F <- MakeTape(f, rep(0, 6))
expect_equal(F(rep(0, 6)), f(rep(0, 6)), info="repeated eigenvalues: value")
expect_true(all(is.nan(F$jacobian(rep(0, 6)))),
            info="repeated eigenvalues: NaN derivatives")
x <- c(0, 0.3, -0.2, 0.1, 0.4, -0.3)
h <- 1e-6
nd <- sapply(1:6, function(i) {
    e <- replace(rep(0, 6), i, h)
    (f(x + e) - f(x - e)) / (2 * h)
})
expect_equal(as.vector(F$jacobian(x)), nd, tol=1e-6,
             info="distinct eigenvalues: derivatives")

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
if (at_home()) {
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
}

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
