## Test interoperability with Matrix package
library(Matrix)
library(RTMB)
tol <- sqrt(.Machine$double.eps)

################################################################################
## Test 1 (sparse matrices)
################################################################################
set.seed(1)
n <- 100
M <- lapply(1:2, function(i)rsparsematrix(n, n, nnz = 10*n))
v1 <- rnorm(n)
f <- function(x) {
    ## Test linear combination of adsparse
    Y <- x[1]*M[[1]] + x[2]*M[[2]]
    ## Test diag, diag<-
    diag(Y) <- diag(Y) + 1.2
    ## Test [, [<-
    Y[1,] <- Y[1,] + Y[2,]
    ## Test adsparse %*% vector
    v <- Y %*% v1
    ## Test matrix transpose
    Y <- t(Y)
    ## Test opposite order
    v <- v + as.vector(t(v1) %*% Y)
    ## Test sparseMatrix %*% advector
    v <- M[[1]] %*% v
    ## return sum square
    sum(v*v)
}
x <- 1:2
F <- RTMB::MakeTape(f,x)
expect_true(abs(F(x)-f(x)) < tol)
expect_true(abs(F(x+1)-f(x+1)) < tol)
expect_true(abs(F(x-1)-f(x-1)) < tol)

################################################################################
## Test 2 (triplet sparse matrices)
################################################################################
M <- lapply(M, as, "TsparseMatrix")
F <- RTMB::MakeTape(f,x)
expect_true(abs(F(x)-f(x)) < tol)
expect_true(abs(F(x+1)-f(x+1)) < tol)
expect_true(abs(F(x-1)-f(x-1)) < tol)

################################################################################
## Test 3 (dense matrices)
################################################################################
M <- lapply(M, as.matrix)
F <- RTMB::MakeTape(f,x)
expect_true(abs(F(x)-f(x)) < tol)
expect_true(abs(F(x+1)-f(x+1)) < tol)
expect_true(abs(F(x-1)-f(x-1)) < tol)
