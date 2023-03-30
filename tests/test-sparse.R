## Test interoperability with Matrix package
library(Matrix)
library(RTMB)
verbose <- FALSE
print <- function(...)if(verbose)base::print(...)
expect <- function(x) {if(!x)warning("FAIL"); cat(if(x)"OK"else"***NOT*** OK","\n")}
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
    print(class(Y))
    ## Test diag, diag<-
    diag(Y) <- diag(Y) + 1.2
    ## Test [, [<-
    Y[1,] <- Y[1,] + Y[2,]
    ## Test adsparse %*% vector
    v <- Y %*% v1
    print(class(v))
    ## Test matrix transpose
    Y <- t(Y)
    print(class(Y))
    ## Test opposite order
    v <- v + as.vector(t(v1) %*% Y)
    print(class(v))
    ## Test sparseMatrix %*% advector
    v <- M[[1]] %*% v
    print(class(v))
    ## return sum square
    sum(v*v)
}
x <- 1:2
F <- RTMB::MakeTape(f,x)
expect(abs(F(x)-f(x)) < tol)
expect(abs(F(x+1)-f(x+1)) < tol)
expect(abs(F(x-1)-f(x-1)) < tol)

################################################################################
## Test 2 (triplet sparse matrices)
################################################################################
M <- lapply(M, as, "TsparseMatrix")
F <- RTMB::MakeTape(f,x)
expect(abs(F(x)-f(x)) < tol)
expect(abs(F(x+1)-f(x+1)) < tol)
expect(abs(F(x-1)-f(x-1)) < tol)

################################################################################
## Test 3 (dense matrices)
################################################################################
M <- lapply(M, as.matrix)
F <- RTMB::MakeTape(f,x)
expect(abs(F(x)-f(x)) < tol)
expect(abs(F(x+1)-f(x+1)) < tol)
expect(abs(F(x-1)-f(x-1)) < tol)
