library(RTMB)
tol <- sqrt(.Machine$double.eps)

################################################################################
## Test 1 (GH issue 13)
################################################################################

f <- MakeTape(function(x) c(x[1],x[1]+x[2]+3*x[3],x[1]+x[2]+x[3],x[1]*x[2]),numeric(3))
J123 <- structure(c(1, 1, 1, 2, 0, 1, 1, 1, 0, 3, 1, 0), dim = 4:3)
expect_identical(f$jacobian(1:3), J123, info="Jacobian evaluation")
expect_identical(f$jacfun()(1:3), J123, info="Jacobian transformation")

################################################################################
## Test 2 (NA propagation on the tape)
################################################################################

## Arithmetic with NA/NaN constants
## Compiler optimizations (--fast-math) can cause this test to fail:
vec <- c(NA, NaN, 1)
f <- function(x) x + vec
F <- MakeTape(f, numeric(1))
expect_identical(F(0), vec, info="NA propagation on the tape")

## Tapes f1 and f2 should be identical assuming correct propagation of
## NAs from 'obs' to 'nll':
obs <- c(2, NA, NA, NA)
f1 <- MakeTape(function(p){
    nll <- -dnorm(obs, p, 1, TRUE)
    sum(nll, na.rm=TRUE)
}, 0)
f2 <- MakeTape(function(p){
    nll <- -dnorm(obs, p, 1, TRUE)
    sum(nll[!is.na(obs)])
}, 0)
expect_identical(f1(1.234), f2(1.234), info="NA propagation on the tape")

################################################################################
## (GH issue 17)
################################################################################

expect_silent(
MakeTape(function(x) {
    y <- numeric(3)
    f <- MakeTape(sin, 1:10)
    y[1] <- x[1] ## Error here is side-effect of previous line
    y
}, 1:10)
)

################################################################################
## (GH issue 18)
################################################################################

F <- MakeTape(function(x) {
    G <- MakeTape(function(y) y*x, numeric(5))
    G(1:5)
}, numeric(5))
expect_equal(F(1:5), (1:5)^2, info="https://github.com/kaskr/RTMB/issues/18")
