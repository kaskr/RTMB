## Test that selected distributions work as expected
library(RTMB)
tol <- sqrt(.Machine$double.eps)
formals(MakeADFun)$silent <- TRUE

################################################################################
## Test 1 (dmultinom)
################################################################################
set.seed(1) ## For checkConsistency
dat <- list(x=c(1:10,10:1))
f <- function(parms) {
    getAll(dat, parms, warn=FALSE)
    x <- OBS(x)
    prob <- sin( p * (1:10) ) + 1
    prob <- prob / sum(prob)
    ans <- -dmultinom(x[1:10], prob=prob, log=TRUE)
    ans <- ans - dmultinom(x[-(1:10)], prob=prob, log=TRUE)
    ans
}
###############################################
parms <- list(p=.2)
obj <- MakeADFun(f, parms)
expect_true(abs(obj$fn() - f(parms)) < tol)
s <- obj$simulate()
expect_true(length(s$x) == length(dat$x))
res1 <- oneStepPredict(obj,method="cdf",discrete=TRUE,trace=FALSE)
res2 <- oneStepPredict(obj,method="oneStepGeneric",discrete=TRUE,discreteSupport=0:55,trace=FALSE)
expect_true(max(abs(res1$res-res2$res)) < tol)
chk.dmultinom <- checkConsistency(obj)
expect_true(abs(summary(chk.dmultinom)$joint$p.value)>.05)
expect_true(abs(summary(chk.dmultinom)$joint$bias)<.05)

################################################################################
## Test 2 (dpois)
################################################################################
set.seed(1) ## For checkConsistency
dat <- list(x=c(1:10))
f <- function(parms) {
    getAll(parms, dat, warn=FALSE)
    x <- OBS(x)
    -sum(dpois(x, p, log=TRUE))
}
parms <- list(p=.2)
obj <- MakeADFun(f, parms)
expect_true(abs(obj$fn() - f(parms)) < tol)
s <- obj$simulate()
expect_true(length(s$x) == length(dat$x))
chk.dpois <- checkConsistency(obj)
expect_true(abs(summary(chk.dpois)$joint$p.value)>.05)
expect_true(abs(summary(chk.dpois)$joint$bias)<.05)

################################################################################
## Test 3 (dseparable)
################################################################################
parms <- list(s=rep(1,3), u=array(0, c(3, 4, 5)))
C1 <- diag(3)+1
C2 <- diag(4)+2
C3 <- diag(5)+3
f <- function(parms) {
    getAll(parms)
    f1 <- function(x) dmvnorm(x, Sigma = s[1] * C1, log=TRUE)
    f2 <- function(x) dmvnorm(x, Sigma = s[2] * C2, log=TRUE)
    f3 <- function(x) dmvnorm(x, Sigma = s[3] * C3, log=TRUE)
    f123 <- dseparable(f1,f2,f3)
    -f123(u, log=TRUE)
}
C <- C3 %x% C2 %x% C1
expect_equal( f(parms) , -dmvnorm(parms$u, Sigma=C, log=TRUE),
             info="Kronecker covariance")
obj <- MakeADFun(f, parms)
expect_equal(as.double(obj$fn()), 102.11397226694,
             info="Taped kronecker covariance")
obj <- MakeADFun(f, parms, random="u")
expect_equal(as.double(obj$fn()), 0,
             info="Separable density integrates to one")
expect_equal(as.double(obj$gr()), c(0,0,0),
             info="Separable density integral independent of parameters")

################################################################################
## Test 4 (dgamma rate argument)
################################################################################
set.seed(1) ## For checkConsistency
dat <- list(x=1:100)
f <- function(parms) {
    getAll(parms, dat, warn=FALSE)
    x <- OBS(x)
    -sum(dgamma(x, shape, rate, log=TRUE))
}
parms <- list(shape=3, rate=.5)
obj <- MakeADFun(f, parms)
expect_true(abs(obj$fn() - f(parms)) < tol)
s <- obj$simulate()
expect_true(length(s$x) == length(dat$x))
chk.dgamma <- checkConsistency(obj)
expect_true(abs(summary(chk.dgamma)$joint$p.value)>.05)
expect_true(all(abs(summary(chk.dgamma)$joint$bias)<.05))
dat <- s
obj <- MakeADFun(f, parms)
osa <- oneStepPredict(obj, method="cdf", trace=FALSE)
expect_true(ks.test(osa$res, "pnorm")$p.value > .05)
