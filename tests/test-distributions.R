## Test that selected distributions work as expected
library(RTMB)
verbose <- FALSE
print <- function(...)if(verbose)base::print(...)
expect <- function(x) {if(!x)warning("FAIL"); cat(if(x)"OK"else"***NOT*** OK","\n")}
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
expect(abs(obj$fn() - f(parms)) < tol)
s <- obj$simulate()
expect(length(s$x) == length(dat$x))
res1 <- oneStepPredict(obj,method="cdf",discrete=TRUE,trace=FALSE)
res2 <- oneStepPredict(obj,method="oneStepGeneric",discrete=TRUE,discreteSupport=0:55,trace=FALSE)
expect(max(abs(res1$res-res2$res)) < tol)
chk.dmultinom <- checkConsistency(obj)
expect(abs(summary(chk.dmultinom)$joint$p.value)>.05)
expect(abs(summary(chk.dmultinom)$joint$bias)<.05)

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
expect(abs(obj$fn() - f(parms)) < tol)
s <- obj$simulate()
expect(length(s$x) == length(dat$x))
chk.dpois <- checkConsistency(obj)
expect(abs(summary(chk.dpois)$joint$p.value)>.05)
expect(abs(summary(chk.dpois)$joint$bias)<.05)
