library(RTMB)

## Read data
source("spatial_data.R")

## Euclidian distance matrix
dd <- as.matrix(dist(Z))

parameters <- list(
    b = c(0, 0),
    a = 1.428571,
    log_sigma = -0.6931472,
    u = rep(0,n) )

f <- function(b, log_sigma, u, a) {
    eta <- X %*% b + exp(log_sigma) * u
    cov <- exp( - a * dd )
    res <- -dmvnorm(u, 0, cov, log=TRUE)
    res <- res - sum(y * eta - exp(eta))
    res
}

do.call("f", parameters)

obj <- MakeADFun(function(p)do.call("f", p), parameters, random = "u")

system.time(
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  lower=c(-100.0, -100.0, 0.01, -3.0),
                  upper=c( 100.0,  100.0, 3.00,  3.0) )
)

rep <- TMB::sdreport(obj)
rep
