## Simulate data
set.seed(123)
n <- 1000
mu <- 3
sigma <- 1.5
y <- rnorm(n, mu, sigma) ## data

## data and parameters
parameters <- list(s=numeric(n), mu=0, sigma=1)

f <- function(parameters) {
    ## Build inner problem
    s <- parameters$s
    mu <- parameters$mu
    sigma <- parameters$sigma
    K <- sum(mu * s + .5 * s * s * sigma^2) - sum(s * y)
    K
}
library(RTMB)
F <- MakeTape(f, parameters)
F <- F$laplace(random=1:n, SPA=TRUE)

## Tape can be used in other objective functions:
obj <- MakeADFun(function(pl)F(pl[[1]]), list(c(mu=0,sigma=1)) )
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdreport(obj)
