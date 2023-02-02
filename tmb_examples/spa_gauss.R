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
    ## MGF
    K <- sum(mu * s + .5 * s * s * sigma^2)
    ## SPA adjustment
    K - sum(s * y)
}
library(RTMB)
F <- MakeTape(f, parameters)
show(F) ## f
F <- F$laplace(random=1:n, SPA=TRUE)
show(F) ## Laplace is on the tape (verify using F$print())

## Tape can be used directly:
nlminb(c(0, 1), F, F$jacobian)

## Or re-used to build new objective functions
G <- function(x) F(x)

## Feed G into normal TMB interface
obj <- MakeADFun(G, list(mu=0, sigma=1) )
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdreport(obj)
